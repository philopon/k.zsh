import os
import sys
import re

import requests
import requests_cache

from xml.etree import ElementTree


requests_cache.install_cache(os.path.join(
    os.path.dirname(__file__),
    '..',
    '.' + os.path.splitext(os.path.basename(__file__))[0]
))

nan = float('nan')

AtomicSymbols = [
    'H', 'He',

    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',

    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',

    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',

    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',

    'Cs', 'Ba',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',

    'Fr', 'Ra',
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo',
]

AtomicMass = [
        # 1
        1.008, 4.002602,
        # 2
        6.94, 9.012182, 10.81, 12.011, 14.007, 15.999, 18.9984032, 20.1797,
        # 3
        22.98976928, 24.3050, 26.9815386, 28.085, 30.973762, 32.06, 35.45, 39.948,
        # 4
        39.0983, 40.078, 44.955912, 47.867, 50.9415, 51.9961, 54.938045, 55.845, 58.933195,
        58.6934, 63.546, 65.38, 69.723, 72.63, 74.92160, 78.96, 79.904, 83.798,
        # 5
        85.4678, 87.62, 88.90585, 91.224, 92.90638, 95.96, 98, 101.07, 102.90550,
        106.42, 107.8682, 112.411, 114.818, 118.710, 121.760, 127.60, 126.90447, 131.293,
        # 6
        132.9054519, 137.327,
        138.90547, 140.116, 140.90765, 144.242, 145, 150.36, 151.964, 157.25,
        158.92535, 162.500, 164.93032, 167.259, 168.93421, 173.054, 174.9668,
        178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084,
        196.966569, 200.59, 204.38, 207.2, 208.98040, 210, 210, 222,
        # 7
        223, 226,
        227, 232.03806, 231.03588, 238.02891, 237, 244, 243, 247, 247, 251, 252, 257, 258, 259, 262,
        261, 262, 266, 264, 269, 268, 271
]

AtomicMass = {s: m for s, m in zip(AtomicSymbols, AtomicMass)}

def parse_AO(s):
    cnt = {'s': 1, 'p': 3, 'd': 6, 'f': 10}
    return sum(int(o[:-1]) * cnt[o[-1]] for o in s.split(','))

AOs_STO_3G = {s: parse_AO(m) for s, m in zip(
    AtomicSymbols,
    ['1s'] * 2 +
    ['2s,1p'] * 8 +
    ['3s,2p'] * 8 +
    ['4s,3p'] * 2 + ['4s,3p,1d'] * 16 +
    ['5s,4p,1d'] * 2 + ['5s,4p,2d'] * 15
)}

AOs_6_31Gd = {s: parse_AO(m) for s, m in zip(
    AtomicSymbols,
    ['2s'] * 2 +
    ['3s,2p,1d'] * 8 +
    ['4s,3p,1d'] * 8 +
    ['5s,4p,1d'] * 2 + ['5s,4p,2d,1f'] * 10 + ['5s,4p,3d'] * 6
)}

BasisSet = {
    ('STO-3G', False): AOs_STO_3G,
    ('6-31G(d)', False): AOs_6_31Gd,
}


def parse_info(file):
    with file:
        result = {}
        for line in file:
            if ':' not in line:
                continue

            line = line.strip().split(':')
            field = line[0].strip()
            value = ' '.join(line[1:]).strip()
            result[field] = value

        return result


class OutputParser(object):
    def __init__(self, file):
        self.file = file

    def __enter__(self):
        self.iter = iter(self.file)
        return self

    def __exit__(self, *args):
        self.file.close()

    def next(self):
        return self.iter.next()

    def _skip_to(self, p):
        if isinstance(p, str):
            s = p
            p = lambda l: s in l

        while True:
            line = self.next()
            if p(line):
                return line

    def _parse_atoms(self):
        self._skip_to('MOLECULE AND FRAGMENT INFORMATION')
        in_section = False
        self.atoms = {}

        while True:
            line = self.next().split()
            if len(line) != 7:
                if in_section:
                    break
                else:
                    continue

            in_section = True
            self.atoms[int(line[0])] = line[1]

    def _parse_int(self, target):
        line = self._skip_to(target)
        return int(line.split('=')[1])

    def _parse_fragments(self):
        self._skip_to('Frag.   Elec.   ATOM')
        self.fragments = {}
        ifrag = None

        while True:
            line = self.next()
            if not line.strip():
                break

            if line[:21].strip():
                if ifrag is not None:
                    self.fragments[ifrag] = atoms

                ifrag = int(line[:13])
                atoms = set()

            atoms = atoms.union(int(a) for a in line[22:].split())

    def _parse_string(self, target, conv=lambda s: s):
        line = self._skip_to(target)
        return conv(line.split('=')[1].strip())

    def _parse_time(self, target):
        line = self._skip_to(target)
        return float(line.split('=')[1].split()[0])

    @property
    def max_mw(self):
        return max(sum(AtomicMass[self.atoms[i]] for i in f) for _, f in self.fragments.items())

    @property
    def max_ao(self):
        bs = BasisSet.get((self.basis, self.diffuse), None)
        if bs is None:
            return

        return max(sum(bs[self.atoms[i]] for i in f) for _, f in self.fragments.items())

    def parse(self):
        try:
            self.diffuse = self._parse_string('DiffuseOn', lambda l: l in ['ON', 'YES'])
            self._parse_atoms()
            self.nf = self._parse_int('NUMBER OF FRAGMENTS')
            self._parse_fragments()
            self.basis = self._parse_string('BASIS SET =')
            self.monomer_time = self._parse_time('Elapsed time: Monomer (Total)')
            self.dimer_time = self._parse_time('Elapsed time: Dimer (Total)')
            self.node_num = self._parse_int('Number of cores (total)')
            self.total_time = self._parse_time('Total time')
        except StopIteration:
            pass


class JobInfo(object):
    def __init__(self, wg='PPI'):
        self.account = os.environ['USER']
        self.wg = wg

    def fetch_pdb(self, pdbid=None):
        if pdbid is not None:
            self.pdb_id = pdbid

        assert hasattr(self, 'pdb_id')

        resp = requests.get('http://www.rcsb.org/pdb/rest/describeMol?structureId={}'.format(self.pdb_id))
        tree = ElementTree.fromstring(resp.content)
        kinds = [
            e.attrib['description']
            for p in tree[0]
            for e in p
            if e.tag == 'polymerDescription'
        ]

        self.protein_kind = '/'.join(kinds)

    def parse_info(self, path):
        result = parse_info(path)
        def s(attr, field, conv=lambda s: s):
            v = result.get(field)
            if v is not None:
                setattr(self, attr, conv(v))

        s('job_id', 'JOB ID')
        s('account', 'USER')
        s('date', 'ACCEPT DATE', lambda t: t.split()[0])
        s('node_num', 'NODE NUM (REQUIRE)', int)
        s('cpu_num', 'CPU NUM (USE)', int)
        s('total_time', 'CPU TIME (TOTAL)', lambda t: float(t.split()[0]) / self.cpu_num / 1000)

    def parse_output(self, file):
        with OutputParser(file) as o:
            o.parse()

            def s(attr, conv=lambda s: s):
                v = getattr(o, attr, None)
                if v is not None:
                    setattr(self, attr, conv(v))

            s('nf')
            s('monomer_time')
            s('dimer_time')
            s('node_num')
            s('total_time')
            s('max_mw')
            s('max_ao')

    @property
    def node_hour(self):
        return float(self.node_num * self.total_time) / 3600

    def __str__(self):
        def g(name):
            return getattr(self, name, '')

        return '\t'.join([
            g('date'),
            g('job_id'),
            g('account'),
            g('wg'),
            g('script_file'),
            g('protein_kind'),
            g('pdb_id'),
            str(g('node_num')),
            str(g('total_time')),
            str(g('node_hour')),
            str(g('nf')),
            str(g('max_mw')),
            str(g('max_ao')),
            str(g('monomer_time')),
            str(g('dimer_time')),
        ])

    def report(self):
        return '''
date:                     {date}
job id:                   {job_id}

account:                  {account}
warking group:            {wg}

script file:              {script_file}

pdb id:                   {pdb_id}
protein kind:             {protein_kind}

number of fragments:      {nf}
maximum molecular weight: {max_mw}
maximum atomic orbital:   {max_ao}

monomer time:             {monomer_time}
dimer time:               {dimer_time}
total time:               {total_time}

number of nodes:          {node_num}
node hour:                {node_hour}
'''[1:-1].format(node_hour=self.node_hour, **self.__dict__)


def main():
    pdbid = None
    info = None
    output = None
    shell = None

    info_regex = re.compile(r'\.i\d+')

    human = False

    for arg in sys.argv[1:]:
        if arg in ['-h', '--human-readable']:
            human = True
            continue

        base, ext = os.path.splitext(arg)

        if ext in ['.out', '.log'] and os.path.isfile(arg):
            output = arg

        elif info_regex.match(ext) and os.path.isfile(arg):
            info = arg

        elif ext == '.sh' and os.path.isfile(arg):
            shell = arg

        elif len(base) == 4 and len(ext) <= 2:
            pdbid = arg.upper()

    if not any([pdbid, info, output, shell]):
        sys.stderr.write('''
USAGE: jr INPUT1 [INPUT2...]

options:
    -h, --human-readable: human readable output

input types:
    output: output file (*.out, *.log)
    info: job information file (*.i\\d*)
    script: job script file (*.sh)
    PDB: pdb id

informations:
    date:         info
    job_id:       info
    account:      info, $USER
    wg:           hard-coded
    script_file:  script
    protein_kind: PDB
    pdb_id:       PDB
    node_num:     output, info
    total_time:   output, info
    node_hour:    output, info
    nf:           output
    max_mw:       output
    max_ao:       output
    monomer_time: output
    dimer_time:   output
'''[1:])
        sys.exit(1)

    ji = JobInfo()

    if pdbid is not None:
        ji.fetch_pdb(pdbid) 

    if shell is not None:
        ji.script_file = os.path.abspath(shell)

    if info is not None:
        ji.parse_info(open(info))

    if output is not None:
        ji.parse_output(open(output))

    if human:
        print(ji.report())
    else:
        print(ji)


if __name__ == '__main__':
    main()
