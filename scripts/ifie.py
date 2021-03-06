import numpy as np
import sys
import os
from collections import defaultdict
import argparse


class Parser(object):
    res_header = 'Seq. Frag. Residue S-S  N-term.  C-Term. Charge'
    frag_header = 'Frag.   Elec.   ATOM'
    ifie_header = 'IJ-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE'
    hyphens = '-' * 20

    def parse(self, f):
        with f:
            self.it = iter(f)
            return self._parse()

    def __next__(self):
        return next(self.it)

    next = __next__

    def __iter__(self):
        return self.it

    @staticmethod
    def _to_checker(i):
        if isinstance(i, (str, unicode)):
            def check(line):
                return i in line

            return check

        elif hasattr(i, '__call__'):
            return i

        else:
            raise ValueError('unknown check type')

    def _drop_while(self, p):
        check = self._to_checker(p)

        for line in self:
            if check(line):
                return line

    def _feed_until(self, p):
        check = self._to_checker(p)

        for line in self:
            if not check(line):
                return

            yield line

    def _parse_ReadGeom(self):
        line = self._drop_while('ReadGeom')
        self.pdb_file = line.split()[2].strip()

    def _parse_AutoFrag(self):
        line = self._drop_while('AutoFrag')
        self.auto_frag = line.split()[2].strip() == 'ON'

    def _parse_residues_AutoFrag(self):
        self._drop_while(self.res_header)

        residues = []

        for line in self._feed_until(lambda l: l.strip()):
            seq = int(line[1:5])
            res = line[16:19]
            residues.append('{0}{1}'.format(res, seq))

        self.residues = np.array(residues)
        self.N = len(residues)

    def _parse_residues_no_AutoFrag_without_PDB(self):
        line = self._drop_while(' NF ')
        self.N = N = int(line.split()[2])

        self.residues = np.array([
            'XXX{0}'.format(i + 1) for i in xrange(N)
        ])

    def _parse_residues_no_AutoFrag_with_PDB(self):
        with open(self.pdb_file) as pdb_file:
            atoms = self._parse_pdb(pdb_file)

        self._drop_while(self.frag_header)

        current = None
        residues = []

        def append_fragment(frags):
            if not frags:
                return

            reses = list(sorted([
                (n, k) for k, n in frags.items()
            ]))
            if len(reses) == 1:
                return residues.append('{0}{1}'.format(*reses[0][1]))

            larges = [k for (n, k) in reses if n > 2]
            if len(larges) > 0:
                return residues.append(
                    '/'.join(['{0}{1}'.format(*r) for r in larges])
                )

            residues.append('/'.join(['{0}{1}'.format(*r) for r in reses]))

        for line in self._feed_until(lambda l: l.strip()):
            if line[:20].strip():
                append_fragment(current)
                current = defaultdict(int)

            for i in line[24:].split():
                current[atoms[int(i)]] += 1

        append_fragment(current)
        self.N = len(residues)
        self.residues = np.array(residues)

    def _parse_ifie(self):
        dist = np.zeros(shape=(self.N, self.N), dtype='float64')
        ifie = np.zeros(shape=(self.N, self.N, 2), dtype='float64')

        self._drop_while(self.ifie_header)
        self._drop_while(self.hyphens)

        for line in self._feed_until(lambda l: l.strip()):
            i = int(line[:13]) - 1
            j = int(line[13:18]) - 1
            d = float(line[18:30])

            hf = float(line[39:50])
            mp2 = float(line[50:61])

            ifie[i, j, :] = hf, mp2
            ifie[j, i, :] = hf, mp2
            dist[i, j] = d
            dist[j, i] = d

        self.distance = dist
        self.ifie = ifie * 627.509

    def _parse_pieda(self):
        try:
            self._drop_while('## PIEDA')
        except StopIteration:
            self.pieda = None
            return

        self._drop_while(self.hyphens)

        pieda = np.zeros(shape=(self.N, self.N, 4), dtype='float64')

        for line in self._feed_until(lambda l: l.strip()):
            i = int(line[:13]) - 1
            j = int(line[13:18]) - 1

            if self.distance[i, j] == 0:
                continue

            es = float(line[18:33])
            ex = float(line[33:48])
            ct = float(line[48:63])
            di = float(line[63:78])

            pieda[i, j, :] = es, ex, ct, di
            pieda[j, i, :] = es, ex, ct, di

        self.pieda = pieda

    def _parse(self):
        self._parse_ReadGeom()
        self._parse_AutoFrag()

        if self.auto_frag:
            self._parse_residues_AutoFrag()

        elif not os.path.isfile(self.pdb_file):
            self._parse_residues_no_AutoFrag_without_PDB()

        else:
            self._parse_residues_no_AutoFrag_with_PDB()

        self._parse_ifie()
        self._parse_pieda()

    def _parse_pdb(self, f):
        rs = [None]

        for line in f:
            if line[:6] not in ['ATOM  ', 'HETATM']:
                continue

            res = line[17:20]
            resi = int(line[22:26])
            rs.append((res, resi))

        return rs


class IFIE(object):
    def parse(self, f):
        p = Parser()
        p.parse(f)
        self.N = p.N
        self.residues = p.residues
        self.distance = p.distance
        self.ifie = p.ifie
        self.pieda = p.pieda

        self.exclude_mask = np.zeros(p.N, dtype='bool')
        self.ligand_mask = np.zeros(p.N, dtype='bool')

    def frag_mask(self, seqs):
        m = np.zeros(self.N, dtype='bool')
        s = np.array(list(seqs)) - 1
        m[s] = True
        return m

    def add_ligand(self, ligands):
        i = self.frag_mask(ligands)

        self.ligand_mask[i] = True
        self.exclude_mask[i] = True
        self.exclude_mask[(self.distance[i, :] == 0).any(axis=0)] = True

    def exclude_far(self, d=8.0):
        far_mask = self.distance[self.ligand_mask, :].min(axis=0) > d
        self.exclude_mask[far_mask] = True

    def exclude_water(self):
        get_name = np.vectorize(lambda x: x[:3])
        self.exclude_mask[get_name(self.residues) == 'HOH'] = True
        self.exclude_mask[get_name(self.residues) == 'WAT'] = True

    def exclude(self, ex):
        self.exclude_mask[self.frag_mask(ex)] = True

    def only(self, inc):
        self.exclude_mask[~self.frag_mask(inc)] = True

    def important(self, threshold=100.0):
        ene = np.abs(self.total[self.ligand_mask, :, 0].sum(axis=0))
        self.exclude_mask[ene < threshold] = True

    @property
    def total(self):
        return self.ifie.sum(axis=2)[:, :, np.newaxis]

    def get_masked(self, attr):
        return getattr(self, attr)[self.ligand_mask][:, ~self.exclude_mask]

    def _gen_csv(header_str, ene_attr):
        def dec(conv):
            def f(self, header=True, sep=','):
                if header:
                    yield header_str

                r = self.residues[~self.exclude_mask]
                if self.ligand_mask.any():
                    d = self.get_masked('distance').min(axis=0)
                else:
                    d = [0] * len(r)

                e = self.get_masked(ene_attr).sum(axis=0)

                I = np.arange(1, self.N + 1)[~self.exclude_mask]

                for i, res, d, ene in zip(I, r, d, e):
                    yield sep.join(
                        [str(i)] + list(conv(self, res, d, *ene))
                    )

            return f

        return dec

    @staticmethod
    def _ene_fmt(f):
        return '{0:.6f}'.format(f)

    @_gen_csv('i,residue,distance,ES,EX,CT+mix,DI', 'pieda')
    def gen_pieda_csv(self, res, *ene):
        return [res] + [self._ene_fmt(e) for e in ene]

    @_gen_csv('i,residue,distance,HF,MP2', 'ifie')
    def gen_ifie_csv(self, res, *ene):
        return [res] + [self._ene_fmt(e) for e in ene]

    @_gen_csv('i,residue,distance,IFIE', 'total')
    def gen_total_csv(self, res, *ene):
        return [res] + [self._ene_fmt(e) for e in ene]

    def gen_csv(self, mode='ifie', header=True, sep=','):
        if mode == 'ifie':
            f = self.gen_ifie_csv
        elif mode == 'pieda':
            f = self.gen_pieda_csv
        else:
            f = self.gen_total_csv

        return f(header=header, sep=sep)

    def plot(self, mode='ifie', ylim=None, figsize=None, legend=True,
             distance=False,
             title='Interaction energy'):

        if self.exclude_mask.all():
            raise ValueError('all residues are ignored.')

        from matplotlib import pyplot
        from matplotlib.font_manager import FontProperties

        if figsize:
            fig = pyplot.figure(figsize=figsize)
        else:
            fig = pyplot.figure()
        ax = fig.add_subplot(111)

        if title:
            japanese = False
            try:
                title.decode("ascii")
            except UnicodeEncodeError:
                japanese = True

            if japanese:
                fp = FontProperties(
                    fname="/usr/share/fonts/ipa-pgothic/ipagp.ttf",
                    size=16
                )
                ax.set_title(title, fontproperties=fp)
            else:
                ax.set_title(title)

        Ylabel = {
            'pieda': ('ES', 'EX', 'CT+mix', 'DI'),
            'ifie': ('HF', 'MP2'),
            'total': ('IFIE',),
        }[mode]

        Y = self.get_masked(mode).sum(axis=0).T
        X = np.arange(Y.shape[1]) * (len(Ylabel) + 1)

        colors = pyplot.rcParams['axes.prop_cycle']()

        if mode != 'pieda':
            next(colors)
            next(colors)

        for i, l, Y, c in zip(range(4), Ylabel, Y, colors):
            ax.bar(X + i, Y, color=c['color'])

        if legend:
            ax.legend(Ylabel, loc='best')

        ax.set_ylabel('IFIE (kcal/mol)')

        ax.set_xticks(X + len(Ylabel) / 2.0)

        if distance:
            xlabels = map(
                '{0}\n{1:.1f}'.format,
                self.residues[~self.exclude_mask],
                self.get_masked('distance').min(axis=0)
            )
        else:
            xlabels = self.residues[~self.exclude_mask]

        ax.set_xticklabels(xlabels, rotation='vertical')

        ax.set_xticks(X - 0.5, minor=True)

        ax.grid(False, which='major', axis='x')
        ax.grid(True, which='minor')

        ax.tick_params(axis='x', which='major', bottom=False, top=False)
        ax.set_xlim(X.min() - 0.5, X.max() + len(Ylabel) + 0.5)

        if ylim:
            ax.set_ylim(*ylim)

    @staticmethod
    def parse_indices(s):
        i = set()
        for r in s.split(','):
            if '-' in r:
                lo, hi = map(int, r.split('-'))
                i.update(range(lo, hi + 1))

            else:
                i.add(int(r))

        return i


class QueryAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, IFIE.parse_indices(values))


class FuzzyChoiceAction(argparse.Action):
    def __init__(self, option_strings, dest, choices=None, **kwargs):
        super(FuzzyChoiceAction, self).__init__(option_strings, dest, **kwargs)
        self.candidates = choices
        self.metavar = '{{{0}}}'.format(','.join(map(str, choices)))

    def __call__(self, parser, namespace, values, option_string=None):
        L = len(values)
        cand = [c for c in self.candidates if c[:L].lower() == values.lower()]
        if len(cand) != 1:
            raise ValueError('ambiguous option')

        setattr(namespace, self.dest, cand[0])


def main():
    parser = argparse.ArgumentParser(
        prog='ifie',
        epilog='query example: "1", "1-5", "1-8,12"'
    )

    mode_group = parser.add_argument_group('mode')
    mode_group.add_argument(
        '-u', '--use', help='data to use. (default: %(default)s)',
        default='total', choices=('ifie', 'pieda', 'total'),
        action=FuzzyChoiceAction
    )
    mode_group.add_argument(
        '-m', '--mode', default='csv', choices=('csv', 'plot'),
        help='set mode (default: %(default)s)',
        action=FuzzyChoiceAction
    )

    file_group = parser.add_argument_group('files')
    file_group.add_argument(
        'input', type=argparse.FileType('r'),
        metavar='OUT_FILE', help='abinit-mp output file'
    )
    file_group.add_argument(
        '-o', '--output',
        type=str, help='image/csv output path', metavar='PATH'
    )

    selector_group = parser.add_argument_group('selectors')
    selector_group.add_argument(
        '-l', '--ligand', help='ligand query', metavar='QUERY',
        action=QueryAction
    )
    selector_group.add_argument(
        '-e', '--exclude', help='exclude query', metavar='QUERY',
        action=QueryAction
    )
    selector_group.add_argument(
        '-i', '--only', help='use only these fragments', metavar='QUERY',
        action=QueryAction
    )
    selector_group.add_argument(
        '-d', '--exclude_far', help='exclude far fragments', metavar='DIST',
        type=float
    )
    selector_group.add_argument(
        '-H', '--exclude-water', action='store_true',
        help='exclude WAT/HOH fragments'
    )
    selector_group.add_argument(
        '-I', '--important', type=float,
        help='only important residues (default: %(const)s)',
        default=None, const=100.0, nargs='?', metavar='ENERGY'
    )

    csv_group = parser.add_argument_group('csv mode')
    csv_group.add_argument(
        '--delimiter', help='csv delimiter (default: %(default)s)',
        type=str, default=','
    )

    plot_group = parser.add_argument_group('plot mode')
    plot_group.add_argument(
        '-s', '--size', help='output image size', metavar=('WIDTH', 'HEIGHT'),
        nargs=2, type=int
    )
    plot_group.add_argument(
        '-y', '--ylim', help='y axis range', metavar=('MIN', 'MAX'),
        nargs=2, type=float
    )
    plot_group.add_argument(
        '-f', '--font-size', help='base font size', metavar='PX', type=int
    )
    plot_group.add_argument(
        '--theme', help='set matplotlib theme',
        type=str, default='ggplot'
    )
    plot_group.add_argument(
        '-D', '--label-distance', help='add distance from ligand to label',
        action='store_true'
    )
    plot_group.add_argument(
        '-L', '--no-legend', help='plot without legend', action='store_true'
    )
    plot_group.add_argument(
        '-t', '--title', help='override title', type=str,
        default='Interaction energy'
    )

    opts = parser.parse_args()

    try:
        os.chdir(os.path.dirname(opts.input.name))
    except OSError:
        pass

    p = IFIE()
    p.parse(opts.input)

    if opts.ligand is not None:
        p.add_ligand(opts.ligand)

    if opts.exclude is not None:
        p.exclude(opts.exclude)

    if opts.exclude_far is not None and p.ligand_mask.any():
        p.exclude_far(opts.exclude_far)

    if opts.exclude_water:
        p.exclude_water()

    if opts.only is not None:
        p.only(opts.only)

    if opts.important is not None:
        p.important(opts.important)

    if opts.mode == 'csv':
        csv_mode(opts, p)
    else:
        plot_mode(opts, p)


def csv_mode(opts, p):
    if opts.output is None:
        output = sys.stdout

    else:
        output = open(output, 'w')

    with output:
        for line in p.gen_csv(mode=opts.use,
                              sep=opts.delimiter, header=True):
            output.write(line)
            output.write('\n')


def plot_mode(opts, p):
    import matplotlib

    if opts.output is None:
        matplotlib.use('Qt4Agg')

    from matplotlib import pyplot

    try:
        pyplot.style.use(opts.theme)
    except IOError:
        sys.stderr.write(
            '''WARNING: theme {0} is not available.
choose one of {1}
'''.format(opts.theme, ','.join(pyplot.style.available))
        )

    if opts.font_size is not None:
        pyplot.rcParams['font.size'] = opts.font_size

    p.plot(
        mode=opts.use, ylim=opts.ylim, figsize=opts.size,
        legend=not opts.no_legend, title=opts.title.decode("UTF-8"),
        distance=opts.label_distance
    )
    pyplot.tight_layout()

    if opts.output is None:
        pyplot.show()
    else:
        pyplot.savefig(opts.output)


if __name__ == '__main__':
    main()
