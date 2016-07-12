import numpy as np
import click
import sys
import os
from collections import defaultdict


class IFIE(object):
    def parse(self, f):
        with f:
            return self._parse(iter(f))

    def _parse(self, it):
        res_header = 'Seq. Frag. Residue S-S  N-term.  C-Term. Charge'
        frag_header = 'Frag.   Elec.   ATOM'
        ifie_header = 'IJ-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE'
        hyphens = '-' * 20

        while True:
            line = next(it)
            if 'ReadGeom' in line:
                self.pdb_file = line.split()[2].strip()
                break

        while True:
            line = next(it)
            if 'AutoFrag' in line:
                self.auto_frag = line.split()[2] == 'ON'
                break

        if self.auto_frag:
            while res_header not in next(it):
                pass

            self.sequences = []
            self.residues = []

            while True:
                line = next(it)
                if not line.strip():
                    break

                seq = int(line[1:5])
                res = line[16:19]
                self.sequences.append(seq)
                self.residues.append(res)

            self.sequences = np.array(self.sequences)
            self.residues = np.array(self.residues)
            self.N = len(self.sequences)

        elif not os.path.isfile(self.pdb_file):
            while True:
                line = next(it)
                if ' NF ' in line:
                    self.N = int(line.split()[2])
                    break

            self.sequences = np.arange(self.N) + 1
            self.residues = np.array(['XXX'] * self.N)

        else:
            with open(self.pdb_file) as pdb_file:
                atoms = self._parse_pdb(pdb_file)

            while frag_header not in next(it):
                pass

            current_fragment = None
            residues = []
            sequences = []

            while True:
                line = next(it)
                if not line.strip():
                    break

                if line[:20].strip():
                    if current_fragment is not None:
                        res, i = list(sorted([
                            (n, k) for k, n in current_fragment.items()
                        ]))[-1][1]

                        residues.append(res)
                        sequences.append(i)

                    current_fragment = defaultdict(int)

                for i in map(int, line[24:].split()):
                    current_fragment[atoms[i]] += 1

            if current_fragment:
                res, i = list(sorted([
                    (n, k) for k, n in current_fragment.items()
                ]))[-1][1]

                residues.append(res)
                sequences.append(i)

            self.N = len(residues)
            self.residues = np.array(residues)
            self.sequences = np.array(sequences)

        self.distance = np.zeros(shape=(self.N, self.N), dtype='float64')
        self.ifie = np.zeros(shape=(self.N, self.N, 2), dtype='float64')
        self.ligand_mask = np.zeros(shape=self.N, dtype='bool')
        self.exclude_mask = np.zeros(shape=self.N, dtype='bool')

        while ifie_header not in next(it):
            pass

        while hyphens not in next(it):
            pass

        while True:
            line = next(it)
            if not line.strip():
                break

            i = int(line[:13]) - 1
            j = int(line[13:18]) - 1
            d = float(line[18:30])

            hf = float(line[39:50])
            mp2 = float(line[50:61])

            self.ifie[i, j, :] = hf, mp2
            self.ifie[j, i, :] = hf, mp2
            self.distance[i, j] = d
            self.distance[j, i] = d

        self.ifie *= 627.509

        try:
            while '## PIEDA' not in next(it):
                pass

            self.pieda = np.zeros(shape=(self.N, self.N, 4), dtype='float64')
        except StopIteration:
            self.pieda = None

        if self.pieda is None:
            return

        while hyphens not in next(it):
            pass

        while True:
            line = next(it)
            if not line.strip():
                break

            i = int(line[:13]) - 1
            j = int(line[13:18]) - 1

            if self.distance[i, j] == 0:
                continue

            es = float(line[18:33])
            ex = float(line[33:48])
            ct_mix = float(line[48:63])
            di = float(line[63:78])

            self.pieda[i, j, :] = es, ex, ct_mix, di
            self.pieda[j, i, :] = es, ex, ct_mix, di

    def _parse_pdb(self, f):
        rs = [None]

        for line in f:
            if line[:6] not in ['ATOM  ', 'HETATM']:
                continue

            res = line[17:20]
            resi = int(line[22:26])
            rs.append((res, resi))

        return rs

    def frag_mask(self, seqs):
        m = np.zeros(len(self.sequences), dtype='bool')
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
        self.exclude_mask[self.residues == 'HOH'] = True
        self.exclude_mask[self.residues == 'WAT'] = True

    def exclude(self, ex):
        self.exclude_mask[self.frag_mask(ex)] = True

    def only(self, inc):
        self.exclude_mask[~self.frag_mask(inc)] = True

    def get_masked(self, attr):
        return getattr(self, attr)[self.ligand_mask][:, ~self.exclude_mask]

    def _gen_csv(header_str, ene_attr):
        def dec(conv):
            def f(self, header=True, sep=','):
                if header:
                    yield header_str

                i = self.sequences[~self.exclude_mask]
                r = self.residues[~self.exclude_mask]
                if self.ligand_mask.any():
                    d = self.get_masked('distance').min(axis=0)
                else:
                    d = [0] * len(i)

                e = self.get_masked(ene_attr).sum(axis=0)

                I = np.arange(1, self.N + 1)[~self.exclude_mask]

                for i, seq, res, d, ene in zip(I, i, r, d, e):
                    yield sep.join(
                        [str(i)] + list(conv(self, seq, res, d, *ene))
                    )

            return f

        return dec

    @staticmethod
    def _res_fmt(i, r):
        return '{0}{1}'.format(r, i)

    @staticmethod
    def _ene_fmt(f):
        return '{0:.6f}'.format(f)

    @_gen_csv('i,residue,distance,ES,EX,CT+mix,DI', 'pieda')
    def gen_pieda_csv(self, i, res, *ene):
        return [self._res_fmt(i, res)] + [self._ene_fmt(e) for e in ene]

    @_gen_csv('i,residue,distance,HF,MP2', 'ifie')
    def gen_ifie_csv(self, i, res, *ene):
        return [self._res_fmt(i, res)] + [self._ene_fmt(e) for e in ene]

    def gen_csv(self, pieda=False, header=True, sep=','):
        f = self.gen_pieda_csv if pieda else self.gen_ifie_csv
        return f(header=header, sep=sep)

    def plot(self, pieda=False, ylim=None, figsize=None, legend=True,
             title='Interaction energy'):

        from matplotlib import pyplot

        if figsize:
            fig = pyplot.figure(figsize=figsize)
        else:
            fig = pyplot.figure()
        ax = fig.add_subplot(111)

        if title:
            ax.set_title(title)

        Ylabel = ('ES', 'EX', 'CT+mix', 'DI') if pieda else ('HF', 'MP2')

        Y = self.get_masked('pieda' if pieda else 'ifie').sum(axis=0).T
        X = np.arange(Y.shape[1]) * (len(Ylabel) + 1)

        colors = pyplot.rcParams['axes.prop_cycle']()

        if not pieda:
            next(colors)
            next(colors)

        for i, l, Y, c in zip(range(4), Ylabel, Y, colors):
            ax.bar(X + i, Y, color=c['color'])

        if legend:
            ax.legend(Ylabel, loc='best')

        ax.set_ylabel('IFIE (kcal/mol)')

        ax.set_xticks(X + len(Ylabel) / 2)
        ax.set_xticklabels(
            map(self._res_fmt,
                self.sequences[~self.exclude_mask],
                self.residues[~self.exclude_mask]
                ),
            rotation='vertical'
        )

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


@click.command()
# files
@click.argument(
    'input', type=click.File('r'), metavar='OUT_FILE'
)
@click.option(
    '-o', '--output',
    help='png/csv output file',
    default=None, type=click.Path()
)
# selector
@click.option(
    '-l', '--ligand',
    help='ligand indices', metavar='QUERY',
    default=None, type=IFIE.parse_indices
)
@click.option(
    '-e', '--exclude',
    help='exclude indices', metavar='QUERY',
    default=None, type=IFIE.parse_indices
)
@click.option(
    '-i', '--only',
    help='include only these indices', metavar='QUERY',
    default=None, type=IFIE.parse_indices
)
@click.option(
    '-d', '--exclude-far',
    help='exclude far fragments', metavar='DISTANCE',
    default=None, type=float
)
# flags
@click.option(
    '-p', '--pieda',
    help='use pieda',
    is_flag=True
)
@click.option(
    '-H', '--exclude-water',
    help='exclude HOH/WAT fragments',
    is_flag=True
)
# mode
@click.option(
    '-c', '--csv',
    help='csv mode',
    is_flag=True
)
# csv specific
@click.option(
    '-s', '--separator',
    help='[csv] specify separator',
    type=str, default=','
)
# plot specific
@click.option(
    '-s', '--size',
    help='[plot] output image size', metavar='WIDTH HEIGHT',
    nargs=2, type=int
)
@click.option(
    '-y', '--ylim',
    help='[plot] y axis range', metavar='MIN MAX',
    nargs=2, type=float
)
@click.option(
    '-f', '--font-size',
    help='[plot] font size', metavar='SIZE',
    type=int, default=None
)
@click.option(
    '-t', '--theme',
    help='[plot] theme', metavar='NAME',
    type=str, default='ggplot'
)
@click.option(
    '-L', '--no-legend',
    help='[plot] plot without legend',
    is_flag=True, default=False
)
@click.option(
    '-t', '--title',
    help='[plot] override title',
    type=str, default='Interaction energy'
)
def main(ligand, only, exclude, exclude_far, exclude_water, pieda, csv,
         input, separator, output, ylim, size,
         font_size, theme, no_legend, title):

    p = IFIE()
    p.parse(input)

    if ligand:
        p.add_ligand(ligand)

    if exclude:
        p.exclude(exclude)

    if exclude_far:
        p.exclude_far(exclude_far)

    if exclude_water:
        p.exclude_water()

    if only:
        p.only(only)

    if csv:
        if output is None:
            output = sys.stdout
        else:
            output = open(output, 'w')

        with output:
            for line in p.gen_csv(pieda=pieda, sep=separator, header=True):
                output.write(line)
                output.write('\n')

    else:
        import matplotlib

        if output is None:
            matplotlib.use('Qt4Agg')

        from matplotlib import pyplot

        try:
            pyplot.style.use(theme)
        except IOError:
            sys.stderr.write(
                '''WARNING: theme {0} is not available.
choose one of {1}
'''.format(theme, ','.join(pyplot.style.available))
            )

        if font_size is not None:
            pyplot.rcParams['font.size'] = font_size

        p.plot(pieda=pieda, ylim=ylim, figsize=size,
               legend=not no_legend, title=title)

        pyplot.tight_layout()

        if output is None:
            pyplot.show()
        else:
            pyplot.savefig(output)


if __name__ == '__main__':
    main()
