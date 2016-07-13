import numpy as np
import sys
import os
from collections import defaultdict
import argparse


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

            self.residues = []

            while True:
                line = next(it)
                if not line.strip():
                    break

                seq = int(line[1:5])
                res = line[16:19]
                self.residues.append('{0}{1}'.format(res, seq))

            self.residues = np.array(self.residues)
            self.N = len(self.residues)

        elif not os.path.isfile(self.pdb_file):
            while True:
                line = next(it)
                if ' NF ' in line:
                    self.N = int(line.split()[2])
                    break

            self.residues = np.array([
                'XXX{0}'.format(i + 1) for i in range(self.N)
            ])

        else:
            with open(self.pdb_file) as pdb_file:
                atoms = self._parse_pdb(pdb_file)

            while frag_header not in next(it):
                pass

            current_fragment = None
            residues = []

            def get_residues(frags):
                reses = list(sorted([
                    (n, k) for k, n in frags.items()
                ]))
                if len(reses) == 1:
                    return '{0}{1}'.format(*reses[0][1])

                larges = [k for (n, k) in reses if n > 2]
                if len(larges) > 0:
                    return '/'.join(['{0}{1}'.format(*r) for r in larges])

                return '/'.join(['{0}{1}'.format(*r) for r in reses])

            while True:
                line = next(it)
                if not line.strip():
                    break

                if line[:20].strip():
                    if current_fragment is not None:
                        residues.append(get_residues(current_fragment))

                    current_fragment = defaultdict(int)

                for i in map(int, line[24:].split()):
                    current_fragment[atoms[i]] += 1

            if current_fragment:
                residues.append(get_residues(current_fragment))

            self.N = len(residues)
            self.residues = np.array(residues)

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

        if figsize:
            fig = pyplot.figure(figsize=figsize)
        else:
            fig = pyplot.figure()
        ax = fig.add_subplot(111)

        if title:
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


def main():
    parser = argparse.ArgumentParser(
        prog='ifie',
        epilog='query example: "1", "1-5", "1-8,12"'
    )

    mode_group = parser.add_argument_group('mode')
    mode_group.add_argument(
        '-u', '--use', help='data to use.',
        default='ifie', choices=('ifie', 'pieda', 'total')
    )
    mode_group.add_argument(
        '-m', '--mode', default='csv', choices=('csv', 'plot'),
        help='set mode (default: %(default)s)',
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
        '-I', '--important', type=float, help='only important residues',
        default=100.0, nargs='?'
    )

    csv_group = parser.add_argument_group('csv mode')
    csv_group.add_argument(
        '--delimiter', help='csv delimiter', type=str, default=','
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

    os.chdir(os.path.dirname(opts.input.name))

    p = IFIE()
    p.parse(opts.input)

    if opts.ligand:
        p.add_ligand(opts.ligand)

    if opts.exclude:
        p.exclude(opts.exclude)

    if opts.exclude_far and p.ligand_mask.any():
        p.exclude_far(opts.exclude_far)

    if opts.exclude_water:
        p.exclude_water()

    if opts.only:
        p.only(opts.only)

    if opts.important:
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
        legend=not opts.no_legend, title=opts.title,
        distance=opts.label_distance
    )
    pyplot.tight_layout()

    if opts.output is None:
        pyplot.show()
    else:
        pyplot.savefig(opts.output)


if __name__ == '__main__':
    main()
