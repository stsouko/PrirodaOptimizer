from itertools import dropwhile, islice, takewhile, tee
from logging import warning
from multiprocessing.pool import ThreadPool
from os.path import dirname, join
from shutil import rmtree
from subprocess import call
from tempfile import mkdtemp

from sklearn.base import BaseEstimator, TransformerMixin

from .atom_map import atom_map, reverse_atom_map
from .conformer import Conformer


class PrirodaOptimizer(BaseEstimator, TransformerMixin):
    def __init__(self, relativistic=False, basis='L1', memory=200, tmp_dir='/tmp', tmp_ram=10, n_jobs=1, n_process=1,
                 steps=100):
        self.relativistic = relativistic
        self.basis = basis
        self.memory = memory
        self.tmp_dir = tmp_dir
        self.tmp_ram = tmp_ram
        self.n_jobs = n_jobs
        self.n_process = n_process
        self.steps = steps

    def fit(self, x, y=None):
        return self

    def transform(self, x, y=None):
        if not isinstance(x, (list, tuple)):
            raise TypeError
        for i in x:
            if not isinstance(i, Conformer):
                raise TypeError

        jobs_1, jobs_2 = tee(self._prepare_calc(x))
        results = []
        with ThreadPool(self.n_jobs) as p:
            for rc, (dir_, _, log), conf in zip(p.imap(call,
                                                       (('priexec', '-np', str(self.n_process), 'pribin', inp, out) for
                                                        _, inp, out in jobs_1)),
                                                jobs_2, x):
                if not rc:
                    try:
                        a, c, hc, mc, mb, e, h = _parse_output(open(log))
                    except:
                        warning('log file parsing error')
                        results.append(conf)
                    else:
                        results.append(Conformer(a, c, conf.charge, conf.multiplicity,
                                                 _hirshfeld_charges=hc, _mulliken_charges=mc, _mulliken_bonds=mb,
                                                 _energy=e, _hessian=h))
                else:
                    warning('priroda exit with error code')
                    results.append(conf)
                rmtree(dir_)
        return results

    def _prepare_input(self, conf: Conformer, tmp_dir):
        basis = join(dirname(__file__), f'basis{int(self.relativistic)}')
        out = [f'$system memory={self.memory} disk={self.tmp_ram} path={tmp_dir} $end',
               f'$control task=optimize+hessian theory=dft four={int(self.relativistic)} print=+bonds+charges '
               f'basis={basis} $end',
               '$dft functional=pbe $end', f'$optimize steps={self.steps} $end',
               f'$molecule charge={conf.charge} mult={conf.multiplicity}',
               ' cartesian',
               f' set={self.basis}']
        for atom, (x, y, z) in zip(conf.atoms, conf.coords):
            out.append(f' {atom_map[atom]} {x:.4f} {y:.4f} {z:.4f}')
        out.append('$end')
        return '\n'.join(out)

    def _prepare_calc(self, confs):  # generator
        for conf in confs:
            tmp_dir = mkdtemp(dir=self.tmp_dir)
            task = self._prepare_input(conf, tmp_dir)
            inp = join(tmp_dir, 'input')
            out = join(tmp_dir, 'output')
            with open(inp, 'w') as f:
                f.write(task)
            yield tmp_dir, inp, out


def _parse_output(log):
    atoms = []
    coords = []
    hirshfeld_charges = []
    mulliken_charges = []
    mulliken_bonds = {}
    mapping = {}

    for line in takewhile(lambda x: not x.startswith('MOL>$end'),
                          islice(dropwhile(lambda x: not x.startswith('MOL> set'), log), 1, None)):
        _, atom, x, y, z = line.split()
        atoms.append(reverse_atom_map[atom])
        coords.append((float(x), float(y), float(z)))

    energy = float(next(log).split()[-1])

    for n, line in enumerate(islice(dropwhile(lambda x: not x.startswith(' Overlap populations:'), log),
                                    2, 2 + len(atoms))):
        a, c, _ = line.split()
        mapping[a] = n
        mulliken_charges.append(float(c))

    for line in islice(log, 2, 2 + len(atoms)):
        n, _, *ms = line.split()
        mulliken_bonds[mapping[n]] = b = {}
        for m, _, o in zip(*[iter(ms)] * 3):
            b[mapping[m]] = float(o)

    for line in islice(log, 4, 4 + len(atoms)):
        hirshfeld_charges.append(float(line.split()[2]))

    freq = next(islice(dropwhile(lambda x: not x.startswith(' | Mode |'), log), 2, None)).split()[3]
    gibbs = float(next(islice(dropwhile(lambda x: not x.startswith(' translational'), log), 3, None)).split()[-1])

    bonds = []
    seen = set()
    for n, ms in mulliken_bonds.items():
        seen.add(n)
        for m, b in ms.items():
            if m not in seen:
                bonds.append((n, m, b))

    return (tuple(atoms), tuple(coords), tuple(hirshfeld_charges), tuple(mulliken_charges), tuple(bonds),
            energy * 627.51 + gibbs, 'i' not in freq)


__all__ = ['PrirodaOptimizer']
