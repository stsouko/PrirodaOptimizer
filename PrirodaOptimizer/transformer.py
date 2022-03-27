# -*- coding: utf-8 -*-
#
#  Copyright 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of PrirodaOptimizer.
#
#  PrirodaOptimizer is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from copy import deepcopy
from importlib.util import find_spec
from itertools import dropwhile, islice, takewhile, tee
from logging import warning
from multiprocessing.pool import ThreadPool
from os import environ
from os.path import dirname, join
from shutil import rmtree
from subprocess import call
from tempfile import mkdtemp
from .atom_map import atom_map, reverse_atom_map
from .conformer import Conformer


if find_spec('scikit-learn'):
    from sklearn.base import BaseEstimator, TransformerMixin
else:
    class BaseEstimator:
        ...

    class TransformerMixin:
        ...


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

    def fit(self, x, y=None):  # sklearn api
        return self

    def transform(self, x, y=None):  # sklearn api
        return list(self.imap(x))

    def imap(self, x):
        if not isinstance(x, (list, tuple)):
            raise TypeError
        for i in x:
            if not isinstance(i, Conformer):
                raise TypeError

        jobs_1, jobs_2 = tee(self._prepare_calc(x))
        with ThreadPool(self.n_jobs) as p:
            for rc, (dir_, _, log), conf in zip(p.imap(lambda x: call(x, env=environ),
                                                       (('priexec', '-np', str(self.n_process), 'pribin', inp, out) for
                                                        _, inp, out in jobs_1)),
                                                jobs_2, x):
                if not rc:
                    try:
                        c = _parse_output(open(log), self.steps, conf.charge, conf.multiplicity)
                    except:
                        warning('log file parsing error')
                        conf = deepcopy(conf)
                        conf.log = open(log).read()
                        yield conf
                    else:
                        c.log = open(log).read()
                        yield c
                else:
                    warning('priroda exit with error code')
                    yield conf
                rmtree(dir_)

    def _prepare_input(self, conf: Conformer, tmp_dir):
        basis = join(dirname(__file__), f'basis{int(self.relativistic)}')
        out = [f'$system\n memory={self.memory}\n disk={self.tmp_ram}\n path={tmp_dir}\n$end',
               f'$control\n task=optimize+hessian\n theory=dft\n four={int(self.relativistic)}\n print=+bonds+charges',
               f' basis={basis}\n$end', '$dft\n functional=pbe\n$end', f'$optimize\n steps={self.steps}\n$end',
               f'$molecule\n charge={conf.charge}\n mult={conf.multiplicity}\n cartesian\n set={self.basis}']
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


def _parse_output(log, steps, ch, ml):
    a, c, hc, mc, mb, e, g = _parse_step_zero(log)
    confs = [Conformer(a, c, ch, ml, _hirshfeld_charges=hc, _mulliken_charges=mc, _mulliken_bonds=mb, _energy=e,
                       _gradient=g)]

    for _ in range(steps):
        a, c, hc, mc, mb, e, g, f = _parse_step(log)
        confs.append(Conformer(a, c, ch, ml, _hirshfeld_charges=hc, _mulliken_charges=mc, _mulliken_bonds=mb, _energy=e,
                               _gradient=g))
        if f:
            break

    freq = 'i' not in next(islice(dropwhile(_nil('| Mode |'), log), 2, None)).split()[3]
    gibbs = float(next(islice(dropwhile(lambda x: not x.startswith(' translational'), log), 3, None)).split()[-1])

    c = confs.pop(-1)
    c.gibbs = c.energy * 627.51 + gibbs
    c.converged = freq
    c.steps = confs
    return c


def _nil(k):
    def f(line):
        return k not in line
    return f


def _parse_step_zero(log):
    atoms = []
    coords = []

    for line in takewhile(_nil('#'),
                          islice(dropwhile(_nil('Atomic Coordinates'),
                                           takewhile(_nil('Geometry Optimization Step'), log)), 1, None)):
        atom, x, y, z = line.split()
        atoms.append(reverse_atom_map[atom])
        coords.append((float(x), float(y), float(z)))

    step = takewhile(_nil('Geometry Optimization Step'),
                     islice(dropwhile(_nil('Geometry Optimization Step'), log), 1, None))
    hc, mc, mb = _parse_charges_bonds(step, len(atoms))
    e, g = _parse_energy(step, len(atoms))
    return tuple(atoms), tuple(coords), hc, mc, mb, e, g


def _parse_step(log):
    step = takewhile(_nil('Geometry Optimization Step'),
                     islice(dropwhile(_nil('Geometry Optimization Step'), log), 1, None))
    a, c, f = _parse_mol(step)
    hc, mc, mb = _parse_charges_bonds(step, len(a))
    e, g = _parse_energy(step, len(a))
    return a, c, hc, mc, mb, e, g, f


def _parse_charges_bonds(log, atoms_count):
    hirshfeld_charges = []
    mulliken_charges = []
    mulliken_bonds = {}
    mapping = {}

    for n, line in enumerate(islice(dropwhile(lambda x: 'Overlap populations:' not in x, log),
                                    2, 2 + atoms_count)):
        a, c, _ = line.split()
        mapping[a] = n
        mulliken_charges.append(float(c))

    for line in takewhile(lambda x: 'time' not in x, islice(log, 2, None)):
        ls = line.split()
        if ls[0] == '|':
            ms = ls[1:]
        else:
            n, _, *ms = ls
        mulliken_bonds[mapping[n]] = b = {}
        for m, _, o in zip(*[iter(ms)] * 3):
            b[mapping[m]] = float(o)

    assert len(mulliken_bonds) == atoms_count

    for line in islice(log, 3, 3 + atoms_count):
        hirshfeld_charges.append(float(line.split()[2]))

    bonds = []
    seen = set()
    for n, ms in mulliken_bonds.items():
        seen.add(n)
        for m, b in ms.items():
            if m not in seen:
                bonds.append((n, m, b))
    return tuple(hirshfeld_charges), tuple(mulliken_charges), tuple(bonds)


def _parse_energy(log, atoms_count):
    e = float(next(dropwhile(lambda x: 'eng> E=' not in x, log))[7:])
    gradient = []
    for line in islice(dropwhile(lambda x: 'eng> G=' not in x, log), atoms_count):
        x, y, z = line[7:].split()
        gradient.append((float(x), float(y), float(z)))
    return e, tuple(gradient)


def _parse_mol(log):
    atoms = []
    coords = []

    for line in takewhile(_nil('$end'),
                          islice(dropwhile(lambda x: not x.startswith(('mol> set', 'MOL> set')), log), 1, None)):
        p, a, x, y, z = line.split()
        atoms.append(reverse_atom_map[a])
        coords.append((float(x), float(y), float(z)))
    return tuple(atoms), tuple(coords), p == 'MOL>'


__all__ = ['PrirodaOptimizer']
