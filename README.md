# Optimizer for PRIRODA

Example of usage:

    from PrirodaOptimizer import Conformer, PrirodaOptimizer

    c1 = Conformer.from_xyz('data/methane.xyz')

    # rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem

    m = Chem.MolFromSmiles('C')
    m = Chem.AddHs(m)
    AllChem.EmbedMultipleConfs(m, numConfs=10)
    c2 = Conformer.from_rdkit(m, conformer=-1)

    p = PrirodaOptimizer(tmp_ram=-1000, n_jobs=3)
    for n, o in enumerate(p.transform([c1, c2]), 1):
        print(c.energy)
        print(c.to_chython())  # get chython.MoleculeContainer
        with open(f'opt_{n}.xyz', 'w') as f:
            c.to_xyz(f)

## Overview

PRIRODA Optimizer was written as one of master's projects by our student group at the 
Chemoinformatics and Molecular modeling laboratory, KFU located in Kazan, Russia.
The reason why we decided to create it was to prepare input files when the user needs to do some
quantum chemistry calculations at "Priroda".
At first, it needs a lot of time, because all preparation should be done manually.
The next reason is the consequences of the first one.
You need to be very accurate to avoid potential mistakes that could follow to mistakes in the
result of calculations or to increase the amount of calculation time that in quantum chemistry is
already required.
To make it faster, clear and correct we have written "PRIRODA Optimizer".

Dmitry x2, Zhenya x2, Daria, Polina, Il'nura, Roma, Anya