from typing import List

import numpy as np
import tqdm

from _libg96 import G96Mol
from _libgro import GroMol
from _libtxyz import TXYZMol


def txyz_to_g96(
    txyz_mol: TXYZMol,
) -> G96Mol:
    return G96Mol(
        title=txyz_mol.title,
        timestep=[0, 0.0],
        position=txyz_mol.position * 0.1,  # Angstrom -> nm
        velocity=np.zeros_like(txyz_mol.position),
        box_vector=txyz_mol.box_vector * 0.1,  # Angstrom -> nm
    )


def __txyz_connectivity(
    txyz_mol: TXYZMol,
) -> List[List[int]]:
 # generate residue number list
    resi_map = {}

    for i in tqdm.trange(txyz_mol.num_atoms):
        atom_idx = i + 1
        connect = txyz_mol.connectivity[i]

        if len(resi_map.keys()) == 0:
            # add to map
            resi_map[1] = connect
            resi_map[1].append(atom_idx)

            continue

        for (key, val) in resi_map.items():
            # if current atom can be found in saved connectivity
            if atom_idx in val:
                # set resi_num
                resi_num[i] = key

                break

            # if saved atoms in current atom's connectivy
            if any([v in connect for v in val]):
                # add to map
                resi_map[key].append(atom_idx)

                # set resi_num
                resi_num[i] = key

                break

        if resi_num[i] == 0:
            # current atom belongs to a new residue
            resi_idx = max(resi_map.keys()) + 1

            # add to map
            resi_map[resi_idx] = connect
            resi_map[resi_idx].append(atom_idx)

            # set resi_num
            resi_num[i] = resi_idx


def txyz_to_gro(
    txyz_mol: TXYZMol,
) -> GroMol:
    resi_num = [1] * txyz_mol.num_atoms
    resi_name = ["MOL"] * txyz_mol.num_atoms

    return GroMol(
        title=txyz_mol.title,
        num_atoms=txyz_mol.num_atoms,
        resi_num=resi_num,
        resi_name=resi_name,
        atom_name=txyz_mol.atom_name,
        atom_num=txyz_mol.atom_num,
        position=txyz_mol.position * 0.1,  # Angstrom -> nm
        velocity=np.zeros_like(txyz_mol.position),
        box_vectors=txyz_mol.box_vector * 0.1,  # Angstrom -> nm
    )
