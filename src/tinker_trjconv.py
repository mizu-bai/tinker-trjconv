import argparse
import os
import warnings

import numpy as np

from _libpbc import _PBC_ENUM, PBCFixer
from _libtxyz import TXYZMol


def _parse_args():
    # prepare argument parser
    parser = argparse.ArgumentParser(
        prog="tinker-trjconv",
        description="A Tinker Trajectory Converter",
    )

    # options to specify input file
    parser.add_argument(
        "-f",
        type=str,
        help="Trajectory: tinker xyz",
        required=True,
    )

    # options to specify output file
    parser.add_argument(
        "-o",
        type=str,
        help="Trajectory: tinker xyz",
        required=False,
        default="mol.txyz",
    )

    # options for pbc
    parser.add_argument(
        "-pbc",
        type=str,
        help="PBC treatment: mol, res, atom",
        required=True,
    )

    args = parser.parse_args()

    return args


def _check_args(args):
    # check args:
    if not os.path.exists(args.f):
        raise FileNotFoundError(
            f"Input tinker xyz file {args.f} does not exist."
        )

    if args.pbc not in _PBC_ENUM:
        raise ValueError(
            f"Available -pbc options are {_PBC_ENUM}"
        )


if __name__ == "__main__":
    # parse args
    args = _parse_args()

    # check args
    _check_args(args)

    # load input tinker xyz
    txyz_mol = TXYZMol.from_file(args.f)

    # check pbc info
    if np.allclose(txyz_mol.box_vector, np.zeros(3)):
        warnings.warn(
            f"Input tinker xyz file {args.f} contains no pbc info. "
            "No modification will be done."
        )

    # process the position of atoms
    pbc_fixer = PBCFixer(pbc=args.pbc)

    txyz_mol.position = pbc_fixer.fix_pbc(
        txyz_mol.position * 0.1,  # Angstrom -> nm
        txyz_mol.box_vector * 0.1,  # Angstrom -> nm
    ) * 10.0  # nm -> Angstrom

    # write to output tinker xyz file
    with open(args.o, "w") as f:
        f.write(f"{txyz_mol}\n")
