import argparse
import glob
import os
import warnings

import numpy as np
from libmol import G96Mol, GroMol, TXYZMol

_FILE_TYPE = {
    "tinker": ["xyz", "txyz", "arc"],
    "g96": ["g96"],
    "gro": ["gro"],
}


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
        help="Trajectory: tinker xyz or arc, g96, gro",
        required=True,
    )

    parser.add_argument(
        "-s",
        type=str,
        help="Structure: tinker xyz, gro",
        required=False,
    )

    # options to specify output file
    parser.add_argument(
        "-o",
        type=str,
        help="Trajectory: tinker xyz or arc, gro, g96",
        required=False,
        default="traj.g96",
    )

    # other options
    parser.add_argument(
        "-t0",
        type=float,
        help="Start time: defaults to 0.0 ps",
        required=False,
        default=0.0,
    )

    parser.add_argument(
        "-timestep",
        type=float,
        help="Time step between input frames: in ps",
        required=True,
    )

    parser.add_argument(
        "-dt",
        type=float,
        help="Only write frame when t MOD dt = first time (ps)",
        required=False,
    )

    args = parser.parse_args()

    return args


def _determine_type(file: str) -> str:
    for (key, val) in _FILE_TYPE.items():
        for v in val:
            if v in file:
                return key

    raise ValueError(
        f"Cannot determine the type of {file}. Avaliable file suffixes are "
        f"{[v for val in _FILE_TYPE.values() for v in val]}."
    )


def _check_args(args):
    # check args
    f_type = _determine_type(args.f)
    print(f"Input file: {args.f}, type: {f_type}")

    o_type = _determine_type(args.o)
    print(f"Output file: {args.o}, type: {o_type}")

    if args.s:
        s_type = _determine_type(args.s)
        print(f"Structure file: {args.s}, type: {s_type}")

    if f_type != "tinker" and o_type != "tinker":
        print(
            "The tinker trjconv is developed to convert trajectories bewteen "
            "tinker and other formats. At least the input (-f) or output (-o) "
            "file should be in tinker format . For converting trajectories "
            "between gro and g96 formats, use GROMACS please."
        )

        exit()

    if f_type == "tinker" and o_type == "tinker":
        print(
            "Why do you want to convert a trajectory from tinker format to "
            "tinker format?"
        )

        exit()

    # input: tinker
    if f_type == "tinker" and o_type == "gro":
        # output: gro
        if (not args.s) or (args.s and s_type != "gro"):
            raise RuntimeError(
                "To convert to trajectory in gro format, the input file "
                "(-f) or structure file (-s) should be specified to a gro "
                "file."
            )

    if o_type == "tinker":
        if not args.s and f_type == "tinker":
            print(
                "No structure (-s) file specified, use input tinker file as "
                "template."
            )

            args.s = args.f
        elif args.s and s_type != "tinker":
            raise RuntimeError(
                "To convert to trajectory in tinker format, the input file "
                "(-f) or structure file (-s) should be specified to a tinker "
                "file."
            )


def _auto_backup(
    file: str,
) -> None:
    if not os.path.exists(file):
        return

    file_dir = os.path.dirname(file)
    file_name = os.path.basename(file)

    old_file_pattern = f"#{file_name}.*#"
    old_files = sorted(glob.glob(os.path.join(
        file_dir,
        old_file_pattern,
    )))

    num_old_files = len(old_files)

    new_backup_file_name = f"#{file_name}.{num_old_files + 1}#"
    new_backup_file = os.path.join(file_dir, new_backup_file_name)

    os.rename(file, new_backup_file)


if __name__ == "__main__":
    # parse args
    args = _parse_args()

    # check args
    _check_args(args)

    if _determine_type(args.f) == "tinker":
        traj = TXYZMol.traj_from_file(args.f)
    elif _determine_type(args.f) == "gro":
        raise NotImplementedError()
    elif _determine_type(args.f) == "g96":
        raise NotImplementedError()

    if _determine_type(args.o) == "tinker":
        raise NotImplementedError()

    # auto backup
    _auto_backup(args.o)

    with open(args.o, "w") as f:
        for (idx, txyz_mol) in enumerate(traj):
            if _determine_type(args.o) == "g96":
                mol = G96Mol(
                    title=txyz_mol.title,
                    timestep=[idx, (args.t0 + args.timestep * idx)],
                    position=txyz_mol.position * 0.1,  # Angstrom -> gro
                    velocity=np.zeros_like(txyz_mol.position),
                    box_vector=txyz_mol.box_vector,
                )
            elif _determine_type(args.o) == "gro":
                raise NotImplementedError()

            print(mol, file=f, flush=True)

    print("Done!")
