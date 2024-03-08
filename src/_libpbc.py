import os
import shutil
import subprocess
import warnings
from typing import Literal

import numpy as np

from _libg96 import G96Mol

_PBC_ENUM = [
    "mol",
    "res",
    "atom",
]

_BACKEND_ENUM = [
    "gmx",
    "python",
]


class PBCFixer:
    def __init__(
        self,
        pbc: Literal["mol", "res", "atom"] = "atom",
        backend: Literal["gmx", "python"] = "gmx",
        workdir: str = ".pbcfixer",
    ) -> None:
        _pbc = pbc.lower()
        _backend = backend.lower()

        # check pbc option
        if _pbc not in _PBC_ENUM:
            raise ValueError(
                f"Avaliable pbc options are: {_PBC_ENUM}"
            )

        # check backend
        if _backend not in _BACKEND_ENUM:
            raise ValueError(
                f"Avaliable backends are: {_BACKEND_ENUM}"
            )

        # check compatibility between pbc option and backend
        if _backend == "gmx" and not os.environ.get("gmx"):
            if _pbc in ["mol", "res"]:
                raise RuntimeError(
                    f"The option -pbc {pbc}"
                    "Please make sure you have GROMACS installed and loaded. "
                )
            else:
                _backend = "python"

                warnings.warn(
                    "Cannot find gmx in environment variables, backend fall "
                    "back to python. NOTE: only -pbc atom is supported!"
                )

        self.pbc = _pbc
        self.backend = _backend
        self.workdir = os.path.abspath(workdir)

    def _python_pbc(
        self,
        position: np.array,
        box_vector: np.array,
    ) -> np.array:
        new_position = position.copy()

        for i in range(3):
            new_position[:, i] %= box_vector[i]

        return new_position

    def _gmx_pbc(
        self,
        position: np.array,
        box_vector: np.array,
    ) -> np.array:
        mol = G96Mol(
            title="Created by PBCFixer",
            timestep=[0, 0.0],
            position=position,
            velocity=np.zeros_like(position),
            box_vector=box_vector,
        )

        with open(os.path.join(self.workdir, "mol.g96"), "w") as f:
            f.write(f"{mol}\n")

        # call gmx
        subprocess.run(
            f"gmx trjconv -f mol.g96 -o mol_fix.g96 -pbc {self.pbc}"
        )

        mol_fix = G96Mol.from_file(os.path.join(self.workdir, "mol.g96"))

        return mol_fix.position

    def fix_pbc(
        self,
        position: np.array,
        box_vector: np.array,
    ) -> np.array:
        """Fix PBC issue.

        Args:
            position (np.array): Input positions, in nm
            box_vector (np.array): PBC box vector, in nm

        Returns:
            new_position (np.array): Output positions, in nm
        """

        return eval(f"self._{self.backend}_pbc")(position, box_vector)
