import warnings
from dataclasses import dataclass, field
from typing import List

import numpy as np


@dataclass
class TXYZMol:
    num_atoms: int
    title: str
    atom_num: List[int]
    atom_name: List[str]
    position: np.array  # Angstrom
    atom_type: np.array
    connectivity: List[List[int]]
    box_vector: np.array = field(
        default_factory=lambda: np.zeros(3)
    )  # Angstrom

    def __str__(self) -> str:
        contents = []

        # 1st line: number of atoms + title
        contents.append(f" {self.num_atoms:5d}  {self.title}")

        # if with pbc info
        if not np.allclose(self.box_vector, np.zeros(3)):
            contents.append(
                " "
                f"{self.box_vector[0]:12.6f}"
                f"{self.box_vector[1]:12.6f}"
                f"{self.box_vector[2]:12.6f}"
                f"{90.0:12.6f}"
                f"{90.0:12.6f}"
                f"{90.0:12.6f}"
            )

        # atoms
        for i in range(self.num_atoms):
            contents.append(
                " "
                f"{(i + 1):5d}  {self.atom_name[i]:3s}"
                f"{self.position[i, 0]:12.6f}"
                f"{self.position[i, 1]:12.6f}"
                f"{self.position[i, 2]:12.6f}"
                f"{self.atom_type[i]:6d}"
            )

            for connect in self.connectivity[i]:
                contents[-1] += f"{connect:6d}"

        return "\n".join(contents)

    @staticmethod
    def from_file(
        file: str,
    ) -> "TXYZMol":
        """Load molecule from input tinker xyz file.

        Args:
            file (str): Path to input tinker xyz file

        Returns:
            txyz_mol (TXYZMol): A tinker xyz molecule
        """
        with open(file) as f:
            contents = f.readlines()

        num_atoms = 0
        title = None
        atom_num = []
        atom_name = []
        position = []
        atom_type = []
        connectivity = []
        box_vector = None

        for (line_idx, line) in enumerate(contents):
            line = line.rstrip()
            # title line
            if line_idx == 0:
                (num_atoms, title) = line.split(maxsplit=1)
                num_atoms = int(num_atoms)

                continue

            # the 2nd line may contains the pbc box vector
            if line_idx == 1:
                try:
                    box_vector = np.array(
                        [float(x) for x in line.split()[0: 3]]
                    )
                    print(f">>> PBC box info found: {box_vector}")
                    continue
                except ValueError:
                    print(">>> No PBC box info found")
                finally:
                    pass

            # parse position
            arr = line.split()
            atom_num.append(int(arr[0]))
            atom_name.append(arr[1])
            position.append([float(x) for x in arr[2: 5]])
            atom_type.append(int(arr[5]))
            connectivity.append([int(x) for x in arr[6:]])

        position = np.array(position)

        if len(position) != num_atoms:
            warnings.warn(
                f"The number of atoms ({num_atoms}) does not match the size of"
                f" positions {position.shape}, reset the number of atoms to "
                f"{len(position)}."
            )

        num_atoms = len(position)

        return TXYZMol(
            num_atoms=num_atoms,
            title=title,
            atom_num=atom_num,
            atom_name=atom_name,
            position=position,
            atom_type=atom_type,
            connectivity=connectivity,
            box_vector=box_vector,
        )
