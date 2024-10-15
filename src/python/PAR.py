import numpy as np

class PAR():
    def __init__(self):
        pass

    def get_int_(self, f):
        
        return int.from_bytes(f.read(4), byteorder='little', signed=True)

    def read_atom_(self, f):

        residue_name = f.read(8).decode("ascii")
        residue_number = self.get_int_(f)
        atom_name = f.read(8).decode("ascii")
        associated_molecule = f.read(31).decode("ascii")

        return residue_name, residue_number, atom_name, associated_molecule

    def read_dihedral_(self, f):
        
        dihedral = []
        for i in range(4):
            dihedral.append(self.get_int_(f))
        dihedral_type = self.get_int_(f)
        parent_dihedral = self.get_int_(f)

        return dihedral, dihedral_type, parent_dihedral
            

    def read_file(self, par_file_path):
        with open(par_file_path, "rb") as f:        
            self.version = self.get_int_(f)
            assert 2 <= self.version <= 4
            
            self.double_prec = self.get_int_(f)
            assert self.double_prec in [0,1]

            self.n_dihedrals = self.get_int_(f)
            assert 0 <= self.n_dihedrals <= 19997
            self.n_angles = self.n_dihedrals + 1
            self.n_bonds = self.n_angles + 1
            self.n_atoms = self.n_bonds + 1

            self.n_frames = self.get_int_(f)
            assert self.n_frames >= 0

            self.bonds_bins_1D = self.get_int_(f)
            assert self.bonds_bins_1D > 0

            self.angles_bins_1D = self.get_int_(f)
            assert self.angles_bins_1D > 0

            self.dihderals_bins_1D = self.get_int_(f)
            assert self.dihderals_bins_1D > 0

            self.bonds_bins_2D = self.get_int_(f)
            assert self.bonds_bins_2D > 0

            self.angles_bins_2D = self.get_int_(f)
            assert self.angles_bins_2D > 0

            self.dihedrals_bins_2D = self.get_int_(f)
            assert self.dihedrals_bins_2D > 0

            self.residue_names = []
            self.residue_numbers = []
            self.atom_names = []
            self.associated_molecules = []
            if self.version >= 3:
                for _ in range(self.n_atoms):
                    residue_name, residue_number, atom_name, associated_molecule = self.read_atom_(f)
                    self.residue_names.append(residue_name)
                    self.residue_numbers.append(residue_number)
                    self.atom_names.append(atom_name)
                    self.associated_molecules.append(associated_molecule)

            self.dihedrals = []
            self.dihedral_types = []
            self.parent_dihedrals = []
            for _ in range(self.n_dihedrals):
                dihedral, dihedral_type, parent_dihedral = self.read_dihedral_(f)
                self.dihedrals.append(dihedral)
                self.dihedral_types.append(dihedral_type)
                self.parent_dihedrals.append(parent_dihedral)

            self.masses = np.fromfile(f, count=self.n_atoms, dtype=np.dtype("f4"))

            prec = "f4" if self.double_prec == 0 else "f8"
            
            self.bonds_entropy_1D = np.fromfile(f, count=self.n_bonds, dtype=np.dtype(prec))
            self.angles_entropy_1D = np.fromfile(f, count=self.n_angles, dtype=np.dtype(prec))
            self.dihedrals_entropy_1D = np.fromfile(f, count=self.n_dihedrals, dtype=np.dtype(prec))
            self.bonds_bonds_entropy_2D = np.fromfile(f, count=self.n_bonds * (self.n_bonds - 1) // 2, dtype=np.dtype(prec))
            self.bonds_angles_entropy_2D = np.fromfile(f, count=self.n_bonds * self.n_angles, dtype=np.dtype(prec))
            self.bonds_dihedrals_entropy_2D = np.fromfile(f, count=self.n_bonds * self.n_dihedrals, dtype=np.dtype(prec))
            self.angles_angles_entropy_2D = np.fromfile(f, count=self.n_angles * (self.n_angles - 1) // 2, dtype=np.dtype(prec))
            self.angles_dihedrals_entropy_2D = np.fromfile(f, count=self.n_angles * self.n_dihedrals, dtype=np.dtype(prec))
            self.dihedrals_dihedrals_entropy_2D = np.fromfile(f, count=self.n_dihedrals * (self.n_dihedrals - 1) // 2, dtype=np.dtype(prec))
    
    def get_version(self):
        return self.version

    def get_double_precision(self):
        return self.double_prec

    def get_n_dihedrals(self):
        return self.n_dihedrals

    def get_n_frames(self):
        return self.n_frames


    def get_bonds_bins_1D(self):
        return self.bonds_bins_1D

    def get_angles_bins_1D(self):
        return self.angles_bins_1D

    def get_dihedrals_bins_1D(self):
        return self.dihderals_bins_1D

    def get_bonds_bins_2D(self):
        return self.bonds_bins_2D

    def get_angles_bins_2D(self):
        return self.angles_bins_2D

    def get_dihedrals_bins_2D(self):
        return self.dihedrals_bins_2D

    def getMutual(self, type1, type2, index1, index2):
        assert index1 >= 0
        assert index2 >= 0
        assert type1 in ("bond", "angle", "dihedral")
        assert type2 in ("bond", "angle", "dihedral")
        assert not ((type1 == type2) and (index1 == index2))

        if type1 == "bond":
            assert index1 < self.n_bonds
        if type1 == "angle":
            assert index1 < self.n_angles
        if type1 == "dihedral":
            assert index1 < self.n_dihedrals
        if type2 == "bond":
            assert index2 < self.n_bonds
        if type2 == "angle":
            assert index2 < self.n_angles
        if type2 == "dihedral":
            assert index2 < self.n_dihedrals

        if (type1 == type2 == "bond"):
            smaller = index1 if index1 < index2 else index2
            bigger = index2 if index1 < index2 else index1
            index = (self.n_bonds - smaller) * (self.n_bonds - smaller - 1) // 2 + smaller - bigger # the 2D-entropies bonds-bonds (also angles-angles and dihedrals-dihedrals) were stored in reverse order, as documented in "Parent.cpp"
            return self.bonds_entropy_1D[smaller] + self.bonds_entropy_1D[bigger] - self.bonds_bonds_entropy_2D[index]

        if (type1 == "bond") and (type2 == "angle"):
            return self.bonds_entropy_1D[index1] + self.angles_entropy_1D[index2] - self.bonds_angles_entropy_2D[index1 * self.n_angles + index2]

        if (type2 == "bond") and (type1 == "angle"):
            return self.bonds_entropy_1D[index2] + self.angles_entropy_1D[index1] - self.bonds_angles_entropy_2D[index2 * self.n_angles + index1]

        if (type1 == "bond") and (type2 == "dihedral"):
            return self.bonds_entropy_1D[index1] + self.dihedrals_entropy_1D[index2] - self.bonds_dihedrals_entropy_2D[index1 * self.n_dihedrals + index2]

        if (type2 == "bond") and (type1 == "dihedral"):
            return self.bonds_entropy_1D[index2] + self.dihedrals_entropy_1D[index1] - self.bonds_dihedrals_entropy_2D[index2 * self.n_dihedrals + index1]

        if (type1 == type2 == "angle"):
            smaller = index1 if index1 < index2 else index2
            bigger = index2 if index1 < index2 else index1
            index = (self.n_angles - smaller) * (self.n_angles - smaller - 1) // 2 + smaller - bigger
            return self.angles_entropy_1D[smaller] + self.angles_entropy_1D[bigger] - self.angles_angles_entropy_2D[index]

        if (type1 == "angle") and (type2 == "dihedral"):
            return self.angles_entropy_1D[index1] + self.dihedrals_entropy_1D[index2] - self.angles_dihedrals_entropy_2D[index1 * self.n_dihedrals + index2]

        if (type2 == "angle") and (type1 == "dihedral"):
            return self.angles_entropy_1D[index2] + self.dihedrals_entropy_1D[index1] - self.angles_dihedrals_entropy_2D[index2 * self.n_dihedrals + index1]

        if (type1 == type2 == "dihedral"):
            smaller = index1 if index1 < index2 else index2
            bigger = index2 if index1 < index2 else index1
            index = (self.n_dihedrals - smaller) * (self.n_dihedrals - smaller - 1) // 2 + smaller - bigger
            return self.dihedrals_entropy_1D[smaller] + self.dihedrals_entropy_1D[bigger] - self.dihedrals_dihedrals_entropy_2D[index]

    def getEntropy2D(self, type1, type2, index1, index2):
        assert index1 >= 0
        assert index2 >= 0
        assert type1 in ("bond", "angle", "dihedral")
        assert type2 in ("bond", "angle", "dihedral")
        assert not ((type1 == type2) and (index1 == index2))

        if type1 == "bond":
            assert index1 < self.n_bonds
        if type1 == "angle":
            assert index1 < self.n_angles
        if type1 == "dihedral":
            assert index1 < self.n_dihedrals
        if type2 == "bond":
            assert index2 < self.n_bonds
        if type2 == "angle":
            assert index2 < self.n_angles
        if type2 == "dihedral":
            assert index2 < self.n_dihedrals

        if (type1 == type2 == "bond"):
            smaller = index1 if index1 < index2 else index2
            bigger = index2 if index1 < index2 else index1
            index = (self.n_bonds - smaller) * (self.n_bonds - smaller - 1) // 2 + smaller - bigger # the 2D-entropies bonds-bonds (also angles-angles and dihedrals-dihedrals) were stored in reverse order, as documented in "Parent.cpp"
            return self.bonds_bonds_entropy_2D[index]

        if (type1 == "bond") and (type2 == "angle"):
            return self.bonds_angles_entropy_2D[index1 * self.n_angles + index2]

        if (type2 == "bond") and (type1 == "angle"):
            return self.bonds_angles_entropy_2D[index2 * self.n_angles + index1]

        if (type1 == "bond") and (type2 == "dihedral"):
            return self.bonds_dihedrals_entropy_2D[index1 * self.n_dihedrals + index2]

        if (type2 == "bond") and (type1 == "dihedral"):
            return self.bonds_dihedrals_entropy_2D[index2 * self.n_dihedrals + index1]

        if (type1 == type2 == "angle"):
            smaller = index1 if index1 < index2 else index2
            bigger = index2 if index1 < index2 else index1
            index = (self.n_angles - smaller) * (self.n_angles - smaller - 1) // 2 + smaller - bigger
            return self.angles_angles_entropy_2D[index]

        if (type1 == "angle") and (type2 == "dihedral"):
            return self.angles_dihedrals_entropy_2D[index1 * self.n_dihedrals + index2]

        if (type2 == "angle") and (type1 == "dihedral"):
            return self.angles_dihedrals_entropy_2D[index2 * self.n_dihedrals + index1]

        if (type1 == type2 == "dihedral"):
            smaller = index1 if index1 < index2 else index2
            bigger = index2 if index1 < index2 else index1
            index = (self.n_dihedrals - smaller) * (self.n_dihedrals - smaller - 1) // 2 + smaller - bigger
            return self.dihedrals_dihedrals_entropy_2D[index]

