import AbinsModules
import io
import numpy as np
from mantid.kernel import Atom


class LoadVASP(AbinsModules.GeneralAbInitioProgram):
    """
    Class for loading VASP ab initio vibrational or phonon data.
    """

    def __init__(self, input_ab_initio_filename=None):
        """
                :param input_ab_initio_filename: name of a file with vibrational or phonon data (foo.out)
                """
        super(LoadVASP, self).__init__(input_ab_initio_filename=input_ab_initio_filename)

        self._parser = AbinsModules.GeneralAbInitioParser()
        self._ab_initio_program = "VASP"

    def read_vibrational_or_phonon_data(self):
        """
        Reads vibrational or phonon data from VASP output files. Saves frequencies, weights of k-point vectors,
        k-point vectors, amplitudes of atomic displacements, hash of the vibrational or phonon data file (hash) to
        <>.hdf5.
        :return  object of type AbinsData.
        """
        data = {}  # container to store read data

        with io.open(self._clerk.get_input_filename(), "rb", ) as vasp_file:

            # read info about atoms and construct atoms data
            self._read_atomic_coordinates(file_obj=vasp_file, data=data)

            # read lattice vectors
            self._read_lattice_vectors(obj_file=vasp_file, data=data)

            # read frequencies, corresponding atomic displacements and construct k-points data
            self._read_modes(file_obj=vasp_file, data=data)

        # save data to hdf file
        self.save_ab_initio_data(data=data)

        # return AbinsData object
        return self._rearrange_data(data=data)

    def _read_lattice_vectors(self, obj_file=None, data=None):
        """
        Reads lattice vectors from .outmol VASP file.
        :param obj_file: file object from which we read
        :param data: Python dictionary to which found lattice vectors should be added
        """
        self._parser.find_last(file_obj=obj_file, msg="Lattice vectors:", move_next_line=True)

        dim = 3
        vectors = []
        obj_file.readline()  # blank line

        for i in range(dim):
            line = obj_file.readline().replace(",", " ").replace(")", "").split()[3:]
            vector = [float(s) for s in line]
            vectors.append(vector)
        data["unit_cell"] = np.asarray(vectors).astype(dtype=AbinsModules.AbinsConstants.FLOAT_TYPE)

    def _read_atomic_coordinates(self, file_obj=None, data=None):
        """
        Reads atomic coordinates from OUTCAR VASP file.

        :param file_obj: file object from which we read
        :param data: Python dictionary to which atoms data should be added
        """
        # self._parser.find_first(file_obj=file_obj, msg="VRHFIN")  # find  the first symbol of an element
        atoms = {}
        tot_num_atoms = 0  # counts number of atoms in the system
        masses = []  # stores masses read from the file

        # collect symbols for all atoms in the system
        pos = None
        while True:

            symbol = self._parser.find_first(file_obj=file_obj, msg="VRHFIN")

            if symbol is None:
                break
            else:
                symbol = symbol.replace("=", " ").replace(":", " ").split()[1]

            # read symbol of an element
            atoms["atom_{}".format(tot_num_atoms)] = {"symbol": symbol,
                                                      "mass": Atom(symbol=symbol).mass, "sort":tot_num_atoms}
            # read corresponding masses
            masses.append(float(self._parser.find_first(file_obj=file_obj,
                                                  msg="POMASS").replace(";", " ").split()[2]))
            tot_num_atoms += 1
            pos = file_obj.tell()

        # check if isotopes are present and update mass if necessary
        self.check_isotopes_substitution(atoms=atoms, masses=masses)

        # find equilibrium coordinates
        file_obj.seek(pos)  # go back with file pointer to the end of atoms section and search for equilibrium positions
        self._parser.find_last(file_obj=file_obj, msg="position of ions in cartesian coordinates  (Angst):",
                               move_next_line=True)

        for i in range(tot_num_atoms):
            atoms["atom_{}".format(i)].update({"coord": [float(j) for j in file_obj.readline().split()]})

        data["atoms"] = atoms

    def _read_modes(self, file_obj=None, data=None):
        """
        Reads vibrational or phonon modes (frequencies and atomic displacements).
        :param file_obj: file object from which we read
        :param data: Python dictionary to which k-point data should be added
        """
        block_begin = self._parser.find_first(file_obj=file_obj,
                                              msg="Eigenvectors and eigenvalues of the dynamical matrix")
        if block_begin is None:
            raise ValueError("No block with eigenvectors and eigenvalues of the dynamical matrix has been found.")

