import AbinsModules
import io
import numpy as np
import six
from mantid.kernel import Atom


class LoadVASP(AbinsModules.GeneralAbInitioProgram):
    """
    Class for loading VASP ab initio vibrational or phonon data.
    """

    def __init__(self, input_ab_initio_filename=None):
        """
        :param input_ab_initio_filename: name of a file with vibrational or phonon data (OUTCAR*)
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
        Reads lattice vectors from OUTCAR* VASP file.
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
        atoms = {}
        tot_num_atoms = 0  # counts number of atoms in the system
        masses = []  # stores masses read from the file

        # collect symbols for all atoms in the system
        pos = None
        while True:

            # read symbol of an element
            symbol = self._parser.find_first(file_obj=file_obj, msg="VRHFIN")
            if symbol is None:
                break
            else:
                symbol = symbol.replace("=", " ").replace(":", " ").split()[1]

            # fill dictionary with atoms data
            atoms["atom_{}".format(tot_num_atoms)] = {"symbol": symbol,
                                                      "mass": Atom(symbol=symbol).mass, "sort": tot_num_atoms}
            # read corresponding masses
            masses.append(float(self._parser.find_first(file_obj=file_obj, msg="POMASS").replace(";", " ").split()[2]))
            tot_num_atoms += 1
            pos = file_obj.tell()

        # check if isotopes are present and update mass if necessary
        self.check_isotopes_substitution(atoms=atoms, masses=masses)

        # find equilibrium coordinates
        file_obj.seek(pos)  # go back with file pointer to the end of atoms section and search for equilibrium positions
        self._parser.find_last(file_obj=file_obj, msg="position of ions in cartesian coordinates  (Angst):",
                               move_next_line=True)

        for i in range(tot_num_atoms):
            atoms["atom_{}".format(i)].update({"coord": np.asarray([float(j) for j in file_obj.readline().split()])})

        data["atoms"] = atoms

    def _read_modes(self, file_obj=None, data=None):
        """
        Reads vibrational modes (frequencies and atomic displacements). Currently only for Gamma point calculations.
        :param file_obj: file object from which we read
        :param data: Python dictionary to which k-point data should be added
        """
        self._parser.find_last(file_obj=file_obj, msg="Eigenvectors and eigenvalues of the dynamical matrix",
                               move_next_line=True)
        # beginning of a block
        file_obj.readline()  # ----------------------------------------------------

        freq = []
        displacements = []

        # ending of a block
        end_msg = b"----------------------------------------------"

        # parse data from ASCII file
        while not (self._parser.file_end(file_obj=file_obj) or
                   self._parser.block_end(file_obj=file_obj, msg=[end_msg])):

            # In each iteration data for one mode is read
            self._read_freq_block(file_obj=file_obj, freq=freq)
            self._read_coord_block(file_obj=file_obj, disp=displacements)

        # convert parsed data to AbINS data structure
        self._num_k = 1  # Only Gamma point
        freq = [freq]
        weights = [1.0]
        k_coordinates = [[0.0, 0.0, 0.0]]
        self._num_atoms = len(data["atoms"])

        data["frequencies"] = np.asarray(freq).astype(dtype=AbinsModules.AbinsConstants.FLOAT_TYPE, casting="safe")
        data["k_vectors"] = np.asarray(k_coordinates).astype(dtype=AbinsModules.AbinsConstants.FLOAT_TYPE,
                                                             casting="safe")
        data["weights"] = np.asarray(weights).astype(dtype=AbinsModules.AbinsConstants.FLOAT_TYPE, casting="safe")

        # num_freq: number of modes (frequencies; fundamentals)
        # num_atom: number of atoms in the system
        # dim: spacial dimension (atoms vibrate in 3D so dim=3)
        #
        # disp[num_freq, num_atom, dim]
        disp = np.asarray(displacements).astype(dtype=AbinsModules.AbinsConstants.COMPLEX_TYPE, casting="safe")

        # multiply displacements by sqrt of mass
        masses = np.asarray([data["atoms"]["atom_{}".format(i)]["mass"] for i in range(self._num_atoms)])
        disp = np.einsum('ijk,j->ijk', disp, np.sqrt(masses))

        # normalize atomic displacements

        # [num_freq, num_atom, dim] -> [num_freq, num_atom, dim, dim] -> [num_freq, num_atom] -> [num_freq]
        norm = np.sum(np.trace(np.einsum('kin, kim->kinm', disp, disp.conjugate()), axis1=2, axis2=3), axis=1)

        factor = 1.0 / np.sqrt(norm)
        disp = np.einsum('jkl, j-> jkl', disp, factor)

        # [num_freq, num_atom, dim] -> [num_atom, num_freq, dim]
        disp = np.transpose(a=disp, axes=(1, 0, 2))

        # k_num: number of k-points
        # [num_atom, num_freq, dim] -> [k_num, num_atom, num_freq, dim]
        data["atomic_displacements"] = np.asarray([disp])

    def _read_freq_block(self, file_obj=None, freq=None):
        """
        Parses block with frequencies.
        :param file_obj: file object from which we read
        :param freq: list with frequencies which we update
        """
        mode_marker = "f"
        soft_mode_marker = "f/i"

        if not six.PY2:
            mode_marker = bytes(mode_marker, "utf8")
            soft_mode_marker = bytes(soft_mode_marker, "utf8")

        line = self._parser.find_first(file_obj=file_obj, msg="f").replace("=", "").split()

        if mode_marker in line:
            freq.append(float(line[6]))
        elif soft_mode_marker in line:
            freq.append(-float(line[6]))
        else:
            raise ValueError("Invalid frequency.")

    def _read_coord_block(self, file_obj=None, disp=None):
        """
        Parses block with coordinates.
        :param file_obj: file object from which we read
        :param disp: list with x, y, z coordinates which we update
        """
        # beginning of a block with atomic displacements for the given mode
        self._parser.find_first(file_obj=file_obj, msg=" X ")

        mode_disp = []
        modes_end_block = b"----------------------------------------------"  # end of block with modes
        begin_freq_block = "f"  # a new block with data for the next mode

        while not (self._parser.file_end(file_obj=file_obj) or
                   self._parser.block_end(file_obj=file_obj, msg=[modes_end_block, begin_freq_block])):

            line = file_obj.readline()

            if line.strip():
                line = line.strip(b"\n").split()
                mode_disp.append([float(line[3]), float(line[4]), float(line[5])])

        disp.append(mode_disp)
