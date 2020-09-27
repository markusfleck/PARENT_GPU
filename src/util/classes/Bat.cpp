//    A class for handling an EntropyMatrix from a .par file of the PARENT_GPU
//    suite Copyright (C) 2020  Markus Fleck
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License  version 3
//    as published by the Free Software Foundation.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "Bat.h"
#include "My_Error.cpp"

#pragma pack(1)

using namespace std;

Bat::Bat(char const *infile_str) {
  infile.exceptions(std::ifstream::badbit);
  infile.open(infile_str, ios::binary | ios::out); // open the .par file;
  if (!infile.good()) {
    My_Error my_error((string("ERROR OPENING FILE ") + string(infile_str) +
                       string("! ABORTING.\n"))
                          .c_str());
    throw my_error;
  }
  infile.exceptions(std::ifstream::badbit | std::ifstream::failbit |
                    std::ifstream::eofbit);

  try {
    read_BAT_header();
  } catch (My_Error my_error) {
    throw my_error;
  } catch (...) {
    My_Error my_error((string("ERROR WHILE READING THE HEADER OF THE FILE ") +
                       string(infile_str) + string("! ABORTING."))
                          .c_str());
    throw my_error;
  }
}

// to read the header of the binary .bat file
void Bat::read_BAT_header() {
  int version;
  int fail = 0;
  char dummystring[31];

  infile.read((char *)&version,
              sizeof(int)); // first read the .bat version number as an integer
  fail = fail | (infile.rdstate() & std::ifstream::failbit);
  infile.read((char *)&precision,
              sizeof(int)); // then read an integer declaring if the trajectory
                            // is stored in double precision
  fail = fail | (infile.rdstate() & std::ifstream::failbit);
  infile.read(
      (char *)&n_dihedrals,
      sizeof(int)); // read an integer containing the number of dihedrals
  fail = fail | (infile.rdstate() & std::ifstream::failbit);
  infile.read((char *)&n_frames,
              sizeof(int)); // and an integer containing the number of frames
  fail = fail | (infile.rdstate() & std::ifstream::failbit);

  if (version < 0) {
    ostringstream oss;
    oss << "ERROR: FILE HEADER CORRUPTED. VERSION NUMBER (" << version
        << ") < 0!";
    My_Error my_error(oss.str());
    throw my_error;
  }
  if ((precision != 0) && (precision != 1)) {
    ostringstream oss;
    oss << "ERROR: FILE HEADER CORRUPTED. DOUBLE PRECISION VALUE (" << precision
        << ") NEITHER 0 NOR 1!";
    My_Error my_error(oss.str());
    throw my_error;
  }
  if (n_dihedrals > 19997) {
    cerr << "WARNING: " << n_dihedrals + 3
         << " ATOMS DECLARED IN THE FILE HEADER (CORRUPTED?). THIS WILL LEAD "
            "TO LARGE OUTPUT."
         << endl;
  }
  if (n_dihedrals < 0) {
    ostringstream oss;
    oss << "ERROR: FILE HEADER CORRUPTED. NUMBER OF DIHEDRALS (" << n_dihedrals
        << ") < 0!";
    My_Error my_error(oss.str());
    throw my_error;
  }
  if (n_frames < 1) {
    ostringstream oss;
    oss << "ERROR: FILE HEADER CORRUPTED. NUMBER OF FRAMES (" << n_frames
        << ") < 1!";
    My_Error my_error(oss.str());
    throw my_error;
  }

  if (version >= 3) {
    for (unsigned int i = 0; i < n_dihedrals + 3;
         i++) { // for ever atom in the system
      infile.read(
          dummystring,
          8 * sizeof(char)); // read the name of the residue it belongs to
      residues.push_back(dummystring);
      fail = fail | (infile.rdstate() & std::ifstream::failbit);
      residue_numbers.push_back(0);
      infile.read(
          (char *)&(residue_numbers[i]),
          sizeof(float)); // read the number of the residue it belongs to
      fail = fail | (infile.rdstate() & std::ifstream::failbit);
      infile.read(dummystring, 8 * sizeof(char)); // read the name of the atom
      atom_names.push_back(dummystring);
      fail = fail | (infile.rdstate() & std::ifstream::failbit);
      infile.read(
          dummystring,
          31 * sizeof(char)); // read the molecule of the residue it belongs to
      belongs_to_molecule.push_back(dummystring);
      fail = fail | (infile.rdstate() & std::ifstream::failbit);
    }
  }

  vector<int> dummyvec;
  dummyvec.push_back(0);
  dummyvec.push_back(0);
  dummyvec.push_back(0);
  dummyvec.push_back(0);
  dummyvec.push_back(0);
  dummyvec.push_back(0);
  for (unsigned int i = 0; i < n_dihedrals; i++) { // then for all dihedrals
    dihedrals_top.push_back(dummyvec);
    infile.read((char *)&(dihedrals_top[i][0]),
                sizeof(int)); // read the atomnumber of the first atom
    fail = fail | (infile.rdstate() & std::ifstream::failbit);
    infile.read((char *)&(dihedrals_top[i][1]), sizeof(int)); // second atom
    fail = fail | (infile.rdstate() & std::ifstream::failbit);
    infile.read((char *)&(dihedrals_top[i][2]), sizeof(int)); // third atom
    fail = fail | (infile.rdstate() & std::ifstream::failbit);
    infile.read((char *)&(dihedrals_top[i][3]), sizeof(int)); // fourth atom
    fail = fail | (infile.rdstate() & std::ifstream::failbit);
    infile.read((char *)&(dihedrals_top[i][4]),
                sizeof(int)); // an integer containing the type of the
                              // dihedral(physical=0,pseudo=1,improper=-1)
    fail = fail | (infile.rdstate() & std::ifstream::failbit);
    infile.read(
        (char *)&(dihedrals_top[i][5]),
        sizeof(int)); // and an integer containing the "parent"dihedral(for
                      // phaseangles, -1 if no "parent")
    fail = fail | (infile.rdstate() & std::ifstream::failbit);
  }
  for (unsigned int i = 0; i < n_dihedrals + 3;
       i++) { // and read the whole massestor of the atoms in the system in
              // single precision (float)
    masses.push_back(0);
    infile.read((char *)&(masses[i]), sizeof(float));
    fail = fail | (infile.rdstate() & std::ifstream::failbit);
  }

  dofs_begin = infile.tellg();
}

streamoff Bat::get_dofs_begin() { return dofs_begin; }

ifstream *Bat::get_file() { return &infile; }

unsigned int Bat::get_n_frames() { return n_frames; }

unsigned int Bat::get_n_dihedrals() { return n_dihedrals; }

unsigned int Bat::get_precision() { return precision; }

#pragma pack(0)
