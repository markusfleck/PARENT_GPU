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

#ifndef BAT_H
#define BAT_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

class Bat {
public:
  Bat(char const *infileInput);
  std::streamoff get_dofs_begin();
  std::ifstream *get_file();
  unsigned int get_n_frames();
  unsigned int get_n_dihedrals();
  unsigned int get_precision();

private:
  void read_BAT_header();

  std::ifstream infile;
  int precision, n_frames;
  unsigned int n_bonds, n_angles, n_dihedrals;
  std::vector<std::vector<int>> dihedrals_top;
  std::vector<float> masses;
  std::vector<std::string> residues;
  std::vector<int> residue_numbers;
  std::vector<std::string> atom_names;
  std::vector<std::string> belongs_to_molecule;
  std::streamoff dofs_begin;
};

#endif
