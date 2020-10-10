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
#include <string>
#include "../types.h"

class Bat {
public:
  Bat(char const *infile_str);
    ~Bat();  
    std::streamoff get_dofs_begin();
  std::ifstream *get_file();
  int get_n_frames();
  int get_n_dihedrals();
  int get_precision();
    void set_precision(int precision);
    void write_GBAT(char const* outfile_str);
    template <class T>
    void load_dofs(T* type_addr[3], int type_id_start[3], int type_id_end[3]); 
private:
  void read_BAT_header(char const *infile_str);
  void write_BAT_header(char const* outfile_str, int version);
    void read_BAT_frame();
    void write_BAT_frame();
    void write_GBAT_frame(unsigned int frame_number);


  std::ifstream infile;
    std::ofstream outfile;
    std::string infile_string;
    std::string outfile_string;
  int precision, n_frames;
  int n_bonds, n_angles, n_dihedrals;
  std::vector<std::vector<int>> dihedrals_top;
  std::vector<float> masses;
  std::vector<std::string> residues;
  std::vector<int> residue_numbers;
  std::vector<std::string> atom_names;
  std::vector<std::string> belongs_to_molecule;
  std::streamoff dofs_begin;

    float time, xtc_prec;
    float dbox[3][3];
    double root_origin_cartesian[3];
    double root_origin_theta, root_origin_phi, root_origin_dihedral;
    double *bonds_frame, *angles_frame, *dihedrals_frame;
    int version;


};

#endif
