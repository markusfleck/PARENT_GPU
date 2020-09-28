//    The utility library for the PARENT program suite
//    Copyright (C) 2016  Markus Fleck (member of the laboratory of Bojan
//    Zagrovic, University of Vienna)
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

//    A scientific publication about this program has been released in the
//    Journal of Chemical Theory and Computation:
//		"PARENT: A Parallel Software Suite for the Calculation of
//Configurational Entropy in Biomolecular Systems" 		DOI: 10.1021/acs.jctc.5b01217

//   We kindly ask you to include this citation in works that publish
//   results generated using this program or any modifications of it.

#include "types.h"
#include <string>

char *getCmdOption(char **begin, char **end, const std::string &option);
bool cmdOptionExists(char **begin, char **end, const std::string &option);
unsigned char get_dof_type_from_id(unsigned int dof_id,
                                   unsigned int n_dihedrals);
unsigned int get_min_id_for_type(unsigned char type, unsigned int n_dihedrals);
unsigned int get_max_id_for_type(unsigned char type, unsigned int n_dihedrals);