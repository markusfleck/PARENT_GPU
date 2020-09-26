//    The utility library for the PARENT program suite
//    Copyright (C) 2016  Markus Fleck (member of the laboratory of Bojan Zagrovic, University of Vienna)
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





//    A scientific publication about this program has been released in the Journal of Chemical Theory and Computation:
//		"PARENT: A Parallel Software Suite for the Calculation of Configurational Entropy in Biomolecular Systems"
//		DOI: 10.1021/acs.jctc.5b01217

//   We kindly ask you to include this citation in works that publish
//   results generated using this program or any modifications of it.







#include "util.h"
#include <algorithm>
#include <iostream>

using namespace std;


char* getCmdOption(char ** begin, char ** end, const string & option)
{
    char ** itr = find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const string& option)
{
    return find(begin, end, option) != end;
}

unsigned char get_dof_type_from_id(unsigned int dof_id, unsigned int n_dihedrals)
{
	if(dof_id < n_dihedrals + 2) return TYPE_B;
	if(dof_id < 2 * n_dihedrals + 3) return TYPE_A;
	return TYPE_D;
}



unsigned int get_min_id_for_type(unsigned char type, unsigned int n_dihedrals)
{
	switch(type) {
  		case TYPE_B:
    			return 0;
  		case TYPE_A:
    			return n_dihedrals + 2;
		case TYPE_D:
    			return 2 * n_dihedrals + 3;
	}
	return 42;
}



unsigned int get_max_id_for_type(unsigned char type, unsigned int n_dihedrals)
{
	switch(type) {
  		case TYPE_B:
    			return n_dihedrals + 1;
  		case TYPE_A:
    			return 2 * n_dihedrals + 2;
		case TYPE_D:
    			return 3 * n_dihedrals + 2;
	}
	return 42;
}

