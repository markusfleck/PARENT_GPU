// BAT_builder, a program to convert a molecular dynamics trajectory from Cartesian to internal bond-angle-torsion coordinates
// Copyright (C) 2020  Markus Fleck
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 3 as 
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "topology.h"
#include "BAT_topology.h"
#include "BAT_trajectory.h"
#include "../util/io/io.h"

#include <iostream>
#include <cstdlib>
#include <cstring>


using namespace std;

int main(int argc, char* argv[])
{
    
    
    Arg_Parser arg_parser(argc, argv);
    if(argc==9){ //conversion from Cartesians to BAT in double precision
					if(!arg_parser.exists("-t")||!arg_parser.exists("-x")||!arg_parser.exists("-o")||!arg_parser.exists("-bb")){
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
						return 1;
					}
        vector <std::string> backboneAtomNames;

        //check if command line arguments were supplied correctly, else put out help and quit
        if(((strcmp(arg_parser.get_ext(arg_parser.get("-t")),"top"))||(strcmp(arg_parser.get_ext(arg_parser.get("-x")),"xtc"))||(strcmp(arg_parser.get_ext(arg_parser.get("-o")),"bat")))) {
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
            return 1;
        }

        //Read in the requested names of backbone atoms
        char* ptr = strtok (arg_parser.get("-bb")," ");
        while (ptr != NULL)
        {
            backboneAtomNames.push_back(string(ptr));
            ptr = strtok (NULL, " ");
        }

				//Read the topology file and build a table for the molecules bonds as well as a list of backbone atoms
        cout<<"Reading topology."<<endl;
        Topology proteins(string(arg_parser.get("-t")).c_str(), backboneAtomNames);

        //Build the BAT Topology using the molecules bonds (and the backbone atoms to make use of phase angles)
        cout<<"Building BAT topology."<<endl;
        BAT_Topology bat(&proteins.bonds_table,&proteins.backbone,&proteins.roots);

        //Use the list of dihedrals and the original trajectory file to build the BAT binary trajectory in double precision (masses are included in the .bat file since they might come in handy at a later point). Also add information about the atoms.
        cout<<"Writing .bat trajectory."<<endl;
        BAT_Trajectory trj(string(arg_parser.get("-x")).c_str(),string(arg_parser.get("-o")).c_str(), &bat.dihedrals, &proteins.masses,&proteins.residues,&proteins.residueNumbers,&proteins.atomNames,&proteins.belongsToMolecule);
		
		}
    else if(argc==10) { //conversion from Cartesians to BAT in single precision
					if(!arg_parser.exists("-t")||!arg_parser.exists("-x")||!arg_parser.exists("-o")||!arg_parser.exists("-bb")||!arg_parser.exists("--single_precision")){
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
						return 1;
					}
                
        vector <std::string> backboneAtomNames;

        //check if command line arguments were supplied correctly, else put out help and quit
        if(((strcmp(arg_parser.get_ext(arg_parser.get("-t")),"top"))||(strcmp(arg_parser.get_ext(arg_parser.get("-x")),"xtc"))||(strcmp(arg_parser.get_ext(arg_parser.get("-o")),"bat")))) {
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
            return 1;
        }

        //Read in the requested names of backbone atoms
        char* ptr = strtok (arg_parser.get("-bb")," ");
        while (ptr != NULL)
        {
            backboneAtomNames.push_back(string(ptr));
            ptr = strtok (NULL, " ");
        }


        //Read the topology file and build a table for the molecules bonds as well as a list of backbone atoms
        cout<<"Reading topology."<<endl;
        Topology proteins(string(arg_parser.get("-t")).c_str(), backboneAtomNames);

        //Build the BAT Topology using the molecules bonds (and the backbone atoms to make use of phase angles)
        cout<<"Building BAT topology."<<endl;
        BAT_Topology bat(&proteins.bonds_table,&proteins.backbone,&proteins.roots);

        //Use the list of dihedrals and the original trajectory file to build the BAT binary trajectory in single precision (masses are included in the .bat file since they might come in handy at a later point). Also add information about the atoms.
        cout<<"Writing .bat trajectory."<<endl;
        BAT_Trajectory trj(string(arg_parser.get("-x")).c_str(),string(arg_parser.get("-o")).c_str(), &bat.dihedrals,&proteins.masses,&proteins.residues,&proteins.residueNumbers,&proteins.atomNames,&proteins.belongsToMolecule,0);

    }
    else  if (argc==5) { // conversion from BAT to Cartesians
					if(!arg_parser.exists("-b")||!arg_parser.exists("-o")){
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
						return 1;
					}

        //check if command line arguments were supplied correctly, else put out help and quit
        if(((strcmp(arg_parser.get_ext(arg_parser.get("-b")),"bat"))||(strcmp(arg_parser.get_ext(arg_parser.get("-o")),"xtc")))) {
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
            return 1;
        }


        //do the conversion
        cout<<"Converting .bat to .xtc."<<endl;
        BAT_Trajectory trj(string(arg_parser.get("-b")).c_str(),string(arg_parser.get("-o")).c_str());

    }
    else {
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
        return 1;
    }


    cout<<"Finished writing trajectory."<<endl;
    exit(EXIT_SUCCESS);
}















