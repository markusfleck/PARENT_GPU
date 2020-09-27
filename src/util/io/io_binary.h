//    The binary IO library for the PARENT program suite
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

#ifndef IO_BINARY_H
#define IO_BINARY_H

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../classes/Entropy_Matrix.h"
#include "../types.h"

int write_BAT_header(std::ofstream *outfile, int double_prec, int numframes,
                     std::vector<std::vector<int>> *dihedrals_top,
                     std::vector<float> *masses,
                     std::vector<std::string> *residues,
                     std::vector<int> *residueNumbers,
                     std::vector<std::string> *atomNames,
                     std::vector<std::string> *belongsToMolecule);
int read_BAT_header(std::ifstream *infile, int *double_prec, int *numFrames,
                    std::vector<std::vector<int>> *dihedrals_top,
                    std::vector<float> *masses,
                    std::vector<std::string> *residues,
                    std::vector<int> *residueNumbers,
                    std::vector<std::string> *atomNames,
                    std::vector<std::string> *belongsToMolecule);
int write_BAT_frame(std::ofstream *outfile, int double_prec, int nDihedrals,
                    float time, float xtcPrec, float **box,
                    double *root_origin_cartesian, double root_origin_theta,
                    double root_origin_phi, double root_origin_dihedral,
                    double *bonds, double *angles, double *dihedrals);
int read_BAT_frame(std::ifstream *infile, int precision, int nDihedrals,
                   float *time, float *xtcPrec, float **dbox,
                   double *root_origin_cartesian, double *root_origin_theta,
                   double *root_origin_phi, double *root_origin_dihedral,
                   double *my_bonds, double *my_angles, double *my_dihedrals);

int write_PAR_header(std::ofstream *outfile, int nDihedrals, int double_prec,
                     int numFrames,
                     std::vector<std::vector<int>> *dihedrals_top,
                     std::vector<float> *masses, int bDens1D, int aDens1D,
                     int dDens1D, int bDens, int aDens, int dDens,
                     std::vector<std::string> *residues,
                     std::vector<int> *residueNumbers,
                     std::vector<std::string> *atomNames,
                     std::vector<std::string> *belongsToMolecule);
int read_PAR_header(std::ifstream *infile, int *nDihedrals, int *double_prec,
                    int *numFrames,
                    std::vector<std::vector<int>> *dihedrals_top,
                    std::vector<float> *masses, int *version, int *bDens,
                    int *aDens, int *dDens, int *bDens1D, int *aDens1D,
                    int *dDens1D, std::vector<std::string> *residues,
                    std::vector<int> *residueNumbers,
                    std::vector<std::string> *atomNames,
                    std::vector<std::string> *belongsToMolecule);
int write_PAR_body(std::ofstream *par_file, int nDihedrals,
                   double *bondsEntropy1D, double *anglesEntropy1D,
                   double *dihedralsEntropy1D, double *bbEntropy,
                   double *baEntropy, double *bdEntropy, double *aaEntropy,
                   double *adEntropy, double *ddEntropy);
int read_PAR_body(std::ifstream *par_file, int nDihedrals,
                  double **bondsEntropy1D, double **anglesEntropy1D,
                  double **dihedralsEntropy1D, double **bbEntropy,
                  double **baEntropy, double **bdEntropy, double **aaEntropy,
                  double **adEntropy, double **ddEntropy);

#endif
