// The header file for the class for handling an Entropy_Matrix from a .par file of the PARENT
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


#ifndef ENTROPY_MATRIX_H
#define ENTROPY_MATRIX_H

#define PRECISION double

#include "../util.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class Entropy_Matrix {
public:
  Entropy_Matrix(char const *infileInput);
  Entropy_Matrix(unsigned int nAtoms);
  Entropy_Matrix(char const *bat_file, PRECISION *storage, unsigned int n_bins);
  ~Entropy_Matrix();

  PRECISION getEntropy(
      int type,
      unsigned int index); // get the 1D entropy of a BAT degree of freedom: e.
                           // g. getEntropy(TYPE_A, 27) gives the entropy of the
                           // 27th angle (indexing starts at 1)
  PRECISION get2DEntropy(
      int type1, int type2, unsigned int index1,
      unsigned int index2); // get the 2D entropy of 2 degrees of freedom:
                            // get2DEntropy(TYPE_B, TYPE_D, 12, 38) gives the 2D
                            // entropy of the 12th bond and the 38th dihedral
                            // (indexing starts at 1)
  PRECISION getMutual(
      int type1, int type2, unsigned int index1,
      unsigned int index2); // get the mutual information between 2 degrees of
                            // freedom: getMutual(TYPE_A, TYPE_B, 27, 9) gives
                            // the mutual information between the 27th angle and
                            // the 9th bond (indexing starts at 1)

  void
  setEntropy(int type, unsigned int index,
             PRECISION value); // set the 1D entropy of a BAT degree of freedom: e.
                            // g. setEntropy(TYPE_A, 27) sets the entropy of the
                            // 27th angle (indexing starts at 1)
  void setEntropy(unsigned int dof_id,
                  PRECISION value); // the same, but dof_int starts at zero and is
                                 // global, i. e. ranges from 0 to n_atoms - 1
  void
  set2DEntropy(int type1, int type2, unsigned int index1, unsigned int index2,
               PRECISION value); // set the 2D entropy of 2 degrees of freedom:
                              // set2DEntropy(TYPE_B, TYPE_D, 12, 38) sets the 2D
                              // entropy of the 12th bond and the 38th dihedral
                              // (indexing starts at 1)
  void setMutual(int type1, int type2, unsigned int index1, unsigned int index2,
                 PRECISION value); // set the mutual information between 2 degrees
                                // of freedom: setMutual(TYPE_A, TYPE_B, 27, 9)
                                // sets the mutual information between the 27th
                                // angle and the 9th bond (indexing starts at 1)

  void write(char const *outfileName); // writes the Entropy_Matrix using a .par
                                       // file format

  unsigned int getNBonds();     // returns the number of bonds in the system
  unsigned int getNAngles();    // returns the number of angles in the system
  unsigned int getNDihedrals(); // returns the number of dihedrals in the system

  std::string getResidueName(
      unsigned int atomNumber); // get the name of the residue of the atom with
                                // the according number (atomnumbers start at 1)
  int getResidueNumber(
      unsigned int atomNumber); // get the name of the residue of the atom with
                                // the according number (atomnumbers start at 1)
  std::string getAtomName(
      unsigned int atomNumber); // get the name of the atom with
                                // the according number (atomnumbers start at 1)
  std::string getMoleculeName(
      unsigned int atomNumber); // get the name of the residue of the atom with
                                // the according number (atomnumbers start at 1)

  int getBondAtom(
      unsigned int bondNumber,
      unsigned int atom); // returns the atomnumber of an atom in a given bond,
                          // e.g. getBondAtom(334,2) gives the number of the
                          // second atom in bond 334 (all indices start at 1)
  int getAngleAtom(
      unsigned int anglelNumber,
      unsigned int atom); // returns the atomnumber of an atom in a given angle,
                          // e.g. getAngleAtom(575,1) gives the number of the
                          // first atom in angle 575 (all indices start at 1)
  int getDihedralAtom(
      unsigned int dihedralNumber,
      unsigned int
          atom); // returns the atomnumber of an atom in a given dihedral, e.g.
                 // getDihedralAtom(727,3) gives the number of the third atom in
                 // dihedral 727 (all indices start at 1)

  void setPseudoZero(); // set mutual information terms involving pseudo degrees
                        // of freedom zero;
  std::streamoff
  get_bat_file_dofs_begin(); // if a .bat file was read, get the location in the
                             // file where the header has ended and the body
                             // starts

private:
  void write_PAR_header();
  void read_PAR_header();
  void write_PAR_body();
  void read_PAR_body();
  void read_BAT_header();

  int bDens1D, aDens1D, dDens1D, bDens, aDens, dDens;
  int double_prec, numFrames, version;
  unsigned int nBonds, nAngles, nDihedrals;
  PRECISION *bondsEntropy1D, *anglesEntropy1D, *dihedralsEntropy1D;
  PRECISION *bbEntropy, *baEntropy, *bdEntropy, *aaEntropy, *adEntropy, *ddEntropy;
  std::vector< std::vector<int> > dihedrals_top;
  std::vector<float> masses;
  std::vector<std::string> residues;
  std::vector<int> residueNumbers;
  std::vector<std::string> atomNames;
  std::vector<std::string> belongsToMolecule;
  std::string infilestring;
  std::ifstream infile;
  std::ofstream outfile;
  std::streamoff bat_file_dofs_begin;
};

#endif
