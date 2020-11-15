// A program analysing the mutual information between two residues
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

#include <iostream>

#include "../util/classes/Residue_Representation.h"
#include "../util/classes/Arg_Parser.h"

using namespace std;


int main(int argc, char* argv[]){
    int counter;

    Arg_Parser arg_parser(argc, argv);
    unsigned int resid1;
    unsigned int resid2;
  
    if( !( arg_parser.exists( string("-f") ) && arg_parser.exists( string("-resid1") )  && arg_parser.exists( string("-resid2") )&& (argc==7) ) ){
        cerr<<"USAGE:\n"<<argv[0]<<" -f input.par -resid1 #residue1 -resid2 #residue2"<<endl;
        return 1;
    }
    
    if ( strcmp( arg_parser.get_ext( arg_parser.get("-f") ) , "par") ) {
    // check for the extensions of the input file
    cerr<<"USAGE:\n"<<argv[0]<<" -f input.par -resid1 #residue1 -resid2 #residue2"<<endl;
    exit(EXIT_FAILURE);
  }
  if (sscanf(arg_parser.get("-resid1"), "%ud", &resid1) != 1) {
    // read the residue number and check for correctness
    cerr << "ERROR: Could not read the first residue number from command line! Aborting"
         << endl;
    exit(EXIT_FAILURE);
  }
  if (sscanf(arg_parser.get("-resid2"), "%ud", &resid2) != 1) {
    // read the residue number and check for correctness
    cerr << "ERROR: Could not read the second residue number from command line! Aborting"
         << endl;
    exit(EXIT_FAILURE);
  }
  

    Residue_Representation rep(arg_parser.get("-f"), false, MODE_TOTAL);
    Entropy_Matrix* mat = rep.getEntropy_Matrix();
    
    
    cout<<rep.getResidueName(resid1)<<rep.getResidueNumber(resid1)<<":"<<rep.getMoleculeName(resid1)<<endl<<endl;
    
    vector< int > atoms1 = rep.getAtoms(resid1);
    vector <vector< int > >dofs1;
    dofs1.push_back(rep.getBonds(resid1));
    dofs1.push_back(rep.getAngles(resid1));
    dofs1.push_back(rep.getDihedrals(resid1));
  
    struct Max_Ent_Dof{
        string type;
        int id;
        double value;
    } max_ent_dof1 = {string("bond"), 1, mat->getEntropy(TYPE_B, dofs1[TYPE_B][0])};

    cout<<"Atoms:"<<endl;
    for(int i: atoms1){
        cout<<rep.getAtomName(i)<<endl;
    }
    cout<<endl;
    
    counter = 1;
    cout<<"Bonds:"<<endl;
    double b_ent1 = 0;
    for(int i: dofs1[TYPE_B]){
        cout<<counter<<" "<<rep.getAtomName( mat->getBondAtom(i,1) )<<" "<<rep.getAtomName( mat->getBondAtom(i,2) )<<" "<<mat->getEntropy(TYPE_B, i)<<endl;
        b_ent1 += mat->getEntropy(TYPE_B, i);
        if(mat->getEntropy(TYPE_B, i) > max_ent_dof1.value){
            max_ent_dof1 = {string("bond"), counter, mat->getEntropy(TYPE_B, i)};
        }
        counter++;
    }
    cout<<endl;
    
    counter = 1;
    cout<<"Angles:"<<endl;
    double a_ent1 = 0;
    for(int i: dofs1[TYPE_A]){
        cout<<counter<<" "<<rep.getAtomName( mat->getAngleAtom(i,1) )<<" "<<rep.getAtomName( mat->getAngleAtom(i,2) )<<" "<<rep.getAtomName( mat->getAngleAtom(i,3) )<<" "<<mat->getEntropy(TYPE_A, i)<<endl;
        a_ent1 += mat->getEntropy(TYPE_A, i);
        if(mat->getEntropy(TYPE_A, i) > max_ent_dof1.value){
            max_ent_dof1 = {string("angle"), counter, mat->getEntropy(TYPE_A, i)};
        }
        counter++;
    }
    cout<<endl;
    
    counter = 1;
    cout<<"Dihedrals:"<<endl;
    double d_ent1 = 0;
    for(int i: dofs1[TYPE_D]){
        cout<<counter<<" "<<rep.getAtomName( mat->getDihedralAtom(i,1) )<<" "<<rep.getAtomName( mat->getDihedralAtom(i,2) )<<" "<<rep.getAtomName( mat->getDihedralAtom(i,3) )<<" "<<rep.getAtomName( mat->getDihedralAtom(i,4) )<<" "<<mat->getEntropy(TYPE_D, i)<<endl;
        d_ent1 += mat->getEntropy(TYPE_D, i);
        if(mat->getEntropy(TYPE_D, i) > max_ent_dof1.value){
            max_ent_dof1 = {string("dihedral"), counter, mat->getEntropy(TYPE_D, i)};
        }
        counter++;
    }
    cout<<endl<<endl<<endl<<endl;
    
    
    
    
    
    cout<<rep.getResidueName(resid2)<<rep.getResidueNumber(resid2)<<":"<<rep.getMoleculeName(resid2)<<endl<<endl;
    
    vector< int > atoms2 = rep.getAtoms(resid2);
    vector <vector< int > >dofs2;
    dofs2.push_back(rep.getBonds(resid2));
    dofs2.push_back(rep.getAngles(resid2));
    dofs2.push_back(rep.getDihedrals(resid2));
  
    Max_Ent_Dof max_ent_dof2 = {string("bond"), 1, mat->getEntropy(TYPE_B, dofs2[TYPE_B][0])};
  
    cout<<"Atoms:"<<endl;
    for(int i: atoms2){
        cout<<rep.getAtomName(i)<<endl;
    }
    cout<<endl;
    
    counter = 1;
    cout<<"Bonds:"<<endl;
    double b_ent2 = 0;
    for(int i: dofs2[TYPE_B]){
        cout<<counter<<" "<<rep.getAtomName( mat->getBondAtom(i,1) )<<" "<<rep.getAtomName( mat->getBondAtom(i,2) )<<" "<<mat->getEntropy(TYPE_B, i)<<endl;
        b_ent2 += mat->getEntropy(TYPE_B, i);
        if(mat->getEntropy(TYPE_B, i) > max_ent_dof2.value){
            max_ent_dof2 = {string("bond"), counter, mat->getEntropy(TYPE_B, i)};
        }
        counter++;
    }
    cout<<endl;
    
    counter = 1;
    cout<<"Angles:"<<endl;
    double a_ent2 = 0;
    for(int i: dofs2[TYPE_A]){
        cout<<counter<<" "<<rep.getAtomName( mat->getAngleAtom(i,1) )<<" "<<rep.getAtomName( mat->getAngleAtom(i,2) )<<" "<<rep.getAtomName( mat->getAngleAtom(i,3) )<<" "<<mat->getEntropy(TYPE_A, i)<<endl;
        a_ent2 += mat->getEntropy(TYPE_A, i);
        if(mat->getEntropy(TYPE_A, i) > max_ent_dof2.value){
            max_ent_dof2 = {string("angle"), counter, mat->getEntropy(TYPE_A, i)};
        }
        counter++;
    }
    cout<<endl;
    
    counter = 1;
    cout<<"Dihedrals:"<<endl;
    double d_ent2 = 0;
    for(int i: dofs2[TYPE_D]){
        cout<<counter<<" "<<rep.getAtomName( mat->getDihedralAtom(i,1) )<<" "<<rep.getAtomName( mat->getDihedralAtom(i,2) )<<" "<<rep.getAtomName( mat->getDihedralAtom(i,3) )<<" "<<rep.getAtomName( mat->getDihedralAtom(i,4) )<<" "<<mat->getEntropy(TYPE_D, i)<<endl;
        d_ent2 += mat->getEntropy(TYPE_D, i);
        if(mat->getEntropy(TYPE_D, i) > max_ent_dof2.value){
            max_ent_dof2 = {string("dihedral"), counter, mat->getEntropy(TYPE_D, i)};
        }
        counter++;
    }
    cout<<endl<<endl<<endl<<endl;
    
    cout<<"Mutual information:"<<endl;
    double bb_mut = 0;
    double ba_mut = 0;
    double bd_mut = 0;
    double aa_mut = 0;
    double ad_mut = 0;
    double dd_mut = 0;
    
    struct Max_Mut_Dof_Pair{
        string type1;
        string type2;
        int id1;
        int id2;
        double value;
    } max_mut_dof_pair = {string("bond"), string("bond"), 1, 1, mat->getMutual(TYPE_B, TYPE_B, dofs1[TYPE_B][0], dofs2[TYPE_B][0])};
    
    
    for (unsigned char type1 = 0; type1 < 3; type1++){ // for bonds, angles and dihedrals of the first member of the dof pair
        for(unsigned int idx1 = 0; idx1 < dofs1[type1].size(); idx1++){ // and all dofs of the current type of the first member of the dof pair
            for (unsigned char type2 = 0; type2 < 3; type2++){ // and all "later" types for the second member of the dof pair
                for(unsigned int idx2 = 0; idx2 < dofs2[type2].size(); idx2++ ){ // and all "later" dofs for the second member of the dof pair
                    
                    double mutual = mat->getMutual(type1, type2, dofs1[type1][idx1], dofs2[type2][idx2]);
                    
                    string type1_str("");
                    string type2_str("");
                    switch(type1){
                          case TYPE_B:
                          type1_str = string("bond");
                            switch(type2){
                              case TYPE_B:
                               type2_str = string("bond");
                                bb_mut += mutual;
                                break;
                              case TYPE_A:
                                type2_str = string("angle");
                                ba_mut += mutual;
                                break;
                              case TYPE_D:
                                type2_str = string("dihedral");
                                bd_mut += mutual;
                                break;
                            }
                            break;
                          case TYPE_A:
                            type1_str = string("angle");
                            switch(type2){
                              case TYPE_B:
                               type2_str = string("bond");
                                ba_mut += mutual;
                                break;
                              case TYPE_A:
                                type2_str = string("angle");
                                aa_mut += mutual;
                                break;
                              case TYPE_D:
                                type2_str = string("dihedral");
                                ad_mut += mutual;
                                break;
                            }
                            break;
                          case TYPE_D:
                            type1_str = string("dihedral");
                            switch(type2){
                              case TYPE_B:
                               type2_str = string("bond");
                                bd_mut += mutual;
                                break;
                              case TYPE_A:
                                type2_str = string("angle");
                                ad_mut += mutual;
                                break;
                              case TYPE_D:
                                type2_str = string("dihedral");
                                dd_mut += mutual;
                                break;
                            }
                            break;
                    }
                
                    if(mutual > max_mut_dof_pair.value){
                        max_mut_dof_pair.type1 = type1_str;
                        max_mut_dof_pair.type2 = type2_str;
                        max_mut_dof_pair.id1 = idx1 + 1;
                        max_mut_dof_pair.id2 = idx2 + 1;
                        max_mut_dof_pair.value = mutual;
                    }
                    
                    cout<<type1_str<<" "<< idx1 + 1<<" "<<type2_str<<" "<< idx2 + 1<<" "<<mutual<<endl;
                }
            }
        }
    }
    cout<<endl;
    
    
    cout<<"Highest entropy degree of freedom residue "<<resid1<<":"<<endl;
    cout<<max_ent_dof1.type<<" "<<max_ent_dof1.id<<" "<<max_ent_dof1.value<<endl;
    cout<<endl;
    
    cout<<"Highest entropy degree of freedom residue "<<resid2<<":"<<endl;
    cout<<max_ent_dof2.type<<" "<<max_ent_dof2.id<<" "<<max_ent_dof2.value<<endl;
    cout<<endl;
    
    cout<<"Highest mutual information degree of freedom pair (only between residues):"<<endl;
    cout<<max_mut_dof_pair.type1<<" "<<max_mut_dof_pair.id1<<" "<<max_mut_dof_pair.type2<<" "<<max_mut_dof_pair.id2<<" "<<max_mut_dof_pair.value<<endl;
    cout<<endl;
    
    
    cout << "TOTAL 1D BONDS ENTROPY RESIDUE "<<resid1<< " = " << b_ent1 << endl;
    cout << "TOTAL 1D ANGLES ENTROPY RESIDUE "<<resid1<< " = " << a_ent1 << endl;
    cout << "TOTAL 1D DIHEDRALS ENTROPY RESIDUE "<<resid1<< " = " << d_ent1 << endl;
    cout<<endl;
    cout << "TOTAL 1D BONDS ENTROPY RESIDUE "<<resid2<< " = " << b_ent1 << endl;
    cout << "TOTAL 1D ANGLES ENTROPY RESIDUE "<<resid2<< " = " << a_ent1 << endl;
    cout << "TOTAL 1D DIHEDRALS ENTROPY RESIDUE "<<resid2<< " = " << d_ent1 << endl;
    cout<<endl;
    cout << "TOTAL 2D BONDS-BONDS MUTUAL INFORMATION BETWEEN RESIDUES = " << bb_mut<< endl;
    cout << "TOTAL 2D BONDS-ANGLES MUTUAL INFORMATION BETWEEN RESIDUES = " << ba_mut<< endl;
    cout << "TOTAL 2D BONDS-DIHEDRALS MUTUAL INFORMATION BETWEEN RESIDUES = " << bd_mut<< endl;
    cout << "TOTAL 2D ANGLES-ANGLES MUTUAL INFORMATION BETWEEN RESIDUES = " << aa_mut<< endl;
    cout << "TOTAL 2D ANGLES-DIHEDRALS MUTUAL INFORMATION BETWEEN RESIDUES = " << ad_mut<< endl;
    cout << "TOTAL 2D DIHEDRALS-DIHEDRALS MUTUAL INFORMATION BETWEEN RESIDUES = "<< dd_mut << endl;
    cout << "TOTAL MUTUAL INFORMATION BETWEEN RESIDUES = "<< bb_mut + ba_mut + bd_mut + aa_mut + ad_mut + dd_mut << endl;
    

    return 0;
}
