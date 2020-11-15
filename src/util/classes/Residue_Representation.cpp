// A class transforming an Entropy_Matrix class of the PARENT program suite to residue representation
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



#include "Residue_Representation.h"
    #include<iostream>

using namespace std;


Residue_Representation::Residue_Representation(char const * infileInput, bool include_full, int mode){
    this->mode = mode;

    if((mode!=MODE_TOTAL)&&(mode!=MODE_AVER)&&(mode!=MODE_MAX)){
        My_Error my_error(string("ERROR: REQUEST FOR AN UNKNOWN RESIDUE REPRESENTATION MODE(")+to_string(mode)+string(").").c_str());
        throw my_error;
    }
    mat = new Entropy_Matrix(infileInput);//load the .par file into the Entropy_Matrix object
    int nDihedrals=mat->getNDihedrals();
    int found;
  
    calcBonds=1;
    calcAngles=1;
    calcDihedrals=1;

    //partition the atoms into residues (groups)
    vector<int> tmpintvec;
    tmpintvec.clear();
    tmpintvec.push_back(0);


    for(int i=1; i<=nDihedrals+3; i++){//for all atoms
            found=0;
            for(unsigned int j=0; j<residueNames.size(); j++){
                if( (mat->getResidueName(i)==residueNames[j]) && (mat->getResidueNumber(i)==residueNumbers[j]) && (mat->getMoleculeName(i)==moleculeNames[j]) ){//if the according residue has already been created
                    found=1;
                    groups[j].push_back(i);//add the atom to the according residue
                    break;
                }
            }
            if(!found){//else
                residueNames.push_back(mat->getResidueName(i));//add the new residues name to the residueNames array
                residueNumbers.push_back(mat->getResidueNumber(i)); //also add the new residue number to the 
                moleculeNames.push_back(mat->getMoleculeName(i));//and the new molecule name of the new residue to the moleculeNames array
                tmpintvec[0]=i;
                groups.push_back(tmpintvec); //start a group for the new residue containing the atomnumbers (with indices starting at 1) 
            }
    }
    nResidues = groups.size();
    
    //assign bonds to bondindices[k] if both atoms are part of residue group[k]    
    tmpintvec.clear();
    tmpintvec.push_back(-1);//create a dummy vector to start the 2D vector  
    bool lead_flag = false; 
    if(calcBonds){
      for(unsigned int k=0; k<groups.size(); k++){
        bondIndices.push_back(tmpintvec);//and use it as a dummy for bondIndices[k] initialization
        for(int j=1; j<=nDihedrals+2; j++){ //in order to check for all bonds if they are in residue group[k]
          found=0;
          bool lead_flag = false;
          for(unsigned int i=0; i<groups[k].size(); i++){//check if the first atom of bond j matches any atom in groups[k]
            if(groups[k][i] == mat->getBondAtom(j,1)){
              found++;
              break;   
            }
          }
          for(unsigned int i=0; i<groups[k].size(); i++){//check if the second atom of bond j matches any atom in groups[k]
            if(groups[k][i] == mat->getBondAtom(j,2)){
              found++;
              lead_flag = true;
              break;   
            }
          }
          if((!include_full)&&lead_flag){bondIndices[k].push_back(j);} //if the leading atom if the bond is in the residue, add the bond to bondIndices[k] 
          if(include_full&&(found==2)){bondIndices[k].push_back(j);} //if the include_full option was set, only include the degree of freedom if all its atoms are part of the residue
        }
       bondIndices[k].erase(bondIndices[k].begin()); //remove the first dummy tmpintvec which was used for bondIndices[k] initialization
      }
    }
        
        
        
    //assign angles to angleindices[k] if all three atoms are part of residue group[k]    
    tmpintvec.clear();
    tmpintvec.push_back(-1);//create a dummy vector to start the 2D vector 
    if(calcAngles){
      for(unsigned int k=0;k<groups.size();k++){
        angleIndices.push_back(tmpintvec);//and use it as a dummy for angleIndices[k] initialization
        for(int j=1;j<=nDihedrals+1;j++){ //in order to check for all angles if they are in residue group[k]
          found=0;
          lead_flag = false;
          for(unsigned int i=0;i<groups[k].size();i++){//check if the first atom of angle j matches any atom in groups[k]
            if(groups[k][i]==mat->getAngleAtom(j,1)){
              found++;
              break;   
            }
          }
          for(unsigned int i=0;i<groups[k].size();i++){//check if the second atom of angle j matches any atom in groups[k]
            if(groups[k][i]==mat->getAngleAtom(j,2)){
              found++;
              break;   
            }
          }
                    for(unsigned int i=0;i<groups[k].size();i++){//check if the third atom of angle j matches any atom in groups[k]
            if(groups[k][i]==mat->getAngleAtom(j,3)){
              found++;
              lead_flag = true;
              break;   
            }
          }
          if((!include_full)&&lead_flag){angleIndices[k].push_back(j);} //if the leading atom if the angle is in the residue, add the angle to angleIndices[k] 
          if(include_full&&(found==3)){angleIndices[k].push_back(j);} //if the include_full option was set, only include the degree of freedom if all its atoms are part of the residue
        }
       angleIndices[k].erase(angleIndices[k].begin()); //remove the first dummy tmpintvec which was used for angleIndices[k] initialization
      }
    }
        
        
        
    //assign dihedrals to dihedralindices[k] if all four atoms are part of residue group[k]    
    tmpintvec.clear();
    tmpintvec.push_back(-1);//create a dummy vector to start the 2D vector
    if(calcDihedrals){
      for(unsigned int k=0;k<groups.size();k++){
        dihedralIndices.push_back(tmpintvec);//and use it as a dummy for dihedralIndices[k] initialization
        for(int j=1;j<=nDihedrals;j++){ //in order to check for all dihedrals if they are in residue group[k]
          found=0;
          lead_flag = false;
          for(unsigned int i=0;i<groups[k].size();i++){//check if the first atom of dihedral j matches any atom in groups[k]
            if(groups[k][i]==mat->getDihedralAtom(j,1)){
              found++;
              break;   
            }
          }
          for(unsigned int i=0;i<groups[k].size();i++){//check if the second atom of dihedral j matches any atom in groups[k]
            if(groups[k][i]==mat->getDihedralAtom(j,2)){
              found++;
              break;   
            }
          }
            for(unsigned int i=0;i<groups[k].size();i++){//check if the third atom of dihedral j matches any atom in groups[k]
            if(groups[k][i]==mat->getDihedralAtom(j,3)){
              found++;
              break;   
            }
          }
            for(unsigned int i=0;i<groups[k].size();i++){//check if the fourth atom of dihedral j matches any atom in groups[k]
            if(groups[k][i]==mat->getDihedralAtom(j,4)){
              found++;
              lead_flag = true;
              break;   
            }
          }
          if((!include_full)&&lead_flag){dihedralIndices[k].push_back(j);} //if the leading atom if the dihedral is in the residue, add the dihedral to dihedralIndices[k] 
          if(include_full&&(found==4)){dihedralIndices[k].push_back(j);} //if the include_full option was set, only include the degree of freedom if all its atoms are part of the residue
        }
       dihedralIndices[k].erase(dihedralIndices[k].begin()); //remove the first dummy tmpintvec which was used for dihedralIndices[k] initialization
      }
    }
    
    for(unsigned int i=0;i<nResidues;i++){
      nBondsVec.push_back(bondIndices[i].size());
      nAnglesVec.push_back(angleIndices[i].size());
      nDihedralsVec.push_back(dihedralIndices[i].size());
    }
}

void Residue_Representation::calculate_matrix(){
    //create a mutual information matrix between all residues
    mutualArray = new double[nResidues*(nResidues-1)/2];
    int counter=0;
    for(unsigned int i=0;i<nResidues-1;i++){
      for(unsigned int j=i+1;j<nResidues;j++){
        double totalBbMutual=0;
        double totalBaMutual=0;
        double totalBdMutual=0;
        double totalAaMutual=0;
        double totalAdMutual=0;
        double totalDdMutual=0;
        
        int averCounter=0;
        double highestMutual=0;
        double mutual;
        
        //for every found bond-bond pair between the residues calculate the mutual information
        for(unsigned int k=0;k<bondIndices[i].size();k++){
          for(unsigned int l=0;l<bondIndices[j].size();l++){
            mutual=mat->getMutual(TYPE_B,TYPE_B,bondIndices[i][k],bondIndices[j][l]);
            totalBbMutual+=mutual;
            if(mutual>highestMutual){
              highestMutual=mutual;
            }
            averCounter++;   
          }    
        }
        
        //same for all found bond-angle pairs
        for(unsigned int k=0;k<bondIndices[i].size();k++){
          for(unsigned int l=0;l<angleIndices[j].size();l++){
            mutual=mat->getMutual(TYPE_B,TYPE_A,bondIndices[i][k],angleIndices[j][l]); 
            totalBaMutual+=mutual;
            if(mutual>highestMutual){
              highestMutual=mutual;
            }
            averCounter++;            
          }
        }

        //bond-dihedral pairs
        for(unsigned int k=0;k<bondIndices[i].size();k++){
          for(unsigned int l=0;l<dihedralIndices[j].size();l++){
            mutual=mat->getMutual(TYPE_B,TYPE_D,bondIndices[i][k],dihedralIndices[j][l]);
            totalBdMutual+=mutual;
            if(mutual>highestMutual){
              highestMutual=mutual;
            }
            averCounter++;    
          }
        }
        
        
        //angle-bond pairs
        for(unsigned int k=0;k<angleIndices[i].size();k++){
          for(unsigned int l=0;l<bondIndices[j].size();l++){
            mutual=mat->getMutual(TYPE_A,TYPE_B,angleIndices[i][k],bondIndices[j][l]);
            totalBaMutual+=mutual;
            if(mutual>highestMutual){
              highestMutual=mutual;
            }
            averCounter++;      
          }
        }
        
        //angle-angle pairs
        for(unsigned int k=0;k<angleIndices[i].size();k++){
          for(unsigned int l=0;l<angleIndices[j].size();l++){
            mutual=mat->getMutual(TYPE_A,TYPE_A,angleIndices[i][k],angleIndices[j][l]);
            totalAaMutual+=mutual;
            if(mutual>highestMutual){
              highestMutual=mutual;
            }
            averCounter++;                
          }    
        }

        //angle-dihedral pairs
        for(unsigned int k=0;k<angleIndices[i].size();k++){
          for(unsigned int l=0;l<dihedralIndices[j].size();l++){
            mutual=mat->getMutual(TYPE_A,TYPE_D,angleIndices[i][k],dihedralIndices[j][l]);
            totalAdMutual+=mutual;
            if(mutual>highestMutual){
              highestMutual=mutual;
            }
            averCounter++;                
          }
        }
        
        
        //dihedral-bond pairs
        for(unsigned int k=0;k<dihedralIndices[i].size();k++){
          for(unsigned int l=0;l<bondIndices[j].size();l++){
            mutual=mat->getMutual(TYPE_D,TYPE_B,dihedralIndices[i][k],bondIndices[j][l]);
            totalBdMutual+=mutual;
            if(mutual>highestMutual){
              highestMutual=mutual;
            }
            averCounter++;                
          }
        }
        
        //dihedral-angle pairs
        for(unsigned int k=0;k<dihedralIndices[i].size();k++){
          for(unsigned int l=0;l<angleIndices[j].size();l++){
            mutual=mat->getMutual(TYPE_D,TYPE_A,dihedralIndices[i][k],angleIndices[j][l]); 
            totalAdMutual+=mutual;
            if(mutual>highestMutual){
              highestMutual=mutual;
            }
            averCounter++;                
          }
        }

        //dihedral-dihedral pairs
        for(unsigned int k=0;k<dihedralIndices[i].size();k++){
          for(unsigned int l=0;l<dihedralIndices[j].size();l++){
            mutual=mat->getMutual(TYPE_D,TYPE_D,dihedralIndices[i][k],dihedralIndices[j][l]);
            totalDdMutual+=mutual;
            if(mutual>highestMutual){
              highestMutual=mutual;
            }
            averCounter++;                
          }    
        }
        
        if(mode==MODE_TOTAL){
          mutualArray[counter]=totalBbMutual+totalBaMutual+totalBdMutual+totalAaMutual+totalAdMutual+totalDdMutual; // the sum of all terms is is the mutual information between residues 
        }
        if(mode==MODE_AVER){
          mutualArray[counter]=(totalBbMutual+totalBaMutual+totalBdMutual+totalAaMutual+totalAdMutual+totalDdMutual)/averCounter; // the sum of all terms is is the mutual information between residues 
        }
        if(mode==MODE_MAX){
          mutualArray[counter]=highestMutual; // the sum of all terms is is the mutual information between residues 
        }
        counter++;
       }
    }
    matrix_calculated = true;
}


Residue_Representation::~Residue_Representation(){

    delete mat;
    if(matrix_calculated) delete[] mutualArray;
}


//gives the name of the residue with index "residueIndex" (indexing starts at 1) 
std::string Residue_Representation::getResidueName(unsigned int residueIndex){
    if((residueIndex<1)||(residueIndex>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR THE NAME OF A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex)+string(").").c_str());
        throw my_error;
    }
  return residueNames[residueIndex-1];  
}


//gives the number (according to the used topology file) of the residue with index "residueIndex" (indexing starts at 1) 
int Residue_Representation::getResidueNumber(unsigned int residueIndex){
    if((residueIndex<1)||(residueIndex>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR THE NUMBER OF A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex)+string(").").c_str());
        throw my_error;
    }
  return residueNumbers[residueIndex-1];  
}


//gives the name of the molecule the residue with index "residueIndex" belongs to (indexing starts at 1) 
std::string Residue_Representation::getMoleculeName(unsigned int residueIndex){ 
    if((residueIndex<1)||(residueIndex>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR THE NAME OF THE MOLECULE OF A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex)+string(").").c_str());
        throw my_error;
    }
  return moleculeNames[residueIndex-1];  
}


//returns the number of residues
unsigned int Residue_Representation::getNResidues(){
  return nResidues;
}


//returns the number of bonds in the residue with index "residueIndex" (indexing starts at 1) 
unsigned int Residue_Representation::getNBonds(unsigned int residueIndex){
  if((residueIndex<1)||(residueIndex>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR THE NUMBER OF BONDS OF A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex)+string(").").c_str());
        throw my_error;
    }
  return nBondsVec[residueIndex-1];  
}

//returns the number of angles in the residue with index "residueIndex" (indexing starts at 1) 
unsigned int Residue_Representation::getNAngles(unsigned int residueIndex){
  if((residueIndex<1)||(residueIndex>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR THE NUMBER OF ANGLES OF A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex)+string(").").c_str());
        throw my_error;
    }
  return nAnglesVec[residueIndex-1];  
}

//returns the number of dihedrals in the residue with index "residueIndex" (indexing starts at 1) 
unsigned int Residue_Representation::getNDihedrals(unsigned int residueIndex){
  if((residueIndex<1)||(residueIndex>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR THE NUMBER OF DIHEDRALS OF A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex)+string(").").c_str());
        throw my_error;
    }
  return nDihedralsVec[residueIndex-1];  
}

vector< int > Residue_Representation::getAtoms(unsigned int residueIndex){
  if((residueIndex<1)||(residueIndex>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR THE ATOMS OF A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex)+string(").").c_str());
        throw my_error;
    }
  return groups[residueIndex-1];  
}

vector< int > Residue_Representation::getBonds(unsigned int residueIndex){
  if((residueIndex<1)||(residueIndex>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR THE BONDS OF A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex)+string(").").c_str());
        throw my_error;
    }
  return bondIndices[residueIndex-1];  
}

vector< int > Residue_Representation::getAngles(unsigned int residueIndex){
  if((residueIndex<1)||(residueIndex>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR THE ANGLES OF A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex)+string(").").c_str());
        throw my_error;
    }
  return angleIndices[residueIndex-1];  
}

vector< int > Residue_Representation::getDihedrals(unsigned int residueIndex){
  if((residueIndex<1)||(residueIndex>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR THE DIHEDRALS OF A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex)+string(").").c_str());
        throw my_error;
    }
  return dihedralIndices[residueIndex-1];  
}

std::string Residue_Representation::getAtomName(unsigned int atomNumber){
    return mat->getAtomName(atomNumber);
}


//returns the mutual information according to the mode set in the constructor between the residues (indexing starts at 1)
double Residue_Representation::getMutual(unsigned int residueIndex1,unsigned int residueIndex2){
  if((residueIndex1<1)||(residueIndex1>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR A MUTUAL INFORMATION VALUE INVOLVING A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex1)+string(").").c_str());
        throw my_error;
    }
  if((residueIndex2<1)||(residueIndex2>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST FOR A MUTUAL INFORMATION VALUE INVOLVING A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex2)+string(").").c_str());
        throw my_error;
    }
  
  int smaller,bigger;

  if(residueIndex1!=residueIndex2) {//for mutual information between two different residues calculate the index in the half-matrix and return it
      smaller=residueIndex1<residueIndex2?residueIndex1:residueIndex2;
      bigger=residueIndex1<residueIndex2?residueIndex2:residueIndex1;
      return mutualArray[smaller*(2*nResidues-(smaller+1))/2+bigger-1-nResidues];
  }
  
  return 0; //if the mutual information between a residue and itself is requested return 0
  
}



//sets the mutual information between the according residues (indexing starts at 1)
void Residue_Representation::setMutual(unsigned int residueIndex1,unsigned int residueIndex2, double value){
  if((residueIndex1<1)||(residueIndex1>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST TO SET A MUTUAL INFORMATION VALUE INVOLVING A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex1)+string(").").c_str());
        throw my_error;
    }
  if((residueIndex2<1)||(residueIndex2>nResidues)){		
        My_Error my_error(string("ERROR: REQUEST TO SET A MUTUAL INFORMATION VALUE INVOLVING A RESIDUE WITH AN INDEX OUT OF RANGE (")+to_string(residueIndex2)+string(").").c_str());
        throw my_error;
    }
  if(residueIndex1==residueIndex2) {
    My_Error my_error(string("ERROR: REQUEST TO SET A MUTUAL INFORMATION VALUE INVOLVING A RESIDUE WITH ITSELF (")+to_string(residueIndex2)+string(").").c_str());
        throw my_error;
    }
  
  //for mutual information between two different residues calculate the index in the half-matrix and use it for setting
  int smaller,bigger;
  smaller=residueIndex1<residueIndex2?residueIndex1:residueIndex2;
  bigger=residueIndex1<residueIndex2?residueIndex2:residueIndex1;
  mutualArray[smaller*(2*nResidues-(smaller+1))/2+bigger-1-nResidues]=value;
}


Entropy_Matrix* Residue_Representation::getEntropy_Matrix(){
    return mat;
}

