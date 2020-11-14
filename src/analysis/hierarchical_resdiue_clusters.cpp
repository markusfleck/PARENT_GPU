// A program identifying residue interaction networks based on mutual information and hierarchical clustering
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
#include <algorithm>



#include "../util/classes/Residue_Representation.h"
#include "../util/classes/Arg_Parser.h"

using namespace std;


double getClosenessMax(vector <unsigned int> &group1, vector <unsigned int> &group2, Residue_Representation *rep, unsigned int minDist){
  
  minDist*=minDist;
  double maxMut=0;
  for (unsigned int i=0;i<group1.size();i++){
    for (unsigned int j=0;j<group2.size();j++){
      if(((group1[i]-group2[j])*(group1[i]-group2[j])>minDist)||(rep->getMoleculeName(group1[i])!=rep->getMoleculeName(group2[j]))){
        if(rep->getMutual(group1[i],group2[j])>maxMut){
          maxMut=rep->getMutual(group1[i],group2[j]);
        }
      }
    }
  }
  return maxMut; 
}

double getClosenessAver(vector <unsigned int> &group1, vector <unsigned int> &group2, Residue_Representation *rep, unsigned int minDist){
  
  minDist*=minDist;
  double aver=0;
  for (unsigned int i=0;i<group1.size();i++){
    for (unsigned int j=0;j<group2.size();j++){
      if(((group1[i]-group2[j])*(group1[i]-group2[j])>minDist)||(rep->getMoleculeName(group1[i])!=rep->getMoleculeName(group2[j]))){
          aver+=rep->getMutual(group1[i],group2[j]);
        }
    }
  }
  return aver/(group1.size()*group2.size()); 
}


int main(int argc, char* argv[]){
  double perc;
  unsigned int minDist; //minimum residues BETWEEN resdiues for mutual information to be taken into account (eliminating nearest neighbor effects)
  
  vector< vector<unsigned int> >clusters;
  vector <unsigned int> tmpintvec;
  tmpintvec.push_back(0);
  vector <double> MIvec;
  
  ofstream vmdFile;
  
  double cutoff=-1;

    Arg_Parser arg_parser(argc, argv);
  
    if( !( ( arg_parser.exists( string("-f") ) && arg_parser.exists( string("-perc") ) && arg_parser.exists( string("-dist") ) && arg_parser.exists( string("-clustermode") ) && arg_parser.exists( string("-residuemode") ) )
        && ( (argc==11) || ( (argc==15) && arg_parser.exists( string("-vmd") ) && arg_parser.exists( string("-gro") ) ) ) ) ){
        cerr<<"USAGE:\n"<<argv[0]<<" -p input.par -perc CutoffPercenatage -dist MinimumResidueDistance -clustermode MAX|AVER -residuemode MAX|AVER|TOTAL [-vmd output.vmd -gro input.gro]"<<endl;
        return 1;
    }
    char* inputFilename = arg_parser.get("-f");
    char* vmdFilename = arg_parser.get("-vmd");
    char* groFilename = arg_parser.get("-gro");
    string clusterMode( arg_parser.get("-clustermode") );
    string residueMode( arg_parser.get("-residuemode") );
    
    if(((clusterMode!=string("MAX"))&&(clusterMode!=string("AVER")))||((residueMode!=string("MAX"))&&(residueMode!=string("AVER"))&&(residueMode!=string("TOTAL")))){
        cerr<<"USAGE:\n"<<argv[0]<<" -p input.par -perc CutoffPercenatage -dist MinimumResidueDistance -clustermode MAX|AVER -residuemode MAX|AVER|TOTAL [-vmd output.vmd -gro input.gro]"<<endl;
        return 1;
    }
    
    if(sscanf(arg_parser.get("-perc"),"%lf",&perc)!=1) {
        cerr<<"ERROR: COULD NOT READ THE PERCENTAGE FOR CUTOFF CALCULATION!"<<endl;
        return 1;
    }
    if(perc<0){
      cerr<<"ERROR: PERCANTAGE FOR CUTOFF CALCULATION SMALLER ZERO!"<<endl;
      return 1;
    }
    if(perc>100){
      cerr<<"ERROR: PERCANTAGE FOR CUTOFF CALCULATION GREATER 100!"<<endl;
      return 1;
    }
    perc/=100;
    
    if(sscanf(arg_parser.get("-dist"),"%ud",&minDist)!=1) {
      cerr<<"ERROR: COULD NOT READ THE MINIMUM RESIDUE DISTANCE!"<<endl;
      return 1;
    }
    if(minDist<0){
      cerr<<"ERROR: MINIMUM RESIDUE DISTANCE SMALLER ZERO!"<<endl;
      return 1;
    }
  

  if(arg_parser.exists( string("-vmd") ) ){
    vmdFile.open(vmdFilename,ios::out);
  }
  
  
  Residue_Representation* rep;//create a residue representation of the mutual information matrix
  
  if(residueMode==string("MAX")){
    rep = new Residue_Representation(inputFilename,true,MODE_MAX);
  }
  else if(residueMode==string("AVER")){
    rep = new Residue_Representation(inputFilename,true,MODE_AVER);
  }
  else if(residueMode==string("TOTAL")){
    rep = new Residue_Representation(inputFilename,true,MODE_TOTAL);
  }
  else{
    cerr<<"ERROR: UNKNOWN RESIDUE REPRESENTATION MODE. ABORTING!"<<endl;
    return 1;
  }
  
  rep->calculate_matrix();
  unsigned int nRes=rep->getNResidues();


  for (unsigned int i=1;i<=nRes;i++){//initialize the clusters with one element (residue) each
    tmpintvec[0]=i;
    clusters.push_back(tmpintvec);
    for(unsigned int j=i+1;j<=nRes;j++){
        MIvec.push_back(rep->getMutual(i,j));
    }
  }
  
  int index; //find the absolute cutoff value for the relative "perc"
  int count=int(perc*MIvec.size());
  for(int k=0;k<count;k++){//by finding and erasing the highest MI value int(perc*MIvec.size()) times, we get the int(perc*MIvec.size()) highest MI value as a cutoff
    cutoff=0;
    index=0;
    for (unsigned int i=0;i<MIvec.size();i++){
      if(MIvec[i]>cutoff){
        cutoff=MIvec[i];
        index=i;
      }
    }
    MIvec.erase (MIvec.begin()+index);
  }
  if(cutoff==0){cerr<<"WARNING: A CUTOFF VALUE OF ZERO IS SET! THE WHOLE MOLECULAR SYSTEM WILL BE PUT IN ONE BIG CLUSTER!"<<endl;}

  
  
  
  double maxCloseness;
  double closeness;
  unsigned int min1 = -1, min2 = -1;
  
  
  
  while(clusters.size()>1){//perform hierarchical clustering
    maxCloseness=-1;
    for (unsigned int i=0;i<clusters.size();i++){
      for (unsigned int j=i+1;j<clusters.size();j++){
        if(clusterMode==string("MAX")){
          closeness=getClosenessMax(clusters[i],clusters[j],rep,minDist);
        }
        else if(clusterMode==string("AVER")){
          closeness=getClosenessAver(clusters[i],clusters[j],rep,minDist);
        }
        else{
          cerr<<"ERROR: UNKNOWN CLUSTERING MODE. ABORTING!"<<endl;
          return 1;
        }
        if(closeness>maxCloseness){
          maxCloseness=closeness;
          min1=i;
          min2=j;
        }
      }
    }

    if(maxCloseness<cutoff){break;}
    clusters[min1].insert(clusters[min1].end(),clusters[min2].begin(),clusters[min2].end());
    clusters.erase (clusters.begin()+min2);
  }

  for(unsigned int i=0;i< clusters.size();i++){
    if(clusters[i].size()>1){
      sort(clusters[i].begin(),clusters[i].end());
      for(unsigned int j=0;j<clusters[i].size();j++){
        cout<<clusters[i][j]<<" ("<<rep->getResidueName(clusters[i][j])<<rep->getResidueNumber(clusters[i][j])<<":"<<rep->getMoleculeName(clusters[i][j])<<")"<<endl;
      }
      cout<<endl<<endl;
    }
  }
  


  if( arg_parser.exists( string("-vmd") ) ){
    vmdFile<<"mol new "<<groFilename<<endl;
    vmdFile<<endl;
    vmdFile<<"mol delrep 0 top"<<endl;
    vmdFile<<"mol representation NewCartoon"<<endl;
    vmdFile<<"mol addrep top"<<endl;

				int color_counter=1;
    for(unsigned int i=0;i< clusters.size();i++){
      if(clusters[i].size()>1){
        vmdFile<<"set sel [atomselect top \"resid";
        for(unsigned int j=0;j<clusters[i].size();j++){
          vmdFile<<" "<<clusters[i][j];
        }
        vmdFile<<"\"]"<<endl;
        vmdFile<<"mol selection \"[$sel text]\""<<endl;
        vmdFile<<"mol representation VDW"<<endl;
        vmdFile<<"mol color ColorID "<<color_counter<<endl;
        vmdFile<<"mol addrep top"<<endl;
								color_counter++;
							}
    }
    vmdFile<<endl<<endl;
    vmdFile.close();
  }
  cout<<"FINISHED SUCCESSFULLY."<<endl;
  return 0;
  
}






