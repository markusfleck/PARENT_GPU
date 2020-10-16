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

#include "Bat.h"
#include "My_Error.cpp"

#pragma pack(1)

using namespace std;

Bat::Bat(char const *infile_str) {
    read_BAT_header(infile_str);
    n_bonds = n_dihedrals + 2;
    n_angles = n_dihedrals + 1;
    n_dofs = n_bonds + n_angles + n_dihedrals;
    bonds_frame = new double[n_bonds];
    angles_frame = new double[n_angles];
    dihedrals_frame = new double[n_dihedrals];
}

Bat::~Bat(){
    delete[] bonds_frame;
    delete[] angles_frame;
    delete[] dihedrals_frame;
}


// to read the header of the binary .bat file
void Bat::read_BAT_header(char const *infile_str) {
    try{
        this->infile_string = string(infile_str);    
        infile.exceptions(std::ifstream::badbit);
        infile.open(infile_str, ios::binary); // open the .bat file;
        if (!infile.good()) {
            My_Error my_error((string("ERROR OPENING FILE ") + string(infile_str) +
                        string("! ABORTING.\n"))
                            .c_str());
            throw my_error;
        } 
    infile.exceptions(std::ifstream::badbit | std::ifstream::failbit |
                    std::ifstream::eofbit);

      char dummystring[31];

      infile.read((char *)&version,
                  sizeof(int)); // first read the .bat version number as an integer
      
      infile.read((char *)&precision,
                  sizeof(int)); // then read an integer declaring if the trajectory
                                // is stored in double precision
      
      infile.read(
          (char *)&n_dihedrals,
          sizeof(int)); // read an integer containing the number of dihedrals
      
      infile.read((char *)&n_frames,
                  sizeof(int)); // and an integer containing the number of frames
      

      if (version < 0) {
        ostringstream oss;
        oss << "ERROR: FILE HEADER CORRUPTED. VERSION NUMBER (" << version
            << ") < 0!";
        My_Error my_error(oss.str());
        throw my_error;
      }
      if ((precision != 0) && (precision != 1)) {
        ostringstream oss;
        oss << "ERROR: FILE HEADER CORRUPTED. DOUBLE PRECISION VALUE (" << precision
            << ") NEITHER 0 NOR 1!";
        My_Error my_error(oss.str());
        throw my_error;
      }
      if (n_dihedrals > 19997) {
        cerr << "WARNING: " << n_dihedrals + 3
             << " ATOMS DECLARED IN THE FILE HEADER (CORRUPTED?). THIS WILL LEAD "
                "TO LARGE OUTPUT."
             << endl;
      }
      if (n_dihedrals < 0) {
        ostringstream oss;
        oss << "ERROR: FILE HEADER CORRUPTED. NUMBER OF DIHEDRALS (" << n_dihedrals
            << ") < 0!";
        My_Error my_error(oss.str());
        throw my_error;
      }
      if (n_frames < 1) {
        ostringstream oss;
        oss << "ERROR: FILE HEADER CORRUPTED. NUMBER OF FRAMES (" << n_frames
            << ") < 1!";
        My_Error my_error(oss.str());
        throw my_error;
      }

      if (version >= 3) {
        for (int i = 0; i < n_dihedrals + 3;
             i++) { // for ever atom in the system
          infile.read(
              dummystring,
              8 * sizeof(char)); // read the name of the residue it belongs to
          residues.push_back(dummystring);
          
          residue_numbers.push_back(0);
          infile.read(
              (char *)&(residue_numbers[i]),
              sizeof(float)); // read the number of the residue it belongs to
          
          infile.read(dummystring, 8 * sizeof(char)); // read the name of the atom
          atom_names.push_back(dummystring);
          
          infile.read(
              dummystring,
              31 * sizeof(char)); // read the molecule of the residue it belongs to
          belongs_to_molecule.push_back(dummystring);
          
        }
      }

      vector<int> dummyvec;
      dummyvec.push_back(0);
      dummyvec.push_back(0);
      dummyvec.push_back(0);
      dummyvec.push_back(0);
      dummyvec.push_back(0);
      dummyvec.push_back(0);
      for (int i = 0; i < n_dihedrals; i++) { // then for all dihedrals
        dihedrals_top.push_back(dummyvec);
        infile.read((char *)&(dihedrals_top[i][0]),
                    sizeof(int)); // read the atomnumber of the first atom
        
        infile.read((char *)&(dihedrals_top[i][1]), sizeof(int)); // second atom
        
        infile.read((char *)&(dihedrals_top[i][2]), sizeof(int)); // third atom
        
        infile.read((char *)&(dihedrals_top[i][3]), sizeof(int)); // fourth atom
        
        infile.read((char *)&(dihedrals_top[i][4]),
                    sizeof(int)); // an integer containing the type of the
                                  // dihedral(physical=0,pseudo=1,improper=-1)
        
        infile.read(
            (char *)&(dihedrals_top[i][5]),
            sizeof(int)); // and an integer containing the "parent"dihedral(for
                          // phaseangles, -1 if no "parent")
        
      }
      for (int i = 0; i < n_dihedrals + 3;
           i++) { // and read the whole massestor of the atoms in the system in
                  // single precision (float)
        masses.push_back(0);
        infile.read((char *)&(masses[i]), sizeof(float));
        
      }

      dofs_begin = infile.tellg();
    } catch (My_Error my_error) {
        throw my_error;
        } catch (...) {
        My_Error my_error((string("ERROR WHILE READING THE HEADER OF THE FILE ") +
                       infile_string + string("! ABORTING."))
                          .c_str());
        throw my_error;
    }
}

// to write the header of the binary .bat file
void Bat::write_BAT_header(char const* outfile_str, int version_tmp) { //version 4 is first gbat version
    try{
        this->outfile_string = string(outfile_str);
      outfile.exceptions(std::ifstream::badbit);
      outfile.open(outfile_str, ios::binary); // open the .bat file;
      if (!outfile.good()) {
        My_Error my_error((string("ERROR OPENING FILE ") + string(outfile_str) +
                           string("! ABORTING.\n"))
                              .c_str());
        throw my_error;
      }
      outfile.exceptions(std::ofstream::badbit | std::ofstream::failbit |
                        std::ofstream::eofbit);


        int dummy = dihedrals_top.size();
      char dummystring[31];

      outfile.write(
          (char *)&version_tmp,
          sizeof(int)); // first write the .bat version number as an integer
      
      outfile.write((char *)&precision,
                       sizeof(int)); // then write an integer declaring if the
                                     // trajectory is stored in double precision
      
      outfile.write(
          (char *)&dummy,
          sizeof(int)); // write an integer containing the number of dihedrals
      
      outfile.write(
          (char *)&n_frames,
          sizeof(int)); // and an integer containing the number of frames
      

      for (int i = 0; i < dummy + 3; i++) { // for all atoms in the system
        for (int j = 0; j < 8; j++) {
          dummystring[j] = '\0';
        }
        sprintf(dummystring, "%.7s", residues[i].c_str());
        outfile.write(
            dummystring,
            8 * sizeof(char)); // write the name of the residue it belongs to
        
        outfile.write(
            (char *)&(residue_numbers[i]),
            sizeof(int)); // write the number of the residue it belongs to
        
        for (int j = 0; j < 8; j++) {
          dummystring[j] = '\0';
        }
        sprintf(dummystring, "%.7s",
            (atom_names[i]).c_str()); // write the name of the atom
            outfile.write(dummystring, 8 * sizeof(char));
    
        for (int j = 0; j < 31; j++) {
          dummystring[j] = '\0';
        }
        sprintf(dummystring, "%.30s", belongs_to_molecule[i].c_str());
        outfile.write(
            dummystring,
            31 * sizeof(char)); // write the name of the molecule it belongs to
        
      }

      for (int i = 0; i < dummy; i++) { // then for all dihedrals
        outfile.write((char *)&dihedrals_top[i][0],
                         sizeof(int)); // write the atomnumber of the first atom
        
        outfile.write((char *)&dihedrals_top[i][1],
                         sizeof(int)); // second atom
        
        outfile.write((char *)&dihedrals_top[i][2], sizeof(int)); // third
                                                                        // atom
        
        outfile.write((char *)&dihedrals_top[i][3],
                         sizeof(int)); // fourth atom
        
        outfile.write((char *)&dihedrals_top[i][4],
                         sizeof(int)); // an integer containing the type of the
                                       // dihedral(physical=0,pseudo=1,improper=-1)
        
        outfile.write(
            (char *)&dihedrals_top[i][5],
            sizeof(int)); // and an integer containing the "parent"dihedral(for
                          // phaseangles, -1 if no "parent")
        
      }
      for (int i = 0; i < dummy + 3; i++) {
        outfile.write(
            (char *)&(masses[i]),
            sizeof(float)); // and write the whole massestor of the atoms in the
                            // system in single precision (float)
        
      }
    } catch (My_Error my_error) {
        throw my_error;
    } catch (...) {
        My_Error my_error((string("ERROR WHILE WRITING THE HEADER OF THE FILE ") +
                       outfile_string + string("! ABORTING."))
                          .c_str());
        throw my_error;
    }
}



void Bat::read_BAT_frame() {
    try{
      // to read a frame from the .bat trajectory
      float fdummy;
      int b_counter = 2;
      int a_counter = 1;
      int d_counter = 0;

      infile.read((char *)&time, sizeof(float)); // read the time of the current
                                                   // according .xtc frame as a float
      
      infile.read((char *)&xtc_prec,
                     sizeof(float)); // read the precision of the current frame
                                     // according .xtc frame for back conversion
      
      infile.read(
          (char *)dbox,
          9 * sizeof(float)); // read the box vectors of the current according
                              // according .xtc frame for back conversion
      

      if (precision == 1) { // if double precision is used
        infile.read(
            (char *)root_origin_cartesian,
            3 * sizeof(double)); // read the Cartestians of the first root atom
                                 // (external coordinates) in double precision
        
        infile.read(
            (char *)&root_origin_theta,
            sizeof(double)); // read the polar coordinates of the second root atom
                             // relative to the first (external coordinates)
        
        infile.read((char *)&root_origin_phi, sizeof(double));
        
        infile.read((char *)&root_origin_dihedral,
                       sizeof(double)); // and the dihedral the root atoms form with
                                        // the origin (external coordinates)
        

        infile.read(
            (char *)bonds_frame,
            2 * sizeof(double)); // read the lengths of the two bonds connecting the
                                 // root atoms (internal coordinates)
        
        infile.read((char *)angles_frame,
                       sizeof(double)); // and the angle between the two rootbonds
                                        // (internal coordinates)
        
        for (int i = 0; i < n_dihedrals;
             i++) { // then for all dihedrals in the system
          infile.read((char *)&(bonds_frame[b_counter]),
                         sizeof(double)); // read the bondlength between the last
                                          // two atoms in the dihedral
          
          b_counter++;
          infile.read((char *)&(angles_frame[a_counter]),
                         sizeof(double)); // read the angle between the last threee
                                          // atoms of the dihedral#
          
          a_counter++;
          infile.read((char *)&(dihedrals_frame[d_counter]),
                         sizeof(double)); // and the value of the dihedral itself
          
          d_counter++;
        }
      } else if (precision == 0) { // if single precision is used, do the same but
                                   // use float instead of double
        infile.read((char *)&fdummy, sizeof(float));
        root_origin_cartesian[0] = fdummy;
        infile.read((char *)&fdummy, sizeof(float));
        root_origin_cartesian[1] = fdummy;
        infile.read((char *)&fdummy, sizeof(float));
        root_origin_cartesian[2] = fdummy;
        infile.read((char *)&fdummy, sizeof(float));
        root_origin_theta = fdummy;
        infile.read((char *)&fdummy, sizeof(float));
        root_origin_phi = fdummy;
        infile.read((char *)&fdummy, sizeof(float));
        root_origin_dihedral = fdummy;
        infile.read((char *)&fdummy, sizeof(float));
        bonds_frame[0] = fdummy;
        infile.read((char *)&fdummy, sizeof(float));
        bonds_frame[1] = fdummy;
        infile.read((char *)&fdummy, sizeof(float));
        angles_frame[0] = fdummy;
        
        for (int i = 0; i < n_dihedrals; i++) {
          infile.read((char *)&fdummy, sizeof(float));
          bonds_frame[b_counter] = fdummy;
          b_counter++;
          infile.read((char *)&fdummy, sizeof(float));
          angles_frame[a_counter] = fdummy;
          a_counter++;
          infile.read((char *)&fdummy, sizeof(float));
          dihedrals_frame[d_counter] = fdummy;
          d_counter++;
        }
      }
    } catch (My_Error my_error) {
        throw my_error;
    } catch (...) {
        My_Error my_error((string("ERROR WHILE READING A FRAME FROM THE FILE ") +
                       infile_string + string("! ABORTING."))
                          .c_str());
        throw my_error;
    }
}


void Bat::write_BAT_frame() {
    try{
      // to attach a frame to the .bat trajectory
      outfile.write((char *)&time,
                       sizeof(float)); // write the time of the current frame according
                                       // .xtc frame as a float
      
      outfile.write((char *)&xtc_prec,
                       sizeof(float)); // write the precision of the current
                                       // according .xtc frame for back conversion
      
      outfile.write(
          (char *)dbox,
          9 * sizeof(float)); // write the box vectors of the current according .xtc
                              // frame for back conversion
      

      if (precision == 1) { // if double precision is requested
        outfile.write(
            (char *)root_origin_cartesian,
            3 * sizeof(double)); // write the Cartestians of the first root atom
                                 // (external coordinates) in double precision
        
        outfile.write(
            (char *)&root_origin_theta,
            sizeof(double)); // write the polar coordinates of the second root atom
                             // relative to the first (external coordinates)
        
        outfile.write((char *)&root_origin_phi, sizeof(double));
        
        outfile.write((char *)&root_origin_dihedral,
                         sizeof(double)); // and the dihedral the root atoms form
                                          // with the origin (external coordinates)
        
        outfile.write(
            (char *)bonds_frame,
            2 * sizeof(double)); // write the lengths of the two bonds connecting
                                 // the root atoms (internal coordinates)
        
        outfile.write((char *)angles_frame,
                         sizeof(double)); // and the angle between the two rootbonds
                                          // (internal coordinates)
        
        for (int i = 0; i < n_dihedrals;
             i++) { // then for all dihedrals in the system
          outfile.write((char *)&bonds_frame[i + 2],
                           sizeof(double)); // write the bondlength between the last
                                            // two atoms in the dihedral
          
          outfile.write((char *)&angles_frame[i + 1],
                           sizeof(double)); // write the angle between the last
                                            // three atoms of the dihedral
          
          outfile.write((char *)&dihedrals_frame[i],
                           sizeof(double)); // and the value of the dihedral itself
          
        }
      } else if (precision ==
                 0) { // if single precision is requested do the analogue

        float dummy;
        dummy = root_origin_cartesian[0]; // conversion from double to float is done
                                          // by assigning the double value to a float
                                          // variable
        outfile.write((char *)&dummy, sizeof(float));
        
        dummy = root_origin_cartesian[1];
        outfile.write((char *)&dummy, sizeof(float));
        
        dummy = root_origin_cartesian[2];
        outfile.write((char *)&dummy, sizeof(float));
        
        dummy = root_origin_theta;
        outfile.write((char *)&dummy, sizeof(float));
        
        dummy = root_origin_phi;
        outfile.write((char *)&dummy, sizeof(float));
        
        dummy = root_origin_dihedral;
        outfile.write((char *)&dummy, sizeof(float));
        
        dummy = bonds_frame[0];
        outfile.write((char *)&dummy, sizeof(float));
        
        dummy = bonds_frame[1];
        outfile.write((char *)&dummy, sizeof(float));
        
        dummy = angles_frame[0];
        outfile.write((char *)&dummy, sizeof(float));
        
        for (int i = 0; i < n_dihedrals; i++) {
          dummy = bonds_frame[i + 2];
          outfile.write((char *)&dummy, sizeof(float));
          
          dummy = angles_frame[i + 1];
          outfile.write((char *)&dummy, sizeof(float));
          
          dummy = dihedrals_frame[i];
          outfile.write((char *)&dummy, sizeof(float));
          
        }
      }
    } catch (My_Error my_error) {
        throw my_error;
    } catch (...) {
        My_Error my_error((string("ERROR WHILE WRITING A FRAME TO THE FILE ") +
                       outfile_string + string("! ABORTING."))
                          .c_str());
        throw my_error;
    }
}

void Bat::write_GBAT_frame(unsigned int frame_number){
    try{

        unsigned char inc;
        if (precision == 1) {
            inc = sizeof(double); // if trajectory is in double precision
        } else {
            inc = sizeof(float);
        }

        outfile.seekp(dofs_begin + frame_number *sizeof(float));
        outfile.write((char *)&time, sizeof(float)); // write the time of the current according .xtc frame as a float
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * sizeof(float) + frame_number * sizeof(float));
        outfile.write((char *)&xtc_prec, sizeof(float)); // write the precision of the current frame according .xtc frame for back conversion
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * sizeof(float) + frame_number * sizeof(float) * 9);
        outfile.write((char *)dbox, 9 * sizeof(float)); // write the box vectors of the current frame according .xtc frame for back conversion
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * sizeof(float) * 9 + frame_number * inc);
        outfile.write((char *)root_origin_cartesian, inc); // write the Cartestians of the first root atom (external coordinates) in double precision
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
        outfile.write((char *)&root_origin_cartesian[1], inc);
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
        outfile.write((char *)&root_origin_cartesian[2], inc);
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
        outfile.write((char *)&root_origin_theta, inc); // write the polar coordinates of the second root atom relative to the first (external coordinates)
            
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
        outfile.write((char *)&root_origin_phi, inc);
            
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
        outfile.write((char *)&root_origin_dihedral, inc); // and the dihedral the root atoms form with the origin (external coordinates)
        
        for(int i = 0; i < n_bonds; i++){
            outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
            outfile.write((char *)&bonds_frame[i], inc);
        }
        
        for(int i = 0; i < n_angles; i++){
            outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
            outfile.write((char *)&angles_frame[i], inc);
        }
        
        for(int i = 0; i < n_dihedrals; i++){
            outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
            outfile.write((char *)&dihedrals_frame[i], inc);
        }

    
    } catch (My_Error my_error) {
        throw my_error;
    } catch (...) {
        My_Error my_error((string("ERROR WHILE WRITING A FRAME TO THE FILE ") +
                       outfile_string + string("! ABORTING."))
                          .c_str());
        throw my_error;
    }
}

void Bat::write_GBAT(char const* outfile_str){ //slow conversion due to scattered writes. Fast conversion is done in src/Bat_builder/convert_BAT_to_GBAT.cpp
    try{
        write_BAT_header(outfile_str, 4);
    
        for(int i = 0; i< n_frames; i++){
            read_BAT_frame();
            write_GBAT_frame(i);
        }
    } catch (My_Error my_error) {
        throw my_error;
    } catch (...) {
        My_Error my_error((string("ERROR WHILE WRITING GBAT ") +
                       outfile_string + string("! ABORTING."))
                          .c_str());
        throw my_error;
    }
}


template <class T>
void Bat::load_dofs(T* type_addr[3], int type_id_start[3], int type_id_end[3], unsigned int padding) {
    try{
        if( (padding > 1) && (n_frames % padding > 0) ){
            padding = padding - ( n_frames % padding );
        }
        else{ padding = 0;}
        
        size_t inc;
        if (precision == 1) {
            inc = sizeof(double); // if trajectory is in double precision
        } 
        else {
            inc = sizeof(float);
        }
    
        if(version == 3){
            
            T* bonds = type_addr[TYPE_B];
            T* angles = type_addr[TYPE_A];
            T* dihedrals = type_addr[TYPE_D];
        
            infile.seekg(dofs_begin);
            for (int frame = 0; frame < n_frames; frame++) {
              // to read a frame from the .bat trajectory
              T ddummy[6];
              float fdummy[11];
              int a_start_g = n_bonds;
              int d_start_g = n_bonds + n_angles;
              int b_counter_g = 0;
              int a_counter_g = a_start_g;
              int d_counter_g = d_start_g;
              size_t b_counter_lt = 0;
              size_t a_counter_lt = 0;
              size_t d_counter_lt = 0;

              if (frame % 100000 == 0) {
                cout << "Reading frame " << frame
                     << " and the following.\n"; // every 10000 frames issue an
                                                 // information to stdout
                cout.flush();
              }

              infile.read(
                  (char *)fdummy,
                  11 *
                      sizeof(float)); // read time, precision and box vectors to dummies
              
              
              infile.read((char *)ddummy, 6 * inc); // external coordinates to dummies
              

              infile.read((char *)ddummy,
                         inc); // read the lengths of the two bonds connecting the root
                               // atoms (internal coordinates)
              //~ cout<<"H1"<<endl;
              if ((b_counter_g >= type_id_start[TYPE_B]) &&
                  (b_counter_g <= type_id_end[TYPE_B]))
                bonds[b_counter_lt++ * (n_frames + padding) + frame] = ddummy[0];
              b_counter_g++;
                //~ cout<<"H2"<<endl;
              infile.read((char *)ddummy,
                         inc); // read the lengths of the two bonds connecting the root
                               // atoms (internal coordinates)
              
              if ((b_counter_g >= type_id_start[TYPE_B]) &&
                  (b_counter_g <= type_id_end[TYPE_B]))
                bonds[b_counter_lt++ * (n_frames + padding) + frame] = ddummy[0];
              b_counter_g++;
                //~ cout<<"H3"<<endl;
              infile.read((char *)ddummy, inc); // and the angle between the two
                                               // rootbonds (internal coordinates)
              
              if ((a_counter_g >= type_id_start[TYPE_A]) &&
                  (a_counter_g <= type_id_end[TYPE_A]))
                angles[a_counter_lt++ * (n_frames + padding) + frame] = ddummy[0];
              a_counter_g++;
                //~ cout<<"H4"<<endl;
              for (int i = 0; i < n_dihedrals;
                   i++) {                        // then for all dihedrals in the system
                infile.read((char *)ddummy, inc); // read the bondlength between the last
                                                 // two atoms in the dihedral
                //~ cout<<"H40"<<endl;
                
                if ((b_counter_g >= type_id_start[TYPE_B]) &&
                    (b_counter_g <= type_id_end[TYPE_B]))
                        bonds[b_counter_lt++ * (n_frames + padding) + frame] = ddummy[0];
                b_counter_g++;
                //~ cout<<"H41"<<endl;

                infile.read((char *)ddummy, inc); // read the angle between the last
                                                 // threee atoms of the dihedral#
                
                if ((a_counter_g >= type_id_start[TYPE_A]) &&
                    (a_counter_g <= type_id_end[TYPE_A]))
                  angles[a_counter_lt++ * (n_frames + padding) + frame] = ddummy[0];
                a_counter_g++;
                //~ cout<<"H42"<<endl;

                infile.read((char *)ddummy, inc); // and the value of the dihedral itself
                
                if ((d_counter_g >= type_id_start[TYPE_D]) &&
                    (d_counter_g <= type_id_end[TYPE_D]))
                  dihedrals[d_counter_lt++ * (n_frames + padding) + frame] = ddummy[0];
                d_counter_g++;
                //~ cout<<"H43"<<endl;
              }
              //~ cout<<"H5"<<endl;
            }
            if (padding > 0){
                for(unsigned int type = 0; type < 3; type++){
                    for(int i = 0; i < (type_id_end[type] - type_id_start[type] + 1) ; i++) memset(&type_addr[type][ i * (n_frames + padding) + n_frames], 0.0, padding);
                }
            }
        }
        else if(version == 4){
            if(padding>0){
                for(unsigned int type = 0; type < 3; type++){
                    if ( (type_id_start[type] >= 0) && (type_id_end[type] >= 0) ){
                        infile.seekg(size_t(dofs_begin) + n_frames * (11 * sizeof(float) + (6 + type_id_start[type]) * inc) );
                        for(int i = 0; i < (type_id_end[type] - type_id_start[type] + 1) ; i++){
                            infile.read((char *)(type_addr[type] + i * (n_frames + padding)), size_t(n_frames) * inc);
                            memset(&type_addr[type][ i * (n_frames + padding) + n_frames], 0.0, padding);
                        }
                    }
                }
            }
            else{
                for(unsigned int type = 0; type < 3; type++){
                    if ( (type_id_start[type] >= 0) && (type_id_end[type] >= 0) ){
                        infile.seekg(size_t(dofs_begin) + n_frames * (11 * sizeof(float) + (6 + type_id_start[type]) * inc) );
                        infile.read((char *)type_addr[type], size_t(n_frames) * (type_id_end[type] - type_id_start[type] + 1) * inc);
                    }
                }
            }
        
        }
        else{
            My_Error my_error((string("ERROR WHILE READING BAT ") +
                       outfile_string + string("! VERSION NOT SUPPORTED. ABORTING."))
                          .c_str());
            throw my_error;
        }
    } catch (My_Error my_error) {
        throw my_error;
    } catch (...) {
        My_Error my_error((string("ERROR WHILE READING BAT ") +
                       infile_string + string("! ABORTING."))
                          .c_str());
        throw my_error;
    }
  }

template void Bat::load_dofs(float* type_addr[3], int type_id_start[3], int type_id_end[3], unsigned int padding);
template void Bat::load_dofs(double* type_addr[3], int type_id_start[3], int type_id_end[3], unsigned int padding);
  

template <class T>  
void Bat::load_externals(float* tpd, T* externals) {//---------------------
    try{
    
        unsigned char inc;
        if (precision == 1) {
            inc = sizeof(double); // if trajectory is in double precision
        } 
        else {
            inc = sizeof(float);
        }
        if(version == 3){
            for (int frame = 0; frame < n_frames; frame++) {
              // to read a frame from the .bat trajectory
                T ddummy[6];
                float fdummy[11];

                infile.seekg(dofs_begin + frame * ( (n_dofs + 6) * sizeof(T) + 11 * sizeof(float) ) );

              if (frame % 100000 == 0) {
                cout << "Reading frame " << frame
                     << " and the following.\n"; // every 10000 frames issue an
                                                 // information to stdout
                cout.flush();
              }

              infile.read(
                  (char *)fdummy, 11 * sizeof(float)); // read time, precision and box vectors to dummies
              
              for(unsigned int i = 0; i < 2; i++){
                tpd[i * n_frames + frame] = fdummy[i];
              }
              memcpy(&tpd[2 * n_frames + frame * 9], &fdummy[2], 9 *sizeof(float));
            
              infile.read((char *)ddummy, 6 * inc); // external coordinates to dummies
              
              for(unsigned int i = 0; i < 6; i++){
                externals[i * n_frames + frame] = ddummy[i];
              }
            }
        }
        else if(version == 4){
            infile.seekg(dofs_begin);
            infile.read((char *)tpd, n_frames * 11 * sizeof(float)); // read time, precision and box vectors to dummies
            infile.read((char *)externals, n_frames * 6 * inc); // read external coordinates
        }
        else{
            My_Error my_error((string("ERROR WHILE READING BAT ") +
                       infile_string + string("! VERSION NOT SUPPORTED. ABORTING."))
                          .c_str());
            throw my_error;
        }
        
    } catch (My_Error my_error) {
        throw my_error;
    } catch (...) {
        My_Error my_error((string("ERROR WHILE READING BAT ") +
                       infile_string + string("! ABORTING."))
                          .c_str());
        throw my_error;
    }
  }

template void Bat::load_externals(float* tpd, float* externals);
template void Bat::load_externals(float* tpd, double* externals);





streamoff Bat::get_dofs_begin() { return dofs_begin; }

ifstream *Bat::get_infile() { infile.seekg(dofs_begin); return &infile; }

ofstream *Bat::get_outfile() { outfile.seekp(dofs_begin); return &outfile; }

int Bat::get_n_frames() { return n_frames; }
int Bat::get_n_frames_padded(unsigned int padding) { 
    
    unsigned int full = n_frames / padding;
    full *= padding;
    if(n_frames - full > 0){
        return full + padding;
    }
    return n_frames; 
}

int Bat::get_n_bonds() { return n_bonds; }
int Bat::get_n_angles() { return n_angles; }
int Bat::get_n_dihedrals() { return n_dihedrals; }
int Bat::get_n_dofs() { return n_dofs; }

int Bat::get_precision() { return precision; }

void Bat::set_precision(int precision){this->precision = precision; };

#pragma pack(0)
