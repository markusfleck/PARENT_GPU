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

      int version;
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
void Bat::write_BAT_header(char const* outfile_str, int version) { //version 4 is first gbat version
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
          (char *)&version,
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
        
        outfile.seekp(size_t(size_t(outfile.tellp())) + (n_frames - frame_number - 1) * sizeof(float) + frame_number * sizeof(float));
        outfile.write((char *)&xtc_prec, sizeof(float)); // write the precision of the current frame according .xtc frame for back conversion
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * sizeof(float) + frame_number * sizeof(float) * 9);
        outfile.write((char *)dbox, 9 * sizeof(float)); // write the box vectors of the current frame according .xtc frame for back conversion
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * sizeof(float) * 9 + frame_number * inc * 3);
        outfile.write((char *)root_origin_cartesian, 3 * inc); // write the Cartestians of the first root atom (external coordinates) in double precision
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc * 3 + frame_number * inc);
        outfile.write((char *)&root_origin_theta, inc); // write the polar coordinates of the second root atom relative to the first (external coordinates)
            
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
        outfile.write((char *)&root_origin_phi, inc);
            
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
        outfile.write((char *)&root_origin_dihedral, inc); // and the dihedral the root atoms form with the origin (external coordinates)
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
        outfile.write((char *)bonds_frame, inc); // write the lengths of the two bonds connecting the root atoms (internal coordinates)
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
        outfile.write((char *)&bonds_frame[1], inc); 
        
        outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
        outfile.write((char *)angles_frame, inc); // and the angle between the two rootbonds (internal coordinates)
        
        for (int i = 0; i < n_dihedrals; i++) { // then for all dihedrals in the system
        
            outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
            outfile.write((char *)&bonds_frame[i + 2], inc); // write the bondlength between the last two atoms in the dihedral
            
            outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
            outfile.write((char *)&angles_frame[i + 1], inc); // write the angle between the last three atoms of the dihedral
            
            outfile.seekp(size_t(outfile.tellp()) + (n_frames - frame_number - 1) * inc + frame_number * inc);
            outfile.write((char *)&dihedrals_frame[i], inc); // and the value of the dihedral itself
              
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

void Bat::write_GBAT(char const* outfile_str){
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









streamoff Bat::get_dofs_begin() { return dofs_begin; }

ifstream *Bat::get_file() { return &infile; }

int Bat::get_n_frames() { return n_frames; }

int Bat::get_n_dihedrals() { return n_dihedrals; }

int Bat::get_precision() { return precision; }

void Bat::set_precision(int precision){this->precision = precision; };

#pragma pack(0)
