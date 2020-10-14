#include <iostream>
#include <sstream>
#include "../util/classes/Bat.h"
#include "../util/util.h"
#include "../util/io/io.h"
#include "../util/classes/My_Error.cpp"

using namespace std;

int main(int argc, char *argv[]) {
    try{

      if (argc != 7) {
        cerr << "USAGE: " << argv[0] << " -f input.bat -o ouput.gbat --ram #GiB\n";
        exit(EXIT_FAILURE);
      }

      Arg_Parser arg_parser(argc, argv);
      if (!arg_parser.exists("-f") ||
          !arg_parser.exists("-o")) {
        // check for correct command line options
        cerr << "USAGE: " << argv[0] << " -f input.bat -o ouput.gbat --ram #GiB\n";
        exit(EXIT_FAILURE);
      }

      if (strcmp(arg_parser.get_ext(arg_parser.get("-f")),
                 "bat") ||
          strcmp(arg_parser.get_ext(arg_parser.get("-o")),
                 "gbat")) {
        // check for the extensions of the input and output file
        cerr << "USAGE: " << argv[0] << " -f input.bat -o ouput.gbat --ram #GiB\n";
        exit(EXIT_FAILURE);
      }
      
      Bat bat(arg_parser.get("-f"));
      cout<<"Converting "<<arg_parser.get("-f")<< " to GBAT."<<endl;
      //~ bat.write_GBAT(arg_parser.get("-o"));
      
      bat.write_BAT_header(arg_parser.get("-o"), 4);
      
      unsigned int inc = (bat.get_precision() == 0) ? sizeof(float) : sizeof(double);
      size_t n_frames = bat.get_n_frames();
      
      unsigned int n_dofs = bat.get_n_dofs();
      unsigned int n_dihedrals = bat.get_n_dihedrals();
      
      stringstream ram_str(arg_parser.get("--ram"));
      double ram_provided;
      ram_str >> ram_provided;
      size_t cpu_ram_available = static_cast<size_t>(1024) * 1024 * 1024 * ram_provided;
      
      unsigned int n_dofs_load = cpu_ram_available / (n_frames * inc);
      size_t externals_bytes = n_frames * (11 * sizeof(float) + 6 * inc);
      if(externals_bytes > cpu_ram_available){
        My_Error my_error((string("ERROR: OUT OF MEMORY! YOU NEED TO PROVIDE AT LEAST ") +
                       to_string(externals_bytes/1024/1024/1024) + string(" GiB TO STORE THE EXTERNAL DEGREES OF FREEDOM (MORE RAM WILL SIGNIFICANTLY INCREASE PERFORMANCE). ABORTING."))
                          .c_str());
        throw my_error;
      } 
      
      if(n_dofs_load > n_dofs) {
        n_dofs_load = n_dofs;
        cpu_ram_available = n_dofs_load * n_frames * inc;
        if(cpu_ram_available < externals_bytes) cpu_ram_available = externals_bytes;
      }
      char* mem = new char[cpu_ram_available];

      ofstream* outfile = bat.get_outfile();
      
      if(bat.get_precision() == 0){ //TODO: check if the RAM is even sufficient to hold the externals
        bat.load_externals((float*) mem, (float*) &mem[n_frames * 11 * sizeof(float)]);
        outfile->write(mem, n_frames * 17 *sizeof(float));
      }
      else if(bat.get_precision() == 1){
        bat.load_externals((float*) mem, (double*) &mem[n_frames * 11 * sizeof(float)]);
        outfile->write(mem, n_frames * 11 * sizeof(float));
        outfile->write(&mem[n_frames * 11 * sizeof(float)], n_frames * 6 * sizeof(double));
      }
      else{
            My_Error my_error((string("ERROR WHILE READING BAT ") +
                       arg_parser.get("-f") + string("! UNKNOWN PRECISION. ABORTING."))
                          .c_str());
            throw my_error;
      }
      
      unsigned int dof_id_start_g = 0;
      unsigned int dof_id_end_g = n_dofs_load - 1;
      if(dof_id_end_g > n_dofs - 1) dof_id_end_g = n_dofs - 1;
      while( dof_id_start_g < n_dofs ){
          cout<<endl<<"Processing dofs "<<dof_id_start_g + 1<<" to "<<dof_id_end_g + 1<<" from a total of "<<n_dofs<<"."<<endl;
          int type_id_start[3] = {-1, -1, -1};
          int type_id_end[3] = {-1, -1, -1}; // inclusive
          unsigned int type_n_dofs[3] = {0, 0, 0};
          unsigned int dofs_loaded = 0;

          for (unsigned short type = get_dof_type_from_id(dof_id_start_g, n_dihedrals);
               type <= get_dof_type_from_id(dof_id_end_g, n_dihedrals); type++) {
            if (dof_id_start_g < get_min_id_for_type(type, n_dihedrals)) {
              type_id_start[type] = get_min_id_for_type(type, n_dihedrals);
            } else {
              type_id_start[type] = dof_id_start_g;
            }

            if (dof_id_end_g > get_max_id_for_type(type, n_dihedrals)) {
              type_id_end[type] = get_max_id_for_type(type, n_dihedrals);
            } else {
              type_id_end[type] = dof_id_end_g;
            }
            type_n_dofs[type] = type_id_end[type] - type_id_start[type] + 1;
            dofs_loaded += type_n_dofs[type];
          }

          if(bat.get_precision() == 0){
            float *type_addr[3];
            type_addr[TYPE_B] = (float*) mem;
            type_addr[TYPE_A] = (float*) mem + n_frames * type_n_dofs[TYPE_B];
            type_addr[TYPE_D] = (float*) mem + n_frames * (type_n_dofs[TYPE_B] + type_n_dofs[TYPE_A]);
            
            bat.load_dofs(type_addr, type_id_start, type_id_end);
            outfile->write(mem, dofs_loaded * n_frames * sizeof(float));
          }
          else{
            double *type_addr[3];
            type_addr[TYPE_B] = (double*) mem;
            type_addr[TYPE_A] = (double*) mem + n_frames * type_n_dofs[TYPE_B];
            type_addr[TYPE_D] = (double*) mem + n_frames * (type_n_dofs[TYPE_B] + type_n_dofs[TYPE_A]);
            
            bat.load_dofs(type_addr, type_id_start, type_id_end);
            outfile->write(mem, dofs_loaded * n_frames * sizeof(double));
          }
          
          dof_id_start_g += n_dofs_load;
          dof_id_end_g += n_dofs_load;
          if(dof_id_end_g > n_dofs - 1) dof_id_end_g = n_dofs - 1;
          
        }
        cout << "GBAT written successsfully." << endl << endl << endl;
    } catch (My_Error my_error) {
        throw my_error;
    } catch (...) {
        My_Error my_error((string("ERROR WHILE WRITING GBAT! ABORTING."))
                          .c_str());
        throw my_error;
    }

  return 0;
}
