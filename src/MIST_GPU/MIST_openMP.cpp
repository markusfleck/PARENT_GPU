#define PRECISION double

#include <iostream>
#include <sys/time.h>
#include <vector>

#include "../util/io/io.h"
#include "../util/types.h"
#include "../util/classes/Entropy_Matrix.h"
#include "../util/classes/My_Error.cpp"

using namespace std;

struct Dof {
    unsigned char type;
    unsigned int id;

};

struct Dof_Pair{
    Dof dof1;
    Dof dof2;
    PRECISION value;
};

struct Mutual {
    PRECISION value;
    unsigned int id;
};


int main(int argc, char *argv[]){
    // start the stopwatch for the execution time
    timeval tv_start, tv_start_calc, tv_end, tv_end_calc;
    gettimeofday(&tv_start, NULL);

    if (argc != 5) {
        cerr << "USAGE: " << argv[0] << " -f input.par -o output.par\n";
        exit(EXIT_FAILURE);
    }

    Arg_Parser arg_parser(argc, argv);
    if (!arg_parser.cmd_option_exists("-f") ||
        !arg_parser.cmd_option_exists("-o")) {
        // check for correct command line options
        cerr << "USAGE: " << argv[0] << " -f input.par -o output.par\n";
        exit(EXIT_FAILURE);
    }

    if (strcmp(arg_parser.get_extension(arg_parser.get_cmd_option("-f")),
             "par") ||
      strcmp(arg_parser.get_extension(arg_parser.get_cmd_option("-o")),
             "par")) {
    // check for the extensions of the input and output file
    cerr << "USAGE: " << argv[0] << " -f input.par -o output.par\n";
    exit(EXIT_FAILURE);
    }

    try {
        Entropy_Matrix ent_mat(arg_parser.get_cmd_option("-f"));
        unsigned int n_bonds = ent_mat.getNBonds();
        unsigned int n_angles = ent_mat.getNAngles();
        unsigned int n_dihedrals = ent_mat.getNDihedrals();
        unsigned int n_type[3] = {n_bonds, n_angles, n_dihedrals};
        unsigned int n_atoms = n_dihedrals + 3;
        unsigned int n_dofs = 3 * n_atoms - 6;
    
        vector<Dof> processed;
        vector<Dof> unprocessed;
        
        for (unsigned char type = 0; type < 3; type++){
            for(unsigned int id = 0; id < n_type[type]; id++){
                Dof tmp_dof = {type, id + 1};
                unprocessed.push_back(tmp_dof);
            }
        }
        
        processed.push_back(unprocessed[0]);
        unprocessed.erase(unprocessed.begin() + 0);
        
        unsigned int n_mutual_max = n_dofs *n_dofs / 4;
        if(n_dofs % 2 == 1) n_mutual_max += n_dofs / 2;
        Mutual* h_mut_mat = new Mutual[n_mutual_max];
        
        
        PRECISION max_value = 0.0;
        unsigned int max_index = 0;
        Dof_Pair* mi_mist = new Dof_Pair[n_dofs - 1];
        
        
        gettimeofday(&tv_start_calc, NULL);
        #pragma omp parallel
        {
            while(unprocessed.size() > 0){
                
                PRECISION max_value_thread = 0.0;
                unsigned int max_index_thread = 0;
            
                #pragma omp for nowait
                for(unsigned int i = 0; i < processed.size(); i++){
                    for(unsigned int j = 0; j < unprocessed.size(); j++){
                        unsigned int idx = i * unprocessed.size() + j;
                        Mutual tmp_mut = {ent_mat.getMutual(processed[i].type, unprocessed[j].type, processed[i].id, unprocessed[j].id), idx};
                        h_mut_mat[idx] = tmp_mut;
                        
                        if(h_mut_mat[idx].value > max_value_thread){
                            max_value_thread = h_mut_mat[idx].value;
                            max_index_thread = h_mut_mat[idx].id;
                        }
                    }
                }

                #pragma omp critical 
                {
                    if(max_value_thread > max_value) {
                        max_value = max_value_thread;
                        max_index = max_index_thread;
                    }
                }
                #pragma omp barrier
                            
        
                #pragma omp single
                {
                    unsigned int idx_processed = max_index / unprocessed.size();
                    unsigned int idx_unprocessed = max_index % unprocessed.size();
                    
                    Dof_Pair tmp_dof_pair = {processed[idx_processed], unprocessed[idx_unprocessed], max_value};
                    mi_mist[processed.size() - 1] = tmp_dof_pair;
                    
                    processed.push_back(unprocessed[idx_unprocessed]);
                    unprocessed.erase(unprocessed.begin() + idx_unprocessed);
                    
                    max_value = 0.0;
                    max_index = 0;
                }
            }
        }
        gettimeofday(&tv_end_calc, NULL);
        
        
        for (unsigned char type1 = 0; type1 < 3; type1++){
            for(unsigned int idx1 = 0; idx1 < n_type[type1]; idx1++){
                for (unsigned char type2 = type1; type2 < 3; type2++){
                    unsigned int start_idx2 = (type1 == type2) ? idx1 + 1 : 0;     
                    for(unsigned int idx2 = start_idx2; idx2 < n_type[type2]; idx2++){
                        ent_mat.setMutual(type1, type2, idx1 + 1, idx2 + 1, 0.0);
                    }
                }
            }
        }
        
        for (unsigned int i = 0 ; i <  n_dofs - 1; i++){
            ent_mat.setMutual(mi_mist[i].dof1.type, mi_mist[i].dof2.type, mi_mist[i].dof1.id, mi_mist[i].dof2.id, mi_mist[i].value);
        }

        ent_mat.write(arg_parser.get_cmd_option("-o"));
        
        gettimeofday(&tv_end, NULL);
        cout << endl << endl;
        cout << "Calculation time: "<< tv_end_calc.tv_sec + 1e-6 * tv_end_calc.tv_usec - tv_start_calc.tv_sec - 1e-6 * tv_start_calc.tv_usec << endl;
        cout << "Total execution time: " << tv_end.tv_sec + 1e-6 * tv_end.tv_usec - tv_start.tv_sec - 1e-6 * tv_start.tv_usec << endl;
        cout << "PROGRAM FINISHED SUCCESSFULLY." << endl << endl << endl;
        
        
    } catch (My_Error my_error) {
        cerr << my_error.what() << endl;
        cerr << "USAGE:\n" << argv[0] << " -p input.par [--short]" << endl;
        return 1;
    } catch (...) {
        cerr << "AN UNIDENTIFIED ERROR HAS OCCURRED! ABORTING.\n" << endl;
        cerr << "USAGE:\n" << argv[0] << " -p input.par [--short]" << endl;
        return 1;
    }
    
 
    return 0;
}
