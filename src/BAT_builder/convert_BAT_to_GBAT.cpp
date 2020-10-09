#include <iostream>
#include "../util/classes/Bat.h"
#include "../util/io/io.h"
#include "../util/classes/My_Error.cpp"

using namespace std;

int main(int argc, char *argv[]) {
    try{

      if (argc != 5) {
        cerr << "USAGE: " << argv[0] << " -f input.bat -o ouput.gbat\n";
        exit(EXIT_FAILURE);
      }

      Arg_Parser arg_parser(argc, argv);
      if (!arg_parser.cmd_option_exists("-f") ||
          !arg_parser.cmd_option_exists("-o")) {
        // check for correct command line options
        cerr << "USAGE: " << argv[0] << " -f input.bat -o ouput.gbat\n";
        exit(EXIT_FAILURE);
      }

      if (strcmp(arg_parser.get_extension(arg_parser.get_cmd_option("-f")),
                 "bat") ||
          strcmp(arg_parser.get_extension(arg_parser.get_cmd_option("-o")),
                 "gbat")) {
        // check for the extensions of the input and output file
        cerr << "USAGE: " << argv[0] << " -f input.bat -o ouput.gbat\n";
        exit(EXIT_FAILURE);
      }
      
      Bat bat(arg_parser.get_cmd_option("-f"));
      bat.write_GBAT(arg_parser.get_cmd_option("-o"));

      cout << "GBAT WRITTEN SUCCESSFULLY." << endl << endl << endl;
    } catch (My_Error my_error) {
        throw my_error;
    } catch (...) {
        My_Error my_error((string("ERROR WHILE WRITING GBAT! ABORTING."))
                          .c_str());
        throw my_error;
    }

  return 0;
}
