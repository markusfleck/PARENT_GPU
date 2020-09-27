#include "Arg_Parser.h"

using namespace std;


Arg_Parser::Arg_Parser(int argc, char *argv[]){
    this->argc = argc;
    this->argv = argv;
    begin = argv;
    end = argv + argc;

}

char *Arg_Parser::get_cmd_option(const string &option) {
    char **itr = find(begin, end, option);
    if (itr != end && ++itr != end) return *itr;
    return 0;
}

bool Arg_Parser::cmd_option_exists(const string &option) {
  return find(begin, end, option) != end;
}

char* Arg_Parser::get_extension(char* file_str){
    char *ptr, *type;
    char delimiter[] = ".";
    string tmp_str = string(file_str);

    ptr = strtok((char*)tmp_str.c_str(), delimiter);
    while (ptr != NULL) {
        type = ptr;
        ptr = strtok(NULL, delimiter);
    }
    return type;
}






