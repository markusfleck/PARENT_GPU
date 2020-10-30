// Command line parsing argument class for the PARENT_GPU program suite
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

#include "Arg_Parser.h"

using namespace std;


Arg_Parser::Arg_Parser(int argc, char *argv[]){
    this->argc = argc;
    this->argv = argv;
    begin = argv;
    end = argv + argc;

}

char *Arg_Parser::get(const string &option) {
    char **itr = find(begin, end, option);
    if (itr != end && ++itr != end) return *itr;
    return 0;
}

bool Arg_Parser::exists(const string &option) {
  return find(begin, end, option) != end;
}

char* Arg_Parser::get_ext(char* file_str){
    char *ptr, *type;
    char delimiter[] = ".";
    string tmp_str = string(file_str);

    ptr = strtok((char*)tmp_str.c_str(), delimiter);
    type = ptr;
    while (ptr != NULL) {
        type = ptr;
        ptr = strtok(NULL, delimiter);
    }
    return type;
}






