// Header file for the command line parsing argument class for the PARENT_GPU program suite
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

#ifndef ARG_PARSER_H
#define ARG_PARSER_H

#include <string>
#include <algorithm>
#include <cstring>

class Arg_Parser{
    public:
        Arg_Parser(int argc, char *argv[]);
        char *get(const std::string &option);
        bool exists(const std::string &option);
        char* get_ext(char* file_str);
        char **begin;
        char **end;
        int argc;
        char **argv;
};
#endif
