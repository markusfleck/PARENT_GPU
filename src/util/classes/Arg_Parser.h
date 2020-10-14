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
