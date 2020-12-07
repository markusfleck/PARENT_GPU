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
#include "My_Error.cpp"


using namespace std;


Arg_Parser::Arg_Parser(int argc, char* argv[]){
    for(int i = 0; i < argc; i++) args.push_back(argv[i]);
}

bool Arg_Parser::compare(char const * cstr1, char const * cstr2){

    int counter = 0;
    if(cstr1[counter] != cstr2[counter]) return false;

    while( (cstr1[counter] != '\0') ){
        if(cstr1[counter + 1] != cstr2[counter + 1]) return false;
        counter++;
    }

    return true;

}

bool Arg_Parser::exists(char const * option){
    
    unsigned int counter = 0;
    for(unsigned int i = 0; i < args.size(); i++){
        if( compare( option, args[i] ) ) counter++;
    }
    
    if(counter > 1){
        My_Error my_error( (string("ERROR: COMMAND LINE OPTION \"") +
                       string(option) + string("\" SPECIFIED MULTIPLE TIMES! ABORTING.") )
                          .c_str() );
        throw my_error;
    }

    return counter;
}

char const * Arg_Parser::get(char const * option){
    if(!exists(option)){
        My_Error my_error( (string("ERROR: COMMAND LINE OPTION \"") +
                       string(option) + string("\" NOT SPECIFIED! ABORTING.") )
                          .c_str() );
        throw my_error;
    }
    
    unsigned int counter = 0;
    for(unsigned int i = 0; i < args.size(); i++){
        if( compare( (char*)option, args[i]) ) break;
        counter++;
    }
    
    if(counter > args.size() - 2){
        My_Error my_error( (string("ERROR: NO VALUE SPECIFIED FOR COMMAND LINE OPTION \"") +
                       string(option) + string("\"! ABORTING.") )
                          .c_str() );
        throw my_error;
    }
    

return args[counter+1];
}

char const * Arg_Parser::get_ext(char const * option){
    char const * value = get(option);
    char const * ext = "\0";

    unsigned int counter = 0;
    unsigned int dot_counter = 0;
    while ( value[counter] != '\0'){
        if(value[counter] == '.') {
            ext = value + counter + 1;
            dot_counter++;
        }
        counter++;
    }
    
    if( (dot_counter == 0) || ext[0] == '\0' ){
        My_Error my_error( (string("ERROR: FILE SPECIFIED FOR COMMAND LINE OPTION \"") +
                string(option) + string("\" (\"") + string(value) + string("\") DOES NOT HAVE AN EXTENSION! ABORTING.") )
                .c_str() );
        throw my_error;
    }
    
    return ext;
}

bool Arg_Parser::check_ext(char const * option, char const * target_ext){
    char const * ext = get_ext(option);
    return compare(ext, target_ext);
}






