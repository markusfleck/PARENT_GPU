// The header file of the ASCII IO library for the PARENT program suite
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

#ifndef IO_TEXT_H
#define IO_TEXT_H

// #pragma pack(1) causes problems!

std::string strip_line(std::string line);
std::string delete_char(std::string line, char del);
std::string strip_blanks(std::string line);
int read_ndx_file(std::ifstream *ndxfile, std::vector<int> *group1,
                  std::vector<int> *group2, std::string name1,
                  std::string name2);

#endif
