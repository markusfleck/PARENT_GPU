// Custom errror class for the PARENT program suite
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

#ifndef MY_ERROR_CPP
#define MY_ERROR_CPP

#include <stdexcept>
#include <string>

class My_Error : public std::runtime_error {
public:
  My_Error(const std::string &msg = "") : std::runtime_error(msg) {}
};

#endif
