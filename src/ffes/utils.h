/*  
 *  Copyright (C) 2016 Maxence Reberol
 *
 *  This file is part of Fast Finite Element Sampling (FFES).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 *
 *  Contact: Maxence Reberol
 *           maxence.reberol@inria.fr
 *
 *           ALICE Team
 *           LORIA, INRIA Lorraine, 
 *           Campus Scientifique, BP 239
 *           54506 VANDOEUVRE LES NANCY CEDEX 
 *           FRANCE
 */

#include <string>
#include <iomanip>
#include <sstream>
#include <vector>

#define check_gl_error() _check_gl_error(__FILE__,__LINE__)
void _check_gl_error(const char *file, int line);

void error_callback(int error, const char* description);

void print_gpu_info();

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

void texture_to_console(unsigned int texture_id, int width, int height, int field_dim);

// infos[2*i+0] is the key, 
// infos[2*i+1] is the value
void write_json_file(const std::string filename, std::vector<std::string>& infos);


void glTexImage2DXPM(const char** xpm_data);
