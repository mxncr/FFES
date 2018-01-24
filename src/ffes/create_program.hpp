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

/* 
 * Inspired from https://github.com/prideout/blog-source/blob/master/p60/CreateProgram.cpp
 * Webpage: http://prideout.net/blog/?p=11
 */

#include <iostream>
#include <string>
#include "third_party/glsw/glsw.h"

void replaceAll( std::string &s, const std::string &search, const std::string &replace ) {
    for( size_t pos = 0; ; pos += replace.length() ) {
        // Locate the substring to replace
        pos = s.find( search, pos );
        if( pos == std::string::npos ) break;
        // Replace by erasing and inserting
        s.erase( pos, search.length() );
        s.insert( pos, replace );
    }
}

GLuint CreateProgram(
        const std::string& shader_path,
        const char* vsKey, const char* gsKey, const char* fsKey,
        const std::string& directives,
        const std::string& mapping_function,
        const std::string& interpolation_function)
{
    glswInit();
    std::string spath = shader_path + "/";
    glswSetPath(spath.c_str(), ".glsl");
    glswAddDirectiveToken("", "#version 430 core");
    if(directives.size() != 0){
        glswAddDirectiveToken("", directives.c_str());
    }
    
    const char* vsSource = glswGetShader(vsKey);
    const char* gsSource = glswGetShader(gsKey);
    const char* fsSource = glswGetShader(fsKey);
    if(!vsSource){
        std::cerr << "Vertex shader not found: " << vsKey << std::endl;
        exit(EXIT_FAILURE);
    }
    if(!fsSource){
        std::cerr << "Fragment shader not found: " << fsKey << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string vsSource_with_mapping = std::string(vsSource);
    if(mapping_function.size() != 0){
        replaceAll(vsSource_with_mapping, "<mapping function>", mapping_function);
    } 

    // Interpolation function is not used when printing reference coordinates
    std::string fsSource_with_interpolation = std::string(fsSource);
    if(interpolation_function.size() != 0){
        replaceAll(fsSource_with_interpolation, "<interpolation function>", interpolation_function);
    } 
    
    GLint compileSuccess;
    GLchar compilerSpew[256];
    GLuint programHandle = glCreateProgram();

    GLuint vsHandle = glCreateShader(GL_VERTEX_SHADER);
    const char * vsp = vsSource_with_mapping.c_str();
    glShaderSource(vsHandle, 1, &vsp, 0);
    glCompileShader(vsHandle);
    { // check for errors
        GLint compiled;
        glGetShaderiv( vsHandle, GL_COMPILE_STATUS, &compiled );
        if ( !compiled ) {
            GLsizei len;
            glGetShaderiv( vsHandle, GL_INFO_LOG_LENGTH, &len );

            GLchar* log = new GLchar[len+1];
            glGetShaderInfoLog( vsHandle, len, &len, log );
            std::cerr << "Shader compilation failed: " << log << std::endl;
            std::cerr << "(for shader " << fsKey << ")" << std::endl;
            std::cerr << vsp << std::endl;
            delete [] log;
                exit(EXIT_FAILURE);
        }
    }
    glGetShaderiv(vsHandle, GL_COMPILE_STATUS, &compileSuccess);
    glGetShaderInfoLog(vsHandle, sizeof(compilerSpew), 0, compilerSpew);

    glAttachShader(programHandle, vsHandle);

    GLuint gsHandle;
    if (gsKey) {
        gsHandle = glCreateShader(GL_GEOMETRY_SHADER);
        glShaderSource(gsHandle, 1, &gsSource, 0);
        glCompileShader(gsHandle);
        { // check for errors
            GLint compiled;
            glGetShaderiv( gsHandle, GL_COMPILE_STATUS, &compiled );
            if ( !compiled ) {
                GLsizei len;
                glGetShaderiv( gsHandle, GL_INFO_LOG_LENGTH, &len );

                GLchar* log = new GLchar[len+1];
                glGetShaderInfoLog( gsHandle, len, &len, log );
                std::cerr << "Shader compilation failed: " << log << std::endl;
                std::cerr << "(for shader " << gsKey << ")" << std::endl;
                std::cerr << gsSource << std::endl;
                delete [] log;
                exit(EXIT_FAILURE);
            }
        }
        glGetShaderiv(gsHandle, GL_COMPILE_STATUS, &compileSuccess);
        glGetShaderInfoLog(gsHandle, sizeof(compilerSpew), 0, compilerSpew);
        glAttachShader(programHandle, gsHandle);
    }
    
    GLuint fsHandle;
    if (fsKey) {
        fsHandle = glCreateShader(GL_FRAGMENT_SHADER);
        const char* fsp = fsSource_with_interpolation.c_str();
        glShaderSource(fsHandle, 1, &fsp, 0);
        glCompileShader(fsHandle);
        { // check for errors
            GLint compiled;
            glGetShaderiv( fsHandle, GL_COMPILE_STATUS, &compiled );
            if ( !compiled ) {
                GLsizei len;
                glGetShaderiv( fsHandle, GL_INFO_LOG_LENGTH, &len );

                GLchar* log = new GLchar[len+1];
                glGetShaderInfoLog( fsHandle, len, &len, log );
                std::cerr << "Shader compilation failed: " << log << std::endl;
                std::cerr << "(for shader " << fsKey << ")" << std::endl;
                std::cerr << fsp << std::endl;
                delete [] log;
                exit(EXIT_FAILURE);
            }
        }
        glGetShaderiv(fsHandle, GL_COMPILE_STATUS, &compileSuccess);
        glGetShaderInfoLog(fsHandle, sizeof(compilerSpew), 0, compilerSpew);
        glAttachShader(programHandle, fsHandle);
    }

    glLinkProgram(programHandle);
    GLint linkSuccess;
    glGetProgramiv(programHandle, GL_LINK_STATUS, &linkSuccess);
    glGetProgramInfoLog(programHandle, sizeof(compilerSpew), 0, compilerSpew);

    if (!linkSuccess) {
        GLsizei len;
        glGetProgramiv( programHandle, GL_INFO_LOG_LENGTH, &len );

        GLchar* log = new GLchar[len+1];
        glGetProgramInfoLog( programHandle, len, &len, log );
        std::cerr << "Shader linking failed: " << log << std::endl;
        delete [] log;
    }
    
    glswShutdown();
    return programHandle;
}
