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
#include <iostream>
#include <vector>
#include <fstream>
#include <assert.h>
#include <string.h>

#include "utils.h"
#include <GL/gl3w.h>

void _check_gl_error(const char *file, int line) {
    GLenum err (glGetError());

    while(err!=GL_NO_ERROR) {
        std::string error;

        switch(err) {
            case GL_INVALID_OPERATION:      error="INVALID_OPERATION";      break;
            case GL_INVALID_ENUM:           error="INVALID_ENUM";           break;
            case GL_INVALID_VALUE:          error="INVALID_VALUE";          break;
            case GL_OUT_OF_MEMORY:          error="OUT_OF_MEMORY";          break;
            case GL_INVALID_FRAMEBUFFER_OPERATION:  error="INVALID_FRAMEBUFFER_OPERATION";  break;
        }

        std::cerr << "GL_" << error.c_str() <<" - "<<file<<":"<<line<<std::endl;
        exit(EXIT_FAILURE);
        err=glGetError();
    }
}

void error_callback(int error, const char* description) { 
    fprintf(stderr, "Error %d: %s\n", error, description);
}

void print_gpu_info(){
    int limit_tb = 0;
    int limit_ssbo = 0;
    glGetIntegerv(GL_MAX_TEXTURE_BUFFER_SIZE, &limit_tb);
    glGetIntegerv(GL_MAX_SHADER_STORAGE_BLOCK_SIZE, &limit_ssbo);
    std::cout << "Texture buffer size (MB): " << limit_tb / 1048576 << std::endl;
    std::cout << "SSBO size (MB): " << limit_ssbo / 1048576  << std::endl;

    GLint n, i;
    glGetIntegerv(GL_NUM_EXTENSIONS, &n);
    for (i = 0; i < n; i++) {
        printf("%s\n", glGetStringi(GL_EXTENSIONS, i));
    }
}

void field_dim_to(unsigned int field_dim, GLint& internal_format, GLenum& format){
    if (field_dim == 1){
        internal_format = GL_R32F;
        format = GL_RED;
    } else if (field_dim == 2){
        internal_format = GL_RG32F;
        format = GL_RG;
    } else if (field_dim == 3){
        internal_format = GL_RGB32F;
        format = GL_RGB;
    } else if (field_dim == 4){
        internal_format = GL_RGBA32F;
        format = GL_RGBA;
    } else {
        std::cerr << "Dimension not supported" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void texture_to_console(unsigned int texture_id, int width, int height, int field_dim){
    GLint internal_format;
    GLenum format;
    field_dim_to(field_dim, internal_format, format);
    std::vector<float> fA(field_dim*width*height, -1);
    glBindTexture(GL_TEXTURE_RECTANGLE, texture_id);
    glGetTexImage(GL_TEXTURE_RECTANGLE, 0, format, GL_FLOAT, &fA[0]);
    for(int j = 0; j < 100; ++j){ std::cout << "-"; }; std::cout << '\n';
    for(int i = 0; i < height; ++i){
        std::cout << i << ": ";
        for(int j = 0; j < width; ++j){
            printf("(");
            for(int d = 0; d < field_dim; ++d){
                printf("%.2f ", fA[field_dim*(width*i+j)+d]);
            }
            printf(") ");
        }
        std::cout << '\n';
    }
    for(int j = 0; j < 100; ++j){ std::cout << "-"; }; std::cout << '\n';
}

void write_json_file(const std::string filename, std::vector<std::string>& infos){
    assert(infos.size() % 2 == 0);
    std::ofstream file;
    file.open(filename);
    file << "{" << '\n';
    for(size_t i = 0; i < infos.size() / 2; ++i){
        std::string comma = "";
        if(2*i+2 != infos.size()) comma = ",";
        file << " \"" << infos[2*i] << "\": " << infos[2*i+1] << comma << '\n';
    }
    file << "}" << '\n';
    file.close();
}


/* From Geogram GFX / GLUP */
static int htoi(char digit) {
    if(digit >= '0' && digit <= '9') {
        return digit - '0';
    }
    if(digit >= 'a' && digit <= 'f') {
        return digit - 'a' + 10;
    }
    if(digit >= 'A' && digit <= 'F') {
        return digit - 'A' + 10;
    }
    fprintf(stderr, "xpm: unknown digit\n");
    return 0;
}

/* The colormap. */
static unsigned char i2r[1024];
static unsigned char i2g[1024];
static unsigned char i2b[1024];
static unsigned char i2a[1024];

/* Converts a two-digit XPM color code into
 *  a color index.
 */
static int char_to_index[256][256];

void glTexImage2DXPM(const char** xpm_data) {
    int width, height, nb_colors, chars_per_pixel;
    int line = 0;
    int color = 0;
    int key1 = 0, key2 = 0;
    char* colorcode;
    int x, y;
    unsigned char* rgba;
    unsigned char* pixel;
    sscanf(
        xpm_data[line], "%6d%6d%6d%6d",
        &width, &height, &nb_colors, &chars_per_pixel
    );
    line++;
    if(nb_colors > 1024) {
        fprintf(stderr, "xpm with more than 1024 colors\n");
        return;
    }
    if(chars_per_pixel != 1 && chars_per_pixel != 2) {
        fprintf(stderr, "xpm with more than 2 chars per pixel\n");
        return;
    }
    for(color = 0; color < nb_colors; color++) {
        int r, g, b;
        int none ;
        
        key1 = xpm_data[line][0];
        key2 = (chars_per_pixel == 2) ? xpm_data[line][1] : 0;
        colorcode = strstr(const_cast<char*>(xpm_data[line]), "c #");
        none = 0;
        if(colorcode == NULL) {
            colorcode = (char*) "c #000000";
            if(strstr(xpm_data[line], "None") != NULL) {
                none = 1;
            } else {
                fprintf(
                    stderr, "unknown xpm color entry (replaced with black)\n"
                );
            }
        }
        colorcode += 3;

        r = 16 * htoi(colorcode[0]) + htoi(colorcode[1]);
        g = 16 * htoi(colorcode[2]) + htoi(colorcode[3]);
        b = 16 * htoi(colorcode[4]) + htoi(colorcode[5]);

        i2r[color] = (unsigned char) r;
        i2g[color] = (unsigned char) g;
        i2b[color] = (unsigned char) b;
        i2a[color] = none ? 0 : 255;
        char_to_index[key1][key2] = color;
        line++;
    }
    rgba = (unsigned char*) malloc((size_t) (width * height * 4));
    pixel = rgba;
    for(y = 0; y < height; y++) {
        for(x = 0; x < width; x++) {
            if(chars_per_pixel == 2) {
                key1 = xpm_data[line][2 * x];
                key2 = xpm_data[line][2 * x + 1];
            } else {
                key1 = xpm_data[line][x];
                key2 = 0;
            }
            color = char_to_index[key1][key2];
            pixel[0] = i2r[color];
            pixel[1] = i2g[color];
            pixel[2] = i2b[color];
            pixel[3] = i2a[color];
            pixel += 4;
        }
        line++;
    }
    
    glTexImage2D(
        GL_TEXTURE_2D, 0,
        GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgba
    );
#ifndef __EMSCRIPTEN__    
    glGenerateMipmap(GL_TEXTURE_2D);
#endif    
    free(rgba);
}
