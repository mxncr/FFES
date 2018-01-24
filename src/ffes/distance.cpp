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
#include "distance.h"

#include <iomanip>
#include <numeric>
#include <algorithm>

#include <GL/gl3w.h>

#include "ffes/field_sampling.h"
#include "ffes/utils.h"

#include <geogram_basic/basic/logger.h>
#include <geogram_basic/basic/file_system.h>
#include <geogram_basic/basic/process.h>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/transform.hpp>

using GEO::Logger;
using std::endl;

namespace ffes {

    static void field_dim_to(index_t field_dim, GLint& internal_format, GLenum& format){
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
            Logger::err("Textures init") << "field_dim=" << field_dim 
                << " not supported, please use 1, 2, 3 or 4" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    void init_framebuffer_and_texture(
            int field_dim,
            int width,
            int height,
            GLuint& fbo_id,
            GLuint& texture_id){
        glGenFramebuffers(1, &fbo_id);
        glGenTextures(1, &texture_id);

        /* Get texture format adapted to the dimension of the field */
        GLint internal_format;
        GLenum format;
        field_dim_to(field_dim, internal_format, format);

        /* Init the texture properties and attach them to the framebuffer */
        glBindFramebuffer(GL_FRAMEBUFFER, fbo_id);

        glBindTexture(GL_TEXTURE_RECTANGLE, texture_id);
        glTexImage2D(GL_TEXTURE_RECTANGLE, 0, internal_format, width, height, 0, format, GL_FLOAT, 0);
        glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST); 
        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, texture_id, 0); 
        glBindTexture(GL_TEXTURE_RECTANGLE, 0);
        check_gl_error();

        /* Default parameters for the framebuffer object */
        glViewport(0, 0, width, height);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        glDisable(GL_CULL_FACE); 
        glDisable(GL_DEPTH_TEST);

        GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            Logger::err("GPU init") << "Framebuffer creation failed\n" << std::endl;
            exit(EXIT_FAILURE);
        }

        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    /* Compute the texture_diff(i,j) = texture_A(i,j) - texture_B(i,j) for each pixel */
    void compute_texture_difference(
            GLuint fbo,
            const long voxel_grid[3],
            GLuint prog_diff,
            GLuint texture_A,
            GLuint texture_B,
            GLuint texture_diff,
            GLuint dummy_vao){

        glBindFramebuffer(GL_FRAMEBUFFER, fbo);

        /* Field difference */
        glBindTexture(GL_TEXTURE_RECTANGLE, texture_diff);

        glClear(GL_COLOR_BUFFER_BIT);

        glBindVertexArray(dummy_vao);

        glUseProgram(prog_diff);
        check_gl_error();

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_RECTANGLE, texture_A);
        GLint tf0 = glGetUniformLocation(prog_diff, "fieldA");
        glUniform1i(tf0, 0);
        check_gl_error();

        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_RECTANGLE, texture_B);
        GLint tf1 = glGetUniformLocation(prog_diff, "fieldB");
        glUniform1i(tf1, 1);
        check_gl_error();

        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        check_gl_error();

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_RECTANGLE, 0);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_RECTANGLE, 0);
    }

    /* Pixel values in texture are reduced to one column of texture_contributions.
     * See the initialization of texture_contributions for more details. */
    void slice_reduction(
            const long voxel_grid[3],
            GLuint prog_reduction,
            GLuint texture,
            index_t slice_no,
            GLuint texture_contributions,
            GLuint dummy_vao){
        geo_assert(slice_no >= 0 && slice_no < voxel_grid[2]);

        /* The reduced texture size is: (width, height) = (nz, ny) */
        glViewport(0, 0, voxel_grid[2], voxel_grid[1]);

        /* Texture reduction and contribution to distance computation (on the GPU) */
        glUseProgram(prog_reduction);
        GLint tfv = glGetUniformLocation(prog_reduction, "texture_values");
        glUniform1i(tfv, 0);

        GLint lw = glGetUniformLocation(prog_reduction, "width");
        glUniform1i(lw, voxel_grid[0]);

        /* slice_pos is used to set the line position in the vertex shader, to specify
         * which pixels will be written in the color attachment
         * slice_pos in ]-1, 1[ */
        float slice_pos = 2. * (double(slice_no) + 0.5) / double(voxel_grid[2]) - 1.;
        geo_assert(slice_pos > -1. && slice_pos < 1.);


        GLint lsno = glGetUniformLocation(prog_reduction, "slice_no");
        GLint lnbs = glGetUniformLocation(prog_reduction, "nb_slices");
        glUniform1ui(lsno, (GLuint) slice_no);
        glUniform1ui(lnbs, (GLuint) voxel_grid[2]);

        check_gl_error();

        glDrawBuffer(GL_COLOR_ATTACHMENT0); /* target texture */
        glBindVertexArray(dummy_vao);
        if(slice_no == 0){
            glClearColor(0.,0.,0.,0.); 
            glClear(GL_COLOR_BUFFER_BIT);
        }
        glDisable(GL_CULL_FACE);
        glDisable(GL_DEPTH_TEST);
        check_gl_error();
        glActiveTexture(GL_TEXTURE0);
        check_gl_error();
        glBindTexture(GL_TEXTURE_RECTANGLE, texture); /* input texture */
        check_gl_error();
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glBindTexture(GL_TEXTURE_RECTANGLE, 0);
        check_gl_error();
    }

    void metadata_file(const std::string& filename, long dims[3], double vxl_size){
        std::ofstream out;
        out.open(filename + ".mhd");
        out << "ObjectType = Image" << '\n';
        out << "NDims = 3" << '\n';
        out << std::setprecision(6);
        out << "DimSize = " << dims[0] << " " << dims[1] << " " << dims[2] << '\n';
        out << "ElementSize = " << vxl_size << " " << vxl_size << " " << vxl_size << '\n';
        out << "ElementType = MET_FLOAT" << '\n';
        out << "ElementDataFile = " << filename + ".raw" << '\n';
        out.close();
    }

    void prepare_voxel_grid_export(
            long voxel_grid[3],
            double vxl_size,
            std::ofstream& vg_A, 
            std::ofstream& vg_B, 
            std::ofstream& vg_D,
            const std::string& name_vxl_A,
            const std::string& name_vxl_B,
            const std::string& name_vxl_D
            ){
        if (name_vxl_A != "no_basename"){
            vg_A.open(name_vxl_A + ".raw");
            metadata_file(name_vxl_A, voxel_grid, vxl_size);
        }
        if (name_vxl_B != "no_basename"){
            vg_B.open(name_vxl_B + ".raw");
            metadata_file(name_vxl_B, voxel_grid, vxl_size);
        }
        if (name_vxl_D != "no_basename"){
            vg_D.open(name_vxl_D + ".raw");
            metadata_file(name_vxl_D, voxel_grid, vxl_size);
        }
    }
}
