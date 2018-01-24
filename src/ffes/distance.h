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

#pragma once

#include <string>
#include <vector>

#include <GL/gl3w.h>

namespace ffes {

    /* If name_vxl_X != "no_basename", open vg_X ofstream at name_vxl_X.raw
     * and write metadata_file name_vxl_X.mhd */
    void prepare_voxel_grid_export(
            long voxel_grid[3],
            double vxl_size,
            std::ofstream& vg_A, 
            std::ofstream& vg_B, 
            std::ofstream& vg_D,
            const std::string& name_vxl_A,
            const std::string& name_vxl_B,
            const std::string& name_vxl_D
            );

    /* Generate framebuffer object and texture attached to COLOR_ATTACHMENT0 */
    void init_framebuffer_and_texture(
            int field_dim,
            int width,
            int height,
            GLuint& fbo_id,
            GLuint& texture_id);

    /* Compute the texture_diff(i,j) = texture_A(i,j) - texture_B(i,j) for each pixel */
    void compute_texture_difference(
            GLuint fbo,
            const long voxel_grid[3],
            GLuint prog_diff,
            GLuint texture_A,
            GLuint texture_B,
            GLuint texture_diff,
            GLuint dummy_vao);

    /* Pixel values in texture are reduced to one column of texture_contributions.
     * See the initialization of texture_contributions for more details. */
    void slice_reduction(
            const long voxel_grid[3],
            GLuint prog_reduction,
            GLuint texture,
            unsigned int slice_no,
            GLuint texture_contributions,
            GLuint dummy_vao);

}



