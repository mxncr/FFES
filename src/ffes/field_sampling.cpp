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

#include <iomanip>
#include <numeric>
#include <algorithm>

#include <GL/gl3w.h>

#include "ffes/field_sampling.h"
#include "ffes/utils.h"
#include "ffes/decompositions.h"
#include "ffes/create_program.hpp"

#include "third_party/bitmap_image.hpp"

#include <geogram_basic/basic/logger.h>
#include <geogram_basic/basic/file_system.h>
#include <geogram_basic/basic/process.h>

#include <geogram_basic/colormaps/french.xpm>
#include <geogram_basic/colormaps/black_white.xpm>
#include <geogram_basic/colormaps/viridis.xpm>
#include <geogram_basic/colormaps/rainbow.xpm>
#include <geogram_basic/colormaps/cei_60757.xpm>
#include <geogram_basic/colormaps/inferno.xpm>
#include <geogram_basic/colormaps/magma.xpm>
#include <geogram_basic/colormaps/parula.xpm>
#include <geogram_basic/colormaps/plasma.xpm>
#include <geogram_basic/colormaps/blue_red.xpm>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/transform.hpp>

using GEO::Logger;
using std::endl;

namespace ffes {


    /* Modified from Geogram_GFX/GLUP */
GLuint init_colormap(const char** xpm_data) {
    unsigned int txt_id;
    glGenTextures(1, &txt_id);
    glBindTexture(GL_TEXTURE_2D, txt_id);
    glTexImage2DXPM(xpm_data);
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
    return txt_id;
}

void generate_colormaps(std::vector<GLuint>& texture_ids, std::vector<std::string>& names){
    names.resize(0); texture_ids.resize(0);
    names.push_back("french"     ); texture_ids.push_back(init_colormap(french_xpm));
    names.push_back("black_white"); texture_ids.push_back(init_colormap(black_white_xpm));
    names.push_back("viridis"    ); texture_ids.push_back(init_colormap(viridis_xpm));
    names.push_back("rainbow"    ); texture_ids.push_back(init_colormap(rainbow_xpm));
    names.push_back("cei_60757"  ); texture_ids.push_back(init_colormap(cei_60757_xpm));
    names.push_back("inferno"    ); texture_ids.push_back(init_colormap(inferno_xpm));
    names.push_back("magma"      ); texture_ids.push_back(init_colormap(magma_xpm));
    names.push_back("parula"     ); texture_ids.push_back(init_colormap(parula_xpm));
    names.push_back("plasma"     ); texture_ids.push_back(init_colormap(plasma_xpm));
    names.push_back("blue_red"   ); texture_ids.push_back(init_colormap(blue_red_xpm));
    // TODO delete the textures at some point
}

// obsolete with the new colormaps
std::vector<float> generate_colormap(const rgb_t* cm){
    std::vector<float> data(3*1000);
    for(int i = 0; i < 1000; ++i){
        data[3*i+0] = float(cm[i].red  ) / 255.f;
        data[3*i+1] = float(cm[i].green) / 255.f;
        data[3*i+2] = float(cm[i].blue ) / 255.f;
    }
    return data;
}

ElementDecomposition::~ElementDecomposition(){
    /* delete GPU buffers if uploaded on it */
    if(glIsBuffer(vbo_id)) glDeleteBuffers(1, &vbo_id);
    if(glIsBuffer(ibo_id)) glDeleteBuffers(1, &ibo_id);
}


void ElementDecomposition::upload_to_gpu(GLuint vao_id){
    if(vertices.size() == 0 || tet_indices.size() == 0){
        Logger::err("ElementDecomposition") << "object must be filled with data"
            << " before upload to GPU, exit" << endl;
        exit(EXIT_FAILURE);
    }

    glBindVertexArray(vao_id);
    check_gl_error();

    /* Upload vertex coordinates */
    glGenBuffers(1, &vbo_id);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_id);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), 
            &vertices[0], GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void *)(0));
    check_gl_error();

    /* Upload tetrahedron indices */
    glGenBuffers(1, &ibo_id);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_id);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, tet_indices.size() * sizeof(unsigned int), 
            &tet_indices[0], GL_STATIC_DRAW);
    check_gl_error();
}

float zmin_of(int c, vector<float>& control_pts, int nb_cp_per_cell){
    float z_min = FLT_MAX; 
    for(int lv = 0; lv < nb_cp_per_cell; ++lv){
        z_min = GEO::geo_min(z_min, control_pts[3*(nb_cp_per_cell*c+lv)+2]);
    }
    return z_min;
}

void FieldObject::build_cube_at_zero( 
        const double bbox_min[3],
        const double bbox_max[3],
        index_t field_dim){
    double v[8][3] = {
        {bbox_min[0], bbox_min[1], bbox_min[2]}, 
        {bbox_max[0], bbox_min[1], bbox_min[2]}, 
        {bbox_min[0], bbox_max[1], bbox_min[2]}, 
        {bbox_max[0], bbox_max[1], bbox_min[2]}, 
        {bbox_min[0], bbox_min[1], bbox_max[2]}, 
        {bbox_max[0], bbox_min[1], bbox_max[2]}, 
        {bbox_min[0], bbox_max[1], bbox_max[2]}, 
        {bbox_max[0], bbox_max[1], bbox_max[2]}
    };
	const int hex_geo2mfem[] = {0,1,3,2,4,5,7,6};

    elt_group.resize(1);
    ElementGroup& eg = elt_group.back();

    eg.nb_mesh_cp_per_cell = 8;
    eg.nb_field_cp_per_cell = 8;
    eg.type = HEX;
    eg.mesh_dim = 3;
    eg.field_dim = field_dim;
    eg.mesh_ctrl_points.resize(3*8);
    for(int i = 0; i < 8; ++i){
        for(int d = 0; d < 3; ++d){
            eg.mesh_ctrl_points[3*i+d] = v[hex_geo2mfem[i]][d];
        }
    }
    eg.field_ctrl_points.resize(field_dim*8, 0); /* field at 0 */
    eg.mapping = "vec3 element_mapping(const vec3 ref_pos, const vec3 values[8]){return vec3(0,0,0);}";
    if (field_dim == 1) {
        eg.interpolation = "float element_interpolation(const vec3 ref_pos, const float values[8]){return 0;}";
    } else if (field_dim == 2) {
        eg.interpolation = "vec2 element_interpolation(const vec3 ref_pos, const vec2 values[8]){return vec2(0,0);}";
    } else if (field_dim == 3) {
        eg.interpolation = "vec3 element_interpolation(const vec3 ref_pos, const vec3 values[8]){return vec3(0,0,0);}";
    } else if (field_dim == 4) {
        eg.interpolation = "vec4 element_interpolation(const vec3 ref_pos, const vec4 values[8]){return vec3(0,0,0,0);}";
    }
}

void FieldObject::update_element_decomposition_for_each_element_type(index_t subdivision_level){
    vector<index_t> subdivs(elt_group.size(), subdivision_level);
    update_element_decomposition_for_each_element_type(subdivs);
}

void FieldObject::update_element_decomposition_for_each_element_type(const vector<index_t>& subdivision_levels){
    geo_assert(subdivision_levels.size() == elt_group.size());
    ref_elt.resize(elt_group.size());
    for(index_t g = 0; g < elt_group.size(); ++g){
        ElementDecomposition& ed = ref_elt[g];
        ed.subdivision_level = subdivision_levels[g];
        ed.primitive = elt_group[g].type;
        if(elt_group[g].type == TET){
            switch(subdivision_levels[g]){
                case 0:
                    ed.vertices    = tetra_r0_vertices;
                    ed.tet_indices = tetra_r0_tet_indices;
                    break;
                case 1:
                    ed.vertices    = tetra_r1_vertices;
                    ed.tet_indices = tetra_r1_tet_indices;
                    break;
                case 2:
                    ed.vertices    = tetra_r2_vertices;
                    ed.tet_indices = tetra_r2_tet_indices;
                    break;
                case 3:
                    ed.vertices    = tetra_r3_vertices;
                    ed.tet_indices = tetra_r3_tet_indices;
                    break;
                default:
                    Logger::err("EltDecomposition") << "subdivision level " <<
                        subdivision_levels[g] << " not supported" << endl;
                    exit(EXIT_FAILURE);
            }

        } else if(elt_group[g].type == HEX){
            switch(subdivision_levels[g]){
                case 0:
                    ed.vertices    = hex_r0_vertices;
                    ed.tet_indices = hex_r0_tet_indices;
                    break;
                case 1:
                    ed.vertices    = hex_r1_vertices;
                    ed.tet_indices = hex_r1_tet_indices;
                    break;
                case 2:
                    ed.vertices    = hex_r2_vertices;
                    ed.tet_indices = hex_r2_tet_indices;
                    break;
                case 3:
                    ed.vertices    = hex_r3_vertices;
                    ed.tet_indices = hex_r3_tet_indices;
                    break;
                default:
                    Logger::err("EltDecomposition") << "subdivision level " <<
                        subdivision_levels[g] << " not supported" << endl;
                    exit(EXIT_FAILURE);
            }
        }
    }
}

const std::map<int, std::string> ndofs2basis = {
    {   8, "Q1"}, // In 3D: dim(Qk) = (k+1)^3
    {  27, "Q2"},
    {  64, "Q3"},
    { 125, "Q4"},
    { 216, "Q5"},
    { 343, "Q6"},
    { 512, "Q7"},
    { 729, "Q8"},
    {1000, "Q9"},
    {  4, "P1"}, // In 3D: dim(Pk) = (k+3)(k+2)(k+1)/6
    { 10, "P2"},
    { 20, "P3"},
    { 35, "P4"},
    { 56, "P5"},
    { 84, "P6"},
    {120, "P7"},
    {165, "P8"},
    {220, "P9"}
};

void replace_old_by_new(std::string& str,
        const std::string& oldStr,
        const std::string& newStr) {
    std::string::size_type pos = 0u;
    while((pos = str.find(oldStr, pos)) != std::string::npos){
        str.replace(pos, oldStr.length(), newStr);
        pos += newStr.length();
    }
}

std::string float_to_double(const std::string& shader_float) {
    std::string shader_double = shader_float;
    replace_old_by_new(shader_double, "float u", "double u");
    replace_old_by_new(shader_double, "float v", "double v");
    replace_old_by_new(shader_double, "float w", "double w");
    replace_old_by_new(shader_double, "float z", "double z");
    replace_old_by_new(shader_double, "float result", "double result");
    replace_old_by_new(shader_double, "vec2 result = vec2", "dvec2 result = dvec2");
    replace_old_by_new(shader_double, "vec3 result = vec3", "dvec3 result = dvec3");
    replace_old_by_new(shader_double, "vec4 result = vec4", "dvec4 result = dvec4");
    replace_old_by_new(shader_double, "float r = 0.", "double r = 0.");
    return shader_double;
}

void FieldObject::update_shader_program_for_each_element_type(DebugShader& debug,
        bool computation_with_double){
    std::string path = GEO::FileSystem::dir_name(GEO::Process::executable_filename()) + "/shaders";
    for(index_t g = 0; g < elt_group.size(); ++g){
        ElementGroup& eg = elt_group[g];
        std::string directives = "#define NB_CP " + std::to_string(eg.nb_mesh_cp_per_cell) + "\n"
            + "#define NB_CV " + std::to_string(eg.nb_field_cp_per_cell) + "\n"
            + "#define FIELD_DIM " + std::to_string(eg.field_dim) + "\n";
        if(debug.debug){
            directives += "#define DEBUG_SHADER\n";
            if (debug.show_elements) directives += "#define DEBUG_SHOW_ELEMENTS\n";
            if (debug.show_faces) directives += "#define DEBUG_SHOW_FACES\n";
            if (debug.show_group){
                directives += "#define DEBUG_SHOW_GROUP\n";
                directives += "#define GROUP " + to_string_with_precision(g / (elt_group.size()-0.99),3) + " \n";
            }
        }
        if (computation_with_double) {
            directives += "#define ct double\n";
            directives += "#define ct2 dvec2\n";
            directives += "#define ct3 dvec3\n";
            directives += "#define ct4 dvec4\n";
        } else {
            directives += "#define ct float\n";
            directives += "#define ct2 vec2\n";
            directives += "#define ct3 vec3\n";
            directives += "#define ct4 vec4\n";
        }

        std::string map_fct = "";
        if (eg.mapping.size() == 0) {
            Logger::out("Shader") << "Mapping function should be provided in the input file, exit" << endl; 
            exit(EXIT_FAILURE);
        } else {
            map_fct = eg.mapping;
        }

        bool USE_DOUBLE = false;
        std::string fct_fct = "";
        if (eg.interpolation.size() == 0) {
            Logger::out("Shader") << "Interpolation function should be provided in the input file, exit" << endl; 
            exit(EXIT_FAILURE);
        } else {
            fct_fct = eg.interpolation;
            if (USE_DOUBLE) {
                fct_fct = float_to_double(fct_fct);
            }
        }

        eg.program = CreateProgram(path, "template-ssbo.Vertex", 
                "template-ssbo.Geometry", 
                "template-ssbo.Fragment",
                directives,
                map_fct,
                fct_fct);
    }
}

#define BUFFER_OFFSET(offset) ((void *) (offset))

void FieldObject::finalise(){
    for(index_t g = 0; g < elt_group.size(); ++g){
        ElementGroup& eg = elt_group[g];

        /* Generate OpenGL Vertex Array Object */
        glGenVertexArrays(1, &eg.vao_id);
        glBindVertexArray(eg.vao_id);

        /* Upload the element decomposition data */
        ref_elt[g].upload_to_gpu(eg.vao_id);

        /* Upload mesh control points for the element mappings */
        glGenBuffers(1, &eg.ctrl_points_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, eg.ctrl_points_vbo);
        glBufferData(GL_ARRAY_BUFFER, eg.mesh_ctrl_points.size() * sizeof(float), 
                &eg.mesh_ctrl_points[0], GL_STATIC_DRAW);
        check_gl_error();

        /* Set the Vertex Attributes pointers and divisor so each drawn instance
         * can have access to the mapping control points */
        const int first_loc = 1; // 0 is reference vertices
        for(unsigned int i = 0; i < eg.nb_mesh_cp_per_cell; ++i){
            glVertexAttribPointer(first_loc+i, 
                    3, 
                    GL_FLOAT, GL_FALSE, 
                    sizeof(GLfloat)*3*eg.nb_mesh_cp_per_cell, 
                    BUFFER_OFFSET(sizeof(float)*eg.mesh_dim*i));
            glEnableVertexAttribArray(first_loc+i);
            glVertexAttribDivisor(first_loc+i, 1);  // change at each instance call
        }
        check_gl_error();

        /* Upload the field coeffient values */
        glGenBuffers(1, &eg.coefficients_ssbo);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, eg.coefficients_ssbo);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(GLfloat)*eg.field_ctrl_points.size(),
                &eg.field_ctrl_points[0], GL_STATIC_DRAW);
        GLsizeiptr ssbo_total = eg.field_ctrl_points.size() / (eg.nb_field_cp_per_cell * eg.field_dim)
            * eg.nb_field_cp_per_cell * eg.field_dim * sizeof(GLfloat);
        glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 0, eg.coefficients_ssbo, 0, ssbo_total);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        check_gl_error();

        glBindVertexArray(0);
    }
}

void FieldObject::print_data(){
    for(index_t g = 0; g < elt_group.size(); ++g){
        ElementGroup& eg = elt_group[g];
        std::cout << "# ElementGroup " << g << '\n';
        index_t nb_cells = eg.mesh_ctrl_points.size() / (eg.nb_mesh_cp_per_cell * eg.mesh_dim);
        for(index_t e = 0; e < nb_cells; ++e){
            std::cout << " elt " << e << ": \n";
            for(index_t lv = 0; lv < eg.nb_mesh_cp_per_cell; ++lv){
                std::cout << "   ";
                for(index_t d = 0; d < eg.mesh_dim; ++d){
                    std::cout << eg.mesh_ctrl_points[eg.mesh_dim*eg.nb_mesh_cp_per_cell*e
                        + eg.mesh_dim * lv + d] << " ";
                }
                std::cout << "\n";
            }
        }
    }
}

long FieldObject::nb_cells(){
    long nb_elts = 0;
    for(index_t g = 0; g < elt_group.size(); ++g){
        ElementGroup& eg = elt_group[g];
        nb_elts += eg.mesh_ctrl_points.size() / (eg.nb_mesh_cp_per_cell*eg.mesh_dim);
    }
    return nb_elts;
}

float zmin_of_cpoints(const std::vector<float>& cpoints,
        int prim_id,
        int nb_cp_per_primitive){
    float zmin = FLT_MAX;
    for(int i = 0; i < nb_cp_per_primitive; ++i){
        if(cpoints[prim_id*nb_cp_per_primitive*3+i*3+2] < zmin){
            zmin = cpoints[prim_id*nb_cp_per_primitive*3+i*3+2];
        }
    }
    return zmin;
}

float zmax_of_cpoints(const std::vector<float>& cpoints,
        int prim_id,
        int nb_cp_per_primitive){
    float zmax = -FLT_MAX;
    for(int i = 0; i < nb_cp_per_primitive; ++i){
        if(cpoints[prim_id*nb_cp_per_primitive*3+i*3+2] > zmax){
            zmax = cpoints[prim_id*nb_cp_per_primitive*3+i*3+2];
        }
    }
    return zmax;
}

void FieldObject::update_element_rendering_range(double z_slice){
    for(index_t g = 0; g < elt_group.size(); ++g){
        ElementGroup& eg = elt_group[g];
        if(eg.prev_zslice > z_slice){
            eg.rdr_start = 0;
            eg.rdr_end   = 0;
        }

        long nb_elt = eg.mesh_ctrl_points.size() / (3 * eg.nb_mesh_cp_per_cell);
        double eps = 0.;

        while(zmin_of_cpoints(eg.mesh_ctrl_points, eg.rdr_end, eg.nb_mesh_cp_per_cell) 
                <= z_slice + eps && (int) eg.rdr_end < (int) nb_elt - 1){
            eg.rdr_end += 1;
        }
        while(zmax_of_cpoints(eg.mesh_ctrl_points, eg.rdr_start, eg.nb_mesh_cp_per_cell) 
                <= z_slice - eps && (int) eg.rdr_start < (int) eg.rdr_end - 1){
            eg.rdr_start += 1;
        }
        geo_debug_assert(eg.rdr_start <= eg.rdr_end);
    }
}

void FieldObject::update_bounding_box(double bbox_min[3], double bbox_max[3]){
    /* This approach assumes the mesh bbox is inside the control point bbox */
    for(index_t g = 0; g < elt_group.size(); ++g){
        const ElementGroup& eg = elt_group[g];
        geo_assert(eg.mesh_dim == 3);
        for(index_t v = 0; v < eg.mesh_ctrl_points.size() / 3; ++v){
            bbox_min[0] = GEO::geo_min(bbox_min[0], (double) eg.mesh_ctrl_points[3*v+0]);
            bbox_max[0] = GEO::geo_max(bbox_max[0], (double) eg.mesh_ctrl_points[3*v+0]);
            bbox_min[1] = GEO::geo_min(bbox_min[1], (double) eg.mesh_ctrl_points[3*v+1]);
            bbox_max[1] = GEO::geo_max(bbox_max[1], (double) eg.mesh_ctrl_points[3*v+1]);
            bbox_min[2] = GEO::geo_min(bbox_min[2], (double) eg.mesh_ctrl_points[3*v+2]);
            bbox_max[2] = GEO::geo_max(bbox_max[2], (double) eg.mesh_ctrl_points[3*v+2]);
        }
    }
}

void FieldObject::update_value_ranges(
        double values_min[/*field_dim*/], 
        double values_max[/*field_dim*/]){
    const int vdim = field_dimension();
    for(index_t g = 0; g < elt_group.size(); ++g){
        const ElementGroup& eg = elt_group[g];
        for(index_t v = 0; v < eg.field_ctrl_points.size() / vdim; ++v){
            double magnitude = 0.;
            for(int d = 0; d < vdim; ++d){
                values_min[d] = GEO::geo_min(values_min[d], (double)
                        eg.field_ctrl_points[vdim*v+d]); 
                values_max[d] = GEO::geo_max(values_max[d], (double)
                        eg.field_ctrl_points[vdim*v+d]);
                magnitude += eg.field_ctrl_points[vdim*v+d] * eg.field_ctrl_points[vdim*v+d];
            }
            magnitude = std::pow(magnitude, 0.5);
            values_min[vdim] = GEO::geo_min(values_min[vdim], magnitude); 
            values_max[vdim] = GEO::geo_max(values_max[vdim], magnitude);
        }
    }
}

void FieldObject::draw_slice(
        GLuint fbo,
        int width,
        int height,
        double z_slice,
        double MVP[16],
        bool use_limited_elt_range){

    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glClearColor(ffes::cc, ffes::cc, ffes::cc, ffes::cc); 
    glClear(GL_COLOR_BUFFER_BIT);
    glViewport(0, 0, width, height);

    for(size_t g = 0; g < elt_group.size(); ++g){
        ElementGroup& eg = elt_group[g];
        ElementDecomposition& ed = ref_elt[g];

        /* Update [rdr_start, rdr_end] */
        if(use_limited_elt_range){
            update_element_rendering_range(z_slice);
        } else {
            eg.rdr_start = 0;
            eg.rdr_end   = eg.mesh_ctrl_points.size() / (3 * eg.nb_mesh_cp_per_cell) - 1;
        }

        /* Use the VAO and shader associated to the group */
        glBindVertexArray(eg.vao_id);
        glUseProgram(eg.program); 
        check_gl_error();

        /* Unecessary ? */
        glDisable(GL_CULL_FACE); 
        glDisable(GL_DEPTH_TEST);

        /* Tetrahedron indices of the reference element mesh */
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ed.ibo_id); // out of the loop ?
        check_gl_error();

        /* Set the projection matrix */
        GLuint matrix = glGetUniformLocation(eg.program, "MVP");
        glUniformMatrix4dv(matrix, 1, GL_FALSE, &MVP[0]);
        check_gl_error();

        /* Set the slicing plane */
        float z_slice_f = float(z_slice);
        GLint loc = glGetUniformLocation(eg.program, "z_slice");
        glUniform1f(loc, z_slice_f);
        check_gl_error();

        /* Bind the SSBO for field coefficient access in the fragment shader */
        GLuint nb_instances = eg.rdr_end - eg.rdr_start + 1;
        geo_debug_assert(nb_instances != 0);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, eg.coefficients_ssbo);
        // There is a error when using glBindBufferRange, I don't understand why
        // GLintptr ssbo_start = eg.rdr_start * eg.nb_field_cp_per_cell * eg.field_dim * sizeof(GLfloat);
        // GLsizeiptr ssbo_size = nb_instances * eg.nb_field_cp_per_cell * eg.field_dim * sizeof(GLfloat);
        // GLsizeiptr ssbo_total = eg.field_ctrl_points.size() / (eg.nb_field_cp_per_cell * eg.field_dim)
        //     * eg.nb_field_cp_per_cell * eg.field_dim * sizeof(GLfloat);
        // printf("start: %li, end: %li, total: %li \n", ssbo_start, ssbo_size, ssbo_total);
        // glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 0, eg.coefficients_ssbo,
        //         ssbo_start, ssbo_total);
        check_gl_error();

        /* Send the id of the first instance for coefficient access */
        GLint lfi = glGetUniformLocation(eg.program, "first_instance");
        glUniform1ui(lfi, eg.rdr_start);
        check_gl_error();

        /* Draw the elements of the group inside the rendering range */
        glDrawElementsInstancedBaseInstance(GL_LINES_ADJACENCY, 
                ed.tet_indices.size(),
                GL_UNSIGNED_INT, 
                (void*)0,
                nb_instances,
                eg.rdr_start);                  // first instance
        check_gl_error();

        /* unbind */
        glBindVertexArray(0);

        eg.prev_zslice = z_slice;
    }
}

glm::dmat4 compute_model_view_projection_matrix(const double bbox_min[3], 
        const double bbox_max[3]){
    /* camera is at the center of the scene, looking to z_max of the scene */
    glm::dvec3 center = glm::dvec3(
            0.5*(bbox_min[0] + bbox_max[0]),
            0.5*(bbox_min[1] + bbox_max[1]),
            0.5*(bbox_min[2] + bbox_max[2]));
    glm::dvec3 lookat = glm::dvec3(center.x, center.y, bbox_max[2]);
    glm::dmat4 view   = glm::lookAt(center, lookat, glm::dvec3(0., 1., 0.) );
    double lx = bbox_max[0] - bbox_min[0];
    double ly = bbox_max[1] - bbox_min[1];
    double lz = bbox_max[2] - bbox_min[2];
    /* warning: changing the m value affects the computed norm 
     * if m != 0.5, the h_pxl used in the norm computation is no longer in sync
     * with the slice sampling .*/
    double m = 0.50; 
    glm::dmat4 projection = glm::ortho(-m*lx, m*lx, -m*ly, m*ly, -m*lz, m*lz);
    return projection * view;
}

GLFWwindow* init_gui(bool visible, int& window_w, int& window_h){
    // Setup window
    glfwSetErrorCallback(error_callback);
    if (!glfwInit())
        return NULL;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_DEPTH_BITS, 32);
    glfwWindowHint(GLFW_STENCIL_BITS, GL_DONT_CARE);

    GLFWmonitor* monitor = glfwGetPrimaryMonitor();
    const GLFWvidmode* mode = glfwGetVideoMode(monitor);
    glfwWindowHint(GLFW_RED_BITS, mode->redBits);
    glfwWindowHint(GLFW_GREEN_BITS, mode->greenBits);
    glfwWindowHint(GLFW_BLUE_BITS, mode->blueBits);
    glfwWindowHint(GLFW_REFRESH_RATE, mode->refreshRate);
    window_w = mode->width-30;
    window_h = mode->height-30;
    if(!visible) glfwWindowHint(GLFW_VISIBLE, false);

    GLFWwindow* window = glfwCreateWindow(window_w, window_h, "ffes", NULL, NULL);

    glfwMakeContextCurrent(window);

    if (gl3wInit()) {
        Logger::err("OpenGL") << "failed to initialize OpenGL" << endl;
        exit(EXIT_FAILURE);
    }

    if (!gl3wIsSupported(4, 3)) {
        Logger::err("OpenGL")<< "OpenGL 4.3 not supported. This is required for SSBO support." << std::endl;
        exit(EXIT_FAILURE);
    }

    return window;
}

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

void init_reduction_texture(
        GLuint& fbo_id,
        const long voxel_grid[3],
        GLuint& texture_contributions){

    GLsizei rect_max;
    glGetIntegerv(GL_MAX_RECTANGLE_TEXTURE_SIZE, &rect_max);
    if (voxel_grid[1] > rect_max || voxel_grid[2] > rect_max) {
        Logger::err("GPU init") << "maximal size for rectangle texture: " << rect_max << std::endl;
        exit(EXIT_FAILURE);
    }

    GLint dims[2];
    glGetIntegerv(GL_MAX_VIEWPORT_DIMS, &dims[0]);
    if (voxel_grid[2] > dims[0] || voxel_grid[1] > dims[1]) {
        Logger::err("GPU init") << "maximal viewport size: " << dims[0] << "x" << dims[1];
        Logger::err("GPU init") << ", incompatible with reduction size: " << voxel_grid[1] << "x" << voxel_grid[2] << std::endl;
        exit(EXIT_FAILURE);
    }

    glGenFramebuffers(1, &fbo_id);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo_id);

    glGenTextures(1, &texture_contributions);
    glBindTexture(GL_TEXTURE_RECTANGLE, texture_contributions);
    glTexImage2D(GL_TEXTURE_RECTANGLE, 0, GL_RGBA32F, (GLsizei) voxel_grid[2], (GLsizei) voxel_grid[1], 0, GL_RGBA, GL_FLOAT, 0);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST); 
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, texture_contributions, 0);

    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glClearColor(0.,0.,0.,0.);
    glClear(GL_COLOR_BUFFER_BIT);

    glBindTexture(GL_TEXTURE_RECTANGLE, 0);
    check_gl_error();

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void bbox_to_voxel_grid(
        long nb_pixel_on_shortest_axis,
        const double bbox_min[3], 
        const double bbox_max[3], 
        long voxel_grid[3]
        ){
    double lx = bbox_max[0] - bbox_min[0];
    double ly = bbox_max[1] - bbox_min[1];
    double lz = bbox_max[2] - bbox_min[2];
    geo_assert(lx != 0. && ly != 0. && lz != 0 && nb_pixel_on_shortest_axis != 0);
    if(lx <= ly && lx <= lz){
        double h_pxl = lx / nb_pixel_on_shortest_axis;
        voxel_grid[0] = nb_pixel_on_shortest_axis;
        voxel_grid[1] = (int) (ly / h_pxl + 0.5); 
        voxel_grid[2] = (int) (lz / h_pxl + 0.5); 
    } else if (ly <= lz && ly <= lx){
        double h_pxl = ly / nb_pixel_on_shortest_axis;
        voxel_grid[0] = (int) (lx / h_pxl + 0.5); 
        voxel_grid[1] = nb_pixel_on_shortest_axis;
        voxel_grid[2] = (int) (lz / h_pxl + 0.5); 
    } else if (lz <= lx && lz <= ly){
        double h_pxl = lz / nb_pixel_on_shortest_axis;
        voxel_grid[0] = (int) (lx / h_pxl + 0.5); 
        voxel_grid[1] = (int) (ly / h_pxl + 0.5); 
        voxel_grid[2] = nb_pixel_on_shortest_axis;
    } 
}

GLuint create_difference_shader(index_t field_dim){
    std::string path = GEO::FileSystem::dir_name(GEO::Process::executable_filename()) + "/shaders";
    return CreateProgram(path, "diff.Vertex", 
            NULL, 
            "diff.Fragment",
            "#define FIELD_DIM " + std::to_string(field_dim) + "\n"
            "#define NO_DATA " + to_string_with_precision(cc, 9) + "\n", 
            "", "");
}

GLuint create_reduction_shader(index_t field_dim){
    std::string path = GEO::FileSystem::dir_name(GEO::Process::executable_filename()) + "/shaders";
    return CreateProgram(path, "reduction.Vertex", 
                NULL, 
                "reduction.Fragment",
                "#define FIELD_DIM " + std::to_string(field_dim) + "\n"
                "#define NO_DATA " + to_string_with_precision(cc, 9) + "\n", 
                "", "");
}

GLuint init_texture_colormap(){
    std::vector<float> cmap = generate_colormap(&jet_colormap[0]);
    GLuint cm_width = cmap.size() / 3;
    GLuint texture_colormap;
    glGenTextures(1, &texture_colormap);
    glBindTexture(GL_TEXTURE_1D, texture_colormap);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, cm_width, 0, GL_RGB, GL_FLOAT, &cmap[0]);
    glBindTexture(GL_TEXTURE_1D, 0);
    return texture_colormap;
}

GLuint create_display_texture_shader(index_t field_dim, DebugShader& debug){
    std::string path = GEO::FileSystem::dir_name(GEO::Process::executable_filename()) + "/shaders";

    std::string directives = 
            "#define FIELD_DIM " + std::to_string(field_dim) + "\n"
            "#define NO_DATA " + to_string_with_precision(cc, 9) + "\n";
    if(debug.debug){
        directives += "#define DEBUG_SHADER\n";
        if (debug.show_elements){ 
            directives += "#define DEBUG_SHOW_ELEMENTS\n";
            directives += "#define NB_ELTS " + std::to_string(debug.nb_elts) + "\n";
        }
        if (debug.show_group) directives += "#define DEBUG_SHOW_GROUP\n";
    }
    return CreateProgram(path, "texture_display.Vertex", NULL, 
            "texture_display.Fragment",
            directives,
            "", "");
}

const float quad_vertices[4*2] = {
  -0.3, -0.3, 
   0.3, -0.3, 
  -0.3,  0.3, 
   0.3,  0.3 };

void display_textures(
        GLFWwindow* window, 
        GLuint shader_program,
        GLuint texture_left,
        GLuint texture_center,
        GLuint texture_right,
        int texture_width,
        int texture_height,
        GLuint texture_colormap,
        int repeat,
        float range_min[3],
        float range_max[3],
        index_t field_dim,
        index_t field_coord
        ){
    geo_assert(field_coord <= field_dim);

    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.f, 0.f, 0.f, 1.f); /* black background */
    glClear(GL_COLOR_BUFFER_BIT);
    check_gl_error();

    glUseProgram(shader_program);

    /* Set uniform values in the shader */
    GLint lmin = glGetUniformLocation(shader_program, "val_min");
    GLint lmax = glGetUniformLocation(shader_program, "val_max");
    GLint lfc = glGetUniformLocation(shader_program, "field_coord");
    glUniform1i(lfc, field_coord);
    GLint lma = glGetUniformLocation(shader_program, "field_dim");
    glUniform1i(lma, field_dim);
    GLint lfw = glGetUniformLocation(shader_program, "texture_width");
    glUniform1i(lfw, texture_width);
    GLint lfh = glGetUniformLocation(shader_program, "texture_height");
    glUniform1i(lfh, texture_height);
    GLuint matrix = glGetUniformLocation(shader_program, "MVP");

    /* Basic colormap */
    GLint tcm = glGetUniformLocation(shader_program, "colormap");
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, texture_colormap);
    glUniform1i(tcm, 1);
    GLint lrep = glGetUniformLocation(shader_program, "repeat");
    glUniform1i(lrep, repeat);
    check_gl_error();

    float ratio = (float) display_w / (float) display_h;
    glm::mat4 proj = glm::ortho(-1.f, 1.f, float(-1.f)/ratio, float(1.f/ratio));
    int tmax = GEO::geo_max(texture_width, texture_height);
    glm::mat4 scale = glm::scale(glm::mat4(1.0f), 
            glm::vec3((float) texture_width/tmax, (float) texture_height/tmax, 1.0f));

    /* Left panel of the window */
    GLint lft = glGetUniformLocation(shader_program, "val");
    glUniform1i(lft, 0);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_RECTANGLE, texture_left);
    check_gl_error();

    /* Drawing to left */
    glm::mat4 MVPl = scale * glm::translate(glm::vec3(-0.66f,0.f,0.f)) * proj;
    glUniformMatrix4fv(matrix, 1, GL_FALSE, &MVPl[0][0]);
    glUniform1f(lmin, range_min[0]);
    glUniform1f(lmax, range_max[0]);
    check_gl_error();
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    check_gl_error();

    /* Right panel of the window */
    GLint lftr = glGetUniformLocation(shader_program, "val");
    glUniform1i(lftr, 0);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_RECTANGLE, texture_right);
    check_gl_error();

    /* Drawing to right */
    glm::mat4 MVPr = scale * glm::translate(glm::vec3(+0.66f,0.f,0.f)) * proj;
    glUniformMatrix4fv(matrix, 1, GL_FALSE, &MVPr[0][0]);
    glUniform1f(lmin, range_min[1]);
    glUniform1f(lmax, range_max[1]);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    check_gl_error();

    /* Center panel of the window */
    GLint lftc = glGetUniformLocation(shader_program, "val");
    glUniform1i(lftc, 0);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_RECTANGLE, texture_center);
    glUniform1f(lmin, range_min[2]);
    glUniform1f(lmax, range_max[2]);
    check_gl_error();

    /* Drawing to center */
    glm::mat4 MVPc = scale * glm::translate(glm::vec3(0.f,0.f,0.f)) * proj;
    check_gl_error();
    glUniformMatrix4fv(matrix, 1, GL_FALSE, &MVPc[0][0]);
    check_gl_error();
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    check_gl_error();

    glBindTexture(GL_TEXTURE_RECTANGLE, 0);
    glBindVertexArray(0);

}


void texture_to_vector(GLuint texture_id, int width, int height, int field_dim,
        std::vector<float>& values){
    GLint internal_format;
    GLenum format;
    field_dim_to(field_dim, internal_format, format);
    values.resize(field_dim*width*height, -1);
    glBindTexture(GL_TEXTURE_RECTANGLE, texture_id);
    glGetTexImage(GL_TEXTURE_RECTANGLE, 0, format, GL_FLOAT, &values[0]);
}

void texture_colormap_to_vector(GLuint texture_id, 
        std::vector<unsigned char>& values){

    glBindTexture(GL_TEXTURE_2D, texture_id);
    int width;
    int height;
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &height);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &width);
    values.resize(4*width*height, 0);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, &values[0]); 
    glBindTexture(GL_TEXTURE_2D, 0);
}

float clamp(float val, float min, float max, int repeat = 0){
    float cval = (val - min) / (max-min);
    if(cval < 0){
        cval =  0;
    } else if(cval > 1){
        cval = 1;
    } 
    if (repeat > 0) {
        float t = float(repeat)+1.;
        cval = t*cval - int(t*cval);
    }
    return cval;
}

void export_slice_in_bmp(const std::string& filename, 
        const std::vector<GLuint>& texture_ids,
        int width,
        int height,
        int vdim,
        int c,
        GLuint colormap,
        const std::vector<float>& min,
        const std::vector<float>& max,
        int repeat)
{
    int n = texture_ids.size();
    std::vector<std::vector<float>> vvalues(n);
    for(int i = 0; i < n; ++i){
        texture_to_vector(texture_ids[i], width, height, vdim, vvalues[i]);
    }

    /* Get the colormap */
    std::vector<unsigned char> cm_pixels;
    texture_colormap_to_vector(colormap, cm_pixels);
    int nb_colors = cm_pixels.size() / 4;

    int hspace = int(0.05*width); 
    bitmap_image image(n * width + (n-1)*hspace, height);
    image.set_all_channels(255, 255, 255);

    for(int t = 0; t < n; ++t){
        for(int i = 0; i < height; ++i){
            for(int j = 0; j < width; ++j){
                float v = vvalues[t][vdim*(i*width+j)+c];
                if(v != ffes::cc){
                    float cv = clamp(v, min[t], max[t], repeat);
                    // TODO: linear interpolation of RGB values in the colormap
                    char red = cm_pixels[4*int(cv*(nb_colors-1))+0];
                    char green = cm_pixels[4*int(cv*(nb_colors-1))+1];
                    char blue = cm_pixels[4*int(cv*(nb_colors-1))+2];
                    image.set_pixel(j + t * (width + hspace), height-1-i,
                            red,green,blue);
                }
            }
        }
    }
    image.save_image(filename);
}

void export_values_texture_to_bmp(
        const std::string& filename, 
        const std::vector<long>& values,
        int width,
        int height,
        GLuint colormap,
        float min,
        float max){

    /* Get the colormap */
    std::vector<unsigned char> cm_pixels;
    texture_colormap_to_vector(colormap, cm_pixels);
    int nb_colors = cm_pixels.size() / 4;

    bitmap_image image(width, height);
    image.set_all_channels(255, 255, 255);

    for(int i = 0; i < height; ++i){
        for(int j = 0; j < width; ++j){
            long id = i*width+j;
            float v = values[id];
            if(v != ffes::cc){
                float cv = clamp(v, min, max, 0);
                // TODO: linear interpolation of RGB values in the colormap
                char red = cm_pixels[4*int(cv*(nb_colors-1))+0];
                char green = cm_pixels[4*int(cv*(nb_colors-1))+1];
                char blue = cm_pixels[4*int(cv*(nb_colors-1))+2];
                image.set_pixel(j, height-1-i, red,green,blue);
            }
        }
    }
    image.save_image(filename);
}


void append_texture_to_raw_stream(std::ostream& out,
        GLuint texture_id,
        int width, int height,
        int vdim, int c){

    /* Get the vector field from the GPU */
    std::vector<float> values;
    texture_to_vector(texture_id, width, height, vdim, values);

    /* Extract the scalar field if necessary */
    if (vdim > 1) {
        printf("[Export raw] Vector field export not implemented\n");
        exit(EXIT_FAILURE);
    }

    out.write(reinterpret_cast<const char*>(&values[0]), values.size() * sizeof(float));
}

std::string get_interpolation_space(const ElementGroup& eg){
    if(ndofs2basis.find(eg.nb_field_cp_per_cell) == ndofs2basis.end()){
        return "Space not found";
    } else {
        return ndofs2basis.at(eg.nb_field_cp_per_cell);
    }
}
std::string get_mapping_space(const ElementGroup& eg){
    if(ndofs2basis.find(eg.nb_mesh_cp_per_cell) == ndofs2basis.end()){
        return "Space not found";
    } else {
        return ndofs2basis.at(eg.nb_mesh_cp_per_cell);
    }
}

}
