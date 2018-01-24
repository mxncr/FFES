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

#include <GLFW/glfw3.h>

/* glm is only used for the projection matrix,
 * should be replaced in the future */
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

typedef unsigned int index_t;
using std::vector;

namespace ffes {
    enum ReferenceElementPrimitive { TET, HEX };
    const float cc = 999999; /* NO_DATA value */

    /* 
     * Tetrahedral mesh of a reference element. 
     * This controls the piecewise-linear approximation of cells with
     * curved/polynomial faces.
     */
    struct ElementDecomposition {
        void upload_to_gpu(GLuint vao_id);
        ~ElementDecomposition(); /* delete GPU buffers if uploaded on it */

        ReferenceElementPrimitive primitive;
        index_t subdivision_level;
        vector<float> vertices;
        vector<index_t> tet_indices;

        /* OpenGL bindings */
        GLuint vao_id;
        GLuint ibo_id;
        GLuint vbo_id;
    };

    /**
     * \brief An ElementGroup store all the cells that share the same
     *        mapping and the same interpolation function.
     */
    struct ElementGroup {
        ReferenceElementPrimitive type;
        /* 
         * The mesh_ctrl_points and field_coefficients are defined cell by
         * cell (chunk by chunk), so values are repeated for
         * points/coefficients shared by adjacent cells. This avoids
         * complicated indirections on the GPU.  
         * */
        index_t nb_mesh_cp_per_cell;
        index_t nb_field_cp_per_cell;
        index_t mesh_dim;
        index_t field_dim;
        vector<float> mesh_ctrl_points;
        vector<float> field_ctrl_points;

        /* Each element group can have specific/custom mapping and interpolation
         * functions. These functions must strictly respect the variable names 
         * used in the vertex shader for mapping and in the fragment shader for
         * interpolation. */
        std::string mapping;
        std::string interpolation;

        /* Only elements in [rdr_start, rdr_end] are rendered when calling the
         * draw command on the ElementGroup. This is useful when elements are
         * sorted and the range is updated at each move of the slicing plane.
         * */
        bool use_rdr_range;
        index_t rdr_start;
        index_t rdr_end;
        double prev_zslice;

        /* OpenGL bindings */
        GLuint vao_id;
        GLuint ctrl_points_vbo;
        GLuint coefficients_ssbo;
        GLuint program;
    };


    /* Options for debugging via shaders */
    struct DebugShader {
        bool debug         = false;
        bool show_elements = false;
        long nb_elts       = 0;
        bool show_group    = false;
        bool show_faces    = false;
    };

    /* A FieldObject is the representation of a mesh and its associated function
     * space.
     * A FieldObject contains one ElementGroup per type of element, where a type
     * of element is defined by a geometric primitive (tetrahedron, hexahedron, etc), 
     * its mapping type (linear, quadratic, trilinear, triquadratic, etc) and the
     * function space defined on the geometric primitive (P1, P2, Q1, etc)
     * */
    class FieldObject {
        public:
            FieldObject(){}
            
            void build_cube_at_zero( 
                    const double bbox_min[3],
                    const double bbox_max[3],
                    index_t field_dim);

            /**
             * \brief  Set the subdivision level for each ElementDecomposition associated to
             *         each ElementGroup.
             *
             * \param subdivision_level Subdivision level
             *      0 - 1 tet per tet or 6 tet per hex
             *      1 - 8 tets per tet or 48 tets per hex
             *      2 - 64 tets per tet or 384 tets per hex
             */
            void update_element_decomposition_for_each_element_type(index_t subdivision_level);

            /**
             * \brief  Set the subdivision level for each ElementDecomposition associated to
             *         each ElementGroup.
             *
             * \param subdivision_levels subdivision level for each ElementGroup
             *      0 - 1 tet per tet or 6 tet per hex
             *      1 - 8 tets per tet or 48 tets per hex
             *      2 - 64 tets per tet or 384 tets per hex
             */
            void update_element_decomposition_for_each_element_type(const vector<index_t>& subdivision_levels);

            /**
             * \brief Build the shaders used to slice each ElementGroup. It retrieves the mapping
             *        and interpolation functions by calling functions in src/shader_functions.hpp.
             *
             * \param debug activate the debug mode in the shaders (different outputs, etc)
             * \param computation_with_double define the types inside functions as double, only
             *          available if the type used is ct/ct2/ct3/ct4
             */
            void update_shader_program_for_each_element_type(
                    DebugShader& debug,
                    bool computation_with_double = false);

            /**
             * \brief  Upload all the content of the ElementGroup's and ElementDecomposition's to
             *         the GPU.
             *         Important: should never be called before the FieldObject is filled (via 
             *         upload_element_decomposition_for_each_element_type() and load_from_...())
             */
            void finalise(); 

            /**
             * \brief  Sample the field at voxel centers. 
             *
             * \param[in] z_slice Position of the slice along the z-axis
             * \param[in] MVP[16] Projection matrix (computed from the bounding box)
             * \param[in] use_limited_elt_range Will update the restricted rendering range based
             *          on z_slice. Cells have to be ordered by z-min in each ElementGroup for
             *          this option to work. It is not the case, the output will be garbage. 
             */
            void draw_slice(
                    GLuint fbo,
                    int width,
                    int height,
                    double z_slice,
                    double MVP[16],
                    bool use_limited_elt_range = true);

            /**
             * \brief Increase the size of the bounding box to include the current FieldObject
             *
             * \param[in,out] bbox_min[3]
             * \param[in,out] bbox_max[3]
             */
            void update_bounding_box(double bbox_min[3], double bbox_max[3]);

            /**
             * \brief Increase the value range (for each dimension of the field) to include
             *        the current FieldObject
             *
             * \param[in,out] values_min[]
             * \param[in,out] values_max[]
             */
            void update_value_ranges(double values_min[/*field_dim*/], 
                    double values_max[/*field_dim*/]);

            /**
             * \brief  Returns the dimension of field (1 for scalar field, etc)
             * \return  the field dimension 
             */
            inline int field_dimension() const {
                assert(elt_group.size() > 0);
                return elt_group[0].field_dim;
            }

            /**
             * \brief  Returns the total number of cells (adding all the ElementGroup's)
             * \return  the number of cells in the FieldObject 
             */
            long nb_cells();

            /**
             * \brief  For debugging only (all the coordinates used in the mappings)
             */
            void print_data();

            /**
             * \brief  Reading access to low level data structure
             */
            const vector<ElementGroup>& element_groups() const {
                return elt_group;
            }

            /**
             * \brief  Read/write access to low level data structure
             */
            const vector<ElementDecomposition>& element_decompositions() const {
                return ref_elt;
            }
            vector<ElementGroup>& element_groups() {
                return elt_group;
            }
            vector<ElementDecomposition>& element_decompositions() {
                return ref_elt;
            }

            const vector<std::string>& get_metadata() const {
                return other_metadata;
            }

        public:
            void update_element_rendering_range(double z_slice);
            vector<ElementGroup> elt_group;
            vector<ElementDecomposition> ref_elt;
            vector<std::string> other_metadata; /* size is pair, even is key, odd is value.
                                                   Can be used to store information */
    };


    /**
     * \brief Build the voxel grid dimensions for the bounding box
     *        and the number of voxels along the shortest axis of the
     *        bounding box
     * \param[in] nb_pixel_on_shortest_axis Nb of pixels along the shortest axis, e.g. 100
     * \param[in] bbox_min[3] Coordinates of the bounding box lower corner
     * \param[in] bbox_max[3] Coordinates of the opposite corner of the bounding box
     * \param[out] voxel_grid[3] The voxel grid dimensions
     */
    void bbox_to_voxel_grid(
            long nb_pixel_on_shortest_axis,
            const double bbox_min[3], 
            const double bbox_max[3], 
            long voxel_grid[3]
            );


    /**
     * \brief  Build the projection matrix used for slicing from the bounding box
     *         Changing this method will affect many things ! (see comments in the code)
     *
     * \param bbox_min[3] Coordinates of the bounding box lower corner
     * \param bbox_max[3] Coordinates of the opposite corner of the bounding box
     *
     * \return  the projection matrix 
     */
    glm::dmat4 compute_model_view_projection_matrix(const double bbox_min[3], 
            const double bbox_max[3]);

    /**
     * \brief Initialize the OpenGL context
     *
     * \param visible If false, no window will be created (OpenGL run in background)
     * \param window_w If visible, width of the window
     * \param window_h If visible, height of the window
     *
     * \return   
     */
    GLFWwindow* init_gui(bool visible, int& window_w, int& window_h);

    /**
     * \brief  Initialize another background framebuffer that will be used for texture reduction.
     *
     * \param[out] fbo_id OpenGL id of the framebuffer (different than the slicing FBO)
     * \param[in] voxel_grid[3] Voxel grid dimensions
     * \param[out] texture_contributions OpenGL id of the texture that will store the reductions (see below)
     *
     * The reduction texture, or texture_contributions, is used to store
     * the contributions from each slice comparison (field A vs field B) to
     * the global distance computations.
     * One slice comparison is reduced into one column of texture_contributions:
     * The value (i) of a given column (j) is the reduction of the i-th row of
     * the j-th slice comparison.
     * For a given RGBA value at (i,j):
     * - R is the number of pixel where field A and field B are both defined in row i
     * - G is the sum of the absolute differences in row i (used for L1 distance)
     * - B is the sum of the squared differences in row i (used for L2 distance)
     * - A is the maximum of the absolute differences in row i (used for L_infty distance)
     */
    void init_reduction_texture(
            GLuint& fbo_id,
            const long  voxel_grid[3],
            GLuint& texture_contributions);

    /**
     * \brief Display the input textures in the window
     *
     * \param[in] window Target window for display
     * \param[in] shader_program Program used to display the textures
     * \param[in] texture_left OpenGL id of the texture to display in left
     * \param[in] texture_center OpenGL id of the texture to display in center
     * \param[in] texture_right OpenGL id of the texture to display in right
     * \param[in] texture_width Width of the textures
     * \param[in] texture_height Height of the textures
     * \param[in] texture_colormap OpenGL id of the colormap texture
     * \param[in] value_min Min range for the field_coord-th coordinate of the left and rigth textures
     * \param[in] value_max Max range for the field_coord-th coordinate of the left and rigth textures
     * \param[in] diff_min Min range for the field_coord-th coordinate of the center texture
     * \param[in] diff_max Max range for the field_coord-th coordinate of the center texture
     * \param[in] field_dim Dimension of the field
     * \param[in] field_coord Coordinate to display (0 <= field_coord < field_dim)
     */
    void display_textures(
            GLFWwindow* window, GLuint shader_program,
            GLuint texture_left, GLuint texture_center, GLuint texture_right,
            int texture_width, int texture_height, 
            GLuint texture_colormap,
            int repeat,
            float range_min[3],
            float range_max[3],
            index_t field_dim, index_t field_coord);

    /**
     * \brief  Build the shader used to compute the difference between two textures
     * \param field_dim Dimension of the field
     * \return  OpenGL id of the shader 
     */
    GLuint create_difference_shader(index_t field_dim);

    /**
     * \brief  Build the shader used to reduce the values in a texture to a column
     * \param field_dim Dimension of the field
     * \return  OpenGL id of the shader 
     */
    GLuint create_reduction_shader(index_t field_dim);

    /**
     * \brief  Build the shader used to display a (single) texture
     * \param field_dim Dimension of the field
     * \return  OpenGL id of the shader 
     */
    GLuint create_display_texture_shader(index_t field_dim, DebugShader& debug);

    /*
     * \brief  Generate a 1D texture with a colormap (for display)
     * \return  OpenGL id of the texture
     */
    GLuint init_texture_colormap();

    void generate_colormaps(std::vector<GLuint>& texture_ids, std::vector<std::string>& names);

    void texture_to_vector(GLuint texture_id, int width, int height, int field_dim,
            std::vector<float>& values);

    void export_slice_in_bmp(const std::string& filename, 
            const std::vector<GLuint>& texture_ids,
            int width,
            int height,
            int vdim,
            int c,
            GLuint colormap_id,
            const std::vector<float>& min,
            const std::vector<float>& max,
            int repeat = 0);

    void export_values_texture_to_bmp(
            const std::string& filename, 
            const std::vector<long>& values,
            int width,
            int height,
            GLuint colormap,
            float min,
            float max);

    void append_texture_to_raw_stream(std::ostream& out,
            GLuint texture_id,
            int width, int height,
            int vdim, int c);

    std::string get_interpolation_space(const ElementGroup& eg);
    std::string get_mapping_space(const ElementGroup& eg);
}



