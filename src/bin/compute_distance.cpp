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

#include <cfloat>
#include <iomanip>

#include <geogram_basic/basic/common.h>
#include <geogram_basic/basic/command_line.h>
#include <geogram_basic/basic/command_line_args.h>
#include <geogram_basic/basic/file_system.h>
#include <geogram_basic/basic/logger.h>
#include <geogram_basic/basic/progress.h>
#include <geogram_basic/basic/assert.h>
#include <geogram_basic/basic/stopwatch.h>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/transform.hpp>

#include <GL/gl3w.h>

#include "imgui.h"
#include "imgui_impl_glfw_gl3.h"

#include "ffes/field_sampling.h"
#include "ffes/utils.h"
#include "ffes/field_io.h"
#include "ffes/distance.h"

#include <unistd.h> // for usleep

using std::string;
using GEO::Logger;

using ffes::FieldObject;

/* Only for debugging and visualization */
static ffes::DebugShader debug;

// for problems with known analytical solution
// only sinus bump implemented (analytical_id == 1)
void fill_texture(
        int analytical_id, 
        const long voxel_grid[3], 
        double h_pxl,
        double z_slice, 
        GLuint texture_id){
    std::vector<float> val;
    int vdim = 0;
    if(analytical_id == 1){
        vdim = 1;
        val.resize(voxel_grid[0]*voxel_grid[1]*vdim, 0);
        for(int i = 0; i < voxel_grid[1]; ++i){
            for(int j = 0; j < voxel_grid[0]; ++j){
                double x = (j + 0.5) * h_pxl;
                double y = (i + 0.5) * h_pxl;
                double z = z_slice;
                double v = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
                val[i*voxel_grid[1]+j] = v;
            }
        }
    } else {
        printf("error: analytical_id %i not implemented\n", analytical_id);
        exit(EXIT_FAILURE);
    }
    glBindTexture(GL_TEXTURE_RECTANGLE, texture_id);
    glTexSubImage2D(GL_TEXTURE_RECTANGLE, 0, 0, 0, voxel_grid[0], voxel_grid[1], GL_RED, GL_FLOAT, &val[0]);
    glBindTexture(GL_TEXTURE_RECTANGLE, 0);
}


namespace GUI {

    bool Combo(const char* label, int* currIndex, std::vector<std::string>& values);

    /* some info to display in the GUI */
    void show_fieldobject_information(const FieldObject& obj){
        const vector<ffes::ElementGroup>& egs = obj.element_groups();
        const vector<ffes::ElementDecomposition>& eds = obj.element_decompositions();
        for(size_t i = 0; i < egs.size(); ++i){
            const ffes::ElementGroup& eg = egs[i];
            std::string mapping = ffes::get_mapping_space(eg);
            if (eg.mapping.size() > 0) mapping += " (embedded)";
            std::string interp = ffes::get_interpolation_space(eg);
            if (eg.interpolation.size() > 0) interp += " (embedded)";
            ImGui::Text(" ElementGroup %i", i);
            ImGui::Text("  interpolation: %s", interp.c_str());
            ImGui::Text("  mapping:       %s", mapping.c_str());
            ImGui::Text("  nb elements: %i", 
                    (int) eg.mesh_ctrl_points.size() / (eg.nb_mesh_cp_per_cell*eg.mesh_dim));
            ImGui::Text("  control points per elt: %i", eg.nb_mesh_cp_per_cell);
            ImGui::Text("  rendering range: [%i, %i]", eg.rdr_start, eg.rdr_end);
            ImGui::Text("  field points per elt: %i", eg.nb_field_cp_per_cell);
            ImGui::Text("  field dim: %i", obj.field_dimension());

            const ffes::ElementDecomposition& ed = eds[i];
            ImGui::Text("  nb sub-tetrahedra: %i", (int) ed.tet_indices.size() / 4);
        }
    }
}

int main(int argc, char** argv) {

    GEO::initialize();
    GEO::Logger::initialize();
 	GEO::Logger::instance()->set_quiet(false);

    /* Parse program input */
	GEO::CmdLine::declare_arg("gui", false, "Launch a GUI (interactive slicing)");
	GEO::CmdLine::declare_arg("samples", 100, "Nb of samples along the shortest bounding box axis");
	GEO::CmdLine::declare_arg("subdiv_A", 0, "Subdivision level for cells of mesh A (between 0 and 3)");
	GEO::CmdLine::declare_arg("subdiv_B", 0, "Subdivision level for cells of mesh B (between 0 and 3)");
	GEO::CmdLine::declare_arg("restricted_range", true, "Optimization: only draw the element that are potentially intersected");
	GEO::CmdLine::declare_arg("analytical_test", 0, "No of the test: 1 - Sinus bump");
	GEO::CmdLine::declare_arg("output", "nofile.json", "Path used to output the distance results (JSON format)");
	GEO::CmdLine::declare_arg("output_voxel_A", "no_basename", "Write header and raw data files for paraview (slow, only if != no_basename)");
	GEO::CmdLine::declare_arg("output_voxel_B", "no_basename", "Write header and raw data files (slow)");
	GEO::CmdLine::declare_arg("output_voxel_diff", "no_basename", "Write header and raw data files (slow)");
	GEO::CmdLine::declare_arg("timings", false, "Query GPU timings for each slice");
	GEO::CmdLine::declare_arg("double_in_shaders", false, "Use double in shaders (for embedded functions with ct types)");
	GEO::CmdLine::set_arg("sys:assert", "abort");
 	vector<string> filenames;
	GEO::CmdLine::parse(argc, argv, filenames, "fieldA <fieldB> <parameter=value ..>");

    int analytical_id = GEO::CmdLine::get_arg_int("analytical_test");
	if(filenames.size() == 1){
        if(analytical_id != 1){
            Logger::warn("Input") << "Only one input field, using a field set to zero in field B" << std::endl;
        }
        filenames.push_back("field_at_zero");
    }
	if(filenames.size() != 2){
        Logger::err("Input") << "Expecting 2 input fields" << std::endl;
        exit(EXIT_FAILURE);
	} 
    bool gui = GEO::CmdLine::get_arg_bool("gui");
    bool rr = GEO::CmdLine::get_arg_bool("restricted_range");
    bool use_double = GEO::CmdLine::get_arg_bool("double_in_shaders");

    /* Window and OpenGL context. If argument is false, OpenGL context will run in background
     * and no window will be created. */
    int window_w, window_h;
    GLFWwindow* window = ffes::init_gui(gui, window_w, window_h);
    if(gui){
        ImGui_ImplGlfwGL3_Init(window, true);
        ImGui::GetIO().IniFilename = NULL; // do not save run-time preferences
    }
    if(window == NULL){
        Logger::err("Window") << "failed to initialize window and/or GUI" << std::endl;
        return 1;
    }

    /* Build the FieldObject from the program inputs */
    /* Get the bounding box of meshes A and B. */
	GEO::Stopwatch chrono_io("Files IO");
    FieldObject A, B;
    ffes::field_load(filenames[0], A);

    double bbox_min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double bbox_max[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    A.update_bounding_box(bbox_min, bbox_max);

    if(filenames[1] != "field_at_zero"){
        ffes::field_load(filenames[1], B);
        B.update_bounding_box(bbox_min, bbox_max);
    } else {
        B.build_cube_at_zero(bbox_min, bbox_max, A.field_dimension());
    }
    geo_assert(A.field_dimension() == B.field_dimension());
    double io_time = chrono_io.elapsed_time();

    /* Console logs */
    Logger::out("Fields") << "[A] path: " << filenames[0] << ", " << A.nb_cells() << " cells";
    Logger::out("Fields") << ", field dim: " << A.field_dimension() << std::endl;
    if(analytical_id == 1){
        Logger::out("Fields") << "[B] analytical solution of the sinus bump problem" << std::endl;
    } else if (filenames[1] != "field_at_zero"){
        Logger::out("Fields") << "[B] path: " << filenames[1] << ", " << B.nb_cells() << " cells";
        Logger::out("Fields") << ", field dim: " << B.field_dimension() << std::endl;
    }

    /* The elapsed time starts here (so it does not include the reading file timings) */
	GEO::Stopwatch chrono("Slicing loop");

    Logger::out("Sampling") << "Bounding box: [" << bbox_min[0] << "," << bbox_max[0] << "]"
        << "[" << bbox_min[1] << "," << bbox_max[1] << "]"
        << "[" << bbox_min[2] << "," << bbox_max[2] << "]" << std::endl;


    /* Convert the bounding box into the regular voxel grid used to sample the
     * fields.  
     * - nb_pxl_shortest_axis is an important approximation parameter for the
     *   distance computation. It controls the number of samples used in the bounding box,
     *   so it controls the distance between two adjacent sampling points. */
    long voxel_grid[3] = {0, 0, 0};
    int nb_pxl_shortest_axis = GEO::CmdLine::get_arg_int("samples");
    ffes::bbox_to_voxel_grid(nb_pxl_shortest_axis, bbox_min, bbox_max, voxel_grid);
    double h_pxl = (bbox_max[2] - bbox_min[2]) / voxel_grid[2];
    Logger::out("Sampling") << "Voxel grid:    " << voxel_grid[0] << " x " << voxel_grid[1]
        << " x " << voxel_grid[2] << ", h_pxl: " << h_pxl << std::endl;


    /* Values range for colormap if visualization, last is norm of vector */
    vector<double> values_min(A.field_dimension()+1,  DBL_MAX);
    vector<double> values_max(A.field_dimension()+1, -DBL_MAX);
    if(gui){
        A.update_value_ranges(&values_min[0], &values_max[0]);
        B.update_value_ranges(&values_min[0], &values_max[0]);
        for(int d = 0; d < A.field_dimension(); ++d){
            Logger::out("Field values") << "[" << values_min[d] << "," << values_max[d] << "]";
        }
        Logger::out("Field values") << std::endl;
        if (A.field_dimension() > 1) {
            Logger::out("Field values") << "magnitude: [" << values_min[A.field_dimension()] << "," << values_max[A.field_dimension()] << "]";
            Logger::out("Field values") << std::endl;
        }
    }


    /* Set the tetrahedron decomposition used to approximate the geometry of curved cells in
     * A and B. */
    int sub_A = GEO::CmdLine::get_arg_int("subdiv_A");
    int sub_B = GEO::CmdLine::get_arg_int("subdiv_B");
    A.update_element_decomposition_for_each_element_type(sub_A);
    B.update_element_decomposition_for_each_element_type(sub_B);


    /* Prepare and upload the data to the GPU */
    A.finalise();
    B.finalise();


    /* Set the shader programs for each ElementGroup of the FieldObject's. */
    A.update_shader_program_for_each_element_type(debug, use_double);
    B.update_shader_program_for_each_element_type(debug, use_double);


    /* Prepare framebuffer and textures
     * one texture per field dimension (= 1 if scalar, = 3 for displacement, etc) 
     * Rendering framebuffers fbo_A, fbo_B will get the slices (rasterization of fields)
     * fbo_D will compute the difference |A-B| (fragment shader called on a quad) */
    GLuint fbo_A, fbo_B, fbo_D;
    GLuint texture_A, texture_B, texture_D;
    ffes::init_framebuffer_and_texture(A.field_dimension(), voxel_grid[0], voxel_grid[1],
            fbo_A, texture_A);
    ffes::init_framebuffer_and_texture(B.field_dimension(), voxel_grid[0], voxel_grid[1],
            fbo_B, texture_B);
    ffes::init_framebuffer_and_texture(A.field_dimension(), voxel_grid[0], voxel_grid[1],
            fbo_D, texture_D);
    GLuint prog_diff = ffes::create_difference_shader(A.field_dimension());

    /* Reduction framebuffer reduces |A-B| into a column of texture_contributions
     * warning: reduction_fbo size is voxel_grid[1] lines and voxel_grid[2] columns */
    GLuint fbo_reduction;
    GLuint texture_contributions;
    ffes::init_reduction_texture(fbo_reduction, voxel_grid, texture_contributions);
    GLuint prog_reduction = ffes::create_reduction_shader(A.field_dimension());

    /* Cannot call glDraw without a VBO, so we use a empty one */
    GLuint dummy_vao;
    glGenVertexArrays(1, &dummy_vao);

    /* Projection matrix for rendering of slices A and B */
    glm::dmat4 MVP = ffes::compute_model_view_projection_matrix(bbox_min, bbox_max);

    /* Interactive visualization of the slices */
    if(gui){
        int field_coord = 0;
        int field_coord_prev = 0;
        int cm = 0;
        int repeat = 0;
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        std::vector<std::string> colormap_names;
        std::vector<GLuint> colormaps;
        ffes::generate_colormaps(colormaps, colormap_names);
        GLuint prog_display = ffes::create_display_texture_shader(A.field_dimension(), debug);
        float z_slice_f = 0.5 * (bbox_min[2] + bbox_max[2]);
        float range_min[3] = {(float) values_min[field_coord], (float) values_min[field_coord], 0.};
        float range_max[3] = {(float) values_max[field_coord], (float) values_max[field_coord], (float) values_max[field_coord]};
        int display_w, display_h;
        while (!glfwWindowShouldClose(window))
        {
            glfwPollEvents();
            ImGui_ImplGlfwGL3_NewFrame();
            glfwGetFramebufferSize(window, &display_w, &display_h);
            ImGui::SetNextWindowSize(ImVec2(int(0.40*display_w), int(0.25*display_h)));
            ImGui::SetNextWindowPos(ImVec2(int(0.01*display_w), int(0.01*display_h)));
            ImGui::Begin("Controls");
            ImGui::SliderFloat("Z slice", &z_slice_f, bbox_min[2], bbox_max[2]);
            if (z_slice_f < bbox_min[2]) z_slice_f = bbox_min[2];
            if (z_slice_f > bbox_max[2]) z_slice_f = bbox_max[2];
            if(A.field_dimension() > 1){
                ImGui::Text("Last field coordinate is the magnitude of the vector field");
                ImGui::SliderInt("Field coordinate", &field_coord, 0, A.field_dimension());
            }
            if (field_coord != field_coord_prev) {
                range_min[0] = (float) values_min[field_coord]; range_min[1] = (float) values_min[field_coord]; range_min[2] = 0.;
                range_max[0] =  (float) values_max[field_coord]; range_max[1] =  (float) values_max[field_coord]; range_max[2] =  (float) values_max[field_coord];
            }

            ImGui::Text("Field visualization");
            ImGui::RangeSliderFloat("Field A - range", &range_min[0], &range_max[0], values_min[field_coord], values_max[field_coord], "[%.2f, %.2f]");
            ImGui::RangeSliderFloat("Field B - range", &range_min[1], &range_max[1], values_min[field_coord], values_max[field_coord], "[%.2f, %.2f]");
            ImGui::RangeSliderFloat("|A-B|",           &range_min[2], &range_max[2],                      0., values_max[field_coord], "[%.2f, %.2f]");
            GUI::Combo("Colormap", &cm, colormap_names);
            ImGui::SliderInt("Repeat colormap", &repeat, 0, 6);

            ImGui::End();
            ImGui::SetNextWindowSize(ImVec2(int(0.2*display_w), int(0.2*display_h)));
            ImGui::SetNextWindowPos(ImVec2(int(0.75*display_w), int(0.01*display_h)));
            ImGui::Begin("Debug");
            if(ImGui::Button("Export bitmap")){
                std::vector<float> mins = {(float) range_min[0], range_min[2], (float) range_min[1]};
                std::vector<float> maxs = {(float) range_max[0], range_max[2], (float) range_max[1]};
                std::vector<GLuint> texts = {texture_A, texture_D, texture_B};
                ffes::export_slice_in_bmp("slice.bmp", texts,
                        voxel_grid[0], voxel_grid[1],
                        A.field_dimension(), field_coord, colormaps[cm], mins, maxs, repeat);

            }
            ImGui::Selectable("Debug shader mode (need reload)", &debug.debug);
            ImGui::Selectable("  show elements", &debug.show_elements);
            ImGui::Selectable("  show group",    &debug.show_group);
            ImGui::Selectable("  show faces",    &debug.show_faces);
            if(ImGui::Button("Reload shaders")){
                A.update_shader_program_for_each_element_type(debug, use_double);
                B.update_shader_program_for_each_element_type(debug, use_double);
                if (debug.show_elements) {
                    debug.nb_elts = A.nb_cells();
                }
                prog_display = ffes::create_display_texture_shader(A.field_dimension(), debug);
                prog_diff = ffes::create_difference_shader(A.field_dimension());
                prog_reduction = ffes::create_reduction_shader(A.field_dimension());
            }
            ImGui::Selectable("Restricted range", &rr);
            if(ImGui::Button("Print coordinates A (console)")){
                A.print_data();
            }
            if(ImGui::Button("Print coordinates B (console)")){
                B.print_data();
            }
            ImGui::End();

            ImGui::SetNextWindowSize(ImVec2(int(0.2*display_w), int(0.2*display_h)));
            ImGui::SetNextWindowPos(ImVec2(int(0.05*display_w), int(0.80*display_h)));
            ImGui::Begin("FieldObject A");
            GUI::show_fieldobject_information(A);
            ImGui::End();

            ImGui::SetNextWindowSize(ImVec2(int(0.2*display_w), int(0.2*display_h)));
            ImGui::SetNextWindowPos(ImVec2(int(0.75*display_w), int(0.80*display_h)));
            ImGui::Begin("FieldObject B");
            GUI::show_fieldobject_information(B);
            ImGui::End();

            /* Slice rendering */
            index_t slice_no = 0; // target column in texture_contributions
            slice_no = (int) ((z_slice_f - bbox_min[2])/(bbox_max[2] - bbox_min[2]) * (voxel_grid[2]-1));
            geo_assert(slice_no >= 0 && slice_no < voxel_grid[2]);

            /* OpenGL rendering to the textures.*/
            /* Draw A */
            A.draw_slice(fbo_A, 
                    voxel_grid[0], voxel_grid[1],
                    z_slice_f, glm::value_ptr(MVP), rr);

            /* Draw B */
            B.draw_slice(fbo_B, 
                    voxel_grid[0], voxel_grid[1],
                    z_slice_f, glm::value_ptr(MVP), rr);

            if(analytical_id != 0){
                fill_texture(analytical_id, voxel_grid, h_pxl, z_slice_f, texture_B);
            }

            /* Slice difference */
            glBindFramebuffer(GL_FRAMEBUFFER, fbo_D);
            glDrawBuffer(GL_COLOR_ATTACHMENT0);
            ffes::compute_texture_difference(fbo_D, voxel_grid, prog_diff,
                    texture_A, texture_B, texture_D, dummy_vao);

            /* Slice reduction (useless in GUI mode ?) */
            glBindFramebuffer(GL_FRAMEBUFFER, fbo_reduction);
            ffes::slice_reduction(voxel_grid, prog_reduction, texture_D,
                    slice_no, texture_contributions, dummy_vao);
            
            /* Screen display */
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            glBindVertexArray(dummy_vao);
            ffes::display_textures(
                    window, 
                    prog_display,
                    texture_A,
                    texture_D,
                    texture_B,
                    voxel_grid[0],
                    voxel_grid[1],
                    colormaps[cm],
                    repeat,
                    range_min,
                    range_max,
                    A.field_dimension(),
                    field_coord);

            ImGui::Render();
            field_coord_prev = field_coord;
            glfwSwapBuffers(window);
        }
        ImGui_ImplGlfwGL3_Shutdown();
        glfwTerminate();
        return 0;
    }

    /* Check if export of 3D grids (for visualization in paraview, export in very slow) */
    std::ofstream vg_A, vg_B, vg_D; // ofstream will be opened only if export
    std::string name_vxl_A = GEO::CmdLine::get_arg("output_voxel_A");
    std::string name_vxl_B = GEO::CmdLine::get_arg("output_voxel_B");
    std::string name_vxl_D = GEO::CmdLine::get_arg("output_voxel_diff");

    ffes::prepare_voxel_grid_export(voxel_grid, h_pxl, vg_A, vg_B, vg_D,
            name_vxl_A, name_vxl_B, name_vxl_D);

    index_t progress_divider = (voxel_grid[2] > 10000) ? 100 : 1;
    GEO::ProgressTask task("slicing",voxel_grid[2]/progress_divider);

    bool timings = GEO::CmdLine::get_arg_bool("timings");
    GLuint query_slices[voxel_grid[2]+1];
    GLuint query_draw[voxel_grid[2]];
    GLuint query_diff[voxel_grid[2]];
    GLuint query_red[voxel_grid[2]];
    if (timings) {
        glGenQueries(voxel_grid[2]+1, query_slices);
        glGenQueries(voxel_grid[2], query_draw);
        glGenQueries(voxel_grid[2], query_diff);
        glGenQueries(voxel_grid[2], query_red);
    }

    /* Slicing loop, main part of the program */ 
    if (timings) glQueryCounter(query_slices[0], GL_TIMESTAMP);
    for(long k = 0; k < voxel_grid[2]; ++k){
        const double z_slice = bbox_min[2] + 0.5*h_pxl + k * h_pxl;

        /* Update progress bar in console */
        if(!(k%progress_divider)) {
            task.progress(k/progress_divider);
        }

        /* Choose the FBOs in the cycle */
        if (k % 100 == 0) glFinish();

        /* Draw A */
        A.draw_slice(fbo_A, 
                voxel_grid[0], voxel_grid[1],
                z_slice, glm::value_ptr(MVP), rr);

        /* Draw B */
        B.draw_slice(fbo_B, 
                voxel_grid[0], voxel_grid[1],
                z_slice, glm::value_ptr(MVP), rr);

        if (timings) glQueryCounter(query_draw[k], GL_TIMESTAMP);

        /* Analytical solution for validation, fill texture B */
        if(analytical_id != 0){
            fill_texture(analytical_id, voxel_grid, h_pxl, z_slice, texture_B);
        }

        /* D = |A-B| */
        /* Doing |A-B| in a different FBO as A, B are attached to fbo_A, fbo_B
         * so they cannot be sampled directly */
        ffes::compute_texture_difference(fbo_D,
                voxel_grid, prog_diff,
                texture_A, texture_B, texture_D, dummy_vao);

        if (timings) glQueryCounter(query_diff[k], GL_TIMESTAMP);

        /* Optional export (avoid, super slow) */
        if (vg_A.is_open()) ffes::append_texture_to_raw_stream(vg_A, texture_A, voxel_grid[0],
                voxel_grid[1], A.field_dimension(), 0);
        if (vg_B.is_open()) ffes::append_texture_to_raw_stream(vg_B, texture_B, voxel_grid[0],
                voxel_grid[1], A.field_dimension(), 0);
        if (vg_D.is_open()) ffes::append_texture_to_raw_stream(vg_D, texture_D, voxel_grid[0],
                voxel_grid[1], A.field_dimension(), 0);


        /* Store the slice contribution, 
         * reduction from D = |A-B| to column of texture_contributions */
        glBindFramebuffer(GL_FRAMEBUFFER, fbo_reduction);
        ffes::slice_reduction(voxel_grid, prog_reduction, texture_D,
                k, texture_contributions, dummy_vao);

        if (timings) glQueryCounter(query_red[k], GL_TIMESTAMP);

        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glFlush();
        if (timings) glQueryCounter(query_slices[1+k], GL_TIMESTAMP);
    }
    glFinish();


    /* Timings of OpenGL commands */
    if (timings) {
        for (int i = 0; i < voxel_grid[2] + 1; ++i) {
            int available = 0;
            while(!available) {
                glGetQueryObjectiv(query_slices[i], GL_QUERY_RESULT_AVAILABLE, &available);
            }
            int available_d = 0;
            while(!available_d) {
                glGetQueryObjectiv(query_slices[i], GL_QUERY_RESULT_AVAILABLE, &available_d);
            }
            int available_di = 0;
            while(!available_di) {
                glGetQueryObjectiv(query_slices[i], GL_QUERY_RESULT_AVAILABLE, &available_di);
            }
            int available_r = 0;
            while(!available_r) {
                glGetQueryObjectiv(query_slices[i], GL_QUERY_RESULT_AVAILABLE, &available_r);
            }
        }
        GLuint64 ts_slices[voxel_grid[2]+1];
        GLuint64 ts_draw[voxel_grid[2]];
        GLuint64 ts_diff[voxel_grid[2]];
        GLuint64 ts_red[voxel_grid[2]];
        for (int i = 0; i < voxel_grid[2] + 1; ++i) {
            glGetQueryObjectui64v(query_slices[i], GL_QUERY_RESULT, &ts_slices[i]);
            if (i != voxel_grid[2]) {
                glGetQueryObjectui64v(query_draw[i], GL_QUERY_RESULT, &ts_draw[i]);
                glGetQueryObjectui64v(query_diff[i], GL_QUERY_RESULT, &ts_diff[i]);
                glGetQueryObjectui64v(query_red[i], GL_QUERY_RESULT, &ts_red[i]);
            }
        }
        printf("\n");
        for (int i = 0; i < voxel_grid[2]; ++i) {
            float t_slice = (ts_slices[i+1]-ts_slices[i])/1000000.0;
            float t_draw = (ts_draw[i]-ts_slices[i])/1000000.0;
            float t_diff = (ts_diff[i]-ts_draw[i])/1000000.0;
            float t_red = (ts_red[i]-ts_diff[i])/1000000.0;
            printf("[Timings] slice %i: %f ms (draw: %f, diff: %f, reduction: %f)\n", i, t_slice,
                    t_draw, t_diff, t_red);
        }
    }

    /* Download the reduced texture from the GPU to the CPU */
    long nbvals = 4*voxel_grid[1]*voxel_grid[2];
    std::vector<float> rdc(nbvals, 0);
    std::vector<float> rdc_i(nbvals, 0);
    glBindTexture(GL_TEXTURE_RECTANGLE, texture_contributions);
    glGetTexImage(GL_TEXTURE_RECTANGLE, 0, GL_RGBA, GL_FLOAT, &rdc_i[0]);
    for (size_t j = 0; j < rdc_i.size(); ++j) {
        rdc[j] += rdc_i[j];
    }

    /* Computation of the global distances */
    long nnz = 0.;
    double sd = 0.;
    double ssqd = 0.;
    double Li = 0.;


    // debug only
    vector<long> dbg(voxel_grid[1]*voxel_grid[2], 0);
    long dbg_c = 0;
    float dbg_min =  FLT_MAX;
    float dbg_max = -FLT_MAX;
    const bool DBG_RDC = true;

    for(size_t i = 0; i < rdc.size() / 4; ++i){
        if ((long) rdc[4*i+0] < 0) {
            Logger::err("Reduction") << "negative number of samples: " << rdc[4*i] << " at pixel ";
            Logger::err("Reduction") << i << ", abort" << std::endl;
            exit(EXIT_FAILURE);
        }
        nnz  += (long) rdc[4*i+0];
        sd   += (double) rdc[4*i+1];
        ssqd += (double) rdc[4*i+2];
        if(rdc[4*i+3] > Li) Li = rdc[4*i+3];

        // debug only
        if (DBG_RDC) {
            dbg[i] = rdc[4*i];
            if (dbg[i] < dbg_min) dbg_min = dbg[i];
            if (dbg[i] > dbg_max) dbg_max = dbg[i];
            if (rdc[4*i] == 0){
                dbg_c += 1;
            }
        }
    }
    if (DBG_RDC) {
        long nb_slices_at_0 = 0;
        const long max = voxel_grid[0] * voxel_grid[1];
        for(long k = 0; k < voxel_grid[2]; ++k){
            long samples_k = 0; 
            for(long i = 0; i < voxel_grid[1]; ++i){
                samples_k += (long) dbg[voxel_grid[2]*i+k];
            }
            if (samples_k < 0.01* max) {
                double pp = (double) samples_k / max * 100;
                printf("[warning] slice %li, got %li samples (max: %li, prop: %.3f%%)\n", k, samples_k, max, pp);
                nb_slices_at_0 += 1;
            }
        }
        if (nb_slices_at_0 > 0) printf("[warning] %li slices with 0 samples inside\n", nb_slices_at_0);

        /* export */
        if (false) {
            std::vector<GLuint> colormaps;
            std::vector<std::string> colormap_names;
            ffes::generate_colormaps(colormaps, colormap_names);
            ffes::export_values_texture_to_bmp("contrib.bmp", dbg, voxel_grid[2], voxel_grid[1],
                    colormaps[0], dbg_min, dbg_max);
        }
    }
    
    double L1 = sd*std::pow(h_pxl, 3);
    double L2 = std::pow(ssqd*std::pow(h_pxl, 3), 0.5);
    double elapsed_time = chrono.elapsed_time();

    long samples_bbox = long(voxel_grid[0])*long(voxel_grid[1])*long(voxel_grid[2]);
    double prop = double(nnz) / double(samples_bbox);


    Logger::out("Timings") << "Files I/O: " << io_time << "s, Distance computation: ";
    Logger::out("Timings") << elapsed_time << "s" << std::endl;

    /* Output in console */
    std::string mt = "-distance";
    if(filenames[1] == "field_at_zero" && analytical_id != 1){
        mt = "-norm";
    }

    Logger::out("Results") << "samples inside: " << nnz << ", " << 100.*prop;
    Logger::out("Results") << "% of the bounding box" << '\n';
    Logger::out("Results") << std::setprecision(4);
    Logger::out("Results") << "L1" << mt << ":    " << L1 << '\n';
    Logger::out("Results") << "L2" << mt << ":    " << L2 << '\n';
    Logger::out("Results") << "Li" << mt << ":    " << Li << std::endl;

    /* Output in json file */
    vector<std::string> ostr;
    ostr.push_back("field_A"); ostr.push_back("\"" + filenames[0] + "\"");
    ostr.push_back("field_B"); ostr.push_back("\"" + filenames[1] + "\"");
    ostr.push_back("nb_samples"); ostr.push_back(std::to_string(nnz));
    ostr.push_back("h_pxl"); ostr.push_back(to_string_with_precision(h_pxl,8));
    ostr.push_back("bbox_prop"); ostr.push_back(to_string_with_precision(prop,3));
    ostr.push_back("subdivision_lvl_A"); ostr.push_back(std::to_string(sub_A));
    ostr.push_back("subdivision_lvl_B"); ostr.push_back(std::to_string(sub_B));
    ostr.push_back("nb_cells_A"); ostr.push_back(std::to_string(A.nb_cells()));
    ostr.push_back("nb_cells_B"); ostr.push_back(std::to_string(B.nb_cells()));
    ostr.push_back("field_dim"); ostr.push_back(std::to_string(A.field_dimension()));
    ostr.push_back("L1"); ostr.push_back(to_string_with_precision(L1, 8));
    ostr.push_back("L2"); ostr.push_back(to_string_with_precision(L2, 8));
    ostr.push_back("Linf"); ostr.push_back(to_string_with_precision(Li, 8));
    ostr.push_back("elapsed_time"); ostr.push_back(to_string_with_precision(elapsed_time,5));
    if (A.get_metadata().size() > 0){
        const vector<std::string>& meta = A.get_metadata();
        for (size_t i = 0; i < meta.size() / 2; ++i){
            ostr.push_back(meta[2*i+0] + "_A");
            ostr.push_back(meta[2*i+1]);
        }
    }
    if (B.get_metadata().size() > 0){
        const vector<std::string>& meta = B.get_metadata();
        for (size_t i = 0; i < meta.size() / 2; ++i){
            ostr.push_back(meta[2*i+0] + "_B");
            ostr.push_back(meta[2*i+1]);
        }
    }
    std::string file = GEO::CmdLine::get_arg("output");
    if (file != "nofile.json") {
        write_json_file(file, ostr);
    }
    return 0;
}

 
namespace GUI {
    /* Helper for ImGUI, from: https://eliasdaler.github.io/using-imgui-with-sfml-pt2 */
    static auto vector_getter = [](void* vec, int idx, const char** out_text)
    {
        auto& vector = *static_cast<std::vector<std::string>*>(vec);
        if (idx < 0 || idx >= static_cast<int>(vector.size())) { return false; }
        *out_text = vector.at(idx).c_str();
        return true;
    };
    bool Combo(const char* label, int* currIndex, std::vector<std::string>& values)
    {
        if (values.empty()) { return false; }
        return ImGui::Combo(label, currIndex, vector_getter,
                static_cast<void*>(&values), values.size());
    }
}
