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

#include <fstream>
#include <numeric>
#include <algorithm>

#include "ffes/field_io.h"
#include "ffes/field_sampling.h"

#include "third_party/geogram_basic/basic/file_system.h"
#include "third_party/geogram_basic/basic/logger.h"
#include "third_party/geogram_basic/basic/string.h"
#include "third_party/json/json.hpp"
#include "third_party/base64/base64.hpp"


namespace ffes {

    using nlohmann::json;

    namespace LoadUtils {
        float zmin_of(int c, vector<float>& control_pts, int nb_cp_per_cell){
            float z_min = FLT_MAX; 
            for(int lv = 0; lv < nb_cp_per_cell; ++lv){
                if(control_pts[3*(nb_cp_per_cell*c+lv)+2] < z_min){
                    z_min = control_pts[3*(nb_cp_per_cell*c+lv)+2];
                }
            }
            return z_min;
        }
    }

    void field_load_json(const std::string& path, FieldObject& obj){
        std::ifstream input(path, std::ios::binary | std::ios::in);
        json jin;
        input >> jin;
        if (!jin.is_object()) {
            printf("[FieldObject IO][json error] was expecting a root object, file: %s\n", path.c_str());
            exit(EXIT_FAILURE);
        }
        std::string version = jin["version"].get<std::string>();
        if (version != "0.1") {
            printf("[FieldObject IO][json error] was expecting version '0.1' (as a string value), file: %s\n", path.c_str());
            exit(EXIT_FAILURE);
        }
        json j = jin["groups"];
        if (j.is_array()) {
            obj.elt_group.resize(obj.elt_group.size() + 1);
            ElementGroup& eg = obj.elt_group.back();
            for (json::iterator it = j.begin(); it != j.end(); ++it) {
                auto jeg = *it;
                eg.mesh_dim  =  jeg["mesh_dim"].get<index_t>();
                eg.field_dim =  jeg["field_dim"].get<index_t>();
                eg.nb_mesh_cp_per_cell = jeg["nb_mesh_cp_per_cell"].get<index_t>();
                eg.nb_field_cp_per_cell = jeg["nb_field_cp_per_cell"].get<index_t>();
                std::string primitive = jeg["primitive"];
                if(primitive == "TET"){
                    eg.type = TET;
                } else if (primitive == "HEX"){
                    eg.type = HEX;
                } else {
                    printf("[FieldObject IO] primitive = %s not supported", primitive.c_str());
                    exit(EXIT_FAILURE);
                }
                eg.interpolation = jeg["interpolation"].get<std::string>();
                eg.mapping = jeg["mapping"].get<std::string>();

                /* Data */
                std::string encoded_mesh_ctrl_points_ns = jeg["mesh_ctrl_points"].get<std::string>();
                size_t nb_mesh_chars = Base64::decode_len(encoded_mesh_ctrl_points_ns);
                std::vector<float> mesh_ctrl_points_ns((nb_mesh_chars + 3)/4);
                Base64::decode(encoded_mesh_ctrl_points_ns, (unsigned char*) &mesh_ctrl_points_ns[0]);
                mesh_ctrl_points_ns.resize(nb_mesh_chars / 4);

                std::string encoded_field_ctrl_points_ns = jeg["field_ctrl_points"].get<std::string>();
                size_t nb_field_chars = Base64::decode_len(encoded_field_ctrl_points_ns);
                std::vector<float> field_ctrl_points_ns((nb_field_chars + 3)/4);
                Base64::decode(encoded_field_ctrl_points_ns, (unsigned char*) &field_ctrl_points_ns[0]);
                field_ctrl_points_ns.resize(nb_field_chars / 4);

                /* Sort cells by z-min */
                size_t nb_cells = mesh_ctrl_points_ns.size() / (eg.nb_mesh_cp_per_cell * eg.mesh_dim);

                vector<int> cell_id(nb_cells);
                std::iota(cell_id.begin(), cell_id.end(), 0);

                // sort the cell indices
                std::sort(cell_id.begin(), cell_id.end(),
                        [&](index_t a, index_t b)->bool{
                        return LoadUtils::zmin_of(a, mesh_ctrl_points_ns, eg.nb_mesh_cp_per_cell) 
                        < LoadUtils::zmin_of(b, mesh_ctrl_points_ns, eg.nb_mesh_cp_per_cell); });

                // fill the data structure
                eg.mesh_ctrl_points.resize(mesh_ctrl_points_ns.size());
                eg.field_ctrl_points.resize(field_ctrl_points_ns.size());
                for(size_t e = 0; e < nb_cells; ++e){
                    int oe = cell_id[e]; // index in the input file/vector
                    for(unsigned int lv = 0; lv < eg.nb_mesh_cp_per_cell; ++lv){
                        for(unsigned int d = 0; d < eg.mesh_dim; ++d){
                            eg.mesh_ctrl_points[eg.mesh_dim*(e*eg.nb_mesh_cp_per_cell+lv)+d]
                                = mesh_ctrl_points_ns[eg.mesh_dim*(oe*eg.nb_mesh_cp_per_cell+lv)+d];
                        }
                    }
                    for(unsigned int lv = 0; lv < eg.nb_field_cp_per_cell; ++lv){
                        for(unsigned int d = 0; d < eg.field_dim; ++d){
                            eg.field_ctrl_points[eg.field_dim*(e*eg.nb_field_cp_per_cell+lv)+d]
                                = field_ctrl_points_ns[eg.field_dim*(oe*eg.nb_field_cp_per_cell+lv)+d];
                        }
                    }
                }
            }
        } else {
            printf("[FieldObject IO][json error] was expecting an array named 'groups', file: %s\n", path.c_str());
            exit(EXIT_FAILURE);
        }
        for (json::iterator it = jin.begin(); it != jin.end(); ++it) {
            if (it.value().is_string()) {
                obj.other_metadata.push_back(it.key());
                obj.other_metadata.push_back(it.value().get<std::string>());
            }
            if (it.value().is_number()) {
                obj.other_metadata.push_back(it.key());
                obj.other_metadata.push_back(std::to_string(it.value().get<double>()));
            }
        }
    }

    void field_load(const std::string& path, FieldObject& obj){
        if(!GEO::FileSystem::is_file(path)){
            GEO::Logger::err("FieldObject IO") << "no file found at " << path << std::endl;
            exit(EXIT_FAILURE);
        }
        if (GEO::FileSystem::extension(path) == "json") {
            field_load_json(path, obj);
        } else {
            GEO::Logger::err("FieldObject IO") << "extension not supported for file: " << path << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}
