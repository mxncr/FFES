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
 * Some functions of this file are copied and/or modified from the 
 * MFEM library (license below). Mainly from the file: mfem/fem/fe.cpp
 *
 *  Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
 *  the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
 *  reserved. See file COPYRIGHT for details.
 *
 *  This file is part of the MFEM library. For more information and source code
 *  availability see http://mfem.org.
 *
 *  MFEM is free software; you can redistribute it and/or modify it under the
 *  terms of the GNU Lesser General Public License (as published by the Free
 *  Software Foundation) version 2.1 dated February 1999.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "mfem.hpp"

#include "third_party/json/json.hpp"
#include "third_party/base64/base64.hpp"

namespace Convert {
    using namespace mfem;
    using nlohmann::json;

    template <typename T>
    std::string to_string_with_precision(const T a_value, const int n = 6)
    {
        std::ostringstream out;
        out << std::setprecision(n) << a_value;
        return out.str();
    }

    const std::string vdim2type[5] = {"", "float", "vec2", "vec3", "vec4"};
    const std::string zero_of_type[5] = { "", "0.", "vec2(0.,0.)",
        "vec3(0.,0.,0.)", "vec4(0.,0.,0.,0.)"
    };

    std::string generate_Qp_function(const mfem::FiniteElement* fe,
            const std::string& function_name,
            int vdim){
        int p = fe->GetOrder();
        const H1_HexahedronElement* hex = dynamic_cast<const H1_HexahedronElement*>(fe);
        const TriLinear3DFiniteElement* Q1 = dynamic_cast<const TriLinear3DFiniteElement*>(fe);
        MFEM_VERIFY(Q1 || hex, "should be H1_HexahedronElement");

        const int precision_f = std::numeric_limits<float>::max_digits10;
        // Get the p+1 lobatto points
        const double *cp = poly1d.ClosedPoints(p); // TODO Update for MFEM 3.3
        std::string lpts = "const float lobatto_pts[" + std::to_string(p+1) + "] = float[](";
        for(int i = 0; i < p + 1; ++i){
            std::string end = "";
            if(i != p + 1 - 1) end = ", ";
            if(i != 0 && i % 12 == 0) end  += "\n    ";
            lpts += to_string_with_precision(cp[i], precision_f) + end;
        }
        lpts += ");\n";

        // Get the w coefficients
        float w[p+1];
        for (int i = 0; i < p+1; ++i){w[i] = 1.;}
        for (int i = 0; i <= p; i++){
            for (int j = 0; j < i; j++){
                float xij = cp[i] - cp[j];
                w[i] *=  xij;
                w[j] *= -xij;
            }
        }
        for (int i = 0; i <= p; i++){
            w[i] = 1.0/w[i];
        }

        std::string wcoefs = "const float w[" + std::to_string(p+1) + "] = float[](";

        for(int i = 0; i < p + 1; ++i){
            std::string end = "";
            if(i != p + 1 - 1) end = ", ";
            if(i != 0 && i % 12 == 0) end  += "\n    ";
            wcoefs += to_string_with_precision(w[i], precision_f) + end;
        }
        wcoefs += ");\n";


        // Get the dof_map
        int nb = 8;
        if (hex) nb = hex->GetDofMap().Size();
        std::string map = "const int dof_map[" + std::to_string(nb) + "] = int[](\n    ";
        if (hex) {
            for(int i = 0; i < nb;  ++i){
                std::string end = "";
                if(i != nb - 1) end = ", ";
                if(i != 0 && i % 12 == 0) end  += "\n    ";
                map += std::to_string(hex->GetDofMap()[i]) + end;
            }
        } else {
            map += "0, 1, 3, 2, 4, 5, 7, 6 ";
        }
        map += ");\n";

        std::string sp = std::to_string(p);
        std::string sp1 = std::to_string(p+1);
        /* Next function is copied from MFEM (mfem/fem/fe.cpp) */
        std::string fct =
            " void poly1d(in float y, out float u["+sp1+"]){            \n";
        fct += 
            "     int i = "+sp+";                                       \n"
            "     int k = "+sp+";                                       \n";
        fct += 
            "     float lk = 1.;                                        \n";
        fct += 
            "     for (k = 0; k < "+sp+"; k++){                         \n"
            "         if (y >= (lobatto_pts[k] + lobatto_pts[k+1])/2){  \n"
            "             lk *= y - lobatto_pts[k];                     \n"
            "         } else {                                          \n"
            "             for (i = k+1; i <= "+sp+"; i++) {             \n"
            "                 lk *= y - lobatto_pts[i];                 \n"
            "             }                                             \n"
            "             break;                                        \n"
            "         }                                                 \n"
            "     }                                                     \n";
        fct +=
            "     float l = lk * (y - lobatto_pts[k]);                  \n";
        fct +=
            "     for (i = 0; i < k; i++) {                             \n"
            "         u[i] = l * w[i] / (y - lobatto_pts[i]);           \n"
            "     }                                                     \n"
            "     u[k] = lk * w[k];                                     \n"
            "     for (i++; i <= "+sp+"; i++) {                         \n"
            "         u[i] = l * w[i] / (y - lobatto_pts[i]);           \n"
            "     }                                                     \n"
            " };                                                        \n";
        std::string header = lpts + wcoefs + map + fct;

        /* The evaluation loop */
        std::string zero = zero_of_type[vdim];
        std::string type = vdim2type[vdim];
        /* Next function is copied from MFEM (mfem/fem/fe.cpp) */
        std::string content =
            " float shape_x["+sp+"+1];                                                      \n"
            " float shape_y["+sp+"+1];                                                      \n"
            " float shape_z["+sp+"+1];                                                      \n"
            " poly1d(ref_pos.x, shape_x);                                                   \n"
            " poly1d(ref_pos.y, shape_y);                                                   \n"
            " poly1d(ref_pos.z, shape_z);                                                   \n"
            " " + type + " result = " + zero + ";                                           \n"
            " int n = 0;                                                                    \n"
            " for (int k = 0; k <= "+sp+"; k++){                                            \n"
            "     for (int j = 0; j <= "+sp+"; j++){                                        \n"
            "         for (int i = 0; i <= "+sp+"; i++){                                    \n"
            "             result += shape_x[i]*shape_y[j]*shape_z[k]*values[dof_map[n]];    \n"
            "             n += 1;                                                           \n"
            "         }                                                                     \n"
            "     }                                                                         \n"
            " }                                                                             \n";

        int nb_values = (p+1)*(p+1)*(p+1);
        std::string eval_str = header + "\n" +
            type + " " + function_name + 
            "(in vec3 ref_pos, in " + type + " values["
            + std::to_string(nb_values) + "]){ \n"
            + content;                                   // declare and assign variable "result"
        eval_str += "    return result;\n}\n";
        return eval_str;

    }

    std::string generate_function_content(int vdim, int p, const std::vector<float>& coefficients){
        const std::string vdim2type[5] = {"", "float", "vec2", "vec3", "vec4"};
        const std::string zero_of_type[5] = { "", "0.", "vec2(0.,0.)",
            "vec3(0.,0.,0.)", "vec4(0.,0.,0.,0.)"
        };
        const std::string type = vdim2type[vdim];
        const std::string zero = zero_of_type[vdim];
        const std::string sp = std::to_string(p);
        int npts = (p+1)*(p+2)*(p+3)/6;
        const std::string npts_str = std::to_string(npts);

        std::string content = 
            "    " + type + " result = " + zero + ";\n"
            "    const float u1 = 1. - ref_pos.x - ref_pos.y - ref_pos.z;                                 \n"    
            "    const float v1 = ref_pos.x;                                                              \n"
            "    const float w1 = ref_pos.y;                                                              \n"
            "    const float z1 = ref_pos.z;                                                              \n";
        std::string vars[4] = {"u", "v", "w", "z"};
        for (int c = 0; c < 4; ++c){
            for (int i = 2; i <= p; ++i) {
                content += "    const float "+vars[c]+std::to_string(i)+" = "+vars[c]+"1";
                for (int j = 2; j <= i; ++j) {
                    content += "*"+vars[c]+"1";
                }
                content += ";\n";
            }
        }
        content += "    float r = 0.; \n";
        for (int c = 0; c < npts; ++c){
            content += "    r = 0.";
            int lc = 0;
            for (int k = 0; k <= p; k++){
                for (int j = 0; j + k <= p; j++){
                    for (int i = 0; i + j + k <= p; i++){
                        int l = p - i - j - k;
                        if (std::abs(coefficients[c*npts+lc]) > 1e-12){
                            content += " + " + to_string_with_precision(coefficients[c*npts+lc], 6);
                            if (i > 0) content += "*u"+std::to_string(i);
                            if (j > 0) content += "*v"+std::to_string(j);
                            if (k > 0) content += "*w"+std::to_string(k);
                            if (l > 0) content += "*z"+std::to_string(l);
                        }
                        lc += 1;
                    }
                }
            }
            content += ";\n";
            content += "    result += r * values["+std::to_string(c)+"];\n";
        }
        return content;
    }

    std::string generate_Pp_function(const mfem::FiniteElement* fe,
            const std::string& function_name,
            int vdim){
        MFEM_VERIFY(fe->GetGeomType() == Geometry::TETRAHEDRON, "only implemented for tets");
        std::string fct;
        /* This function expresses the interpolation function associated to fe
         * as a combination of monomials lambda_0^i lambda_1^j lambda_2^k lambda_3^l.
         * The matrix A is the linear system that enforces phi_i(x_j) = delta_ij
         * The i-th row of the matrix/table coefficients are the coefficients of
         * phi_i expressed in the monomials basis.
         * The ordering is determined by the ordering of the nodes in fe.
         */
        const mfem::IntegrationRule& ir = fe->GetNodes();
        int p = fe->GetOrder();
        int npts = ir.Size();
        mfem::DenseMatrix A(ir.Size());
        A = 0.;
        for (int c = 0; c < ir.Size(); ++c){
            mfem::IntegrationPoint ip = ir.IntPoint(c);
            double u = 1. - ip.x - ip.y - ip.z;
            double v = ip.x;
            double w = ip.y;
            double z = ip.z;
            int lc = 0;
            for (int k = 0; k <= p; k++){
                for (int j = 0; j + k <= p; j++){
                    for (int i = 0; i + j + k <= p; i++){
                        int l = p - i - j - k;
                        // printf(" (%.2f^%i %.2f^%i %.2f^%i %.2f^%i=)", u,i,v,j,w,k,z,l);
                        A(c, lc) = std::pow(u, i)*std::pow(v, j)*std::pow(w, k)*std::pow(z, l);
                        lc += 1;
                    }
                }
            }
        }
        A.Invert();
        std::vector<float> coefficients;
        for (int c = 0; c < ir.Size(); ++c){
            mfem::Vector V(ir.Size());
            V = 0.; V[c] = 1.;
            mfem::Vector R(ir.Size());
            R = 0.;
            A.Mult(V, R);
            for(int i = 0; i < R.Size(); ++i){
                std::string end = ", ";
                std::string newline = "";
                if(i == R.Size() - 1) newline = "\n    ";
                if(i == R.Size() - 1 && c == npts - 1) end = "";
                coefficients.push_back(R[i]);
            }
        }

        const std::string type = vdim2type[vdim];
        const std::string zero = zero_of_type[vdim];
        const std::string sp = std::to_string(p);
        const std::string npts_str = std::to_string(npts);

        std::string content = generate_function_content(vdim, p, coefficients);

        MFEM_VERIFY(function_name == "element_interpolation"
                || function_name == "element_mapping", "Expect element_interpolation or element_mapping");

        std::string eval_str =
            type + " " + function_name + 
            "(in vec3 ref_pos, in " + type + " values["
            + std::to_string(npts) + "]){ \n"
            + content                                   // declare and assign variable "result"
            + "    return result;\n}\n";

        fct = eval_str;
        return fct;
    }

    std::string generate_function(const mfem::FiniteElement* fe,
            const std::string& function_name,
            int vdim){
        const H1_HexahedronElement* Qp = dynamic_cast<const H1_HexahedronElement*>(fe);
        const H1_TetrahedronElement* Pp = dynamic_cast<const H1_TetrahedronElement*>(fe);
        const Quadratic3DFiniteElement* P2 = dynamic_cast<const Quadratic3DFiniteElement*>(fe);
        const Linear3DFiniteElement* P1 = dynamic_cast<const Linear3DFiniteElement*>(fe);
        const TriLinear3DFiniteElement* Q1 = dynamic_cast<const TriLinear3DFiniteElement*>(fe);
        if ((Q1 || Qp) && !(Pp || P2 || P1)) {
            return generate_Qp_function(fe, function_name, vdim);
        } else if ((Pp || P2 || P1) && !Qp) {
            return generate_Pp_function(fe, function_name, vdim);
        } 
        return "failed";
    }

    double largest_edge(const mfem::Mesh* mesh){
        Array<int> v;
        double largest = 0.;
        for (int i = 0; i < mesh->GetNEdges(); ++i){
            mesh->GetEdgeVertices(i, v);
            const double* v0 = mesh->GetVertex(v[0]);
            const double* v1 = mesh->GetVertex(v[1]);
            double length = 0.;
            for (int j = 0; j < 3; ++j) {
                length += std::pow(v0[j]-v1[j],2);
            }
            length = std::pow(length, 0.5);
            if (length > largest) largest = length;
        }
        return largest;
    }

    /* Export per cell data (mesh + field) (+ interpolation functions for tets) */
    void ExportToFieldObjectFormat(
            const std::string& filename,
            mfem::GridFunction* gf){

        /* Initialize the json array containing the ElementGroup */
        json ja = json::array();

        FiniteElementSpace* fes = gf->FESpace();
        const mfem::Mesh* mesh = fes->GetMesh();
        mfem::Array<int> dofs;
        mfem::Array<int> vertices;
        uint32_t nb_elt_groups = 0;
        int first_elt_of_group[3] = {-1, -1, -1};
        int nb_elt_in_group[3] = {0, 0, 0};

        uint32_t mesh_dim = fes->GetMesh()->SpaceDimension();
        uint32_t field_dim = fes->GetVDim();

        /* Assuming only one ElementGroup
         * (default in MFEM library) */
        first_elt_of_group[0] = 0;
        nb_elt_in_group[0] = mesh->GetNE();
        nb_elt_groups = 1;

        /* Loop over the element groups */
        for(unsigned int g = 0; g < nb_elt_groups; ++g){
            /* Test the first element of the group */
            ElementTransformation* elttrans_fe = fes->GetElementTransformation(first_elt_of_group[g]);
            IsoparametricTransformation* isotrans_fe = dynamic_cast<IsoparametricTransformation*>(elttrans_fe);
            uint32_t nb_mesh_cp_per_cell = 0;
            if (isotrans_fe) {
                nb_mesh_cp_per_cell = isotrans_fe->GetFE()->GetDof();
            } else {
                nb_mesh_cp_per_cell = mesh->GetElement(first_elt_of_group[g])->GetNVertices();
            }
            MFEM_ASSERT(nb_mesh_cp_per_cell != 0, "should not happen");
            fes->GetElementDofs(first_elt_of_group[g], dofs);
            uint32_t nb_field_cp_per_cell = dofs.Size();

            /* Get the mesh and field data */
            uint64_t cp_size = mesh_dim * nb_mesh_cp_per_cell * nb_elt_in_group[g];
            uint64_t fc_size = field_dim * nb_field_cp_per_cell * nb_elt_in_group[g];
            std::vector<float> mesh_cp_values(cp_size);
            std::vector<float> field_cp_values(fc_size);
            for(long e = 0; e < mesh->GetNE(); ++e){
                /* Mesh control points */
                float* elt_coords = &mesh_cp_values[e*mesh_dim*nb_mesh_cp_per_cell];
                IsoparametricTransformation* isotrans = dynamic_cast<IsoparametricTransformation*>(fes->GetElementTransformation(e));
                if (!isotrans) {
                    mesh->GetElementVertices(e, vertices);
                    for(int lv = 0; lv < vertices.Size(); ++lv){
                        elt_coords[3*lv+0] = mesh->GetVertex(vertices[lv])[0];
                        elt_coords[3*lv+1] = mesh->GetVertex(vertices[lv])[1];
                        elt_coords[3*lv+2] = mesh->GetVertex(vertices[lv])[2];
                    }
                } else {
                    DenseMatrix& pm = isotrans->GetPointMat();
                    MFEM_ASSERT(pm.Width() == nb_mesh_cp_per_cell, "problem");
                    MFEM_ASSERT(pm.Height() == mesh_dim, "problem");
                    for(int lv = 0; lv < pm.Width(); ++lv){
                        elt_coords[3*lv+0] = pm(0, lv);
                        elt_coords[3*lv+1] = pm(1, lv);
                        elt_coords[3*lv+2] = pm(2, lv);
                    }
                }

                /* Field coefficients */
                fes->GetElementDofs(e, dofs);
                float* field_coefs = &field_cp_values[e*field_dim*nb_field_cp_per_cell];
                for(int lv = 0; lv < dofs.Size(); ++lv){
                    mfem::Array<int> vdofs;
                    vdofs.SetSize(1);
                    vdofs[0] = dofs[lv];
                    fes->DofsToVDofs(vdofs);
                    for(int d = 0; d < vdofs.Size(); ++d){
                        field_coefs[field_dim*lv+d] = (*gf)[vdofs[d]];
                    }
                }
            }

            /* Encode the vectors of floats into strings */
            std::string encoded_mesh = Base64::encode((unsigned char*) &mesh_cp_values[0],mesh_cp_values.size()*sizeof(float));
            std::string encoded_field = Base64::encode((unsigned char*) &field_cp_values[0],field_cp_values.size()*sizeof(float));

            /* Json object associated to this ElementGroup */
            json jg = json::object();

            jg["mesh_dim"] = mesh_dim;
            jg["primitive"] = mesh->GetElementType(first_elt_of_group[g]) == Geometry::TETRAHEDRON ? "TET" : "HEX";
            jg["nb_mesh_cp_per_cell"] = nb_mesh_cp_per_cell;
            jg["mesh_ctrl_points"] = encoded_mesh;
            ElementTransformation* elttrans = fes->GetElementTransformation(first_elt_of_group[g]);
            IsoparametricTransformation* isotrans = dynamic_cast<IsoparametricTransformation*>(elttrans);
            jg["mapping"] = generate_function(isotrans->GetFE(), "element_mapping", 3);

            jg["field_dim"] = field_dim;
            jg["nb_field_cp_per_cell"] = nb_field_cp_per_cell;
            jg["field_ctrl_points"] = encoded_field;
            jg["interpolation"] = generate_function(fes->GetFE(first_elt_of_group[g]), "element_interpolation", fes->GetVDim());

            ja.push_back(jg);
        }

        json j = json::object();
        j["version"] = "0.1";
        j["groups"] = ja;
        j["nb_conforming_dofs"] = fes->GetNConformingDofs();
        j["largest_edge"] = largest_edge(mesh);

        std::ofstream out(filename);
        out << j.dump(2);
        out.close();
        printf("FieldObject saved to %s\n", filename.c_str());
    }
}

int main(int argc, char* argv[]){

    if(argc != 3 && argc != 4){
        printf("Usage: ./convert-from-mfem path_to_mesh.mesh path_to_gridfunction.gf <field_object.json>\n");
        return 1;
    }
    std::string mesh_file(argv[1]);
    std::string gridfunction_file(argv[2]);
    std::string output("field_object.json");
    if(argc == 4){
        output = std::string(argv[3]);
    }

    mfem::Mesh mesh(mesh_file.c_str(), 0, 0, false);
    std::ifstream solfile(gridfunction_file);
    mfem::GridFunction gf(&mesh, solfile);

    Convert::ExportToFieldObjectFormat(output, &gf);
    return 0;
}
