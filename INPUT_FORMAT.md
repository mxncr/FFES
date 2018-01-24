As of version 0.1, the ffes software only accepts inputs in a .json format.

# Informal description

IMPORTANT: The mapping and interpolation functions are defined in the input file ! This
allows easy support for arbitrary functions without modifications of the
program source code.  More importantly, it avoids mismatches between the finite
element code export and this program import (e.g. different DOF numbering or
positioning for high-orders FEM).

The json object contains an entry `groups` which is an array of ElementGroup.
Each ElementGroup is a chunk of finite elements with the same
mapping and interpolation functions. The mesh and field coefficients are stored
in strings using Base64 encoding (see below). The mapping and the interpolation functions are
also stored in strings.

For instance:
- A tetrahedral mesh with P1 elements has only one ElementGroup.
- A hex-dominant meshes with tetrahedra and hexahedra has two ElementGroup.
- A tetrahedral mesh with P1, P2, P3, P4 and P5 finite element has five ElementGroup.
- A tetrahedral mesh with affine geometry (P1), quadratic geometry (P2) and cubic interpolation (P3) has two ElementGroup (P1-P3 and P2-P3).


# FieldObject .json version v0.1

The input file is loaded with the function field_load() defined in the file src/ffes/field_io.cpp.
A short version of the load function is given at the end of this file for reference
A example of writer function is available in src/bin/convert-from-mfem.cpp.

Specification (key ordering has no importance):

    {
        "version": "0.1",
        "groups": [
            {
                "mesh_dim": <int: dimension of the mesh, should be 3>,
                "primitive": <string: TET or HEX, other decompositions not implemented by default>,
                "nb_mesh_cp_per_cell": <int: nb of control points per cell in the mesh>,
                "mapping": <string: mapping function>,
                "mesh_ctrl_points": <string: Base64 encoding of mesh coefficients>,
                "field_dim": <int: dimension of the field, e.g. 1 for scalar field or 3 for displacement>,
                "nb_field_cp_per_cell": <int: nb of control points per cell in the field, e.g. 4 for linear tetrahedra (even if field_dim > 1)>,
                "interpolation": <string: interpolation function>,
                "field_ctrl_points": <string: Base64 encoding of field coefficients>
            },
            {
                // another group definition if multiple ElementGroup
            }
        ],
        "arbitrary_key": <string: other metadata can be stored in this json, they will be copied to the output>
    }

Example:

    {
        "version": "0.1",
        "groups": [
            {
                "primitive": "TET",
                "field_ctrl_points": "AAA..not included because too long..weqeqwezdfar",
                "field_dim": 1,
                "interpolation": "float element_interpolation(in vec3 ref_pos, in float values[4]){ \n    ct result = 0.;\n    const ct u1 = 1. - ref_pos.x - ref_pos.y - ref_pos.z;                                 \n    const ct v1 = ref_pos.x;                                                              \n    const ct w1 = ref_pos.y;                                                              \n    const ct z1 = ref_pos.z;                                                              \n    ct r = 0.; \n    r = 0. + 1*u1;\n    result += r * values[0];\n    r = 0. + 1*v1;\n    result += r * values[1];\n    r = 0. + 1*w1;\n    result += r * values[2];\n    r = 0. + 1*z1;\n    result += r * values[3];\n    return float(result);\n}\n",
                "mapping": "vec3 element_mapping(in vec3 ref_pos, in vec3 values[4]){ \n    ct3 result = ct3(0.,0.,0.);\n    const ct u1 = 1. - ref_pos.x - ref_pos.y - ref_pos.z;                                 \n    const ct v1 = ref_pos.x;                                                              \n    const ct w1 = ref_pos.y;                                                              \n    const ct z1 = ref_pos.z;                                                              \n    ct r = 0.; \n    r = 0. + 1*u1;\n    result += r * values[0];\n    r = 0. + 1*v1;\n    result += r * values[1];\n    r = 0. + 1*w1;\n    result += r * values[2];\n    r = 0. + 1*z1;\n    result += r * values[3];\n    return vec3(result);\n}\n",
                "mesh_ctrl_points": "lnc5Pyjzf...",
                "mesh_dim": 3,
                "nb_field_cp_per_cell": 4,
                "nb_mesh_cp_per_cell": 4
            }
        ],
        "largest_edge": 0.485146494845629,
        "nb_conforming_dofs": 118,
        "time_assembly": "0.003329",
        "time_solver": "0.000402"
    }



# Base64 encoding of coefficients

Base64 is a technique used to store binary data into ASCII strings (with a small size loss).
There are various ways to achieve this. Please use the encoding/decoding implementation which
is defined in the very simple file `./third_party/base64/base64.hpp`. If you prefer to use your
Base64 encoding functions, you have to verify they are compatible with our decoding.

Encoded data is expected to be **32bits floats** !

The mesh and field coefficients are defined cell per cell, as if the fields
were always discontinuous.  This allows easy implementation of mapping and
interpolation and fast processing by the CPU (contiguous data). For a mesh with
*nb_elems* elements, there are *nb_elems * nb_mesh_cp_per_cell * mesh_dim*
coefficients in the vector of floats associated to *mesh_ctrl_points*. Same
for the fields.

- Examples:

Encoding of a vector of floats to a string:

    std::vector<float> values(nb_values);
    /* fill the vector */
    std::string encoded_values = Base64::encode((unsigned char*) &values[0],values.size()*sizeof(float));

Decoding of a string to a vector of floats:

    std::string encoded = /* string value extracted from the json file */;
    size_t nb_chars = Base64::decode_len(encoded);
    std::vector<float> data((nb_chars + 3)/4);
    Base64::decode(encoded, (unsigned char*) &data[0]);
    data.resize(nb_chars / 4);


# Mapping and interpolation functions

### Mappings

The string value associated to the key "mapping" should be a GLSL function implementing
a function with the signature:

    vec3 element_mapping(const vec3 ref_pos, const vec3 values[ <nb_mesh_cp_per_cell> ]){ .. }

The string value can contains "\n" characters for easier debugging / printing.

Examples:

- affine tetrahedron mapping:


    vec3 element_mapping(const vec3 ref_pos, const vec3 values[4]){
        vec3 result =
            (1.0f-ref_pos.x-ref_pos.y-ref_pos.z) * values[0]
            + ref_pos.x                          * values[1]
            + ref_pos.y                          * values[2]
            + ref_pos.z                          * values[3];
        return result;
    }

resulting in a json entry:

    {
        "groups" : [
            {
                "mapping" = "vec3 element_mapping(const vec3 ref_pos, const vec3 values[4]){ vec3 result = (1.0f-ref_pos.x-ref_pos.y-ref_pos.z) * values[0] + ref_pos.x * values[1] + ref_pos.y * values[2] + ref_pos.z * values[3]; return result; }",
                ..
            }
        ]
        ..
    }


- quadratic tetrahedron mapping:


    vec3 element_mapping(const vec3 ref_pos, const vec3 values[10]){
        float L0 = 1. - ref_pos.x - ref_pos.y - ref_pos.z;
        float L1 = ref_pos.x;
        float L2 = ref_pos.y;
        float L3 = ref_pos.z;
        vec3 result =
            L0 * ( 2.0 * L0 - 1.0 )   * values[0]
            + L1 * ( 2.0 * L1 - 1.0 ) * values[1]
            + L2 * ( 2.0 * L2 - 1.0 ) * values[2]
            + L3 * ( 2.0 * L3 - 1.0 ) * values[3]
            + 4.0 * L0 * L1           * values[4]
            + 4.0 * L0 * L2           * values[5]
            + 4.0 * L0 * L3           * values[6]
            + 4.0 * L1 * L2           * values[7]
            + 4.0 * L1 * L3           * values[8]
            + 4.0 * L2 * L3           * values[9];
        return result;
    }

- trilinear hexahedron mapping:


    vec3 element_mapping(const vec3 ref_pos, const vec3 values[8]){
        float x = ref_pos.x, y = ref_pos.y, z = ref_pos.z;
        float ox = 1.-x, oy = 1.-y, oz = 1.-z;
        vec3 result =
            (ox * oy * oz * values[0]
             +   x * oy * oz * values[1]
             +   x *  y * oz * values[2]
             +  ox *  y * oz * values[3]
             +  ox * oy *  z * values[4]
             +   x * oy *  z * values[5]
             +   x *  y *  z * values[6]
             +  ox *  y *  z * values[7]);
        return result;
    }


### Interpolation

The string value associated to the key "interpolation" should be a GLSL function implementing
a function with the signature:

    <field_type> element_interpolation(const vec3 ref_pos, const <field_type> values[ <nb_field_cp_per_cell> ]){ .. }

with <field_type> = **float** for scalar fields,**vec3** for vector fields of dimension 3, **vec4** for vector fields of dimension 4

The string value can contains "\n" characters for easier debugging / printing.

Examples:

- scalar linear interpolation in tetrahedron:


    float element_interpolation(const vec3 ref_pos, const float values[4]){
        float result =
            (1.0f-ref_pos.x-ref_pos.y-ref_pos.z) * values[0]
            + ref_pos.x                          * values[1]
            + ref_pos.y                          * values[2]
            + ref_pos.z                          * values[3];
        return result;
    }

- vector (e.g. displacement) trilinear interpolation in hexahedron:


    vec3 element_interpolation(const vec3 ref_pos, const vec3 values[8]){
        float x = ref_pos.x, y = ref_pos.y, z = ref_pos.z;
        float ox = 1.-x, oy = 1.-y, oz = 1.-z;
        vec3 result =
            (ox * oy * oz * values[0]
             +   x * oy * oz * values[1]
             +   x *  y * oz * values[2]
             +  ox *  y * oz * values[3]
             +  ox * oy *  z * values[4]
             +   x * oy *  z * values[5]
             +   x *  y *  z * values[6]
             +  ox *  y *  z * values[7]);
        return result;
    }


### Remarks

The value associated to the key "mapping" is paste at the beginning of the vertex shader. 
Any GLSL valid code can be used, including multiple function definitions, variables... The
only requirement is that the **element_mapping** signature is implemented.

The value associated to the key "interpolation" is paste at the begging of the fragment shader.
Any valid GLSL code can be used, the only requirement is that the **element_interpolation** is
implemented and that its <field_type> is coherent with the <field_dim> value.

For more advanced functions, see the code of the executable
src/bin/convert-from-mfem.cpp which generates the functions for Lagrange finite
elements of arbitrary order (hex or tet).

# Appendix

### field_load(): short version of the code


    void field_load_json(const std::string& path, FieldObject& obj){
        std::ifstream input(path, std::ios::binary | std::ios::in);
        json jin;
        input >> jin;
        ...
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
        ...
    }

