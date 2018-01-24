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

-- Vertex
precision highp float;
precision highp int;

layout(location = 0) in vec3 reference_coordinates;
layout(location = 1) in vec3 control_points[NB_CP];

uniform dmat4 MVP;

out VertexData
{
    vec3 refCoords;
    vec3 world_pos;
    int instance_id; // this is the ID relative to the first instance in the draw call
} outVS;

// This function will be replaced at the program creation
// function declaration: vec3 element_mapping(const vec3 ref_pos, const vec3 pts[NB_CP]);
// Exemple:
//    vec3 element_mapping(const vec3 ref_pos, const vec3 pts[4]){
//        return (1.0f-ref_pos.x-ref_pos.y-ref_pos.z)*pts[0] 
//            + ref_pos.x * pts[1] + ref_pos.y * pts[2]      
//            + ref_pos.z * pts[3];                              }
<mapping function>



void main(void)
{
    // Mapping of the reference element to the world space
    gl_Position = vec4(0.,0.,0.,1.); // TODO: useless line ?
    outVS.refCoords = reference_coordinates;
    outVS.world_pos = element_mapping(reference_coordinates, control_points);
    outVS.instance_id = gl_InstanceID;
}

-- Geometry

layout (lines_adjacency) in;
layout (triangle_strip, max_vertices = 4) out;

in VertexData
{
    vec3 refCoords;
    vec3 world_pos;
    int instance_id;
} inGS[];

out GeometryData
{
    vec3 refCoords;
    flat int instance_id;
} outGS;

uniform float z_slice;
uniform dmat4 MVP;

// Marching tetra configurations
const int edge_intersection[16][7] = {
    //nb, e01, e02, e03, e12, e13, e23     p3, p2, p1, p0
    {0, 0, 0, 0, 0, 0, 0},             //   0   0   0   0 
    {3, 1, 1, 1, 0, 0, 0},             //   0   0   0   1
    {3, 1, 0, 0, 1, 1, 0},             //   0   0   1   0
    {4, 0, 1, 1, 1, 1, 0},             //   0   0   1   1
    {3, 0, 1, 0, 1, 0, 1},             //   0   1   0   0
    {4, 1, 0, 1, 1, 0, 1},             //   0   1   0   1
    {4, 1, 1, 0, 0, 1, 1},             //   0   1   1   0
    {3, 0, 0, 1, 0, 1, 1},             //   0   1   1   1
    {3, 0, 0, 1, 0, 1, 1},             //   1   0   0   0
    {4, 1, 1, 0, 0, 1, 1},             //   1   0   0   1
    {4, 1, 0, 1, 1, 0, 1},             //   1   0   1   0
    {3, 0, 1, 0, 1, 0, 1},             //   1   0   1   1
    {4, 0, 1, 1, 1, 1, 0},             //   1   1   0   0
    {3, 1, 0, 0, 1, 1, 0},             //   1   1   0   1
    {3, 1, 1, 1, 0, 0, 0},             //   1   1   1   0
    {0, 0, 0, 0, 0, 0, 0}              //   1   1   1   1 
};
const int edge_to_vertex[6][2] = {
    // input : e01, e02, e03, e12, e13, e23  
    // output : {p_i, p_j}
    // see edge_intersection
    {0, 1},
    {0, 2},
    {0, 3},
    {1, 2},
    {1, 3},
    {2, 3}
};

void main(void)
{
    // 1. Check if vertices are over z_slice
    //    Determine the configuration 
    int state[4] = {0,0,0,0};
    if(inGS[0].world_pos.z > z_slice) state[0] = 1;
    if(inGS[1].world_pos.z > z_slice) state[1] = 1;
    if(inGS[2].world_pos.z > z_slice) state[2] = 1;
    if(inGS[3].world_pos.z > z_slice) state[3] = 1;
    int config_no = state[0] * 1 + state[1] * 2 + state[2] * 4 + state[3] * 8;

    if(edge_intersection[config_no][0] == 0) return;

    // 2. Look for vertex indices of intersected edges
    int edge_vertex_indices[6*2]; // worst case: 2 triangles (6 edges)
    int intersection_no = 0;
    for(int j = 1; j <= 6; ++j){ // fixed version
        if(edge_intersection[config_no][j] != 1) continue; // this edge is not intersected

        edge_vertex_indices[2*intersection_no+0] = edge_to_vertex[j-1][0];
        edge_vertex_indices[2*intersection_no+1] = edge_to_vertex[j-1][1];
        intersection_no += 1;
    }

    // 3. Compute triangle coordinates and interpolate parametric coordinates
    for(int e = 0; e < intersection_no; ++e){
        int v1 = edge_vertex_indices[2*e+0];
        int v2 = edge_vertex_indices[2*e+1];
        float z1 = inGS[v1].world_pos.z; 
        float z2 = inGS[v2].world_pos.z; 
        float lambda = (z1 != z2) ? (z_slice-z1)/(z2-z1) : 1.f;

        outGS.refCoords = mix(inGS[v1].refCoords, inGS[v2].refCoords, lambda);
        vec3 wp = mix(inGS[v1].world_pos, inGS[v2].world_pos, lambda);
        gl_Position = vec4(MVP*dvec4(wp, 1.));
        outGS.instance_id = inGS[0].instance_id;
        EmitVertex();
    }
    EndPrimitive();
}

-- Fragment

// This function will be replaced at the program creation
// function declaration: float/vec2/vec3/vec4 element_interpolation(const vec3 ref_pos, 
//                                         const float/vec2/vec3/vec4 values[NB_CV]){ }
// the type has to be chosen accordingly to FIELD_DIM
<interpolation function>

in GeometryData
{
    vec3 refCoords;
    flat int instance_id;
} inFS;

uniform uint first_instance;

layout(std430, binding=0) buffer Coefficients { float coefs[]; };
#if FIELD_DIM == 1
layout(location = 0) out float field_values;
#elif FIELD_DIM == 2
layout(location = 0) out vec2 field_values;
#elif FIELD_DIM == 3
layout(location = 0) out vec3 field_values;
#elif FIELD_DIM == 4
layout(location = 0) out vec4 field_values;
#endif

void main(void)
{
#if FIELD_DIM == 1
    float elt_coefs[NB_CV];
#elif FIELD_DIM == 2
    vec2 elt_coefs[NB_CV];
#elif FIELD_DIM == 3
    vec3 elt_coefs[NB_CV];
#elif FIELD_DIM == 4
    vec4 elt_coefs[NB_CV];
#endif

#if FIELD_DIM == 1
    for(int i = 0; i < NB_CV; ++i){
        elt_coefs[i] = coefs[NB_CV*(first_instance+inFS.instance_id)+i];
    }
#else
    for(int i = 0; i < NB_CV; ++i){
        for(int d = 0; d < FIELD_DIM; ++d){
            elt_coefs[i][d] = coefs[FIELD_DIM*(NB_CV*(first_instance+inFS.instance_id)+i)+d];
        }
    }
#endif
    field_values = element_interpolation(inFS.refCoords, elt_coefs);

#ifdef DEBUG_SHADER
    float dbg_val = -999.;
 #ifdef DEBUG_SHOW_ELEMENTS
    dbg_val = first_instance+inFS.instance_id;
 #endif
 #ifdef DEBUG_SHOW_GROUP
    dbg_val = GROUP;
 #endif
 #ifdef DEBUG_SHOW_FACES
     float seuil = 0.005;
     if (NB_CP == 8) { // Hexahedron
         vec3 v1 = inFS.refCoords.xyz;
         vec3 v2 = vec3(1.-inFS.refCoords.x, 1.-inFS.refCoords.y, 1.-inFS.refCoords.z); 
         if(any(lessThan(v1, vec3(seuil)))){
             dbg_val = 999999;
         }
         if(any(lessThan(v2, vec3(seuil)))){
             dbg_val = 999999;
         }
     } else { // Tetrahedron
         vec4 vBC = vec4(1. - inFS.refCoords.x - inFS.refCoords.y - inFS.refCoords.z, inFS.refCoords.xyz);
         if(any(lessThan(vBC, vec4(seuil)))){
             dbg_val = 999999;
         }
     }
 #endif
 #if FIELD_DIM == 1
     if (dbg_val != -999.){
         field_values = dbg_val;
     }
 #else
     if (dbg_val != -999.){
         field_values.x = dbg_val;
     }
 #endif
#endif
}
