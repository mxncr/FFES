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


// uniform float slice_pos;
// uniform float pxl_width;

uniform uint slice_no;
uniform uint nb_slices;

/* Larger than the screen */
// const vec2 lineVertices[2] = 
// vec2[](
//     vec2( 0.0f,-1.1f), 
//     vec2( 0.0f, 1.1f)
//     );
const vec2 lineVertices[4] = 
vec2[](
    vec2(-1.0f, 1.1f), 
    vec2(-1.0f,-1.1f),
    vec2( 1.0f, 1.1f), 
    vec2( 1.0f,-1.1f)
    );

void main()
{

    // slice_pos range:  ]-1, 1[
    double slice_pos = 2. * ((double(slice_no) + 0.5)/double(nb_slices)) - 1.;
    double pxl_width = 1.99 / double(nb_slices);

    float left_pos  = float(slice_pos - pxl_width / 2.);
    float right_pos = float(slice_pos + pxl_width / 2.);
    // left_pos  = 0.99;
    // right_pos = 0.999;

    vec4 out_pos;
    if (gl_VertexID == 0 || gl_VertexID == 1) {
        out_pos = vec4(
                left_pos,
                lineVertices[gl_VertexID].y, 
                0.0, 
                1.0);
    } else {
        out_pos = vec4(
                right_pos,
                lineVertices[gl_VertexID].y, 
                0.0, 
                1.0);
    }
    // out_pos = vec4(
    //         lineVertices[gl_VertexID].x,
    //         lineVertices[gl_VertexID].y, 
    //         0.0, 
    //         1.0);
    gl_Position = out_pos;

    // These positions determine which fragments will be
    // processed.
}

-- Fragment
uniform sampler2DRect texture_values;
uniform int width;

layout(location = 0) out vec4 reduced;

void main(void)
{
    float sum_L1 = 0;
    float sum_L2 = 0;
    float max_Li = 0;
    uint nnz = 0;
    for(int c = 0; c < width; ++c){
        ivec2 coord = ivec2(c, gl_FragCoord.y);
        vec4 txl = texelFetch(texture_values, coord);
        if(txl.r != NO_DATA){
#if FIELD_DIM == 1
            float d = txl.r;
#elif FIELD_DIM == 2
            float d = length(txl.rg);
#elif FIELD_DIM == 3
            float d = length(txl.rgb);
#elif FIELD_DIM == 4
            float d = length(txl.rgba);
#endif
            sum_L1 += abs(d);
            sum_L2 += abs(d)*abs(d);
            max_Li = max(abs(d), max_Li); 
            nnz +=1 ;
        }
    }
    reduced = vec4(float(nnz), sum_L1, sum_L2, max_Li);
    // reduced = vec4(1, sum_L1, sum_L2, max_Li); // only for debugging synchronisation issies
}
