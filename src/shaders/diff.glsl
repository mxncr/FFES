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

const vec2 quadVertices[4] = 
vec2[](
    vec2(-1.0f, -1.0f), 
    vec2( 1.0f, -1.0f),
    vec2(-1.0f,  1.0f), 
    vec2( 1.0f,  1.0f)
    );

void main()
{
    gl_Position = vec4(quadVertices[gl_VertexID], 0.0, 1.0);
}

-- Fragment
uniform sampler2DRect fieldA;
uniform sampler2DRect fieldB;

#if FIELD_DIM == 1
layout(location = 0) out float field_diff;
void main(void)
{
    ivec2 coord = ivec2(gl_FragCoord.x, gl_FragCoord.y);
    float val_a = texelFetch(fieldA, coord).r;
    float val_b = texelFetch(fieldB, coord).r;
    if(val_a != NO_DATA && val_b != NO_DATA){
        field_diff = abs(val_a - val_b);
    } else {
        field_diff = NO_DATA;
    }
}
#elif FIELD_DIM == 2
layout(location = 0) out vec2 field_diff;
void main(void)
{
    ivec2 coord = ivec2(gl_FragCoord.x, gl_FragCoord.y);
    vec2 val_a = texelFetch(fieldA, coord).rg;
    vec2 val_b = texelFetch(fieldB, coord).rg;
    if(val_a.x != NO_DATA && val_b.x != NO_DATA){
        field_diff = abs(val_a - val_b);
    } else {
        field_diff = vec2(NO_DATA, NO_DATA);
    }
}
#elif FIELD_DIM == 3
layout(location = 0) out vec3 field_diff;
void main(void)
{
    ivec2 coord = ivec2(gl_FragCoord.x, gl_FragCoord.y);
    vec3 val_a = texelFetch(fieldA, coord).rgb;
    vec3 val_b = texelFetch(fieldB, coord).rgb;
    if(val_a.x != NO_DATA && val_b.x != NO_DATA){
        field_diff = abs(val_a - val_b);
    } else {
        field_diff = vec3(NO_DATA, NO_DATA, NO_DATA);
    }
}
#elif FIELD_DIM == 4
layout(location = 0) out vec4 field_diff;
void main(void)
{
    ivec2 coord = ivec2(gl_FragCoord.x, gl_FragCoord.y);
    vec4 val_a = texelFetch(fieldA, coord).rgba;
    vec4 val_b = texelFetch(fieldB, coord).rgba;
    if(val_a.x != NO_DATA && val_b.x != NO_DATA){
        field_diff = abs(val_a - val_b);
    } else {
        field_diff = vec4(NO_DATA, NO_DATA, NO_DATA, NO_DATA);
    }
}
#endif

