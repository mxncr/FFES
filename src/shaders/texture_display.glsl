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

const vec2 quadUV[4] = 
vec2[](
    vec2(0.0, 0.0), 
    vec2(1.0, 0.0), 
    vec2(0.0, 1.0), 
    vec2(1.0, 1.0) 
    );

const vec2 quadVert[4] = 
vec2[](
    vec2(-0.3, -0.3), 
    vec2( 0.3, -0.3), 
    vec2(-0.3,  0.3), 
    vec2( 0.3,  0.3) 
    );

uniform mat4 MVP;
out vec2 uv;

void main()
{
    gl_Position = MVP*vec4(quadVert[gl_VertexID], 0.0f, 1.0f);
    uv = quadUV[gl_VertexID];
}

-- Fragment

uniform sampler2DRect val;
uniform sampler2D colormap;
uniform int repeat;
uniform float val_min;
uniform float val_max;
uniform int field_coord;
uniform int field_dim;
uniform int texture_width;
uniform int texture_height;
in vec2 uv;
out vec4 color;

void main(void)
{

    ivec2 c = ivec2(int(uv.x * texture_width), int(uv.y * texture_height));
    vec4 v = texelFetch(val, c).rgba;

    /* Display white if texture is not defined */
    if(v.r == NO_DATA){
        color = vec4(1.,1.,1.,1.);
        return;
    } 

    /* Value of the current field */
    float cval = 0.;
    if (field_coord < field_dim) {
        cval = v[field_coord];
    } else {
        for (int i = 0; i < field_dim; ++i){
            cval += v[i]*v[i];
        }
        cval = sqrt(cval);
    }

    /* Normalized value */
    float nval = (cval - val_min) / (val_max - val_min);
#ifdef DEBUG_SHADER
#ifdef DEBUG_SHOW_ELEMENTS
    if (int(cval) % 2 == 0){
        nval = cval / (2.*NB_ELTS);
    } else {
        nval = 1. / 2. + cval / (2.*NB_ELTS);
    }
#endif
#ifdef DEBUG_SHOW_GROUP
    nval = cval;
#endif
#endif
    if(nval < 0.f){
        nval = 0.f;
    } else if (nval > 1.f){
        nval = 1.f;
    }

    /* Normalized value to color */
    if (repeat > 0) {
        float t = float(repeat)+1.;
        nval = t*nval - int(t*nval);
    }
    color = vec4(texture(colormap, vec2(nval, 0.5)).xyz, 1.);

    return;
}
