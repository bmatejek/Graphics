/***********************************************************
*   Simple GLSL demo code for COS426 Computer Graphics     *
*                           ****                           *
*   Most is borrowed from                                  *
*   http://www.lighthouse3d.com/opengl/glsl/               *
*                           ****                           *
*   Dependencies: glew, glut, gl                           *
*                           ****                           *
*   Pulled together by Aleksey Boyko(aboyko@princeton.edu) *
************************************************************/

varying vec3 normal;
varying vec3 vertex_light_position1;
varying vec3 vertex_light_half_vector1;
varying vec3 vertex_light_position2;
varying vec3 vertex_light_half_vector2;
uniform sampler2D color_texture;
uniform int texton;

varying vec4 lightDir;

uniform vec4 _WorldSpaceLightPos0;
uniform vec3 _WorldSpaceCameraPos; // camera position in world space
uniform mat4 _Object2World; // model matrix	
	
float luminosity(vec4 color)
{
    vec4 lumcoeff = vec4(0.299,0.587,0.114,0.0);
    return dot(color, lumcoeff);
}

void main()
{
	
	//float intensity;
	vec4 intensity;
	vec4 color;
//	float color;
	vec3 n = normalize(normal);
		

//vec3 viewDirection = normalize(_WorldSpaceCameraPos - vec3(position));



       

	  // Calculatethe ambient term
	  vec4 ambient_color = gl_FrontMaterial.ambient * gl_LightSource[0].ambient + gl_LightModel.ambient * gl_FrontMaterial.ambient;

    // Calculate the diffuse term
    vec4 diffuse_color = gl_FrontMaterial.diffuse * gl_LightSource[0].diffuse;

    // Calculate the specular value
    vec4 specular_color = gl_FrontMaterial.specular * gl_LightSource[0].specular * pow(max(dot(normal, vertex_light_half_vector1), 0.0) , gl_FrontMaterial.shininess);

    // Set the diffuse value (darkness). This is done with a dot product between the normal and the light
    // and the maths behind it is explained in the maths section of the site.
    float diffuse_value = max(dot(normal, vertex_light_position1), 0.0);



	  ambient_color += gl_FrontMaterial.ambient * gl_LightSource[1].ambient;

    // Calculate the diffuse term
     diffuse_color += gl_FrontMaterial.diffuse * gl_LightSource[1].diffuse;

    // Calculate the specular value
    specular_color += gl_FrontMaterial.specular * gl_LightSource[1].specular * pow(max(dot(normal, vertex_light_half_vector2), 0.0) , gl_FrontMaterial.shininess);

    // Set the diffuse value (darkness). This is done with a dot product between the normal and the light
    // and the maths behind it is explained in the maths section of the site.
    diffuse_value += max(dot(normal, vertex_light_position2), 0.0);

    // Set the output color of our current pixel
//    gl_FragColor = ambient_color + diffuse_color * diffuse_value + specular_color;
      //vec4 total = ambient_color + diffuse_color * diffuse_value + specular_color;

      intensity = ambient_color + diffuse_color * diffuse_value + specular_color;
      //total

      if (texton > 0)
            intensity *= texture2D(color_texture, gl_TexCoord[0].st);

      //intensity = luminosity(total);	

// NEED TO DO IT SEPARATELY ON EACH COLOR CHANNEL

//	intensity = dot(vec3(gl_LightSource[0].position),n);
//	intensity -= dot(vec3(gl_LightSource[1].position),n);	

//	color = vec4(1.0, 1.0, 1.0, 1.0);
	color = vec4(0.0, 0.0, 0.0, 1.0);
	
	for (int i = 0; i < 3; i++){ 

	if (intensity[i] > 0.65)  
		//color = vec4(1.0,1.0,1.0,1.0);
		color[i] = 1.25;
	else if (intensity[i] > 0.3)
		//color = vec4(0.6,0.6,0.6,1.0);
		color[i] = 1.0;
	else if (intensity[i] > 0.15)
		//color = vec4(0.4,0.4,0.4,1.0);
		color[i] = .75;
	else 
		//color = vec4(0.2,0.2,0.2,1.0);
		color[i] = 0.0;
	}
	
//	gl_FragColor = color*gl_Color;
//	gl_FragColor = total*color;
	gl_FragColor = vec4(intensity[0]*color[0], intensity[1]*color[1], intensity[2]*color[2], 1.0);


/*
float _UnlitOutlineThickness = 100000.0;
float _LitOutlineThickness =   50000.0;

gl_FragColor = vec4(1.0,1.0,1.0,1.0);

// higher priority: outline
    if (dot(viewDirection, normal)> 0.0 || dot(viewDirection, normal)<= 0.0)
      gl_FragColor = vec4(0.0,1.0,0.0,1.0);
   
   

        if (dot(viewDirection, normal)
               < mix(_UnlitOutlineThickness, _LitOutlineThickness, 
               max(0.0, dot(normal, normalize(vec3(_WorldSpaceLightPos0))))))
            {
               gl_FragColor = vec4(1.0,0.0,0.0,1.0); 
            }
*/
//texture2D(color_texture, gl_TexCoord[0].st);
} 