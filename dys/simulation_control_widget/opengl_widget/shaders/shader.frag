#version 330
in vec3 Color;

void main()
{
	// set color
	// use vertex color
	gl_FragColor = vec4(Color, 1.0f);
}