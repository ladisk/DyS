#version 330

layout(location = 0) in vec4 position;
in vec3 VertexColor;
uniform mat4 gl_ModelViewMatrix;
uniform mat4 gl_ProjectionMatrix;
uniform vec3 vertex_color;

out vec3 Color;
void main()
{
    gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * position;
    // set color
    Color = vertex_color;
    
}