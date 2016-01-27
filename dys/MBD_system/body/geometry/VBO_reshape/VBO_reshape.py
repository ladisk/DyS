'''
Created on 3. dec. 2013

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
from OpenGL import *
from OpenGL.GL import *
from OpenGL.GLU import *

import numpy as np


def create_VBO_array(file_extension, vertices, normals, colors, GL_primitive_type, interleaved_type):
    """
    method creates interleaved VBO from
    in:
        vertices matrix
        normals matrix
        colors matrix
        GL_primitive_type = triangle or quad
    out:
        interleaved_VBO_array - array that is passed to GPU memory
    """
    #    expand normals matrix, geometry data defines normal per surface-triangle that is defined by 3 points
    if GL_primitive_type == "triangle":
        if file_extension == ".stl":
            if len(vertices) != len(normals):
                normals = np.repeat(normals, repeats=3, axis=0)
                
        elif file_extension == ".obj":
            None
        else:
            None
        
    elif GL_primitive_type == "quad":
        normals = np.repeat(normals, repeats=4, axis=0)
    else:
        raise Error, "GL_primitive_type not correct, check string type"
        None
    

    #    check if array type is interleaved or not
    if (interleaved_type == "true"):
        VBO_array = np.hstack((vertices, normals, colors)).flatten()
    elif (interleaved_type == "false"):
        VBO_array = np.vstack((vertices, normals, colors))
    else:
        None
    
    VBO_array = np.array(VBO_array, dtype='float32') 
    return VBO_array
