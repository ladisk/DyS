"""
Created on 3. dec. 2013

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
from OpenGL import *
from OpenGL.GL import *
from OpenGL.GLU import *

import numpy as np


def create_VBO_array(file_extension, vertices, normals, colors, GL_primitive_type, interleaved=True):
    """
    Method creates interleaved VBO from input data
    :param file_extension:      file extension
    :param vertices:            vertices 3D (numpy array)
    :param normals:             normals 3D (numpy array)
    :param colors:              color if mono-color or color per vertex (numpy array)
    :param GL_primitive_type:   type of GL primitive-geometry (string), supported types:
                                triangle, line, lines(todo), quad(todo)
    :param interleaved:    status of vbo data (boolean)
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
        print "todo"

    elif GL_primitive_type == "line":
        normals = None

    elif GL_primitive_type == "lines":
        print "todo"
    else:
        raise Error, "GL_primitive_type not correct, check string type"
        None

    #    check if array type is interleaved or not
    if interleaved and GL_primitive_type == "triangle":
        VBO_array = np.hstack((vertices, normals, colors)).flatten()
    elif interleaved and GL_primitive_type == "line":
        VBO_array = np.hstack((vertices, colors)).flatten()
    else:
        VBO_array = np.vstack((vertices, normals, colors))

    #   convert data to 32bit float to load to GPU vram
    VBO_array = np.array(VBO_array, dtype='float32') 
    return VBO_array
