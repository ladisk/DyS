'''
Created on 30. nov. 2013

@author: lskrinjar
'''
import os
import ctypes
from pprint import pprint
from OpenGL import *
from OpenGL.GL import *
from OpenGL.GL.ARB.vertex_buffer_object import *
from OpenGL.GL.shaders import *
from OpenGL.GLU import *


import VBO_reshape.VBO_reshape as VBO_reshape
from MBD_system.body.geometry.model_loader.model_loader import ModelLoader
import numpy as np
import vertices_offset.vertices_offset as vertices_offset


class Geometry(object):
    """
    classdocs
    """
    def __init__(self, filename, parent=None):
        """
        Constructor
        Class constructor for OpenGL VBO
        in:
            CAD_LCS_GCS - array of center of body mass in GCS
            geom_file - a body geometry full filename (filename = name + extension)
            geom_file_extension - extension 
            color_GL - a normalized array of RGB color vector for the body
            transparent_GL - a transparency factor of a body between 0 and 1
        out:
            VBO interleaved array of data saved to GPU memory
        """
        self._parent = parent

        self.filename = filename

        self._filename, self._file_extension = os.path.splitext(filename)

        #    read geom data from stl (geometry) geom_file with model_loader
        if os.path.isfile(filename):
            self.geom_data = ModelLoader(filename)

        else:
            raise IOError, "File not found!"

        #   shift vertices that body center of mass is the origin
        self.geom_data.vertices = vertices_offset.offset(self.geom_data.vertices, parent.CM_CAD_LCS)

    def create_VBO(self):
        """
        Method creates a Vertex Buffer Object and writes it to VRAM on GPU
        """
        #    generate a new VBO and get the associated vbo_id
        num_of_VBOs = 1
        
        #    create buffer name
        self.vbo_id = GLuint()
        self.vbo_id = glGenBuffers(num_of_VBOs)

        #    bind name to buffer
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)
        #    test if is buffer
        if glIsBuffer(self.vbo_id) == 1:
            self.VBO_created = True
        else:
            self.VBO_created = False
            raise Error, "VBO not created!"

        #    number of vertices
        self.N_vertices = len(self.geom_data.vertices)

        #    color vector
        #    check if color vector is a uniform color - every vertex has the same color or is a color vector different for every vertex
        if len(self._parent.color_GL) == 3:
            color_GL = np.array([np.append(self._parent.color_GL, self._parent.transparent_GL)], dtype='float32')
            self.vertex_color = np.repeat(color_GL, self.N_vertices, axis=0)
        elif color_GL.shape == color_GL.shape:
            None
        else:
            raise Error, "Color data shape not equal to vertices."
        
        #    reshape VBO_data and create interleaved VBO
        self.VBO_data = VBO_reshape.create_VBO_array(self._file_extension, self.geom_data.get_vertices_3D(), self.geom_data.get_normals_3D(), self.vertex_color, GL_primitive_type="triangle", interleaved_type="true")
        
        #    VBO_data size in bytes
        self.VBO_data_size_in_bytes = arrays.ArrayDatatype.arrayByteCount(self.VBO_data)

        #    allocate space and upload the data from CPU to GPU
        #    bind data to buffer
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)
        #    add VBO data to buffer
        glBufferData(GL_ARRAY_BUFFER, self.VBO_data_size_in_bytes, self.VBO_data, GL_DYNAMIC_DRAW)
        glBindBuffer(GL_ARRAY_BUFFER, 0)

    def __del__(self):
        """
        Delete opengl object VBO from VRAM
        """
        glDeleteBuffers(1, int(self.vbo_id))

    def _update_VBO_color(self, color, transparent):
        """
        Update color of VBO object
        """
        #   color array of each vertex
        _color = np.array([np.append(color, transparent)], dtype='float32').flatten()

        #   offset in bits
        _offset = 24

        #   size in bits, to be replaces
        _size = arrays.ArrayDatatype.arrayByteCount(_color)

        glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)
        #   loop through nodes
        for i in range(0, self.N_vertices):
            glBufferSubData(GL_ARRAY_BUFFER, _offset+i*(_offset+16), _size, _color)
        
        glBindBuffer(GL_ARRAY_BUFFER, 0)