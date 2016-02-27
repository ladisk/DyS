"""
Created on 30. nov. 2013

@author: lskrinjar
"""
import os

import numpy as np
from OpenGL.GL import *
from OpenGL.GL.ARB.vertex_buffer_object import *
from OpenGL.GLU import *

import VBO_reshape.VBO_reshape as VBO_reshape
import vertices_offset.vertices_offset as vertices_offset
from model_loader.model_loader import ModelLoader


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

        #   created status
        self.VBO_created = False

        #   shift vertices that body center of mass is the origin
        self.geom_data.vertices = vertices_offset.offset(self.geom_data.vertices, parent.CM_CAD_LCS)
        
#         print "self.geom_data.vertices ="
#         print self.geom_data.vertices
#         
#         print len(self.geom_data.vertices)
#         print np.shape(self.geom_data.vertices)
#         triangles = np.arange(0, len(self.geom_data.vertices)).reshape(3, len(self.geom_data.vertices)/3)+0
#         np.savetxt(filename+"triagnles.txt", triangles, "%4i", delimiter=",", newline=";\n")
        
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
        self.VBO_data = VBO_reshape.create_VBO_array(self._file_extension, self.geom_data.get_vertices_3D(), self.geom_data.get_normals_3D(), self.vertex_color, GL_primitive_type="triangle")
        
        #    VBO_data size in bytes
        self.VBO_data_size_in_bytes = arrays.ArrayDatatype.arrayByteCount(self.VBO_data)

        #    allocate space and upload the data from CPU to GPU
        #    bind data to buffer
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)
        #    add VBO data to buffer
        glBufferData(GL_ARRAY_BUFFER, self.VBO_data_size_in_bytes, self.VBO_data, GL_DYNAMIC_DRAW)
        glBindBuffer(GL_ARRAY_BUFFER, 0)

        if glGetError() == GL_NO_ERROR:
            self.VBO_created = True

    def __del__(self):
        """
        Delete opengl object VBO from VRAM
        """
        try:
            glDeleteBuffers(1, int(self.vbo_id))
        except:
            pass

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

    def _paint_VBO(self):
        """

        :return:
        """
        if self.VBO_created:
            #   bind to buffer
            glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)

            #    pointer - offset in bits
            v_pointer = ctypes.c_void_p(0)
            n_pointer = ctypes.c_void_p(12)
            c_pointer = ctypes.c_void_p(24)

            #    stride in bits (1 float = 4 bits)
            stride_in_bits = 40

            glVertexPointer(3, GL_FLOAT, stride_in_bits, v_pointer)
            glColorPointer(4, GL_FLOAT, stride_in_bits, c_pointer)
            glEnableClientState(GL_NORMAL_ARRAY)
            glNormalPointer(GL_FLOAT, stride_in_bits, n_pointer)

            #    draw by type of visualization
            if self._parent.display_style == "wireframe":
                glDisable(GL_LIGHTING)
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)  # GL_FRONT, GL_FRONT_AND_BACK
                glDrawArrays(GL_TRIANGLES, 0, self.N_vertices)
                glEnable(GL_LIGHTING)
                
            elif self._parent.display_style == "filled":
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
                glDrawArrays(GL_TRIANGLES, 0, self.N_vertices)
                
            elif self._parent.display_style == "points":
                glDisable(GL_LIGHTING)
                glPointSize(0.2)
                glPolygonMode(GL_FRONT_AND_BACK, GL_POINT)
                glDrawArrays(GL_TRIANGLES, 0, self.N_vertices)
                glPointSize(1)
                glEnable(GL_LIGHTING)
            else:
                None

            glBindBuffer(GL_ARRAY_BUFFER, 0)
            