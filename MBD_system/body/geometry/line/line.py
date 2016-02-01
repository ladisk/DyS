"""
Created on 23. feb. 2015

@author: lskrinjar
"""
from pprint import pprint
from OpenGL import *
from OpenGL.GL import *
from OpenGL.GL.ARB.vertex_buffer_object import *
from OpenGL.GL.shaders import *
from OpenGL.GLU import *
import numpy as np


from MBD_system.extract_from_dictionary_by_string_in_key import extract_from_dictionary_by_string_in_key
from MBD_system.body.geometry.VBO_reshape.VBO_reshape import create_VBO_array

class Line(object):
    """
    classdocs
    """
    def __init__(self, node_i=None, node_j=None, normal=None, length=None, parent=None):
        """
        Constructor
        """
        #   parent
        self._parent = parent

        #   node coordinates
        self.node_i = node_i
        self.node_j = node_j
        self.nodes = [self.node_i, self.node_j]
        self.normal_2D = normal

        #   length
        self.length = length

        #   override parent transparent options
        self._parent.transparent_GL = 1

    def add_attributes_from_dict(self, dict):
        """

        :return:
        """
        for key in dict:
            setattr(self, key, dict[key])

    def _create_line_nodes(self):
        """

        :return:
        """
        if self.node_i is not None and self.node_j is not None:
            nodes = np.array([np.append(self.node_i, self._parent.z_dim),
                                np.append(self.node_j, self._parent.z_dim)], dtype='float32')
        else:
            nodes = np.array([[-self.length/2., 0, self._parent.z_dim],
                                [self.length/2., 0, self._parent.z_dim]], dtype='float32')
        return nodes

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
        self.vertices = self._create_line_nodes()
        #   number of vertices
        self.N_vertices = len(self.vertices)

        #    color vector
        #    check if color vector is a uniform color - every vertex has the same color or is a color vector different for every vertex
        if len(self._parent.color_GL) == 3:
            color_GL = np.array([np.append(self._parent.color_GL, self._parent.transparent_GL)], dtype='float32')
            self.color = np.repeat(color_GL, self.N_vertices, axis=0)
        else:
            raise Error, "Color data shape not equal to vertices."

        #    reshape VBO_data and create interleaved VBO
        self.VBO_data = create_VBO_array(None, self.vertices, None, self.color, GL_primitive_type="line")

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
        self.vbo_id = 1
        glDeleteBuffers(1, int(self.vbo_id))

    def _update_VBO_color(self, color, transparent):
        """
        Update color of VBO object
        """
        #   color array of each vertex
        _color = np.array([np.append(color, transparent)], dtype='float32').flatten()

        #   offset in bits
        _offset = 12

        #   size in bits, to be replaces
        _size = arrays.ArrayDatatype.arrayByteCount(_color)

        glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)
        #   loop through nodes
        for i in range(0, self.N_vertices):
            glBufferSubData(GL_ARRAY_BUFFER, _offset+i*(_offset+16), _size, _color)#16

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
            c_pointer = ctypes.c_void_p(12)
            # c_pointer = ctypes.c_void_p(24)

            #    stride in bits (1 float = 4 bits)
            stride_in_bits = 28

            glDisableClientState(GL_NORMAL_ARRAY)

            glVertexPointer(3, GL_FLOAT, stride_in_bits, v_pointer)
            glColorPointer(4, GL_FLOAT, stride_in_bits, c_pointer)
            glNormalPointer(GL_FLOAT, stride_in_bits, 0)

            glDisable(GL_LIGHTING)
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
            glDrawArrays(GL_LINES, 0, self.N_vertices)
            glEnable(GL_LIGHTING)

            #    draw (fill or wireframe)
            # if self._parent.display_style == "wireframe":
            #     glDisable(GL_LIGHTING)
            #     glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)  # GL_FRONT, GL_FRONT_AND_BACK
            #     glDrawArrays(GL_TRIANGLES, 0, self.N_vertices)
            #     glEnable(GL_LIGHTING)
            # elif self._parent.display_style == "filled":
            #     glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
            #     glDrawArrays(GL_TRIANGLES, 0, self.N_vertices)
            # elif self._parent.display_style == "points":
            #     glDisable(GL_LIGHTING)
            #     glPointSize(0.2)
            #     glPolygonMode(GL_FRONT_AND_BACK, GL_POINT)
            #     glDrawArrays(GL_TRIANGLES, 0, self.N_vertices)
            #     glPointSize(1)
            #     glEnable(GL_LIGHTING)
            # else:
            #     None

            glBindBuffer(GL_ARRAY_BUFFER, 0)

