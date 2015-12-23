'''
Created on 9. jul. 2015

@author: reti-luka
'''
import sys
import os
from pprint import pprint
import itertools
import logging

import numpy as np
from PyQt4 import QtCore, QtGui
from OpenGL import *
from OpenGL.GL import *
from OpenGL.GLU import *


class CoordinateSystem(object):
    """
    classdocs
    """
    __id = itertools.count(-1)

    def __init__(self, parent=None):
        """
        Class constructor of VBO 
        """
        self._parent = parent

        self.scale = 1E-2
        
        #    points
        self.origin = self.scale*np.array([0, 0, 0])
        self.x_axis = self.scale*np.array([1, 0, 0])
        self.y_axis = self.scale*np.array([0, 1, 0])
        self.z_axis = self.scale*np.array([0, 0, 1])

        #    colors
        #    red
        self.x_color = np.array([[1, 0, 0]])
        #    green
        self.y_color = np.array([[0, 1, 0]])
        #    blue
        self.z_color = np.array([[0, 0, 1]])

        #   colors vector for ewery node
        x_colors = np.repeat(self.x_color, repeats=3, axis=0)
        y_colors = np.repeat(self.y_color, repeats=3, axis=0)
        z_colors = np.repeat(self.z_color, repeats=3, axis=0)

        #   axis vector for ever axes
        _x_axis = np.array([self.origin, self.x_axis, self.origin])
        _y_axis = np.array([self.origin, self.y_axis, self.origin])
        _z_axis = np.array([self.origin, self.z_axis, self.origin])

        #   axis matrix
        axis = np.vstack((_x_axis, _y_axis, _z_axis))
        #   number of elements - nodes axis attribute
        self.N_axis = len(axis)

        #   colors matrix
        colors = np.vstack((x_colors, y_colors, z_colors))

        #   flatten vbo array
        _vbo_array = np.hstack((axis, colors)).flatten()

        #   vbo array to use with opengl library
        self.vbo_array = np.array(_vbo_array, dtype='float32')

    def _create_VBO(self):
        """
        Function creates VBO of coordinate system
        """
        self.__v_pointer = ctypes.c_void_p(0)
        self.__c_pointer = ctypes.c_void_p(12)
        
        #    generate a new VBO and get the associated vbo_id
        num_of_VBOs = 1
        
        
        #    create buffer name
        self.vbo_id = GLuint()
        self.vbo_id = glGenBuffers(num_of_VBOs)
        
        #    bind name to buffer
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)
        
        
        #    VBO_data size in bytes
        self.VBO_data_size_in_bytes = arrays.ArrayDatatype.arrayByteCount(self.vbo_array)
        
        
        #    allocate space and upload the data from CPU to GPU
        #    bind data to buffer
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)
        
        
        #    add VBO data to buffer
        glBufferData(GL_ARRAY_BUFFER, self.VBO_data_size_in_bytes, self.vbo_array, GL_STATIC_DRAW)
        glBindBuffer(GL_ARRAY_BUFFER, 0)

    def _paintGL_VBO_CS(self):
        """
        Paint local coordinate system VBO
        """
        #   bindy buffer to id
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)

        #    stride in bits (1 float = 4 bits)
        stride_in_bits = 24

        #   pointers
        glVertexPointer(3, GL_FLOAT, stride_in_bits, self.__v_pointer)
        glColorPointer(3, GL_FLOAT, stride_in_bits, self.__c_pointer)

        glDisable(GL_LIGHTING)

        #   draw lines - opengl
        glDrawArrays(GL_LINES, 0, self.N_axis)

        glEnable(GL_LIGHTING)