"""
Created on 7. jan. 2014

@author: lskrinjar
"""

from OpenGL.GL import *


class GLMatrix(object):
    """
    classdocs
    """

    def __init__(self):
        """
        Constructor of object to save opengl modelview and projection matrices for:
        main viewport and
        CS viewport
        in:
            glGetFloatv(GL_MODELVIEW_MATRIX) of main viewport
            glGetDoublev(GL_PROJECTION_MATRIX) of main viewport
            glGetFloatv(GL_MODELVIEW_MATRIX) of CS viewport viewport
            glGetDoublev(GL_PROJECTION_MATRIX) of CS viewport viewport
        """
        self.modelview_matrix = glLoadIdentity()
        self.projection_matrix = glLoadIdentity()
        self.modelview_matrix_CS = glLoadIdentity()
        self.projection_matrix_CS = glLoadIdentity()
        
    def update_modelview_matrix(self, m):
        """
        Function updates attribute value
        :param m: matrix
        :return:
        """
        self.modelview_matrix = m
        
    def update_projection_matrix(self, m):
        """
        Function updates attribute value
        :param m: matrix
        :return:
        """
        self.projection_matrix = m
        
    def update_modelview_matrix_CS(self, m):
        """
        Function updates attribute value
        :param m: matrix
        :return:
        """
        self.modelview_matrix_CS = m
        
    def update_projection_matrix_CS(self, m):
        """
        Function updates attribute value
        :param m: matrix
        :return:
        """
        self.projection_matrix_CS = m
