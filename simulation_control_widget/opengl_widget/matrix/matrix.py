"""
Created on 7. jan. 2014

@author: lskrinjar
"""

from OpenGL import *
from OpenGL.GL import *
from OpenGL.GLU import *


class Matrix(object):
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
        
    def update_modelview_matrix(self, modelview_matrix):
        self.modelview_matrix = modelview_matrix
        
    def update_projection_matrix(self, projection_matrix):
        self.projection_matrix = projection_matrix
        
    def update_modelview_matrix_CS(self, modelview_matrix_CS):
        self.modelview_matrix_CS = modelview_matrix_CS
        
    def update_projection_matrix_CS(self, projection_matrix_CS):
        self.projection_matrix_CS = projection_matrix_CS