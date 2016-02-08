'''
Created on 6. jan. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

from OpenGL import *
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtOpenGL import *

from simulation_control_widget.opengl_widget.paint_text.paint_text import update_text_color
# from ...paint_text_color_update import update_text_color


#    packages
def display(self, gl_matrix, CS_window_width, CS_window_height):
    """
    
    """
    #    update GL_PROJECTION_MATRIX
    #    in display CS viewport opengl projection matrix is always glLoadIdentity()
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()

    #    update GL_MODELVIEW_MATRIX
    glMatrixMode(GL_MODELVIEW)
    glLoadMatrixd(gl_matrix.modelview_matrix_CS)
    
    glViewport(0, 0, int(CS_window_width), int(CS_window_height))
    
    glClear(GL_DEPTH_BUFFER_BIT)
    glDisable(GL_LIGHTING)
    
    #    line width
    glLineWidth(1.2)
    axis_length = 0.8

    
    glBegin(GL_LINES)
    # X axis - RED
    glColor3f(1, 0, 0)
    glVertex3f(0, 0, 0)
    glVertex3f(axis_length, 0, 0)
    
    # Y axis - GREEN
    glColor3f(0, 1, 0)
    glVertex3f(0, 0, 0)
    glVertex3f(0, axis_length, 0)

    # Z axis - BLUE
    glColor3f(0, 0, 1)
    
    glVertex3f(0, 0, 0)
    glVertex3f(0, 0, axis_length)
    glEnd()
    
    #    CS text
    #    change color
    update_text_color.update(self)
        
    position = axis_length + 0.05
    self.renderText(position, 0, 0, str("x"), font=QFont("Consolas", 10))
    self.renderText(0, position, 0, str("y"), font=QFont("Consolas", 10))
    self.renderText(0, 0, position, str("z"), font=QFont("Consolas", 10))
    
    glEnable(GL_LIGHTING)