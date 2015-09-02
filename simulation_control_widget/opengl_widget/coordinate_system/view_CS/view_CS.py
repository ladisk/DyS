'''
Created on 4. nov. 2013

@author: lskrinjar
'''
from OpenGL import *
from OpenGL.GL import *
from OpenGL.GLU import *


# from camera_view.views.views import Views
class View(object):
    '''
    Connection between mouse motion and transformation matrix
    '''
    def __init__(self):
        '''
        Constructor
        '''
        glMatrixMode(GL_MODELVIEW)
        self.current_MODELVIEW_matrix = glLoadIdentity()
        self.reset()

    def reset(self):
        glPushMatrix()
        glLoadIdentity()
        self.current_MODELVIEW_matrix = glGetDoublev(GL_MODELVIEW_MATRIX)
        glPopMatrix()
    
    def rotate(self, rx, ry, rz):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(ry, 1, 0, 0)
        glRotatef(rx, 0, 1, 0)
        glMultMatrixf( self.current_MODELVIEW_matrix )
        self.current_MODELVIEW_matrix = glGetFloatv( GL_MODELVIEW_MATRIX )
        glPopMatrix()
        
    def front(self):
        glPushMatrix()
        glLoadIdentity()
        self.current_MODELVIEW_matrix = glGetFloatv( GL_MODELVIEW_MATRIX )
        glPopMatrix()  
        
    def back(self):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(180, 0, 1, 0)
        self.current_MODELVIEW_matrix = glGetFloatv( GL_MODELVIEW_MATRIX )
        glPopMatrix()

    def left(self):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(-90, 0, 1, 0)
        self.current_MODELVIEW_matrix = glGetFloatv( GL_MODELVIEW_MATRIX )
        glPopMatrix()
        
    def right(self):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(+90, 0, 1, 0)
        self.current_MODELVIEW_matrix = glGetFloatv( GL_MODELVIEW_MATRIX )
        glPopMatrix()
        
    def bottom(self):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(-90, 1, 0, 0)
        self.current_MODELVIEW_matrix = glGetFloatv( GL_MODELVIEW_MATRIX )
        glPopMatrix()
        
    def top(self):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(+90, 1, 0, 0)
        self.current_MODELVIEW_matrix = glGetFloatv( GL_MODELVIEW_MATRIX )
        glPopMatrix()