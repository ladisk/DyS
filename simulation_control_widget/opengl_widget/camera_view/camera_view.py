'''
Created on 4. nov. 2013

@author: lskrinjar
'''
from OpenGL import *
from OpenGL.GL import *
from OpenGL.GLU import *


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
        self.current_PROJECTION_matrix = glLoadIdentity()
        self.reset()
        
        
    def reset(self):
        glPushMatrix()
        glLoadIdentity()
        self.current_MODELVIEW_matrix = glGetDoublev(GL_MODELVIEW_MATRIX)
        glPopMatrix()
              
    def translate(self, tx, ty, tz):
        glPushMatrix()
        glLoadIdentity()
        glTranslatef(tx, ty, tz)
        glMultMatrixf(self.current_MODELVIEW_matrix)
        self.current_MODELVIEW_matrix = glGetFloatv(GL_MODELVIEW_MATRIX)
        glPopMatrix()        


    def rotate(self, rx, ry, rz):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(ry, 1, 0, 0)
        glRotatef(rx, 0, 1, 0)
        glMultMatrixf(self.current_MODELVIEW_matrix)
        self.current_MODELVIEW_matrix = glGetFloatv(GL_MODELVIEW_MATRIX)
        glPopMatrix()


    def zoom(self, aspect, left_GL, right_GL, bottom_GL, top_GL, near_GL, far_GL, zoom):
        """

        :param aspect:
        :param left_GL:
        :param right_GL:
        :param bottom_GL:
        :param top_GL:
        :param near_GL:
        :param far_GL:
        :param zoom:
        :return:
        """
        self._left_GL = -aspect * abs(left_GL * zoom)
        self._right_GL = +aspect * abs(right_GL * zoom)
        
        self._bottom_GL = -abs(bottom_GL * zoom)
        self._top_GL = +abs(top_GL * zoom)

        self._near_GL = near_GL
        self._far_GL = far_GL
        
        glViewport(0,
                   0,
                   2 * int(self._right_GL),
                   2 * int(self._top_GL))
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        self._orthoGL()
        self.current_PROJECTION_matrix = glGetDoublev(GL_PROJECTION_MATRIX)

        glMatrixMode(GL_MODELVIEW)


    def _orthoGL(self):

        glOrtho(self._left_GL,
                self._right_GL,
                self._bottom_GL,
                self._top_GL,
                self._near_GL,
                self._far_GL)


    def front(self):
        glPushMatrix()
        glLoadIdentity()
        self.current_MODELVIEW_matrix = glGetFloatv(GL_MODELVIEW_MATRIX)
        glPopMatrix()  
        
        
    def back(self):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(180, 0, 1, 0)
        self.current_MODELVIEW_matrix = glGetFloatv(GL_MODELVIEW_MATRIX)
        glPopMatrix()


    def left(self):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(-90, 0, 1, 0)
        self.current_MODELVIEW_matrix = glGetFloatv(GL_MODELVIEW_MATRIX)
        glPopMatrix()


    def right(self):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(+90, 0, 1, 0)
        self.current_MODELVIEW_matrix = glGetFloatv(GL_MODELVIEW_MATRIX)
        glPopMatrix()


    def bottom(self):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(-90, 1, 0, 0)
        self.current_MODELVIEW_matrix = glGetFloatv(GL_MODELVIEW_MATRIX)
        glPopMatrix()


    def top(self):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(+90, 1, 0, 0)
        self.current_MODELVIEW_matrix = glGetFloatv(GL_MODELVIEW_MATRIX)
        glPopMatrix()


    def isometric(self):
        glPushMatrix()
        glLoadIdentity()
        glRotatef(+45, 1, -0.783, -0.217)
        self.current_MODELVIEW_matrix = glGetFloatv(GL_MODELVIEW_MATRIX)
        glPopMatrix()

