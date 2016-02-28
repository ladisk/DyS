"""

created by: lskrinjar
date of creation: 28/02/2016
time of creation: 18:34
"""
from OpenGL.GL import *

from coordinate_system import CoordinateSystem


class CoordinateSystemView(CoordinateSystem):
    """
    classdocs
    """

    def __init__(self, width, height, parent=None):
        self.scale = 200
        super(CoordinateSystemView, self).__init__(parent=parent)

        #   parent
        self._parent = parent
        #   coordinate system view
        #   width
        self.width = width
        #   height
        self.height = height

        #   visible
        self._visible = True

    def paintGL(self, matrix):
        """

        :return:
        """
        #    in display CS viewport opengl projection matrix is always glLoadIdentity()
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()

        #    update GL_MODELVIEW_MATRIX
        glMatrixMode(GL_MODELVIEW)
        glLoadMatrixd(matrix.modelview_matrix_CS)

        #   viewport
        glViewport(0, 0, int(self.width), int(self.height))

        glClear(GL_DEPTH_BUFFER_BIT)
        glDisable(GL_LIGHTING)

        # self._paintGL_VBO()
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

        self._parent.qglColor(self._parent.text_color)

        axis_length = 0.7
        position = axis_length + 0.1

        self._parent.renderText(position, 0, 0, str("x"), font=self._parent.font)
        self._parent.renderText(0, position, 0, str("y"), font=self._parent.font)
        self._parent.renderText(0, 0, position, str("z"), font=self._parent.font)

        glEnable(GL_LIGHTING)