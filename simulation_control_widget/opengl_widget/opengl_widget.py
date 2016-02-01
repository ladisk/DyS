"""
Created on 5. feb. 2014

@author: lskrinjar
"""
import os
from pprint import pprint
import time
import numpy as np


from OpenGL import *
from OpenGL.GL import *
from OpenGL.GL.ARB.vertex_buffer_object import *
from OpenGL.GL.shaders import *
from OpenGL.GLU import *
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtOpenGL import *


import coordinate_system.display_CS.display_CS as display_CS
import paint_text.paint_text as paint_text
from simulation_control_widget.opengl_widget.view.view import View
from coordinate_system.coordinate_system import CoordinateSystem
from simulation_control_widget.opengl_widget.matrix.matrix import Matrix
from simulation_control_widget.opengl_widget.window.window import Window


class OpenGLWidget(QtOpenGL.QGLWidget):
    """
    classdocs
    Angle units in opengl are always deg
    """

    def __init__(self, MBD_system=None, parent=None):
        """
        Constructor
        """
        self._parent = parent
        # QtOpenGL.QGLWidget.__init__(self, parent)
        super(OpenGLWidget, self).__init__(parent)
        # super(SimulationControlWidget, self).__init__(parent)
        #    main data file
        self.MBD_system = MBD_system
        self.step = None
        
        #    opengl properties for glOrtho()
        #    the window corner OpenGL coordinates are (-+1, -+1)
        xy_window_values = 1.E-1
        
        self.left_GL = -1 * xy_window_values
        self.right_GL = +1 * xy_window_values
        self.bottom_GL = -1 * xy_window_values
        self.top_GL = +1 * xy_window_values
        self.near_GL = -1000
        self.far_GL = +1000
        self.zoom = 1
        self.delta_zoom = 0
        self.CS_window = 80
        
        self.__scale_factor_Translation = self.scale_factor_Translation = 0.001
        self.__scale_factor_Rotation = self.scale_factor_Rotation = 1
        
        #    window dimensions
        self.initial_window_width = self._parent.geometry().width() #500
        self.initial_window_height = self._parent.geometry().height() #500
        self.resize_factor_width = 1
        self.resize_factor_height = 1

        #    coordinate system
        self.GCS = CoordinateSystem()#parent=self

        # core profile
        glformat = QtOpenGL.QGLFormat()
        glformat.setVersion(4, 1)
        glformat.setProfile(QtOpenGL.QGLFormat.CoreProfile)
        glformat.setSampleBuffers(True)

        #    context menu properties
        self.setContextMenuPolicy(QtCore.Qt.NoContextMenu)  # NoContextMenu, CustomContextMenu

    def update_data(self, bodies):
        """
        
        """
        #    main data file
        self.MBD_system.bodies = bodies

    def initializeOverlayGL(self):
        """
        
        """

    def initializeGL(self):
        """
        Initialize OpenGL, VBOs, upload data on the GPU, etc.
        """
        #    background color - default color
        glClearColor(0, 0, 0, 1)
        
        #    enable functions
        #    z-buffer
        glEnable(GL_DEPTH_TEST)
        #    shading flat or smooth GL_FLAT or GL_SMOOTH
        glShadeModel(GL_FLAT)
        #    normal vectors are normalised to unit length after transformation and before lighting
        glEnable(GL_NORMALIZE)
#         glEnable(GL_CULL_FACE)
        glEnable(GL_COLOR_MATERIAL)
        glEnable(GL_LIGHTING)
        #    for transparent models
        glEnable(GL_BLEND)
        #    activate alpha factor for blending
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        
        glDisable(GL_CULL_FACE)
        #    get OpenGL info
#         print "GL Version:", glGetString(GL_VERSION)
#         print "GL Vendor:", glGetString(GL_VENDOR)
#         print "GL Renderer:", glGetString(GL_RENDERER)
        #         print "GL Extensions:", glGetString(GL_EXTENSIONS)
#         print "---------------------------------------------------------"

        #    lights
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0)
        glEnable(GL_LIGHT0)
        #    gl matrix object
        self.gl_matrix = Matrix()
        #    camera view
        self.view = View()
        #    coordinate system view
        self.view_CS = View()

        #    viewing
        self.geometry()

        #   create VBO of GCS
        self.GCS._create_VBO()

        #    display window
        self.window = Window(width=self.initial_window_width, height=self.initial_window_height, left_GL=self.left_GL, right_GL=self.right_GL, bottom_GL=self.bottom_GL, top_GL=self.top_GL, near_GL=self.near_GL, far_GL=self.far_GL, zoom=self.zoom, scale_factor_Translation=self.scale_factor_Translation, scale_factor_Rotation=self.scale_factor_Rotation)
        self.scale_factor_Translation = self.scale_factor_Translation * float(self.window.width) / float(self._parent.geometry().width())
        self.scale_factor_Rotation = self.scale_factor_Rotation * float(self.window.width) / float(self._parent.geometry().width())

        #    shaders
        self.shader_program = QGLShaderProgram()
        self.shader_program.addShaderFromSourceFile(QGLShader.Vertex, os.path.join(os.path.dirname(os.path.abspath(__file__)), "shaders", "shader.vert"))
        self.shader_program.addShaderFromSourceFile(QGLShader.Fragment, os.path.join(os.path.dirname(os.path.abspath(__file__)), "shaders", "shader.frag"))
        self.shader_program.log()
        self.shader_program.link()

    def paintGL(self):
        """
        Display geometry
        """
        # Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
                
        #    setup camera
        #    projection matrix
        glMatrixMode(GL_PROJECTION)
        glLoadMatrixf(self.view.current_PROJECTION_matrix)
        #    modelview matrix
        glMatrixMode(GL_MODELVIEW)

        #   display GCS in main window
        # if self.MBD_system.GCS_visible:
            #   plot markers in ground
            # pprint(vars(self.MBD_system.ground))
        # self.MBD_system.ground.paintGL_VBO()
        # if self.GCS._visible:
        #     self.GCS._paintGL_VBO()
                # try:
                #     self.GCS._paintGL_VBO()
                # except:
                #     pass
            # display_CS_main.display(self, self.view.current_MODELVIEW_matrix, int(self.window.width), int(self.window.height))
        
        glViewport(0, 0, int(self.window.width), int(self.window.height))

        glLoadMatrixf(self.view.current_MODELVIEW_matrix)

        glPointSize(10)
        glBegin(GL_POINTS)
        for contact in self.MBD_system.contacts:
            contact._paing_GL_GCS(self.step)
        glEnd()
        glLineWidth(1.2)

        # self.MBD_system.ground.paintGL_VBO()

        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)

        #   display main CS (GCS) as coordinate system VBO does not have normals it has
        #   to be displayed before normals are enabled in OpenGL
        if self.GCS._visible:
            self.GCS._paintGL_VBO()

        glEnableClientState(GL_NORMAL_ARRAY)

        #    display main GCS
#         if self.GCS._visible:
#             self.GCS._paintGL_VBO()

        self.MBD_system.bodies.sort(key=lambda body: body.transparent_GL, reverse=True)

        for body in self.MBD_system.bodies:
            glUseProgram(0)
            # if body._visible:# and body.VBO_created:
            #    we don't want each object to move the camera
            glPushMatrix()

            #   paint VBO of global coordinate system
            # try:
            # self.GCS._paintGL_VBO()
            # except:
            #     pass

            #   paint body geometry
            body.paintGL_VBO(self.step, self.shader_program)
            #   paint body AABB tree
            # body._paintGL_VBO_AABBtree(self.shader_program)
            # glUseProgram(0)

#                 #    this has to be done in body.paintGL_VBO() - after transformations matrices
#                 if body.AABBtree != None and body.AABBtree.visibility and body.AABBtree.VBO_created:
#                     self.shader_program.bind()
#                     color_location = self.shader_program.uniformLocation("vertex_color")
#                     self.shader_program.setUniformValue(color_location, QVector3D(body.AABBtree.color_GL[0], body.AABBtree.color_GL[1], body.AABBtree.color_GL[2]))  # cg.cgtypes.vec3(body.AABBtree.color_GL) ctypes.c_floatQVector3D(1, 0, 0)
#                     body.AABBtree.paintGL_VBO_tree()
#                     glUseProgram(0)
            glPopMatrix()

        self.MBD_system.bodies.sort(key=lambda body: body.body_id, reverse=False)

        glDisableClientState(GL_COLOR_ARRAY)
        glDisableClientState(GL_NORMAL_ARRAY)
        glDisableClientState(GL_VERTEX_ARRAY)

        #    display CS orientation
        #    save matrices to gl_matrix object
        self.gl_matrix.update_projection_matrix(glGetDoublev(GL_PROJECTION_MATRIX)) 
        self.gl_matrix.update_modelview_matrix(glGetDoublev(GL_MODELVIEW_MATRIX))
        self.gl_matrix.update_modelview_matrix_CS(self.view_CS.current_MODELVIEW_matrix)
        
        #    display CS symbol in lower left corner
        try:
            display_CS.display(self, self.gl_matrix, self.CS_window, self.CS_window)
        except:
            pass

        #    reload saved matrices from gl_matrix object
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glViewport(0, 0, self.window.width, self.window.height)

        glMatrixMode(GL_MODELVIEW)
        
        #    display text
        paint_text.text(self, self.resize_factor_width, self.resize_factor_height, filename=self.MBD_system._name, simulation_time=self.MBD_system.time, simulation_step_number=self.MBD_system.step_num)
    
    def repaintGL(self, step = None):
        """
        Repaint - update opengl
        """
        self.step = step
        self.updateGL()
    
    def geometry(self):
        """
        Create MBD System (MBD_system) geometry
        """
        if glGetError() == GL_NO_ERROR:
            if not self.GCS._VBO_created:
                self.GCS._create_VBO()

            for body in self.MBD_system.bodies:
                #    create VBO of each body
                if not body.VBO_created and body._name != "ground":
                    body.create_VBO()

            for force in self.MBD_system.forces:
                for marker in force.markers:
                    if not marker._VBO_created:
                        marker._create_VBO()

    def refit(self):
        """
        Function finds max coordinates of a point on a body in GCS to resize window on all system geometries


        :return:
        """
        _uP_i_max = []
        for body in self.MBD_system.bodies:
            _uP_i_max_i = body.get_uP_i_max() + abs(body.R)

            _uP_i_max.append(_uP_i_max_i)

        #   convert to numpy array
        _uP_i_max = np.array(_uP_i_max)

        #   max element in every axis
        [x, y, z] = np.amax(_uP_i_max, axis=0)
        #   max dim is 50% larger than max coordinate
        self._max_dim = 1.5 * np.amax(np.array([x, y]))

        delta = self.view._right_GL / self._max_dim


        self.left_GL = - self._max_dim
        self.right_GL = self._max_dim

        self.bottom_GL = -self._max_dim
        self.top_GL = self._max_dim

        self.view.zoom(self.aspect, self.left_GL, self.right_GL, self.bottom_GL, self.top_GL, self.near_GL, self.far_GL, self.zoom)

        self.scale_factor_Rotation = self.__scale_factor_Rotation / delta
        self.scale_factor_Translation = .1 * self.__scale_factor_Translation / delta

    def get_aspect(self):

        return self.window.aspect

    def takeSnapShot(self):
        """
        Create a snapshot
        """
        screen_shot = self.grabFrameBuffer(withAlpha=True)
        pprint(vars(screen_shot))
        return screen_shot

    def resizeGL(self, width, height):
        """
        Resize

        :param width:
        :param height:

        """
        #    paint within the whole window
        self.width = width
        self.height = height

        #    resize factor
        self.resize_factor_width = self.initial_window_width / self.width
        self.resize_factor_height = self.initial_window_height / self.height

        self.window = Window(width, height, self.left_GL, self.right_GL, self.bottom_GL, self.top_GL, self.near_GL, self.far_GL, self.zoom, self.scale_factor_Translation, self.scale_factor_Rotation)

        glViewport(0, 0, self.window.width, self.window.height)
        
        #    aspect ratio of width by height
        self.aspect = float(width) / height

        #    set orthographic projection (2D only) 
        self.view.zoom(self.aspect, self.left_GL, self.right_GL, self.bottom_GL, self.top_GL, self.near_GL, self.far_GL, self.zoom)

        glMatrixMode(GL_MODELVIEW)

        self.refit()

    def mouseReleaseEvent(self, event):
        """
        
        """
        self.last_pos = QtCore.QPoint(event.pos())

        if event.button() == QtCore.Qt.LeftButton:
            self.setCursor(Qt.CursorShape(Qt.ArrowCursor))
         
        elif event.button() == QtCore.Qt.RightButton:
            self.setCursor(Qt.CursorShape(Qt.ArrowCursor))
            self.contextMenuEvent(event)
        else:
            self.setCursor(Qt.CursorShape(Qt.ArrowCursor))

    def mousePressEvent(self, event):
        """
        Runs when mouse button is pressed
        """
        self.last_pos = QtCore.QPoint(event.pos())
        self.setCursor(Qt.CursorShape(Qt.ArrowCursor))

        # print "event =", event
        # pprint(vars(event))
        if event.button() == QtCore.Qt.LeftButton:
            print event.pos()

        # print "-------------------------------------------"

    def _repaintGL(self):
        """
        
        """
        self.repaintGL()

    def contextMenuEvent(self, event):
        """

        """
        popMenu = QtGui.QMenu(parent=self)

        _updateGLAction = popMenu.addAction("Refresh")
        
        _updateGLAction.triggered.connect(self._repaintGL)
        popMenu.addSeparator()

        popMenu.addMenu(self.viewSubMenu())

        popMenu.addSeparator()
        _show_GCSAction = QtGui.QAction("Show GCS", self, checkable=True, checked=self.GCS._visible)
        _show_GCSAction.triggered.connect(self.GCS._show)
        popMenu.addAction(_show_GCSAction)

        popMenu.exec_(event.globalPos())

        self.repaintGL()

    def mouseMoveEvent(self, event):
        """
        Rotate and translate
        """
        #    relative position coordinates of a mouse
        dxy = event.posF() - self.last_pos
        
        dx = dxy.x()
        dy = dxy.y()

        if event.buttons() & Qt.MidButton:
            #    translate
            if event.modifiers() & Qt.ControlModifier:
                self.setCursor(Qt.CursorShape(Qt.SizeAllCursor))

                tX = dx * self.scale_factor_Translation
                tY = -dy * self.scale_factor_Translation
                self.view.translate(tX, tY, 0)
            #    rotate
            else:
                self.setCursor(Qt.CursorShape(Qt.OpenHandCursor))
                rX = dx * self.scale_factor_Rotation
                rY = dy * self.scale_factor_Rotation
                self.view.rotate(rX, rY, 0)
                self.view_CS.rotate(rX, rY, 0)             

            self.updateGL()
            self.last_pos = event.posF()
    
    def wheelEvent(self, event):
        """
        Zoom
        """
        self.delta_zoom = event.delta() / 120 #120
        
        self.zoom = self.zoom * (1 + 0.1 * self.delta_zoom)
        self.scale_factor_Translation = self.scale_factor_Translation * (1 + 0.1 * self.delta_zoom)
        self.scale_factor_Rotation = self.scale_factor_Rotation * (1 + 0.01 * self.delta_zoom)
        self.aspect = self.window.aspect
        self.view.zoom(self.aspect, self.left_GL, self.right_GL, self.bottom_GL, self.top_GL, self.near_GL, self.far_GL, self.zoom)

        self.updateGL()

    def viewSubMenu(self):
        viewDirectionSubmenu = QtGui.QMenu(parent=self)
        
        viewDirectionSubmenu.setTitle("View")
        
        _viewFrontAction = viewDirectionSubmenu.addAction("Front")
        _viewFrontAction.triggered.connect(self.viewFront)

        _viewBackAction = viewDirectionSubmenu.addAction("Back")
        _viewBackAction.triggered.connect(self.viewBack)

        _viewBottomAction = viewDirectionSubmenu.addAction("Bottom")
        _viewBottomAction.triggered.connect(self.viewBottom)

        _viewTopAction = viewDirectionSubmenu.addAction("Top")
        _viewTopAction.triggered.connect(self.viewTop)

        _viewLeftAction = viewDirectionSubmenu.addAction("Left")
        _viewLeftAction.triggered.connect(self.viewLeft)

        _viewRightAction = viewDirectionSubmenu.addAction("Right")
        _viewRightAction.triggered.connect(self.viewRight)

        _viewIsometricAction = viewDirectionSubmenu.addAction("Isometric")
        _viewIsometricAction.triggered.connect(self.viewIsometric)

        return viewDirectionSubmenu

    @QtCore.pyqtSlot()
    def viewFront(self):
        self.view.front()
        self.view_CS.front()
        self.updateGL()

    @QtCore.pyqtSlot()
    def viewBack(self):
        self.view.back()
        self.view_CS.back()
        self.updateGL()

    @QtCore.pyqtSlot()
    def viewBottom(self):
        self.view.bottom()
        self.view_CS.bottom()
        self.updateGL()

    @QtCore.pyqtSlot()
    def viewTop(self):
        self.view.top()
        self.view_CS.top()
        self.updateGL()

    @QtCore.pyqtSlot()
    def viewRight(self):
        self.view.left()
        self.view_CS.left()
        self.updateGL()

    @QtCore.pyqtSlot()
    def viewLeft(self):
        self.view.right()
        self.view_CS.right()
        self.updateGL()

    @QtCore.pyqtSlot()
    def viewIsometric(self):
        self.view.isometric()
        self.view_CS.isometric()
        self.updateGL()
