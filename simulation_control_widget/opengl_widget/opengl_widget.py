'''
Created on 5. feb. 2014

@author: lskrinjar
'''
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


import camera_view.camera_view as camera_view
import coordinate_system.display_CS.display_CS as display_CS
import coordinate_system.display_CS_main.display_CS_main as display_CS_main
import display_window.display_window as display_window
import gl_matrix.gl_matrix as gl_matrix
import paint_text.paint_text as paint_text
from coordinate_system.coordinate_system import CoordinateSystem


class OpenGLWidget(QtOpenGL.QGLWidget):
    '''
    classdocs
    Angle units in opengl are always deg
    '''

    def __init__(self, MBD_system=None, parent=None):
        '''
        Constructor
        '''
        self.parent = parent
        QtOpenGL.QGLWidget.__init__(self, parent)
        #    main data file
        self.MBD_system = MBD_system
        self.step = None
        
        #    opengl properties for glOrtho()
        #    the window corner OpenGL coordinates are (-+1, -+1)
        xy_window_values = 0.5E-1
        
        self.left_GL = -1 * xy_window_values
        self.right_GL = +1 * xy_window_values
        self.bottom_GL = -1 * xy_window_values
        self.top_GL = +1 * xy_window_values
        self.near_GL = -1000
        self.far_GL = +1000
        self.zoom = 1
        self.delta_zoom = 0
        self.CS_window = 80
        
        self.scale_factor_Translation = 0.005 * xy_window_values
        self.scale_factor_Rotation = 10. * xy_window_values
        
        #    window dimensions
        self.initial_window_width = 500.
        self.initial_window_height = 500.
        self.resize_factor_width = 1
        self.resize_factor_height = 1
    
    
        #    coordinate system
        self.LCS = CoordinateSystem()
        
    
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
        self.gl_matrix = gl_matrix.Matrix()
        #    camera view
        self.camera = camera_view.View()
        #    coordinate system view
        self.camera_CS = camera_view.View()

        #    viewing
        self.geometry()

        #   create VBO of GCS
        self.LCS._create_VBO()

        #    display window
        self.display_window = display_window.Window(width=self.initial_window_width, height=self.initial_window_height, left_GL=self.left_GL, right_GL=self.right_GL, bottom_GL=self.bottom_GL, top_GL=self.top_GL, near_GL=self.near_GL, far_GL=self.far_GL, zoom=self.zoom, scale_factor_Translation=self.scale_factor_Translation, scale_factor_Rotation=self.scale_factor_Rotation)
        
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
        glLoadMatrixf(self.camera.current_PROJECTION_matrix)
        #    modelview matrix
        glMatrixMode(GL_MODELVIEW)


        if self.LCS._visible:
            try:
                self.LCS._paintGL_VBO()
            except:
                pass
            # display_CS_main.display(self, self.camera.current_MODELVIEW_matrix, int(self.display_window.width), int(self.display_window.height))
        
        glViewport(0, 0, int(self.display_window.width), int(self.display_window.height))

        glLoadMatrixf(self.camera.current_MODELVIEW_matrix)

        glPointSize(10)
        glBegin(GL_POINTS)
        for contact in self.MBD_system.contacts:
            contact._paing_GL_GCS(self.step)
        glEnd()
        glLineWidth(1.2)

        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_NORMAL_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)

        #    display main GCS
#         if self.LCS._visible:
#             self.LCS._paintGL_VBO()

        #   plot markers in ground
        self.MBD_system.ground.paintGL_VBO()


        self.MBD_system.bodies.sort(key=lambda body: body.transparent_GL, reverse=True)

        
        for body in self.MBD_system.bodies:
            glUseProgram(0)
            # if body._visible:# and body.VBO_created:
            #    we don't want each object to move the camera
            glPushMatrix()

            #   paint VBO of global coordinate system
            try:
                self.LCS._paintGL_VBO()
            except:
                pass

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
        self.gl_matrix.update_modelview_matrix_CS(self.camera_CS.current_MODELVIEW_matrix)
        
        
        #    display CS symbol in lower left corner
        try:
            display_CS.display(self, self.gl_matrix, self.CS_window, self.CS_window)
        except:
            pass

        #    reload saved matrices from gl_matrix object
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glViewport(0, 0, self.display_window.width, self.display_window.height)

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
            if not self.LCS._VBO_created:
                self.LCS._create_VBO()


            for body in self.MBD_system.bodies:
                #    create VBO of each body
                if body.VBO_created == False:
                    body.create_VBO()

            # for force in self.MBD_system.forces:
            #     for marker in force.markers:
            #         if not marker._VBO_created:
            #             marker._create_VBO()

   
    def takeSnapShot(self):
        """
        
        """
        screen_shot = self.grabFrameBuffer(withAlpha=True)
        return screen_shot
        
        
    def resizeGL(self, width, height):
        """
        Resize
        """
        self.scale_factor_Translation = self.scale_factor_Translation * float(self.display_window.width) / float(self.parent.geometry().width())
        self.scale_factor_Rotation = self.scale_factor_Rotation * float(self.display_window.width) / float(self.parent.geometry().width())
        #    paint within the whole window
        self.width = width
        self.height = height

        #    resize factor
        self.resize_factor_width = self.initial_window_width / self.width
        self.resize_factor_height = self.initial_window_height / self.height

        self.display_window = display_window.Window(width, height, self.left_GL, self.right_GL, self.bottom_GL, self.top_GL, self.near_GL, self.far_GL, self.zoom, self.scale_factor_Translation, self.scale_factor_Rotation)

        glViewport(0, 0, self.display_window.width, self.display_window.height)
        
        #    aspect ratio
        self.aspect = float(width) / height

        #    set orthographic projection (2D only) 
        self.camera.zoom(self.aspect, self.left_GL, self.right_GL, self.bottom_GL, self.top_GL, self.near_GL, self.far_GL, self.zoom)

        glMatrixMode(GL_MODELVIEW)
        
    
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
        
        """
        self.last_pos = QtCore.QPoint(event.pos())
        self.setCursor(Qt.CursorShape(Qt.ArrowCursor))


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
        _show_GCSAction = QtGui.QAction("Show GCS", self, checkable=True, checked=self.LCS._visible)
        _show_GCSAction.triggered.connect(self.LCS._show)
        popMenu.addAction(_show_GCSAction)


        popMenu.exec_(event.globalPos())


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
                self.camera.translate(tX, tY, 0)
            #    rotate
            else:
                self.setCursor(Qt.CursorShape(Qt.OpenHandCursor))
                rX = dx * self.scale_factor_Rotation
                rY = dy * self.scale_factor_Rotation
                self.camera.rotate(rX, rY, 0)
                self.camera_CS.rotate(rX, rY, 0)             

            self.updateGL()
            self.last_pos = event.posF()


    
    def wheelEvent(self, event):
        """
        Zoom
        """
        self.delta_zoom = event.delta() / 120
        
        self.zoom = self.zoom * (1 + 0.1 * self.delta_zoom)
        self.scale_factor_Translation = self.scale_factor_Translation * (1 + 0.1 * self.delta_zoom)
        self.scale_factor_Rotation = self.scale_factor_Rotation * (1 + 0.01 * self.delta_zoom)
        self.aspect = self.display_window.aspect
        
        self.camera.zoom(self.aspect, self.left_GL, self.right_GL, self.bottom_GL, self.top_GL, self.near_GL, self.far_GL, self.zoom)

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
        self.camera.front()
        self.camera_CS.front()
        self.updateGL()


    @QtCore.pyqtSlot()
    def viewBack(self):
        self.camera.back()
        self.camera_CS.back()
        self.updateGL()


    @QtCore.pyqtSlot()
    def viewBottom(self):
        self.camera.bottom()
        self.camera_CS.bottom()
        self.updateGL()


    @QtCore.pyqtSlot()
    def viewTop(self):
        self.camera.top()
        self.camera_CS.top()
        self.updateGL()


    @QtCore.pyqtSlot()
    def viewRight(self):
        self.camera.left()
        self.camera_CS.left()
        self.updateGL()


    @QtCore.pyqtSlot()
    def viewLeft(self):
        self.camera.right()
        self.camera_CS.right()
        self.updateGL()


    @QtCore.pyqtSlot()
    def viewIsometric(self):
        self.camera.isometric()
        self.camera_CS.isometric()
        self.updateGL()
