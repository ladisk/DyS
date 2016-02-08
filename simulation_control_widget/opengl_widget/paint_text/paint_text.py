"""
Created on 13. jan. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import os
import time

from OpenGL import *
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtOpenGL import *

import numpy as np
from simulation_control_widget.opengl_widget.paint_text_color_update import update_text_color


def text(self, resize_factor_width, resize_factor_height, filename, simulation_time=0, simulation_step_number=0):
    """
    in:
        filename as string
    """
    dx_left_pos = 0.02 * resize_factor_width
    dx_right_pos = 1 - (resize_factor_width * (0.01 + 0.7))
    dy_top_pos = 1 - resize_factor_height * (0.01 + 0.07)
    
    glLoadIdentity()

    #    info text
    #    change color
    update_text_color.update(self)
    #    filename
    self.renderText(dx_left_pos - 1, dy_top_pos, 0, QString("Filename: ") + QString(filename), font=QFont("Consolas", 10))
    #    time and date
    self.renderText(dx_right_pos, dy_top_pos, 0, str(time.strftime("%b %d %Y %H:%M:%S")), font=QFont("Consolas", 10))
    #    time
    self.renderText(dx_left_pos - 1, dy_top_pos - resize_factor_height * 0.1, 0, QString("Simulation time: ") + str(simulation_time), font=QFont("Consolas", 10))
    #    step num
    self.renderText(dx_left_pos - 1, dy_top_pos - resize_factor_height * 0.2, 0, QString("Step number: ") + str(simulation_step_number), font=QFont("Consolas", 9))

