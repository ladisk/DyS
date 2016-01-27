'''
Created on 3. feb. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

from OpenGL import *
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtOpenGL import *

import numpy as np


def update(self):
    #    change color
    if np.linalg.norm([np.array(glGetFloat(GL_COLOR_CLEAR_VALUE))], 2)  > 1.2:
        self.qglColor(QtCore.Qt.black)
    else:
        self.qglColor(QtCore.Qt.white)
        
    return None