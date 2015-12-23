"""
Created on 10. maj 2013

@author: lskrinjar
"""
import numpy as np
import platform
import sys
import scipy

#    computer
print "computer"
import socket
print socket.gethostname()


#    OS props
print "OS props"
print platform.platform(aliased=0, terse=0)

#    python props
print "python props"
print platform.python_build()

#    python version
print "python version"
print platform.python_version()

#    OS bit and OS type
print "OS bit and OS type"
print platform.architecture()

#    check if system is 32bit or 64bit
print "check if system is 32bit or 64bit"
print("%x" % sys.maxsize, sys.maxsize < 2 ** 64)

#    check the version of numpy
print "numpy version"
print np.__version__

#    check matplotlib version
import matplotlib
print "matplotlib version"
print matplotlib.__version__

#    check scipy version
print "scipy version"
print scipy.__version__
# 
#     PyOpenGL version
import pkg_resources
print "PyOpenGL version"
print pkg_resources.get_distribution("PyOpenGL").version


print "PyQt4 version"
from PyQt4.QtCore import QT_VERSION_STR
# from PyQt4.pyqtconfig import Configuration

print QT_VERSION_STR

# #   itertools version
# import itertools
from OpenGL import __version__ as OpenGL_version
print "OpenGL_version =", OpenGL_version

from OpenGL_accelerate import __version__
print "OpenGL_accelerate_version =", __version__

# from OpenGL_accelerate import __version__ as OpenGLaccl_version
# print "OpenGLaccl_version =", OpenGLaccl_version


# from OpenGL import *
# from OpenGL.GL import *
# from OpenGL.GLE import *
# from OpenGL.GLU import *
# from OpenGL.GLUT import *
#  
# # wglMakeCurrent
# glutInit()
# print glGetString(GL_VERSION)

import ctypes
print ctypes.sizeof(ctypes.c_voidp)