"""
Created on 30. nov. 2013

@author: lskrinjar
"""
import os

import numpy as np


from MBD_system.MBD_system_items import TreeItem
from model_loader.model_loader import ModelLoader

class Geometry(TreeItem):
    """
    classdocs
    """
    def __init__(self, filename, parent=None):
        """
        Constructor
        Class constructor for OpenGL VBO
        in:
            CAD_LCS_GCS - array of center of body mass in GCS
            geom_file - a body geometry full filename (filename = name + extension)
            geom_file_extension - extension
        out:
            VBO interleaved array of data saved to GPU memory
        """
        super(Geometry, self).__init__(filename, parent)
        #   parent
        self._parent = parent

        #    type
        self._typeInfo = "geometry"

        #    filename
        self.filename = filename
        if self.filename is not None:
            self._filename, self._file_extension = os.path.splitext(filename)

            #    read geom data from stl (geometry) geom_file with model_loader
            if os.path.isfile(filename):
                self.geom_data = ModelLoader(filename)

        else:
            print "File not found!"

        #   visualization properties
        self._visible = True
        self.vtk_actor = None

    def set_vtk_data(self):
        """

        :return:
        """

    def __del__(self):
        """
        Delete vtk object
        """


