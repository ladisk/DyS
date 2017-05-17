"""

created by: lskrinjar
date of creation: 29/12/2016
time of creation: 09:26
"""
import os
import shutil
import threading


import vtk


class MovieMaker(threading.Thread):
    """
    classdocs
    """

    def __init__(self, folder=None, fps=24, parent=None):
        """

        :return:
        """
        threading.Thread.__init__(self)

        #   parent
        self._parent = parent

        #    name
        self._name = "_animation"

        #   set filter
        self.windowToImageFilter = None

        #   movie writer
        self.movieWriter = None

        #   folder abs path
        #   absolute path to store animation images to use them to create video file
        if folder is None:
            self.folder = os.path.join(os.getcwd(), self._name)
        else:
            self.folder = folder

        #   information
        self.completed = False
        self.delete_folder_when_completed = False

    def setRenderWindow(self, renWin):
        """

        :param renWin:
        :return:
        """
        self.windowToImageFilter = vtk.vtkWindowToImageFilter()
        self.windowToImageFilter.SetInput(renWin)
        self.windowToImageFilter.SetInputBufferTypeToRGB()
        self.windowToImageFilter.ReadFrontBufferOn()#Off
        self.windowToImageFilter.Update()

    def setInputConnection(self):
        """

        :return:
        """
        if self.windowToImageFilter is not None:
            self.movieWriter = vtk.vtkAVIWriter()
            # self.movieWriter = vtk.vtkOggTheoraWriter()
            self.movieWriter.SetInputConnection(self.windowToImageFilter.GetOutputPort())

            filename = self._name + ".avi"
            # filename = self._name + ".ogv"
            self.movieWriter.SetFileName(filename)

            # self.movieWriter.Start()