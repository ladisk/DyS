"""

created by: lskrinjar
date of creation: 25/07/2016
time of creation: 16:12
"""
import vtk


class VTKCoordinateSystem(vtk.vtkAxesActor):
    """

    """

    def __init__(self, parent=None):
        """

        :return:
        """
        super(VTKCoordinateSystem, self).__init__(parent)

        self.SetConeRadius(0.)
        self.AxisLabelsOff()

        self.setSize()

    def setSize(self, size=2E-3):
        """
        Set size of coordinate system
        :param size:
        :return:
        """
        self._GCS_size = size
        self.SetTotalLength(self._GCS_size, self._GCS_size, self._GCS_size)

