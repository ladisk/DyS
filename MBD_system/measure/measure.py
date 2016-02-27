"""

created by: lskrinjar
date of creation: 27/01/2016
time of creation: 12:41
"""

from MBD_system.MBD_system_items import MeasureItem

class Measure(MeasureItem):
    """
    classdocs
    """
    def __init__(self, ID, name, measurement_type, parent=None):
        """
        Constructor
        """
        super(Measure, self).__init__(name, parent)

        #   parent
        self._parent = parent
        
        #    id of item to measure its property value
        self.id = ID
        
        #    variable to measure
        self.measurement_type = measurement_type
        self.measurement_types = ["Rx",
                                  "Ry",
                                  "theta",
                                  "dRx",
                                  "dRy",
                                  "dtheta",
                                  "Fx",
                                  "Fy",
                                  "F",
                                  "Mz"]

    def _meaure(self, t=None):
        """

        :param t:
        :return:
        """
