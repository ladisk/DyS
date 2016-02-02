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
    def __init__(self, ID, name, parent=None):
        """
        Constructor
        """
        super(Measure, self).__init__(name, parent)

        #   parent
        self._parent = parent


    def _meaure(self, t=None):
        """

        :param t:
        :return:
        """
