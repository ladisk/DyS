"""

created by: lskrinjar
date of creation: 04/11/2016
time of creation: 12:03
"""

class MeshOptions(object):
    """
    Mesh object constructor
    """

    def __init__(self, element_size=None, number_of_elements=None, parent=None):
        """

        :return:
        """
        self._parent = parent

        self.element_size = element_size
        self.number_of_elements = number_of_elements