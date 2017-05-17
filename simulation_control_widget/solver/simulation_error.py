"""

created by: lskrinjar
date of creation: 15/02/2017
time of creation: 14:05
"""

class SimulationError(object):
    """
    classdocs
    """

    def __init__(self, parent=None):
        """

        :param parent:
        """
        self._parent = parent

        self.error = False
        self.info = ""

    def setError(self, info, q=[]):
        """

        :param info:
        :return:
        """
        self.info = info
        self.error = True

        print "q =", q
        print self.info
        raise ValueError, self.info

    def setWarning(self, info):
        """

        :param info:
        :return:
        """
        self.info = info


if __name__ == "__main__":
    sol_error = SimulationError()

