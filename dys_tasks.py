"""
Created on 6. jun. 2016

@author: skrinjar.luka@gmail.com

"""
import sys
import time


import numpy as np
from PyQt4 import QtCore, QtGui


from dys import DySMainWindow


class Tasks(QtCore.QObject):
    def __init__(self):
        super(Tasks, self).__init__()

        self.pool = QtCore.QThreadPool()
        
        self.c_r = np.arange(0, 1, 0.1)
        self.pool.setMaxThreadCount(1)

    def process_result(self, task):
        print 'Receiving', task

    def start(self):
        app = QtGui.QApplication(sys.argv)


        for i in range(0, len(self.c_r)):
            print "i =", i

            dys_worker = DySMainWindow()
            dys_worker.show()
            dys_worker.MBD_system.contacts[0].contact_model.c_r = self.c_r[i]
            # time.sleep(.5)
            dys_worker.simulation_control_widget.simulationStart()

            dys_worker.simulation_control_widget.simulationReset()

        self.pool.waitForDone()
        sys.exit(app.exec_())

if __name__ == "__main__":
    main = Tasks()
    main.start()
