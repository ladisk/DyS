'''
Created on 13. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''


from pprint import pprint
import sys

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from job_list_widget_ui import Ui_Form


class JobListWidget(QWidget):
    '''
    classdocs
    '''

    def __init__(self, job_list_, flags = 0, parent = None):
        '''
        Constructor
        '''
        QtGui.QWidget.__init__(self,parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.ui.tableWidget.setColumnCount(3)
        self.ui.tableWidget.setRowCount(10)
        self.ui.tableWidget.setHorizontalHeaderLabels(["Name", "Started", "Finished"]) 

        self.job_list_ = job_list_
        
        self.create_job_list(self.job_list_)

    def create_job_list(self, job_list_):
        for i in xrange(0, len(job_list_)):
            self.ui.tableWidget.setItem(i, 0, QtGui.QTableWidgetItem(QString(self.job_list_[i])))
            




