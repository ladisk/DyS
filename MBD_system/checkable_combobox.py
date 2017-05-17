"""

created by: lskrinjar
date of creation: 16/05/2016
time of creation: 19:45
"""
from PyQt4 import QtGui, QtCore


class CheckableComboBox(QtGui.QComboBox):
    def __init__(self, parent):
        super(CheckableComboBox, self).__init__(parent)
        self.view().pressed.connect(self.handleItemPressed)
        self.setModel(QtGui.QStandardItemModel(self))

    def handleItemPressed(self, index):
        item = self.model().itemFromIndex(index)
        if item.checkState() == QtCore.Qt.Checked:
            item.setCheckState(QtCore.Qt.Unchecked)
        else:
            item.setCheckState(QtCore.Qt.Checked)
