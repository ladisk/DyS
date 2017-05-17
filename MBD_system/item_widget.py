__author__ = 'lskrinjar'

from PyQt4 import QtCore, QtGui


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


class ItemWidget(QtGui.QDialog):

    def __init__(self, parent=None):
        super(ItemWidget, self).__init__(parent=parent)



        self._parent = parent


    def closeEvent(self, event):
        """

        :param event:
        :return:
        """
        self._parent._parent.SimulationControlWidget.OpenGLWidget._repaintGL()


    def _move(self):
        """
        Move widget (dialog window) to screen center
        """
        frameGm = self.frameGeometry()
        screen = QtGui.QApplication.desktop().screenNumber(QtGui.QApplication.desktop().cursor().pos())
        centerPoint = QtGui.QApplication.desktop().screenGeometry(screen).center()
        frameGm.moveCenter(centerPoint)
        self.move(frameGm.topLeft())


    def _cancel(self):
        """

        """
        self.close()