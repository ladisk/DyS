'''
Created on 15. mar. 2014

@author: Ales
'''
import sys

from PyQt4 import QtCore, QtGui

from output_widget_ui import Ui_Form


class TextOutputSignal(QtCore.QObject):

    textWritten = QtCore.pyqtSignal(str)

    def write(self, text):
        self.textWritten.emit(str(text))

class OutputWidget(QtGui.QWidget):

    def __init__(self, parent=None, flags=0):
        super(OutputWidget, self).__init__(parent)
        
        #    Install the custom output stream
        sys.stdout = TextOutputSignal(textWritten=self.normalOutputWritten)
        sys.stderr = TextOutputSignal(textWritten=self.errorOutputWritten)
        
        self.ui = Ui_Form()
        self.ui.setupUi(self)

    def __del__(self):
        #    Restore sys.stdout
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

    def setWindowFlags(self, flags):
        super(OutputWidget, self).setWindowFlags(flags)

    def normalOutputWritten(self, text):
        cursor = self.ui.textOutputWidget.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        cursor.insertText(text)
        self.ui.textOutputWidget.ensureCursorVisible()
        
    def errorOutputWritten(self, text):
        self.normalOutputWritten("*** ERROR: " + text)
#         cursor = self.ui.textOutputWidget.textCursor()
#         cursor.movePosition(QtGui.QTextCursor.End)
#         cursor.insertText(text)
#         self.ui.textOutputWidget.ensureCursorVisible()
        
        
    
    