"""
Created on 18. mar. 2016

@author: skrinjar.luka@gmail.com

"""
from IPython.qt.console.rich_ipython_widget import RichIPythonWidget
from IPython.terminal.interactiveshell import TerminalInteractiveShell
from IPython.qt.inprocess import QtInProcessKernelManager
from PyQt4 import QtGui, QtCore


class PythonConsole(RichIPythonWidget):

    def __init__(self, **kwarg):
        super(RichIPythonWidget, self).__init__()
        self.kernel_manager = QtInProcessKernelManager()
        self.kernel_manager.start_kernel()
        self.kernel = self.kernel_manager.kernel
        self.kernel.gui = 'qt4'
        self.kernel.shell.push(kwarg)
        self.kernel_client = self.kernel_manager.client()
        self.kernel_client.start_channels()