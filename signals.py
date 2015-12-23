__author__ = 'lskrinjar'
from PyQt4 import QtCore


class StepSignal(QtCore.QObject):
    signal_step = QtCore.pyqtSignal(int, name='signal_step')

class RepaintGLSignal(QtCore.QObject):
    signal_repaintGL = QtCore.pyqtSignal()

class StatusSignal(QtCore.QObject):
    signal_status = QtCore.pyqtSignal(str, name='status')

class SaveScreenshot(QtCore.QObject):
    signal_saveScreenshot = QtCore.pyqtSignal()

class SolutionFilenameSignal(QtCore.QObject):
    signal_filename = QtCore.pyqtSignal(str, name='filename')

class EnergySignal(QtCore.QObject):
    signal_energy = QtCore.pyqtSignal(float, float, name='energy')

class SaveSignal(QtCore.QObject):
    signal_save = QtCore.pyqtSignal()

class ErrorTimeIntegrationSignal(QtCore.QObject):
    signal_time_integration_error = QtCore.pyqtSignal()

class SolutionSignal(QtCore.QObject):
    """
    :param int(long):   ID of solution data object
    :param str:         solution data filename
    """
    solution_data = QtCore.pyqtSignal(int, str)