__author__ = 'lskrinjar'
from PyQt4 import QtCore


class StepSignal(QtCore.QObject):
    signal_step = QtCore.pyqtSignal(int, name='signal_step')

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

class SignalSimulationStatus(QtCore.QObject):
    signal_simulation_status = QtCore.pyqtSignal(str)

class SignalUpdateGUI(QtCore.QObject):
    signal_simulation_status = QtCore.pyqtSignal()

class RefreshSignal(QtCore.QObject):
    signal_refresh = QtCore.pyqtSignal()

class FinishedSignal(QtCore.QObject):
    signal_finished = QtCore.pyqtSignal()

class CreateAnimationFile(QtCore.QObject):
    signal_createAnimationFile = QtCore.pyqtSignal()

class LoadSolutionFile(QtCore.QObject):
    signal_loadSolutionFile = QtCore.pyqtSignal()