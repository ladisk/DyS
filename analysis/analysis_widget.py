"""
Created on 27. jan. 2014

@author: lskrinjar
"""
import inspect
import os
import subprocess
from pprint import pprint

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from MBD_system.body.body_widget import BodyWidget
from MBD_system.contact.contact_widget import ContactWidget
from MBD_system.force.force_widget import ForceWidget
from MBD_system.joint.joint_widget import JointWidget
from MBD_system.solution_data.solution_data import SolutionData
from MBD_system.spring.spring_widget import SpringWidget
from analysis.analysis_table_model import AnalysisTableModel
from analysis.analysis_tree_model import AnalysisTreeModel
from analysis_ui import Ui_Form


class solutionFilenameSignal(QtCore.QObject):
    signal_filename = QtCore.pyqtSignal(str, name='')

class loadSolutionFile(QtCore.QObject):
    signal_loadSolutionFile = QtCore.pyqtSignal()

class AnalysisWidget(QtGui.QWidget):  # QMainWindow#, QAbstractItemView, QWidget
    """
    classdocs
    """
    def __init__(self, MBD_system_=None, parent=None):
        """
        Constructor
        """
        super(AnalysisWidget, self).__init__(parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        self.filename_signal = solutionFilenameSignal()
        self.load_loadSolutionFile = loadSolutionFile()
        
        self._parent = parent

        #   progress bar to zero
        self.ui.progressBar.setValue(0)

        #    set validators to limit user input
        __validator_dbl = QtGui.QDoubleValidator()
        __validator_int = QtGui.QIntValidator()

        #    number of simulations
        self.ui.number_of_simulations.setValidator(__validator_int)
        self.ui.number_of_simulations.setText(str("1E+1"))

        self.MBD_system = MBD_system_

        self._data = self._build_solution_data_list()



        # _model = MBD_system.TreeModel(self.MBD_system)
        # _model = TreeModel(self.MBD_system)
        
        
        #    run method list_parameters() of MBD_system item object to create a tree structure with item properties
        self.MBD_system._children[0].parameters_construct_tree_view()
        #    tree model view
        #    set model
        _model = AnalysisTreeModel(self.MBD_system._children[0].parameters, parent=self)
        #    set model of tree view
        self.ui.treeView.setModel(_model)
        
        #    tree view settings
        self.ui.treeView.setRootIsDecorated(True)
        self.ui.treeView.setHeaderHidden(False)
        self.ui.treeView.expandAll()
        self.ui.treeView.setColumnWidth(0, 300)
        # self.ui.treeView.resizeColumnsToContents()
        
        #    list model view
        #    set model
        # data = ["luka", "peter", "miha", "eva"]
        _model = AnalysisTableModel(self._data, parent=self)#self.MBD_system._children[0].parameters_selected
        #    set model of table view
        self.ui.tableView.setModel(_model)
        for row in xrange(0, max(len(self._data), int(float(self.ui.number_of_simulations.text())))):
            self.ui.tableView.setRowHeight(row, 17)
#         self.ui.listView.resizeRowsToContents()

        # we want our listview to have a context menu taken from the actions on this widget
        # those actions will be to delete an item :)
        # self.ui.treeView.customContextMenuRequested.connect(self.contextMenuEvent)
        QtCore.QObject.connect(self.ui.treeView, QtCore.SIGNAL("clicked (QModelIndex)"),self.row_clicked)

        #    connections
        self.ui.closeButton.clicked.connect(self.close)

        #    show widget
        self._setWindowFlags()
        # self.showMaximized()
        # self.show()

    def _build_solution_data_list(self):
        """

        """
        _data = []
        _N = int(float(self.ui.number_of_simulations.text()))
        __places = len(str(_N))
        for i in xrange(0, _N):
            num = str(i).zfill(__places)
            _name = self.MBD_system._children[0]._name+"_"+num
            _sol_data = SolutionData(_name)
            _data.append(_sol_data)

        return _data

    def row_clicked(self, index):
        '''
        when a row is clicked... show the name
        '''
        None

    def _setWindowFlags(self):
        """
        Set window flags
        """
        flags = self.windowFlags()
        # flags = Qt.Window
        # self.setWindowFlags(QtCore.Qt.WindowSystemMenuHint |
        #                       QtCore.Qt.WindowMinMaxButtonsHint)
        self.setWindowFlags(self.windowFlags())
        # self.setWindowFlags(self.windowFlags() | Qt.WindowMaximizeButtonHint | Qt.WindowSystemMenuHint)
    #     # self.setWindowFlags(QtCore.Qt.Window | Qt.WindowMaximizeButtonHint)
    #     self.setWindowFlags(Qt.WindowMinimizeButtonHint | ~Qt.WindowMaximizeButtonHint)
    #     self.setWindowFlags(self.windowFlags() |
    #                         QtCore.Qt.WindowSystemMenuHint |
    #                         QtCore.Qt.WindowMinMaxButtonsHint)
        # super(TreeViewWidget, self).setWindowFlags(flags)

    
    @pyqtSlot()
    def onTriggered(self, event):
        """
        
        """
        # tell our model to remove the selected row.
#         self.ui.treeView.model()
  
#         self.model().removeRows(self.currentIndex().row(), 1)
    def create_action(self, event):
        print event.pos()

    def setSelection(self, current, old):

        current = self._proxyModel.mapToSource(current)

    def closeEvent(self, event):
        """

        """
        pass
        # self._parent.SimulationControlWidget.OpenGLWidget._repaintGL()

    def _list_parameters(self):
        """
        
        """
        self._item._list_parameters(item=self.MBD_system)

    
    def _load_group_items(self):
        """
        
        """
        print "Under construction!", os.path.realpath(__file__)


    def _load_solution_data(self):
        """

        """
        #   this has to be put in a thread
        data = self._item._load_dat_file()
        self._parent.SimulationControlWidget.load_solution_file(filename = self._item._name, solution_data=data)


    def _load_solution_file_to_project(self):
        """
        Load solution file to opened project if there is any solution file.
        TODO - check if right solution data is loaded to project
        """
        #   filetype filter
        _filter = "Data File (*.dat);;Excel (*.xlsx);;CSV (*.csv)"

        #   open dialog
        _open_file = QtGui.QFileDialog()
        _open_file.setFileMode(QFileDialog.ExistingFiles)
        
        #   directory to open
        _directory = self._item._parent._parent._children[0].MBD_folder_abs_path_
        #    set directory from where to open
        _open_file.setDirectory(_directory)

        #    file abs path and filetype(extension)
        file, filetype = _open_file.getOpenFileNameAndFilter(self, "Load Solution File", _directory, _filter)
        
        #    get filename from abspath
        filename = os.path.basename(str(file))
        
        #   filename without extension
        name, _extension = filename.split(".")

        #   create solution data object
        _solution = SolutionData(name, file=file)

        #   add solution data object item to treeview
        self.ui.treeView.model().insertRow(len(self._item._children), _solution, self._item)


    def _open(self):
        """
        
        """
        subprocess.call(['notepad.exe', str(self.__filename)])

    def add_solution_data(self, filename):
        """
        
        """
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        #   root index
        root_index = self.ui.treeView.rootIndex()

        #   get root node with root index
        self._item = self.ui.treeView.model().getNode(root_index)
        MBD_system_item = self._item._children[0]
        
        #    solution item
        _solution_group_item = MBD_system_item._children[6]

        #    add child (solution data) to solution item
        pos = len(_solution_group_item._children)
        self._childNode = SolutionData(filename)
        self.ui.treeView.model().insertRow(pos, self._childNode, _solution_group_item)


    def _create(self, index):
        """
        
        """
        self._parentNode = self._item
        self._childNode = None
        if self._item._typeInfo == "group":
            
            if self._item._name.lower() == "bodies":
                self._widget = BodyWidget(self._parentNode, parent=self)

            if self._item._name.lower() == "forces":
                self._widget = ForceWidget(self._parentNode, parent=self)

            if self._item._name.lower() == "joints":
                self._widget = JointWidget(self._parentNode, parent=self)

            if self._item._name.lower() == "contacts":
                self._widget = ContactWidget(self._parentNode, parent=self)

            if self._item._name.lower() == "springs":
                self._widget = SpringWidget(self._parentNode, parent=self)
                # self._childNode = MBD_system_items.SpringItem("Spring_" + str(len(self._parentNode._children)))



            # if self._childNode is not None:
            #     print "created =", self._childNode
            #     if len(self._parentNode._children) == 0:
            #         pos = 1
            #     else:
            #         pos = len(self._parentNode._children)
            #     self.ui.treeView.model().insertRow(pos, self._childNode, self._parentNode)
    
    
    def _insertTreeItem(self):
        """
        
        """
        print "????"
        pos = len(self._parentNode._children)
        self._childNode = self._widget.item
        self.ui.treeView.model().insertRow(pos, self._childNode, self._parentNode)
    

    def __getDOF(self):
        """
        
        """
        C_q, C_qT = self._parent.SimulationControlWidget.solver.solveODE.ode_fun.getDOF()


    def _get_q(self):
        """

        """
        q = self.MBD_system._children[0].get_q()
        print "q ="
        print q


    def _evaluate_C(self):
        """

        """
        t = 0
        q_ = self.MBD_system._children[0].get_q()
        print "C(q, t)"
        print self._parent.SimulationControlWidget.solver.solveODE.ode_fun.create_C(t, q_)


    def _edit(self, index):
        """

        """
        self._parentNode = self._item._parent

        if self._item._typeInfo.lower() == "body":
            self._widget = BodyWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "contact":
            self._widget = ContactWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "force":
            self._widget = ForceWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "joint":
            self._widget = JointWidget(self._parentNode, parent=self)

        if self._item._typeInfo.lower() == "spring":
            self._widget = SpringWidget(self._parentNode, parent=self)

        try:
            self._widget._edit(self._item)
        except:
            pass


    def _delete(self):
        """

        """
        pos = self._item._parent._children.index(self._item)
        self._item._parent.removeChild(pos)
        del(self._item)


    def _properties(self, index):
        print "_properties()"
        pprint(vars(self._item))
        
    
    def _save_contact_solution(self, index):
        """
        
        """
        self._item.save_solution_data()