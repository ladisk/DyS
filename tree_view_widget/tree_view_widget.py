"""
Created on 27. jan. 2014

@author: lskrinjar
"""
import inspect
import os
import subprocess
from pprint import pprint

import numpy as np
from matplotlib import pyplot as plt
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *


from MBD_system import read_and_write
from MBD_system.MBD_system_widget import MBDSystemWidget
from MBD_system.body.body_widget import BodyWidget
from MBD_system.contact.contact_widget import ContactWidget
from MBD_system.force.force_widget import ForceWidget
from MBD_system.item_widget import ItemWidget
from MBD_system.joint.joint_widget import JointWidget
from MBD_system.motion.motion_widget import MotionWidget
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.solution_data.solution_data import SolutionData
from MBD_system.spring.spring_widget import SpringWidget
from analysis.analysis_widget import AnalysisWidget
from tree_model import TreeModel
from tree_view_widget_ui import Ui_Form
from global_variables import GlobalVariables
from signals import SolutionFilenameSignal
from signals import CreateAnimationFile
from signals import LoadSolutionFile
    

class TreeViewWidget(QWidget):  # QMainWindow#, QAbstractItemView 
    """
    classdocs
    """

    def __init__(self, project_item, MBD_system, parent=None, flags=0):
        """
        Constructor
        """
        super(TreeViewWidget, self).__init__(parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)

        #   signals
        self.filename_signal = SolutionFilenameSignal()
        self.load_LoadSolutionFile = LoadSolutionFile()
        self.create_animation_file = CreateAnimationFile()

        #   parent
        self._parent = parent

        self.project_item = project_item
        # _model = MBD_system.TreeModel(self.MBD_system)
        _model = TreeModel(self.project_item)

        self.MBD_system = MBD_system

        self.ui.treeView.setModel(_model)

        #   predefined attributes
        self._item = None
        self._selected_item = None

        #   settings
        self.ui.treeView.setRootIsDecorated(True)
        self.ui.treeView.setHeaderHidden(True)
        self.ui.treeView.expandAll()
        self.highlightSelectedItem = True

        # we want our listview to have a context self.menu taken from the actions on this widget
        # those actions will be to delete an item :)
        self.ui.treeView.customContextMenuRequested.connect(self.contextMenuEvent)
        self.ui.treeView.clicked.connect(self.clickedItem)

        # self._create_animation_file.connect(self._parent.simulation_control_widget._create_animation_file)
 
        self._widget = ItemWidget(parent=self)

        # print "self.screenGeometry() =", self.screenGeometry()

    def setProjectItem(self, project_item, MBD_system):
        """

        :param project:
        :return:
        """
        self.MBD_system = MBD_system

        self.ui.treeView.selectionModel()
        self.ui.treeView.clearSelection()

        _model = TreeModel(project_item)
        self.ui.treeView.setModel(_model)
        # self.ui.treeView.setRootIndex(_model._rootNode)
        self.ui.treeView.expandAll()

    def setWindowFlags(self, flags):
        """

        """
        super(TreeViewWidget, self).setWindowFlags(flags)

    @pyqtSlot()
    def onTriggered(self, event):
        """
        
        """
        # tell our model to remove the selected row.
#         self.ui.treeView.model()
  
#         self.model().removeRows(self.currentIndex().row(), 1)

    def create_action(self, event):
        """

        """
        print event.pos()

    def setSelection(self, current, old):
        """

        """
        current = self._proxyModel.mapToSource(current)

    def clickedItem(self, event):
        """
        Method run when item in tree view is clicked (selected)
        :param event:
        :return:
        """
        self.menu = QtGui.QMenu(self)

        indexes = self.ui.treeView.selectedIndexes()
        #   get item
        self._item = self.ui.treeView.model().getNode(indexes[0])

        #   highlight selected item in visualization widget
        if self.highlightSelectedItem:
            #   highlight selected item
            if hasattr(self._item, "vtk_actor"):
                if self._item.vtk_actor is not None:
                    if hasattr(self._item, "highlightSelected"):
                        self._item.highlightSelected()

            elif hasattr(self._item, "GetVisibility"):
                self._item.highlightSelected()

            #   unhighlight
            if self._selected_item is not None and self._selected_item is not self._item:
                if hasattr(self._selected_item, "unHighlightSelected"):
                    self._selected_item.unHighlightSelected()

            if hasattr(self._parent.simulation_control_widget.vtkWidget, "refresh"):
                self._parent.simulation_control_widget.vtkWidget.refresh()

                self._selected_item = self._item

    def contextMenuEvent(self, event):
        """
        Context menu
        """
        self.menu = QtGui.QMenu(self)

        indexes = self.ui.treeView.selectedIndexes()

        self._item = self.ui.treeView.model().getNode(indexes[0])

        if len(indexes) > 0:
            level = 0
            index = indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1

        if self._item._typeInfo == "group":
            createAction = self.menu.addAction("Create")
            createAction.triggered.connect(self._create)
            
            self.menu.addSeparator()
            
            loadGroupItemsAction = self.menu.addAction("Load items from file")
            loadGroupItemsAction.triggered.connect(self._load_group_items)

            printAlltemsAction = self.menu.addAction("List of " + self._item._name)
            printAlltemsAction.triggered.connect(self._print_group_items)
        
            self.menu.addSeparator()

            if self._item._name == "Motions":
                pprint(vars(self._item))
            
        if self._item._typeInfo != "group" or self._item._typeInfo != "solution":# or self._item._name.lower() != "solution":
            editAction = self.menu.addAction("Edit")
            editAction.triggered.connect(self._edit)

            deleteAction = self.menu.addAction("Delete")
            deleteAction.setEnabled(True)
            deleteAction.triggered.connect(self._delete)

            # saveAction = self.menu.addAction("Save")
            # saveAction.triggered.connect(self._save)

            propertiesAction = self.menu.addAction("Properties")
            propertiesAction.triggered.connect(self.properties)

            saveAction = self.menu.addAction("Save")
            saveAction.triggered.connect(self._save)

            self.menu.addSeparator()

        if self._item._typeInfo == "MBDsystem" or self._item._typeInfo == "settings":
            self.menu.addSeparator()
            if hasattr(self._parent.simulation_control_widget.vtkWidget, "GCS"):
                show_GCS_Action = QtGui.QAction("Show GCS ", self, checkable=True, checked=self._parent.simulation_control_widget.opengl_widget.GCS._visible)
                show_GCS_Action.triggered.connect(self._parent.simulation_control_widget.opengl_widget.GCS._show)
                self.menu.addAction(show_GCS_Action)

            propertiesAction = self.menu.addAction("Properties")
            propertiesAction.triggered.connect(self._properties)

            globalVariablesAction = self.menu.addAction("Global Variables")
            globalVariablesAction.triggered.connect(self._globalVariables)
            
            self.menu.addSeparator()

            profileAction = self.menu.addAction("Profile")
            profileAction.triggered.connect(self._parent.simulation_control_widget.profile_functions)
            
            getDOFsAction = self.menu.addAction("DOF")
            getDOFsAction.triggered.connect(self.__getDOF)

            if self._item._typeInfo == "MBDsystem":
                evaluate_q_Action = self.menu.addAction("q")
                evaluate_q_Action.triggered.connect(self._evaluate_q)

                evaluate_q_Action = self.menu.addAction("q_2")
                evaluate_q_Action.triggered.connect(self._evaluate_q_2)
                
                save_q_Action = self.menu.addAction("Save q")
                save_q_Action.triggered.connect(self._save_q)
                
                load_q_Action = self.menu.addAction("Load q")
                load_q_Action.triggered.connect(self._load_q)
                
                self.menu.addSeparator()

                evaluate_M_Action = self.menu.addAction("M")
                evaluate_M_Action.triggered.connect(self.print_M)

                evaluate_C_Action = self.menu.addAction("C(q, t)")
                evaluate_C_Action.triggered.connect(self._evaluate_C)

                evaluate_C_q_Action = self.menu.addAction("C_q(q)")
                evaluate_C_q_Action.triggered.connect(self.print_C_q)

                evaluate_C_t_Action = self.menu.addAction("C_t(t)")
                evaluate_C_t_Action.triggered.connect(self._evaluate_C_t)

                evaluate_M_Action = self.menu.addAction("Q_g")
                evaluate_M_Action.triggered.connect(self.print_Q_g)

                evaluate_M_Action = self.menu.addAction("Q_e")
                evaluate_M_Action.triggered.connect(self.print_Q_e)

                evaluate_M_Action = self.menu.addAction("Q_s")
                evaluate_M_Action.triggered.connect(self.print_Q_s)

                evaluate_M_Action = self.menu.addAction("Q_q")
                evaluate_M_Action.triggered.connect(self.print_Q_q)

                evaluate_M_Action = self.menu.addAction("Q_dq")
                evaluate_M_Action.triggered.connect(self.print_Q_dq)

                evaluate_M_Action = self.menu.addAction("C_q_L_q")
                evaluate_M_Action.triggered.connect(self.print_C_q_L_q)

                evaluate_M_Action = self.menu.addAction("J")
                evaluate_M_Action.triggered.connect(self.print_J)

                self.menu.addSeparator()
                list_parameters_Action = self.menu.addAction("List parameters")
                list_parameters_Action.triggered.connect(self._list_parameters)

                analysis_Action = self.menu.addAction("Analysis")
                analysis_Action.triggered.connect(self._analysis)
                analysis_Action.setShortcut(QtGui.QKeySequence("Alt+A"))

                self.menu.addSeparator()
                
                #    forces
                viewForces_Action = self.menu.addAction("Forces")
                viewForces_Action.triggered.connect(self._viewForces)
                self.menu.addSeparator()

        elif self._item._typeInfo == "body":
            if self._item._name == "Ground" or self._item.body_id == -1:
                pass
            else:
                if hasattr(self._item.vtk_actor, "GetVisibility"):
                    if self._item.vtk_actor.GetVisibility():
                        show_body_Action = QtGui.QAction("Hide Body", self)
                        show_body_Action.triggered.connect(self._item.vtk_actor.VisibilityOff)
                    else:
                        show_body_Action = QtGui.QAction("Show Body", self)
                        show_body_Action.triggered.connect(self._item.vtk_actor.VisibilityOn)
                    self.menu.addAction(show_body_Action)

                #   show/hide LCS
                if hasattr(self._item.LCS_actor, "GetVisibility"):
                    if self._item.LCS_actor.GetVisibility():
                        show_LCS_Action = QtGui.QAction("Hide LCS", self)
                        show_LCS_Action.triggered.connect(self._item.LCS_actor.VisibilityOff)
                    else:
                        show_LCS_Action = QtGui.QAction("Show LCS", self)
                        show_LCS_Action.triggered.connect(self._item.LCS_actor.VisibilityOn)
                    self.menu.addAction(show_LCS_Action)

                #   show/hide geometry CS
                if hasattr(self._item.geometry_CS_actor, "GetVisibility"):
                    if self._item.geometry_CS_actor.GetVisibility():
                        show_geometry_CS_Action = QtGui.QAction("Hide geometry CS", self)
                        show_geometry_CS_Action.triggered.connect(self._item.geometry_CS_actor.VisibilityOff)
                    else:
                        show_geometry_CS_Action = QtGui.QAction("Show geometry CS", self)
                        show_geometry_CS_Action.triggered.connect(self._item.geometry_CS_actor.VisibilityOn)
                    self.menu.addAction(show_geometry_CS_Action)

                if self._item.body_type == "flexible body":
                    if self._item._visible_elements:
                        show_mesh_elements_Action = QtGui.QAction("Hide Elements", self)
                    else:
                        show_mesh_elements_Action = QtGui.QAction("Show Elements", self)

                    show_mesh_elements_Action.triggered.connect(self._item._element_visibility)
                    self.menu.addAction(show_mesh_elements_Action)

                    if hasattr(self._item.mesh.vtk_actor, "GetVisibility"):
                        if self._item.mesh.vtk_actor.GetVisibility():
                            show_mesh_nodes_Action = QtGui.QAction("Hide Nodes", self)
                            show_mesh_nodes_Action.triggered.connect(self._item.mesh.vtk_actor.VisibilityOff)
                        else:
                            show_mesh_nodes_Action = QtGui.QAction("Show Nodes", self)
                            show_mesh_nodes_Action.triggered.connect(self._item.mesh.vtk_actor.VisibilityOn)

                        self.menu.addAction(show_mesh_nodes_Action)

                    #   print nodes
                    self.menu.addSeparator()
                    print_mesh_nodes_Action = QtGui.QAction("Print Nodes", self)
                    print_mesh_nodes_Action.triggered.connect(self.print_nodes)
                    self.menu.addAction(print_mesh_nodes_Action)
                    self.menu.addSeparator()

                if hasattr(self._item, "CAD_CS"):
                    show_CAD_CS_Action = QtGui.QAction("Show CAD CS", self, checkable=True, checked=self._item.CAD_CS._visible)
                    show_CAD_CS_Action.triggered.connect(self._item.CAD_CS._show)
                    self.menu.addAction(show_CAD_CS_Action)
                    if self._item.CAD_CS._VBO_created:
                        show_CAD_CS_Action.setEnabled(True)
                    else:
                        show_CAD_CS_Action.setEnabled(False)

                show_ID_Action = QtGui.QAction("Show ID", self, checkable=True, checked=self._item._idVisible)
                show_ID_Action.triggered.connect(self._item._showID)
                self.menu.addAction(show_ID_Action)

                self.menu.addSeparator()

                getBody_q_Action = self.menu.addAction("q")
                getBody_q_Action.triggered.connect(self._item.print_q)

                getBody_dq_Action = self.menu.addAction("dq")
                getBody_dq_Action.triggered.connect(self._item.print_dq)

                getBody_qdq_Action = self.menu.addAction("[q dq]")
                getBody_qdq_Action.triggered.connect(self._item.get_qdq)
                self.menu.addSeparator()

                printBody_M_Action = self.menu.addAction("M")
                printBody_M_Action.triggered.connect(self._item.print_M)

                if self._item.body_type == "flexible body":
                    printK_Action = self.menu.addAction("K")
                    printK_Action.triggered.connect(self._item.print_K)

                    print_Q_g_Action = self.menu.addAction("Q_g")
                    print_Q_g_Action.triggered.connect(self._item.print_Q_g)

                self.menu.addSeparator()

        elif self._item._typeInfo == "solution":
            loadAction = self.menu.addAction("Load from File")
            loadAction.triggered.connect(self._load_solution_file_to_project)

        #   solution data item context menu
        elif self._item._typeInfo == "solutiondata":
            loadAction = self.menu.addAction("Load solution")
            loadAction.triggered.connect(self._load_solution_data)

            saveSolutionAction = self.menu.addAction("Save solution to file")
            saveSolutionAction.triggered.connect(self._save_solution)


            propertiesAction = self.menu.addAction("Properties")
            propertiesAction.triggered.connect(self.solution_properties)

            if self._item.loaded:
                saveAnimationFileAction = self.menu.addAction("Save animation to file")
                saveAnimationFileAction.triggered.connect(self._create_animation_file)

        elif self._item._typeInfo == "group":
            propertiesAction = self.menu.addAction("Properties")
            propertiesAction.triggered.connect(self._group_properties)

#         elif self._item._typeInfo == "solutiondata":
#             
#             openAction = self.menu.addAction("Open")
#             openAction.triggered.connect(self._open)
#             
#             self.menu.addSeparator()
# 
#             deleteAction = self.menu.addAction("Delete")
#             deleteAction.triggered.connect(self._delete)
            
        elif self._item._typeInfo == "force":
            name = self.MBD_system.bodies[self._item.body_id]._name
            show_Action = QtGui.QAction("Show force on Body i id = " + str(self._item.body_id) + " ("+name+")", self, checkable=True, checked=self._item._visible)
            show_Action.triggered.connect(self._item._show)
            self.menu.addAction(show_Action)

            self.menu.addSeparator()

            if hasattr(self._item, "vtk_actor"):
                if self._item.vtk_actor is not None:
                    if self._item.vtk_actor.GetVisibility():
                        show_LCS_Action = QtGui.QAction("Hide Force", self)
                        show_LCS_Action.triggered.connect(self._item.vtk_actor.VisibilityOff)
                    else:
                        show_LCS_Action = QtGui.QAction("Show Force", self)
                        show_LCS_Action.triggered.connect(self._item.vtk_actor.VisibilityOff)

                    self.menu.addSeparator()

            evaluate_Q_e_Action = QtGui.QAction("Evaluate Q_e", self)
            evaluate_Q_e_Action.triggered.connect(self.evaluate_Q_e)
            self.menu.addAction(evaluate_Q_e_Action)

        elif self._item._typeInfo == "joint" or self._item._typeInfo == "spring":
            if self._item._typeInfo == "spring":
                if self._item.vtk_actor is not None:
                    if self._item.vtk_actor.GetVisibility():
                        show_spring_Action = QtGui.QAction("Hide Spring", self)
                        show_spring_Action.triggered.connect(self._item.vtk_actor.VisibilityOff)
                    else:
                        show_spring_Action = QtGui.QAction("Show Spring", self)
                        show_spring_Action.triggered.connect(self._item.vtk_actor.VisibilityOn)

                    self.menu.addAction(show_spring_Action)

            for force, ID, body_id in zip(self._item._Fn_list, ["i", "j"], self._item.body_id_list):
                if isinstance(body_id, int):
                    name = self.MBD_system.bodies[body_id]._name
                else:
                    name = "ground"
                show_Action = QtGui.QAction("Show force on Body " + ID + " id = "+str(body_id) + " ("+name+")", self, checkable=True, checked=force._visible)
                show_Action.triggered.connect(force._show)
                self.menu.addAction(show_Action)

            self.menu.addSeparator()
            for marker, ID, body_id in zip(self._item.markers, ["i", "j"], self._item.body_id_list):
                if isinstance(body_id, int):
                    name = self.MBD_system.bodies[body_id]._name
                else:
                    name = "ground"

                if marker.GetVisibility():
                    show_marker_Action = QtGui.QAction("Hide u_P on body " + ID + " id = "+str(body_id) + " ("+name+")", self)
                    show_marker_Action.triggered.connect(marker.VisibilityOff)
                else:
                    show_marker_Action = QtGui.QAction("Show u_P on body " + ID + " id = " + str(body_id) + " (" + name + ")", self)
                    show_marker_Action.triggered.connect(marker.VisibilityOn)
                self.menu.addAction(show_marker_Action)

            if self._item._typeInfo == "joint":
                if self._item.joint_type == "prismatic":
                    marker = self._item.markers[-1]
                    show_marker_uQ_Action = QtGui.QAction("Show u_Q on body "+"i", self, checkable=True, checked=marker._visible)
                    show_marker_uQ_Action.triggered.connect(marker._show)
                    self.menu.addAction(show_marker_uQ_Action)

            self.menu.addSeparator()
            if self._item._typeInfo == "joint":
                evaluate_C_Action = QtGui.QAction("Evaluate C(q, t)", self)
                evaluate_C_Action.triggered.connect(self.evaluate_C)
                self.menu.addAction(evaluate_C_Action)

                evaluate_C0_Action = QtGui.QAction("Evaluate C0(q, t)", self)
                evaluate_C0_Action.triggered.connect(self.evaluate_C0)
                self.menu.addAction(evaluate_C0_Action)

                evaluate_C_q_Action = QtGui.QAction("Evaluate C_q(q)", self)
                evaluate_C_q_Action.triggered.connect(self.evaluate_C_q)
                self.menu.addAction(evaluate_C_q_Action)

                evaluate_Q_d_Action = QtGui.QAction("Evaluate Q_d", self)
                evaluate_Q_d_Action.triggered.connect(self.evaluate_Q_d)
                self.menu.addAction(evaluate_Q_d_Action)

            evaluate_rijP_Action = QtGui.QAction("Evaluate d", self)
            evaluate_rijP_Action.triggered.connect(self.evaluate_d)
            self.menu.addAction(evaluate_rijP_Action)

            if self._item._typeInfo == "spring":
                evaluate_Q_e_Action = QtGui.QAction("Evaluate Q_e", self)
                evaluate_Q_e_Action.triggered.connect(self.evaluate_Q_e)
                self.menu.addAction(evaluate_Q_e_Action)

                evaluate_l0l_Action =  QtGui.QAction("Evaluate l0l", self)
                evaluate_l0l_Action.triggered.connect(self.evaluate_l0l)
                self.menu.addAction(evaluate_l0l_Action)

                if self._item.spring_type == "rotational":
                    evaluate_theta0_action = QtGui.QAction("Evaluate theta0", self)
                    evaluate_theta0_action.triggered.connect(self._evaluate_theta0)
                    self.menu.addAction(evaluate_theta0_action)

                    evaluate_theta_action = QtGui.QAction("Evaluate theta", self)
                    evaluate_theta_action.triggered.connect(self._evaluate_theta)
                    self.menu.addAction(evaluate_theta_action)

                self.menu.addSeparator()
                for force, ID, body_id in zip(self._item._Fn_list, ["i", "j"], self._item.body_id_list):
                    getData_Action = QtGui.QAction("Force solution data on body " + ID, self)
                    getData_Action.triggered.connect(force.get_data)
                    self.menu.addAction(getData_Action)

        elif self._item._typeInfo == "contact":
            for visible, ID, body_id in zip(self._item._visible_F_list, ["i", "j"], self._item.body_id_list):
                if isinstance(body_id, int):
                    name = self.MBD_system.bodies[body_id]._name
                else:
                    name = "ground"

                if visible:
                    show_Action = QtGui.QAction("Hide F on body "+ID+" id = "+str(body_id)+" ("+name+")", self)
                else:
                    show_Action = QtGui.QAction("Show F on body "+ID+" id = "+str(body_id)+" ("+name+")", self)

                if ID == "i":
                    show_Action.triggered.connect(self._item.change_visibility_Fi)
                if ID == "j":
                    show_Action.triggered.connect(self._item.change_visibility_Fj)

                self.menu.addAction(show_Action)

            # self.menu.addSeparator()
            # for force, ID, body_id in zip(self._item._Ft_list, ["i", "j"], self._item.body_id_list):
            #     if isinstance(body_id, int):
            #         name = self.MBD_system.bodies[body_id]._name
            #     else:
            #         name = "ground"
            #
            #     show_Action = QtGui.QAction("Show Ft on body "+ID+" id = "+str(body_id)+" ("+name+")", self, checkable=True, checked=force._visible)
            #     show_Action.triggered.connect(force._show)
            #     self.menu.addAction(show_Action)

            self.menu.addSeparator()
            if self._item._contact_type in ["general", "roughness profile"]:
                # try:
                for body_id, ID in zip(self._item.body_id_list, ["i", "j"]):
                    _body = self._item._bodies[body_id]
                    _AABB = _body.AABBtree

                    # show_AABB_action = QtGui.QAction("Show AABB on body "+ID, self, checkable=True, checked=_AABB._visible)
                    # show_AABB_action.triggered.connect(_AABB._show_AABB)
                    # self.menu.addAction(show_AABB_action)

                self.menu.addSeparator()

                for body_id, ID in zip(self._item.body_id_list, ["i", "j"]):
                    _body = self._item._bodies[body_id]
                    _AABB = _body.AABBtree

                    properties_of_AABB_action = self.menu.addAction("Properties of AABB on body "+ID)
                    properties_of_AABB_action.triggered.connect(_AABB._get_properties)

                self.menu.addSeparator()

                if self._item.AABB_list!=[]:
                    for body_id, AABB in zip(self._item.body_id_list, self._item.AABB_list):
                        print_AABB_nodes_LCS_action = self.menu.addAction("AABB nodes LCS "+ID)
                        print_AABB_nodes_LCS_action.triggered.connect(AABB.print_nodes_LCS)

                    for body_id, AABB in zip(self._item.body_id_list, self._item.AABB_list):
                        print_AABB_nodes_LCS_action = self.menu.addAction("AABB nodes GCS "+ID)
                        print_AABB_nodes_LCS_action.triggered.connect(AABB.print_nodes_GCS)

                    self.menu.addSeparator()

                    for AABB, ID in zip(self._item.AABB_list, ["i", "j"]):
                        if hasattr(AABB.vtk_actor, "GetVisibility"):
                            if AABB.vtk_actor.GetVisibility():
                                show_AABB_Action = QtGui.QAction("Hide AABB of Body " + ID, self)
                                show_AABB_Action.triggered.connect(AABB.vtk_actor.VisibilityOff)

                            else:
                                show_AABB_Action = QtGui.QAction("Show AABB of Body " + ID, self)
                                show_AABB_Action.triggered.connect(AABB.vtk_actor.VisibilityOn)

                            self.menu.addAction(show_AABB_Action)

            try:
                for contact_geometry, ID in zip(self._item.contact_geometry_list, ["i", "j"]):
                    # get_nodes_action = QtGui.QAction("Get contact nodes "+ID, self, checkable=True, checked=_AABB._visible)
                    get_nodes_action = self.menu.addAction("Get contact nodes "+ID)
                    get_nodes_action.triggered.connect(contact_geometry.get_nodes)
                    # self.menu.addAction(get_nodes_action)
            except:
                pass

            contactPointAction = QtGui.QAction("Contact point in GCS", self)
            contactPointAction.triggered.connect(self._contact_point)
            self.menu.addAction(contactPointAction)
            self.menu.addSeparator()

            #   evaluate contact forces
            evaluate_Q_e_Action = QtGui.QAction("Evaluate Q_e", self)
            evaluate_Q_e_Action.triggered.connect(self.evaluate_Q_e)
            self.menu.addAction(evaluate_Q_e_Action)

            #   evaluate contact
            evaluate_contact_Action = QtGui.QAction("Evaluate Contact", self)
            evaluate_contact_Action.triggered.connect(self.evaluate_contact)
            self.menu.addAction(evaluate_contact_Action)

            #   show properties of contact points
            print_contact_point_object_properties_Action = QtGui.QAction("Properties of Contact Point", self)
            print_contact_point_object_properties_Action.triggered.connect(self._item.contact_point_object_properties)
            self.menu.addAction(print_contact_point_object_properties_Action)

            if self._item._contact_type.lower() == "pin-slot clearance joint linear":
                evaluate_pin_in_slot_position_Action = QtGui.QAction("Evaluate pin in slot position", self)
                evaluate_pin_in_slot_position_Action.triggered.connect(self.evaluate_pin_in_slot_position)
                self.menu.addAction(evaluate_pin_in_slot_position_Action)

                evaluate_delta_Action = QtGui.QAction("Evaluate delta", self)
                evaluate_delta_Action.triggered.connect(self.evaluate_delta)
                self.menu.addAction(evaluate_delta_Action)

            self.menu.addSeparator()

            if self._item._contact_type.lower() == "revolute clearance joint":
                for marker, ID in zip(self._item.markers, ["i", "j"]):
                    show_marker_Action = QtGui.QAction("Show u_P on Body " + ID, self, checkable=True, checked=marker._visible)
                    show_marker_Action.triggered.connect(marker._show)
                    self.menu.addAction(show_marker_Action)
                    
                self.menu.addSeparator()
                evaluate_rijP_Action = QtGui.QAction("Evaluate rijP", self)
                evaluate_rijP_Action.triggered.connect(self.evaluate_d)
                self.menu.addAction(evaluate_rijP_Action)
                
                evaluate_delta_Action = QtGui.QAction("Evaluate delta", self)
                evaluate_delta_Action.triggered.connect(self.evaluate_delta)
                self.menu.addAction(evaluate_delta_Action)

            if self._item._contact_type.lower() == "pin-slot clearance joint linear":
                print_r_CP_list_Action = QtGui.QAction("Print r CP list (Center Points)", self)
                print_r_CP_list_Action.triggered.connect(self.print_r_CP_list)
                self.menu.addAction(print_r_CP_list_Action)

                print_r_slot_frame_r_Action = QtGui.QAction("Print r slot frame", self)
                print_r_slot_frame_r_Action.triggered.connect(self.print_r_slot_frame)
                self.menu.addAction(print_r_slot_frame_r_Action)

                self.menu.addSeparator()

                for marker, ID in zip(self._item.markers, ["i", "i", "j"]):
                    if hasattr(marker, "GetVisibility"):
                        if marker.GetVisibility():
                            show_marker_Action = QtGui.QAction("Hide Body "+ID+" joint center", self)
                            show_marker_Action.triggered.connect(marker.VisibilityOff)
                        else:
                            show_marker_Action = QtGui.QAction("Show Body " + ID + " joint center", self)
                            show_marker_Action.triggered.connect(marker.VisibilityOn)
                    self.menu.addAction(show_marker_Action)

            if self._item._contact_type.lower() == "contact point-line":
                for marker, ID in zip(self._item.markers, ["i", "j", "j"]):
                    show_marker_Action = QtGui.QAction("Show u_P on Body " + ID + " id = " + str(marker._parent.body_id) + " (" + marker._parent._name + ")", self, checkable=True, checked=marker._visible)
                    show_marker_Action.triggered.connect(marker._show)
                    self.menu.addAction(show_marker_Action)

            if self._item._contact_type.lower() == "roughness profile":
                for profile, ID in zip(self._item.roughness_profile_list, ["i", "j"]):
                    if hasattr(profile.geometry.vtk_actor, "GetVisibility"):
                        if profile.geometry.vtk_actor.GetVisibility():
                            show_profile_Action = QtGui.QAction("Hide profile of Body " + ID + " id = " + str(profile.body_id) + " (" + profile.body._name + ")", self)
                            show_profile_Action.triggered.connect(profile.geometry.vtk_actor.VisibilityOff)

                        else:
                            show_profile_Action = QtGui.QAction("Show profile of Body " + ID + " id = " + str(profile.body_id) + " (" + profile.body._name + ")", self)
                            show_profile_Action.triggered.connect(profile.geometry.vtk_actor.VisibilityOn)
                        self.menu.addAction(show_profile_Action)

            if hasattr(self._item, "AABB_list"):
                if self._item.AABB_list != []:
                    self.menu.addSeparator()
                    self.menu.addMenu(self.contextSubMenuEventPlotAABB())
                    self.menu.addSeparator()

            contactForceAction = self.menu.addAction("Force solution data")
            contactForceAction.triggered.connect(self._item.get_contact_force)

        elif self._item._typeInfo == "motion":
            evaluate_C_Action = QtGui.QAction("Evaluate C(q, t)", self)
            evaluate_C_Action.triggered.connect(self.evaluate_C)
            self.menu.addAction(evaluate_C_Action)

            evaluate_C_q_Action = QtGui.QAction("Evaluate C_q(q)", self)
            evaluate_C_q_Action.triggered.connect(self.evaluate_C_q)
            self.menu.addAction(evaluate_C_q_Action)

            evaluate_C_t_Action = QtGui.QAction("Evaluate C_t(t)", self)
            evaluate_C_t_Action.triggered.connect(self.evaluate_C_t)
            self.menu.addAction(evaluate_C_t_Action)
        else:
            pass

        if self._item._typeInfo in ["force", "spring"]:
            evaluate_F_Action = QtGui.QAction("Evaluate F", self)
            evaluate_F_Action.triggered.connect(self._item.print_F)

            self.menu.addAction(evaluate_F_Action)

        if self._item._typeInfo in "marker":
            if self._item.GetVisibility():
                show_marker_Action = QtGui.QAction("Hide marker", self)
                show_marker_Action.triggered.connect(self._item.VisibilityOff)

            else:
                show_marker_Action = QtGui.QAction("Show marker", self)
                show_marker_Action.triggered.connect(self._item.VisibilityOn)
            self.menu.addAction(show_marker_Action)

        if hasattr(self._item, "solution_data"):
            print_solution_data_Action = QtGui.QAction("Print solution data", self)
            print_solution_data_Action.triggered.connect(self._item.solution_data.print_solution_containers)
            self.menu.addAction(print_solution_data_Action)

        testingAction = self.menu.addAction("TESTING")
        testingAction.triggered.connect(self.testing)

        self.menu.exec_(event.globalPos())

        if hasattr(self._parent.simulation_control_widget.vtkWidget, "refresh"):
            self._parent.simulation_control_widget.vtkWidget.refresh()

    def contextSubMenuEventPlotAABB(self):
        """
        Subcontext menu for plot of AABB
        :param event:
        :return:
        """
        sub_menu_plot_AABB = QtGui.QMenu(self)
        sub_menu_plot_AABB.setTitle("Plot")

        sub_menu_plot_AABB.addSeparator()

        #   body i
        plot_AABBi_Action = QtGui.QAction("Plot 2D AABB of Body i", self)
        plot_AABBi_Action.triggered.connect(self.plot_AABBi_2D)
        sub_menu_plot_AABB.addAction(plot_AABBi_Action)

        plot_AABBi_tree_Action = QtGui.QAction("Plot 2D AABB Tree of Body i", self)
        plot_AABBi_tree_Action.triggered.connect(self.plot_AABBi_tree_2D)
        sub_menu_plot_AABB.addAction(plot_AABBi_tree_Action)

        sub_menu_plot_AABB.addSeparator()

        #   body j
        plot_AABBj_Action = QtGui.QAction("Plot 2D AABB of Body j", self)
        plot_AABBj_Action.triggered.connect(self.plot_AABBj_2D)
        sub_menu_plot_AABB.addAction(plot_AABBj_Action)

        plot_AABBj_tree_Action = QtGui.QAction("Plot 2D AABB Tree of Body j", self)
        plot_AABBj_tree_Action.triggered.connect(self.plot_AABBj_tree_2D)
        sub_menu_plot_AABB.addAction(plot_AABBj_tree_Action)

        sub_menu_plot_AABB.addSeparator()

        return sub_menu_plot_AABB

    def testing(self):
        self._item.testing()

    def closeEvent(self, event):
        """

        """
        self._parent.simulation_control_widget.vtkWidget.refresh()

    def _list_parameters(self):
        """
        
        """
        self._item._list_parameters(item=self.MBD_system)

    def _analysis(self):
        """
        
        """
        self._widget = AnalysisWidget(self.MBD_system)
        self._widget.show()

    def _load_group_items(self):
        """
        
        """
        print "Under construction!", os.path.realpath(__file__)

    def _print_group_items(self):
        """

        :return:
        """
        if hasattr(self.MBD_system, self._item._name.lower()):
            print "List of", self._item._name
            items = getattr(self.MBD_system, self._item._name.lower())

            for i, item in enumerate(items):
                print i, item._name, item

    def _load_solution_data(self):
        """

        """
        #   this has to be put in a thread
        if hasattr(self._item, "filename"):
            _filename = self._item.filename
        else:
            _filename = None

        _id = id(self._item)

        self._parent.simulation_control_widget.load_solution_file(solution_object_id=_id, filename=_filename)#_name

    def _load_solution_file_to_project(self):
        """
        Load solution file to opened project if there is any solution file.
        TODO - check if right solution data is loaded to project
        """
        #   filetype filter
        # _filter = "Data File (*.dat);;Excel (*.xlsx);;CSV (*.csv)"
        _filter = SolutionData.qfiletypes()

        #   open dialog
        _open_file = QtGui.QFileDialog()
        _open_file.setFileMode(QFileDialog.ExistingFiles)
        
        #   directory to open
        _directory = self._item._parent._parent._children[0].MBD_folder_abs_path
        #    set directory from where to open
        _open_file.setDirectory(_directory)

        #    file abs path and filetype(extension)
        _file, filetype = _open_file.getOpenFileNameAndFilter(self, "Load Solution File", _directory, _filter)

        if os.path.isfile(_file):
            self.setCursor(Qt.CursorShape(Qt.WaitCursor))
            #    get filename from abspath
            filename = os.path.basename(str(_file))

            #   filename without extension
            name, _extension = filename.split(".")

            #   create solution data object
            _solution = SolutionData(name=name, _file=str(_file))
            _solution.read_file()
            # pprint(vars(_solution))

            #   add solution data object item to treeview
            self.ui.treeView.model().insertRow(len(self._item._children), _solution, self._item)

            #   add solution to list of solutions
            self.MBD_system.solutions.append(_solution)
            self.MBD_system.loaded_solution = _solution
            self.setCursor(Qt.CursorShape(Qt.ArrowCursor))

    def _open(self):
        """
        
        """
        subprocess.call(['notepad.exe', str(self.__filename)])

    def add_solution_data(self, solution_data_object_ID=None, filename=None):
        """
        Function adds created solution object item in tree view in folder solutions
        """
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)

        #   root index
        root_index = self.ui.treeView.rootIndex()

        #   get root node with root index
        self._item = self.ui.treeView.model().getNode(root_index)
        MBD_system_item = self._item._children[0]
        # pprint(vars(MBD_system_item))
        #    solutions (group) item
        # _solution_group_item = MBD_system_item._children[6]


        _solution_group_item = [child for child in MBD_system_item._children if child._name == "Solution"][0]

        #    add child (solution data) to solution item
        pos = len(_solution_group_item._children)
        _childNode = MBD_system_item.solutions[-1]
        self.ui.treeView.model().insertRow(pos, _childNode, _solution_group_item)

    def solution_properties(self):
        pprint(vars(self._item))

    def _create(self, index):
        """
        
        """
        self._parentNode = self._item
        self._childNode = None

        if self._item._typeInfo == "group":
            if self._item._name.lower() == "bodies":
                self._widget = BodyWidget(self._parentNode, parent=self)#self, self._parent

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
        pos = len(self._parentNode._children)
        self._childNode = self._widget.item
        self.ui.treeView.model().insertRow(pos, self._childNode, self._parentNode)

    def __getDOF(self):
        """
        
        """
        C_q, C_qT = self._parent.simulation_control_widget.solver.solveODE.ode_fun.getDOF()

    def _evaluate_q_2(self):
        """

        """
        print self.MBD_system.evaluate_q_2()

    def _evaluate_q(self):
        """

        :return:
        """
        q = self.MBD_system.evaluate_q()
        print "q of MBD system"
        print "size =", len(q)
        print q
    
    def _load_q(self):
        """
        Function loads vector q from solution file to MBD system object and its bodies
        """
        _path = self.MBD_system.MBD_folder_abs_path_
        _load_file = QtGui.QFileDialog(self)
        _load_file.setDirectory(_path)

        _supported_types = self.__supported_types_solution()

        #   get file path from open file dialog
        file_path = _load_file.getOpenFileName(parent=self._parent, caption=QString("Load q"), directory=QString(_path), filter=QString(_supported_types))
        file_path.replace("/", "\\")

        #   load file data if file is selected
        if file_path:
            #   read data file and store data to solution data object
            solution_data = SolutionData(_file=file_path, MBD_system=self.MBD_system)

            q = self.MBD_system.q0 = solution_data._q_sol_container

            for body in self.MBD_system.bodies:
                #   q to body
                body.R[0:2] = q2R_i(q, body.body_id)
                #   dq to body
                body.theta[2] = q2theta_i(q, body.body_id) 

    def _save_q(self):
        """
        Function saves curent vector q to solution file
        """
        _path = self.MBD_system.MBD_folder_abs_path_
        _save_file = QtGui.QFileDialog(self)
        _save_file.setDirectory(_path)
        
        #    get supported types as string
        _supported_types = self.__supported_types_solution()
        
        #    filepath and filetype
        file_path, filetype = _save_file.getSaveFileNameAndFilter(parent=self._parent, caption=QString("Save q"), directory=QString(_path), filter=QString(_supported_types))
        
        #    get current vector q of MBD system
        q = self.MBD_system.evaluate_q()
        
        #    create solution data object to store vector q to file
        solution_data = SolutionData(name=None, _file=str(file_path), MBD_system=self.MBD_system)
       
        _data = np.array([np.concatenate((np.zeros(5), q), axis=0)])
        solution_data.add_data(_data)
        
        solution_data.write_to_file(_file_abs_path=file_path)
        
    def __supported_types_solution(self):
        """
        
        """
        #   predefine empty string
        _supported_types = ""
        for solution_type, type_name in zip(self.MBD_system._solution_filetypes, self.MBD_system._soltuion_filetype_software):
            _supported_type = type_name + " (*" + solution_type + ")" +";;"

            _supported_types = _supported_types + _supported_type

        #   remove last ;;
        _supported_types = _supported_types[0:-2]
        
        return _supported_types

    def _evaluate_C(self):
        """

        """
        t = 0
        q = self.MBD_system.evaluate_q()
        print "C(q, t):"
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_C(q, t)

    def print_C_q(self):
        """

        :return:
        """
        q = self.MBD_system.evaluate_q()
        print "C_q(q, t):"
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_C_q(q, 0.)

    def _evaluate_C_t(self):
        """

        :return:
        """
        t = 0
        q = self.MBD_system.evaluate_q()
        print "C_t(t):"
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_C_t(q, 0.)

    def print_Q_q(self):
        """

        :return:
        """
        q = self.MBD_system.evaluate_q()
        if self._parent.simulation_control_widget.solver.analysis.DAE_fun is None:
            self._parent.simulation_control_widget.solver.analysis.create_DAE()

        print "Q_q:"
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_Q_q(q)

    def print_Q_dq(self):
        """

        :return:
        """
        q = self.MBD_system.evaluate_q()
        if self._parent.simulation_control_widget.solver.analysis.DAE_fun is None:
            self._parent.simulation_control_widget.solver.analysis.create_DAE()
        print "Q_dq:"
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_Q_dq(q)

    def print_C_q_L_q(self):
        """

        :return:
        """
        q = self.MBD_system.evaluate_q()
        if self._parent.simulation_control_widget.solver.analysis.DAE_fun is None:
            self._parent.simulation_control_widget.solver.analysis.create_DAE()
        print "C_q_L_q:"
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_C_q_L_q(q, 0.)

    def print_J(self):
        """

        :return:
        """
        q = self.MBD_system.evaluate_q()
        if self._parent.simulation_control_widget.solver.analysis.DAE_fun is None:
            self._parent.simulation_control_widget.solver.analysis.create_DAE()
        print "J:"
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_Jacobian(q, 0.)

    def _edit(self, index):
        """

        """
        self._parentNode = self._item._parent

        print "self._item._typeInfo =", self._item._typeInfo

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

        if self._item._typeInfo.lower() == "motion":
            self._widget = MotionWidget(self._parentNode, parent=self)

        # if self._item._typeInfo.lower() == "MBDsystem":
        #     self._widget = MBDSystemWidget(self._item)

        self._widget.edit(self._item)
        # except:
        #     raise ValueError, "edit() not completed"

    def _delete(self):
        """

        """
        pos = self._item._parent._children.index(self._item)
        self._item._parent.removeChild(pos)
        del(self._item)

    def properties(self):
        """

        :return:
        """
        pprint(vars(self._item))
        if self._item._type == "contact":
            print self._item.contact_model

    def _save(self):
        """

        :return:
        """
        obj = self._item
        name = self._item._name+".dprj"
        read_and_write.write(obj, name)

    def _save_solution(self):
        """
        Function saves solution to file
        :return:
        """
        save_dialog = QtGui.QFileDialog()
        file_abs_path, file_type = save_dialog.getSaveFileNameAndFilter(self, 'Save solution', self.MBD_system.MBD_folder_abs_path, self._item._filetypes)

        self._item.write_to_file(_file_abs_path=file_abs_path)

    def _properties(self, index):
        """
        Open properties window
        """
        _widget = MBDSystemWidget(self._item, parent=self)
        _widget.edit(self._item)

    def _globalVariables(self):
        """

        :return:
        """
        pprint(vars(GlobalVariables))
    
    def _group_properties(self):
        """
        Function shows properties of group objects
        """
        if self._item._name.lower() == "bodies":
            for i, (child, body_MBD) in enumerate(zip(self._item._children, self.MBD_system.bodies)):
                print i, child, body_MBD._name
        
        if self._item._name.lower() == "forces":
            for i, (child, force_MBD) in enumerate(zip(self._item._children, self.MBD_system.forces)):
                print i, child, force_MBD._name
        
        if self._item._name.lower() == "joints":
            for i, (child, joint_MBD) in enumerate(zip(self._item._children, self.MBD_system.joints)):
                print i, child, joint_MBD._name

        if self._item._name.lower() == "contacts":
            for i, (child, contact_MBD) in enumerate(zip(self._item._children, self.MBD_system.contacts)):
                print i, child, contact_MBD._name

        if self._item._name.lower() == "springs":
            for i, (child, spring_MBD) in enumerate(zip(self._item._children, self.MBD_system.springs)):
                print i, child, spring_MBD._name

        if self._item._name.lower() == "motions":
            for i, (child, motion_MBD) in enumerate(zip(self._item._children, self.MBD_system.motions)):
                print i, child, motion_MBD._name

        if self._item._name.lower() == "measures":
            for i, (child, motion_MBD) in enumerate(zip(self._item._children, self.MBD_system.motions)):
                print i, child, motion_MBD._name

    def _save_contact_solution(self, index):
        """
        
        """
        self._item.save_solution_data()

    def _create_animation_file(self):
        """

        :return:
        """
        print "exe in tree view widget"
        self.create_animation_file.signal_createAnimationFile.emit()
        print "signal emited!!!"

    def evaluate_C(self):
        """

        :return:
        """
        #   get current q vector of MBD system
        q = self.MBD_system.evaluate_q()

        print self._item.evaluate_C(q)

    def evaluate_C0(self):
        """

        :return:
        """
        #   get current q vector of MBD system
        self.MBD_system.evaluate_q()

        q0 = self.MBD_system.q0
        print "C0 ="
        print self._item.evaluate_C(q0)

    def evaluate_C_q(self):
        """

        :return:
        """
        #   get current q vector of MBD system
        q = self.MBD_system.evaluate_q()

        C_q = self._item.evaluate_C_q(q)
        print "self._item =", self._item
        if self._item._typeInfo == "joint":
            for C_q_i, id in zip(C_q, ["i", "j"]):
                print "C_q_%s="%id
                print C_q_i
        else:
            print "C_q ="
            print C_q

    def evaluate_C_t(self):
        """

        :return:
        """
        #   get current q vector of MBD system
        q = self.MBD_system._children[0].evaluate_q()

        print self._item.evaluate_C_t(q=q, t=0)

    def evaluate_Q_d(self):
        """

        :return:
        """
        #   get current q vector of MBD system
        q = self.MBD_system.evaluate_q()
        if self._item._typeInfo == "joint":
            print self._item.evaluate_Q_d(q)
        else:
            print self.MBD_system.evaluate_Q_d(q)

    def evaluate_d(self):
        """

        :return:
        """
        #   get current q vector of MBD system
        q = self.MBD_system.evaluate_q()
        if self._item.typeInfo() == "joint" or self._item.typeInfo() == "spring":
            print self._item.evaluate_rijP(q)

        if self._item.typeInfo().lower() == "contact":
            print self._item.evaluate_rijP(q)
        
    def evaluate_delta(self):
        """
        
        """
        #   get current q vector of MBD system
        q = self.MBD_system.evaluate_q()
        print "q =", q
        if self._item.typeInfo().lower() == "contact":
            delta = self._item.evaluate_delta(q)
            print "delta =", delta
    
    def print_M(self):
        """
        
        :return:
        """
        self.check_DAE_fun()

        self._parent.simulation_control_widget.solver.analysis.DAE_fun.preprocessing()
        M = self._parent.simulation_control_widget.solver.analysis.DAE_fun.M
        print "M ="
        print M

        print "M_inv ="
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.M_inv

    def print_Q_g(self):
        """

        :return:
        """
        print "Q_g ="
        self.check_DAE_fun()

        self._parent.simulation_control_widget.solver.analysis.DAE_fun.preprocessing()
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_Q_g()

    def check_DAE_fun(self):
        """

        :return:
        """
        if self._parent.simulation_control_widget.solver.analysis.DAE_fun is None:
            self._parent.simulation_control_widget.solver.analysis.create_DAE()

    def print_Q_e(self):
        """

        :return:
        """
        self.check_DAE_fun()
        print "Q_e ="
        self._parent.simulation_control_widget.solver.analysis.DAE_fun.preprocessing()
        q = self.MBD_system.evaluate_q()
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_Q_e(0., q)

    def print_Q_s(self):
        """

        :return:
        """
        print "Q_s ="
        self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_M()
        q = self.MBD_system.evaluate_q()
        print self._parent.simulation_control_widget.solver.analysis.DAE_fun.evaluate_Q_s(q)

    def print_nodes(self):
        """
        Method prints nodes in a mesh
        :return:
        """
        if self._item.mesh.nodes:
            for i, node in enumerate(self._item.mesh.nodes):
                print i, node

        else:
            print Warning, "Nodes list empty!"

    def print_r_CP_list(self):
        """

        :return:
        """
        q = self.MBD_system.evaluate_q()
        self._item.print_r_CP_list(q)

    def print_r_slot_frame(self):
        """

        :return:
        """
        q = self.MBD_system.evaluate_q()

        self._item.print_r_slot_frame(q)

    def evaluate_pin_in_slot_position(self):
        """

        :return:
        """
        q = self.MBD_system.evaluate_q()

        self._item.evaluate_pin_position(q)

    def evaluate_Q_e(self):
        """

        :return:
        """
        #   evaluate current vector q of MBD system
        q = self.MBD_system.evaluate_q()
        if self._item._typeInfo == "force":
            self._item.print_Q_e(q)
        
        if self._item._typeInfo in ["spring", "contact"]:
            Q_e_list = []
            if self._item._typeInfo == "spring":
                Q_e_list = self._item.evaluate_Q_e(q)

            elif self._item._typeInfo == "contact":
                Q_e_list = self._item.evaluate_Q_e(0., q)

            for Q_e in Q_e_list:
                print "Q_e =", Q_e

    def evaluate_contact(self):
        """

        :return:
        """
        #   evaluate current vector q of MBD system
        q = self.MBD_system.evaluate_q()

        #   check for contact
        self._item.check_for_contact(0, 0., q)

        #   evaluate contact
        if self._item.contact_detected and self._item.status == 1:
            self._item.evaluate_contact(0., q)
    
    def evaluate_l0l(self):
        """
        Function evaluates deformation of spring length
        """
        q = self.MBD_system.evaluate_q()
        
        self._item.evaluate_l0l(q)

    def _evaluate_theta0(self):
        """
        Function evaluates initial angle of rotational spring
        :return:
        """
        q = self.MBD_system.evaluate_q0()

        print "theta0[deg] ="
        print np.rad2deg(self._item._evaluate_theta(q))

    def _evaluate_theta(self):
        """
        Function evaluates initial angle of rotational spring
        :return:
        """
        q = self.MBD_system.evaluate_q()

        print "theta[deg] ="
        print np.rad2deg(self._item._evaluate_theta(q))

    def _contact_point(self):
        """
        Method gets information of point of contact
        :return:
        """
        if self._parent.simulation_control_widget._status == "animation":
            self._item.get_contact_point(step=self._parent.simulation_control_widget._step)
    
    def _viewForces(self):
        """
        
        """
        print "_viewForces()"
        for force in self.MBD_system.forces:
            print force._name, force

    def plot_AABBi_tree_2D(self):
        """

        :return:
        """
        AABB = self._item.AABB_list[0]
        #   set plot
        fig, ax = self._set_plot()
        AABB.plot_2D(_ax=ax)

        plt.show()

    def plot_AABBj_tree_2D(self):
        """

        :return:
        """
        AABB = self._item.AABB_list[1]
        #   set plot
        fig, ax = self._set_plot()
        AABB.plot_2D(_ax=ax)

        plt.show()

    def plot_AABBi_2D(self):
        """

        :return:
        """
        AABB = self._item.AABB_list[0]
        AABB.plot_2D()

    def plot_AABBj_2D(self):
        """

        :return:
        """
        AABB = self._item.AABB_list[1]
        AABB.plot_2D()

    def _set_plot(self):
        """

        :return:
        """
        fig = plt.figure(num=1, figsize=(10, 8), dpi=300, facecolor='w', edgecolor='k')
        ax = plt.subplot(111, aspect="auto")
        ax.ticklabel_format(style='sci', axis='both')

        return fig, ax

