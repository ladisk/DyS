"""
Created on 21. feb. 2014

@author: lskrinjar
"""
import time
import itertools
import os
import sys
import inspect


import numpy as np
from matplotlib import pyplot as plt
import xlsxwriter
import vtk
from sympy import Symbol
from sympy.parsing.sympy_parser import parse_expr
import numpy as np
import scipy as sp
from scipy import integrate


from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.MBD_system_items import ForceItem
from MBD_system.force.force_matrix import Force_Q_e_matrix
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.transform_cs import u_P_cad2cm_lcs
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from global_variables import GlobalVariables
from MBD_system.q2q_body import q2q_body


class Force(ForceItem):
    """
    classdocs
    """
    __id = itertools.count(-1)

    def __init__(self, body_id, force_name=None, Fx=0., Fy=0., Mz=0., u_iP_f=np.array([0., 0., 0.]), r_P_GCS=np.zeros(2), node_id=None, element_id=None, element_ksi=None, L_i=None, filename=None, data=None, visible=False, scale=1E-3, dict={}, parent=None):
        """
        Constructor of force object
        :param body_id:       number of body on which the force is applied as int
        :param force_name:    force name as string
        :param F_x:           component of force
        :param F_y:           component of force
        :param M_z:           component of moment
        :param u_iP_f:        vector of acting force in CAD LCS of a body
        """
        #   body id
        self.body_id = body_id

        #   parent
        self._element = None
        self.created_automatically = True
        self.remove = False

        #   initial conditions
        self.q0 = []

        if hasattr(parent, "_type") and isinstance(self.body_id, int):
            #   force object is created from input file before simulation
            if parent._name.lower() == "forces":
                self._parent = parent
                self.created_automatically = False

                #   body handle
                self._body = self._parent._parent.bodies[self.body_id]

                self._parent._parent.evaluate_q()
                self.q0 = self._parent._parent.q0

            elif parent._typeInfo in ["spring", "joint"]:
                self._parent = parent
                if hasattr(self._parent._parent._parent, "bodies"):
                    self._body = self._parent._parent._parent.bodies[self.body_id]
                elif self._parent._parent._typeInfo in ["joint"]:
                    # print "TESTING =", self._parent._parent, self._parent._parent._typeInfo
                    self._body = self._parent._parent._parent._parent.bodies[self.body_id]
                # self._body = self._parent._parent._parent._name, self._parent._parent._parent.bodies[self.body_id]
                # time.sleep(100)

                # if parent._parent._typeInfo in ["joint"]:
                #     print "parent._parent =", parent._parent._parent._typeInfo
                # else:
                #     self._body = self._parent._parent.bodies[self.body_id]

            elif hasattr(parent, "_parent"):
                if parent._parent._name.lower() == "contacts":
                    self._parent = parent._parent._parent.forces
                    # self._parent = None

                if parent._parent._name.lower() == "springs":
                    self._parent = parent
                    # self._body = parent._parent._parent.ground

                #   body handle
                if isinstance(self.body_id, int):
                    if hasattr(parent._parent._parent, "bodies"):
                        self._body = parent._parent._parent.bodies[self.body_id]

                else:
                    self._body = parent._parent._parent.ground

            else:
                #   body handle
                if isinstance(self.body_id, int):
                    self._body = parent._parent._parent.bodies[self.body_id]
                else:
                    self._body = parent._parent._parent.ground

                # self._parent = parent._parent._parent.forces
                # print "created automatically at contact!"
                # print "parent =", parent

        else:
            self._parent = None
            self._body = None
        # print "self._parent =", self._parent
        # if hasattr(self._parent, "_name"):
        #     print "name =", self._parent._name
        # print "--------------------------------------------"
        # if not inspect.isclass(self._parent):
        #     self._parent = None
        super(Force, self).__init__(force_name, parent=self._parent)
        #    number
        self.force_id = self.__id.next()
        #    additional user input comments
        self._comments = ""

        #    name as string
        if force_name is None:
            self._name = "Force_" + str(self.force_id)
        else:
            self._name = force_name

        #   additional force object property
        self._body_assigned = False

        #   force components Fx, Fy, Mz
        self.Fx = Fx
        self.Fy = Fy
        self.F = np.array([self.Fx, self.Fy])
        self.F_e = np.linalg.norm(self.F, ord=2)
        self.Mz = Mz

        #    generalized force vector
        self.Q_e = np.zeros(3)

        #   current step of integration
        self.step = 0
        
        #    user sub-routine as filename attribute
        self.filename = filename

        #   function of time
        self._Fx_t = 0.
        self._Fy_t = 0.
        self._Mz_t = 0.
        
        #    coordinate system, options: 
        #    LCS - force is applied in body LCS coordinates
        #    GCS - force is applied in GCS coordinates
        self.CS = "GCS"

        #   rigid body
        #    position of acting force in CAD LCS of a body
        self.u_iP_f = u_iP_f[0:2]
        #   position of force vector in CS
        #   LCS
        self.u_P_LCS = np.zeros(2)
        #   GCS
        self.r_P_GCS = r_P_GCS

        #   flexible body
        self.node_id = node_id
        self.element_id = element_id
        self.element_ksi = element_ksi
        self.L_i = L_i

        #    z dimension
        if len(u_iP_f) == 2:
            self.z_dim_lcs = 0.
        else:
            self.z_dim_lcs = u_iP_f[2]

        #   force data from file - matrix F(t)
        self.data = data

        #   visualization properties
        self.vtk_mapper = None
        self.vtk_actor = None
        self.force = None
        self.renderer = None
        self.scale = scale
        self.color = np.array([1, 0, 0], dtype="float32")
        self._visible = visible
        self.P1 = None
        self.P2 = None

        #   force position in body LCS (center of mass)
        if self._body is not None:
            self.u_P_LCS = u_P_cad2cm_lcs(self.body_id, self._body, self.u_iP_f)
        else:
            self.u_P_LCS = self.u_iP_f
        
        #   init solution containers
        self._solution_containers()

        #   add extra attributes from dictionary
        if dict is not None:
            self.add_attributes_from_dict(dict)

        if self.filename:
            self.file, extension = os.path.split(self.filename)

        self.filepath = None
        if self.filename is not None:
            if hasattr(GlobalVariables, "MBDsystem_folder"):
                self.filepath = os.path.join(GlobalVariables.MBDsystem_folder, self.filename)

        #    add force object to body list of forces
        # if hasattr(self._parent, "_typeInfo"):
        #     if self._parent._typeInfo == "group":
        #         self._parent._parent.bodies[self.body_id].forces.append(self)
        
        # self.markers = self._create_markers()

        #   create vtk data for visualization automatically during simulation when force object is created atomatically
        #   e.g. when contact is detected and contact forces are present
        # if self._parent is None:
        #     if self._parent.typeInfo() == "group":
        #         self.set_vtk_data()

        self._update_count = 0

        # print "CREATED FORCE OBJECT =", self, self._name

    def set_element(self):
        """

        :param element:
        :return:
        """
        _element = None
        if self.node_id == -1:
            #   last element
            _element = self._body.mesh.elements[-1]
            self.element_ksi = 1

        elif isinstance(self.node_id, int):
            _element = next((element for element in self._body.mesh.elements if self.node_id in element.node_id_list), None)
            self.element_ksi = 0

        elif self.L_i is not None and self.node_id is None:
            L = 0.
            for i, element in enumerate(self._body.mesh.elements):
                print i, L, self.L_i, element.L
                if L <= self.L_i < L + element.L:
                    _element = element
                    self.element_ksi = (self.L_i - L) / element.L
                else:
                    L += element.L
                    _element = None

        else:
            raise ValueError, "Node id out of range!"

        if _element is None:
            raise ValueError, "Element not found!"

        # print "self.element_ksi =", self.element_ksi, "name =", self._name

        return _element

    def set_vtk_data(self, renderer=None):
        """

        :return:
        """
        if renderer is not None:
            self.renderer = renderer

        #   force vector
        Q_e = self.evaluate_Q_e(0, self._body._parent._parent.q0)

        Q_e_x = self.F[0]
        Q_e_y = self.F[1]

        #   create force vector as line
        self.force = vtk.vtkLineSource()
        if len(self.r_P_GCS) == 2:
            self.P1 = np.append(self.r_P_GCS, 0.)
        elif len(self.r_P_GCS) == 3:
            self.P1 = self.r_P_GCS
        else:
            pass

        self.force.SetPoint1(self.P1[0], self.P1[1], self.P1[2])
        self.P2 = self.P1 + self.scale * np.array([Q_e_x, Q_e_y, 0.])
        self.force.SetPoint2(self.P2[0], self.P2[1], self.P2[2])

        # self.force.SetPoint1(0, 0, 0)
        # self.force.SetPoint2(1, -1, 1)

        #   create mapper
        self.vtk_mapper = vtk.vtkPolyDataMapper()
        self.vtk_mapper.SetInputConnection(self.force.GetOutputPort())

        #   create actors
        self.vtk_actor = vtk.vtkActor()
        self.vtk_actor.SetMapper(self.vtk_mapper)
        self.vtk_actor.GetProperty().SetColor(self.color)
        self.vtk_actor.GetProperty().SetLineWidth(2)

        #   set visibility for object actor
        # if not self._visible:
        #     self.vtk_actor.VisibilityOff
        self.vtk_actor.VisibilityOn()

        if self.renderer is not None:
            self.renderer.AddActor(self.vtk_actor)

    def update_vtk_data(self, t=0., q=None):
        """

        :return:
        """
        if self.active and not self.remove and self.vtk_actor is not None:
            #   force vector
            Fx = self.F[0]
            Fy = self.F[1]

            #   update point of applied force
            # print "self.r_P_GCS =", self.r_P_GCS, self._name
            self.P1 = np.append(self.r_P_GCS, 0.)
            self.force.SetPoint1(self.P1[0], self.P1[1], self.P1[2])

            #   update force vector
            try:
                self.P2 = self.P1 + self.scale * np.array([Fx, Fy, 0.])
            except:
                print "name =", self._name
            self.force.SetPoint2(self.P2[0], self.P2[1], self.P2[2])

            self.force.Modified()
            self.vtk_mapper.Update()

        else:
            self._clear_vtk_data()
            # if list:
            #     list.remove(self)

    def clear_vtk_data(self):
        """

        :return:
        """
        self._clear_vtk_data()

    def _clear_vtk_data(self):
        """

        :return:
        """
        # print "vtk_actor is NOT NONE!@_clear_vtk_data", self._name, self, "parent =", self._parent
        if self.renderer is not None and self.vtk_actor is not None:
            # print "self.vtk_actor ="
            # print self.vtk_actor
            # print self.renderer.GetActors()
            self.renderer.RemoveViewProp(self.vtk_actor)
            self.renderer.RemoveActor(self.vtk_actor)
            # time.sleep(0.1)
            self.vtk_mapper = None
            self.vtk_actor = None
            self.force = None
            self.vtk_actor = None

            # self.force.Delete()
            # self.vtk_mapper.Delete()
            self.renderer = None

    def deactivate(self):
        """

        :return:
        """
        # print "deactivate()", self, self._name
        self.active = False
        self.remove = True
        self.vtk_actor.VisibilityOff()
        # print "active =", self.active, "remove =", self.remove

    def set_deactivate(self):
        """

        :return:
        """
        self.active = False

    def set_remove(self):
        """

        :return:
        """
        self.remove = True

    def _show(self):
        """

        :return:
        """
        if self.vtk_actor is not None:
            if self._visible:
                self._visible = False
                self.vtk_actor.VisibilityOff()
            else:
                self._visible = True
                self.vtk_actor.VisibilityOn()

    def load_data_file(self, filename):
        """

        :return:
        """

    def _load_csv(self):
        """

        :return:
        """

    def _solution_containers(self):
        """
        
        """
        self._solution_container_names = ["_step_num_solution_container",
                                     "_t_solution_container",
                                     "_u_iP_f_solution_container",
                                     "_Q_e_solution_container",
                                     "_Fx_solution_container",
                                     "_Fy_solution_container",
                                     "_Mz_solution_container"]

        self._step_num_solution_container = []
        self._t_solution_container = []

        self._u_iP_f_solution_container = []
        self._Q_e_solution_container = []
        self._Fx_solution_container = []
        self._Fy_solution_container = []
        self._Mz_solution_container = []

    def _evaluate_F(self, t):
        """
        Function evaluates force function if is function of time. String symbol is evaluated using sympy package with
        function evalf()
        :param t:
        :return:
        """
        time = Symbol('time')

        if isinstance(self.Fx, basestring):
            _exp = parse_expr(self.Fx)
            _Fx = _exp.evalf(subs={time:t})
        else:
            _Fx = self.Fx

        if isinstance(self.Fy, basestring):
            _exp = parse_expr(self.Fy)
            _Fy = _exp.evalf(subs={time:t})
        else:
            _Fy = self.Fy

        if self.filepath is not None:
            if os.path.isfile(self.filepath):
                #   set arguments to pass in function variable
                sys.argv = [self.filename, t]
                #   execute file
                execfile(self.filepath)
                #   get local variables as dict
                __locals = locals()

                #   evaluated force
                _Fx = __locals["Fx"]
                _Fy = __locals["Fy"]

        self.Fx = _Fx
        self.Fy = _Fy
        F = np.array([_Fx, _Fy])
        return F

    def print_F(self):
        """

        :return:
        """
        print "F ="
        print self._evaluate_F(0.)

    def print_Q_e(self, q):
        print self.evaluate_Q_e(0, q)

    def evaluate_Q_e(self, t, q):
        """
        Evaluates generalized external force vector from force object properties and vector q (only angle theta is used from q vector)
        in:
            q - MBS system vector
        out:
            vector of force acting on a body (Fx, Fy, Mz)
        """
        #   construct force vector (can be function of time)
        self.F = self._evaluate_F(t)

        if self.active and not self.remove:
            if self._parent is not None or self._body is not None:#this is NOK for force at spring level

                if self._body.body_type == "rigid body":
                    self.Q_e = self._evaluate_Q_e_rigid(t, q, self.F)

                elif self._body.body_type in ["flexible body", "finite element"]:
                    self.Q_e = self._evaluate_Q_e_flexible_F(t, q, self.F) + self._evaluate_Q_e_flexible_M(t, q, self.Mz)
                    # if (self.F != np.zeros(2)).any():
                    #     # print "_evaluate_Q_e_flexible_F()"
                    #     self.Q_e =
                    #
                    # elif self.Mz != 0.:
                    #     # print "_evaluate_Q_e_flexible_M()"
                    #     self.Q_e = self._evaluate_Q_e_flexible_M(t, q, self.Mz)
                    #
                    # else:
                    #
                    #     print "Force Attribute F is not defined properly"
                else:
                    raise ValueError, "Body type not correct! Body type is %s"%self._parent._parent.bodies[self.body_id].body_type

            else:
                self.Q_e = None

        else:
            self.Q_e = np.zeros(3)

        return self.Q_e

    def _evaluate_Q_e_rigid(self, t, q, F):
        """
        Evaluate generalized external force for rigid body
        :param t:
        :param q:
        :return:
        """
        #   position of force in GCS
        self.r_P_GCS = u_P_lcs2gcs(self.u_P_LCS, q, self.body_id)

        #   get theta of body from q vector
        theta = q2theta_i(q, self.body_id)

        #   construct a matrix (function of vector q)
        matrix = Force_Q_e_matrix(self.body_id, self.u_iP_f, theta)

        #   form a vector
        Q_e = np.dot(matrix.matrix, F)

        return Q_e

    def _evaluate_Q_e_point_mass(self, t, q, F):
        """
        Evaluate generalized external force for point mass type of body
        :return:
        """
        #   position of force in GCS
        self.r_P_GCS = q2R_i(q, self.body_id)

        #   generalized external force vector
        Q_e = F

        return Q_e

    def _evaluate_Q_e_flexible_M(self, t, q, M):
        """
        Evaluate generalized external force for ANCF flexible body
        :param t:
        :param q:
        :param M:
        :return:
        """
        if self._element is None:
            self._element = self.set_element()

        #   predefine empty vector
        Q_e = np.zeros(self._element.mesh.n_NC)

        #   get body generalized coordinates
        e_b = q2q_body(q, self.body_id)

        for i, element in enumerate(self._body.mesh.elements):
            if self.node_id in element.node_id_list or self.node_id == -1:

                if self.node_id != -1:
                    ksi = element.node_id_list.index(self.node_id)
                else:
                    ksi = 1

                e_i = element.evaluate_e_i(e_b=e_b)

                Q_e_i = element.evaluate_Q_e_M(e_i, M, ksi)

                #   transformed force vector
                if element.B is None:
                    element.evaluate_B()

                if element.T is None:
                    element.evaluate_T()

                _Q_e_i = reduce(np.dot, [element.B.T, element.T.T, Q_e_i])

                #   sum
                Q_e += _Q_e_i

        return Q_e

    def _evaluate_Q_e_flexible_F(self, t, q, F):
        """

        :param t:
        :param q:
        :return:
        """
        if self._element is None:
            self._element = self.set_element()

        #   predefine empty vector
        Q_e = np.zeros(self._element.mesh.n_NC)

        for i, element in enumerate(self._body.mesh.elements):
            if self.node_id in element.node_id_list or self.node_id == -1 and self.element_ksi is None:
                #   undeformed position
                # self.r_P_GCS = self._element.mesh.nodes[self.node_id]

                #   deformed position
                if self.node_id == -1:
                    node = self._element.geometry_nodes[-1, :]
                else:
                    if self._element.node_id_list.index(self.node_id) == 0:
                        node = self._element.geometry_nodes[0, :]

                    elif self._element.node_id_list.index(self.node_id) == 1:
                        node = self._element.geometry_nodes[-1, :]

                    else:
                        node = np.zeros(2)

                self.r_P_GCS = node

                if self.node_id != -1:
                    ksi = element.node_id_list.index(self.node_id)

                elif self.element_ksi is not None:
                    ksi = self.element_ksi

                else:
                    ksi = 1

                S = element._evaluate_S(ksi)

                if (self.node_id == 0) or (self.node_id == len(self._body.mesh.nodes)):
                    pass
                else:
                    F = 0.5 * F

                Q_e_i = np.dot(S.T, F)

            else:
                self.r_P_GCS = self._element.evaluate_r(self.element_ksi)

                S = self._element._evaluate_S(self.element_ksi)

                Q_e_i = np.dot(S.T, F)

            #   transformed force vector
            if element.B is None:
                element.evaluate_B()

            if element.T is None:
                element.evaluate_T()

            _Q_e_i = reduce(np.dot, [element.B.T, element.T.T, Q_e_i])

            #   sum
            Q_e += _Q_e_i


        return Q_e

    def set_F_u_P(self, F=np.zeros(2), u_P=np.zeros(2)):
        """

        :param F:
        :param u_P:
        :return:
        """
        #   define force vector
        self._update_F(F)
        #   define position vector of point of application
        self._update_u_P(u_P)

    def update_F(self, q, step=None, F=np.zeros(2), u_P=np.zeros(2), r_P=np.zeros(2)):
        """

        :param step:
        :param Fx:
        :param Fy:
        :param Mz:
        :return:
        """
        if step is None:
            pass
        else:
            self._step_num_solution_container.append(step)

        #    update applied force vector
        self._update_F(F)

        #    update position vector of applied force
        #   LCS
        if (u_P == np.zeros(2)).all():
            u_P = self.u_iP_f

        self._update_u_P(u_P)

        #   GCS
        if (r_P == np.zeros(2)).all() and hasattr(q, "__len__"):
            r_P = u_P_lcs2gcs(u_P, q, self.body_id)

        self._update_r_P(r_P)

    def update_Q_e(self, step=None, Q_e=np.zeros(3), u_P=np.zeros(2)):
        """
        Function updates a generalized external force vector Q_e of applied force
        :param step:
        :param Q_e:
        :return:
        """
        #   update generalized force vector
        self.Q_e = Q_e

        #   update force position vector in body LCS
        if (u_P == np.zeros(2)).all:
            u_P = self.u_iP_f
        self._update_u_P(u_P)

    def update_step(self, step):
        """

        :param step:
        :return:
        """
        self.step = step

        self._step_num_solution_container.append(step)

    def _update_F(self, F):
        """
        Function updates force vector
        :return: None
        """
        [self.Fx, self.Fy] = F
        self.F = F

    def _update_Mz(self, Mz):
        """

        :param Mz:
        :return:
        """
        self.Mz = Mz

    def _update_u_P(self, u_P):
        """
        
        """
        self.u_P_LCS = self.u_iP_f = u_P
        # self._u_iP_f_solution_container.append(u_P)

    def _update_r_P(self, r_P):
        """

        :param r_P:
        :return:
        """
        self.r_P_GCS = r_P

    def reset(self, q0):
        """

        :return:
        """
        #   set q0
        self.q0 = q0

        self.Fx = 0.
        self.Fy = 0.
        self.Mz = 0.

        if self._body.body_type == "rigid body":
            self.r_P_GCS = u_P_lcs2gcs(self.u_P_LCS, self.q0, self.body_id)

        if self._body.body_type == "flexible body":
            self.r_P_GCS = self.q0

        self.F = self._evaluate_F(0.)

        self.update_vtk_data()

    def _track_data(self, step, h, t):
        """

        :param step:
        :param h:
        :param t:
        :param q:
        :return:
        """
        self._step_num_solution_container.append(step)

        self._t_solution_container.append(t)

        self._u_iP_f_solution_container.append(self.u_iP_f)
        
        self._Fx_solution_container.append(self.Fx)
        self._Fy_solution_container.append(self.Fy)

        self._Q_e_solution_container.append(self.Q_e)

    def get_data(self):
        """

        :return:
        """
        print "Force solution data:"
        print "t =", len(self._t_solution_container)
        print "uP =", len(self._u_iP_f_solution_container)
        print "Fx =", len(self._Fx_solution_container)
        print "Fy =", len(self._Fy_solution_container)
        print "Mz =", len(self._Mz_solution_container)
        for i, (t, uP, Fx, Fy, Q_e) in enumerate(zip(self._t_solution_container, self._u_iP_f_solution_container, self._Fx_solution_container, self._Fy_solution_container, self._Q_e_solution_container)):
            print i, t, uP, Fx, Fy, Q_e

    def save_solution_data(self, filename=None):
        """

        :return:
        """
        self._comments = "Contact force on body id: "+str(self.body_id)
        if filename is None:
            self._solution_filename = "_" + self._name + "_sol"+".xlsx"
        else:
            self._solution_filename = filename

        #   column headers
        _header = ["i-th step", "u_Px", "u_Py", "Fx", "Fy", "Mz"]

        # Create an new Excel file and add a worksheet.
        workbook = xlsxwriter.Workbook(self._solution_filename, {'nan_inf_to_errors': True})
        format_1 = workbook.add_format({'num_format': '#0.000000000000000'})
        format_2 = workbook.add_format({'num_format': '#0.00000000'})
        worksheet = workbook.add_worksheet("data")

        #   write comments
        worksheet.write(0, 0, self._comments)
        #   write header
        worksheet.write_row(1,0,_header)

        worksheet.write_column(2, 0, np.array(self._step_num_solution_container, dtype="float32"))
        worksheet.write_column(2, 1, np.array(self._u_iP_f_solution_container, dtype="float32")[:,0], format_1)
        worksheet.write_column(2, 2, np.array(self._u_iP_f_solution_container, dtype="float32")[:,1], format_1)
        worksheet.write_column(2, 3, np.array(self._Fx_solution_container, dtype="float32"), format_1)
        worksheet.write_column(2, 4, np.array(self._Fy_solution_container, dtype="float32"), format_1)
        worksheet.write_column(2, 5, np.array(self._Mz_solution_container, dtype="float32"), format_1)

        #   column width
        worksheet.set_column('B:E', 20)

        #   close file
        workbook.close()

    def evaluate_F(self, t, t0=None, t02=None):
        """

        :param t:
        :return:
        """
        F0 = 1.626
        f0 = 38
        f = np.pi*f0
        t0 = 0.076
        t_F = t0/2

        if t < t_F:
            F = F0 * np.sin((2*np.pi/t0)*t)
        else:
            F = 0

        return F

    def evaluate_Q_q(self):
        """

        :return:
        """

    def evaluate_Q_dq(self):
        """

        :return:
        """

    def testing(self):
        """

        :return:
        """
        #   this is working
        self._parent._parent._parent.dys.simulation_control_widget.vtkWidget.renderer.RemoveActor(self.vtk_actor)
        # print " self._parent._parent._parent._parent =",  self._parent._parent._parent._parent
        # self._parent._parent._parent._parent.dys.simulation_control_widget.vtkWidget.renderer.AddActor(self.vtk_actor)

if __name__ == "__main__":
    filename = "input_force.py"
    f = Force(0, force_name="ForceTesting", Fx=1, Fy=0, Mz=0, u_iP_f=np.array([0, 0, 0]), filename=filename)
    f.filepath = os.path.join(os.getcwd(), filename)
    # print "f =", f._evaluate_F(0.)

    t = np.arange(0, 2+1E-2, 1E-2)
    Fx = Fy = np.zeros(len(t))
    for i in range(0, len(t)):
        Fx[i], Fy[i] = f._evaluate_F(t[i])

    fig = plt.figure(num=1, figsize=(6, 4), dpi=72, facecolor='w', edgecolor='k')
    ax = plt.subplot(111, aspect="auto")  # auto, normal

    plt.plot(t, abs(Fy), color="black")

    plt.ylabel("F [N]")
    plt.xlabel("t [s]")

    labels = [item.get_text() for item in ax.get_yticklabels()]
    # print "labels =", labels
    labels[-2] = 'F0'
    # print "labels =", labels
    ax.set_yticklabels(labels)

    plt.grid()
    plt.show()



