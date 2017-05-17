"""
Created on 21. feb. 2014

@author: lskrinjar
"""
import itertools
from pprint import pprint
import numpy as np
import vtk


from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.MBD_system_items import SpringItem
from MBD_system.dr_ij_P_dq import dr_ij_P_dq
from MBD_system.force.force import Force
from MBD_system.force.force_matrix import Force_Q_e_matrix
from MBD_system.force.force_vector import Force_Q_e_vector
from MBD_system.q2R_i import q2R_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.r_ij_P import r_ij_P
from MBD_system.u_P_cad2cm_lcs import u_P_cad2cm_lcs
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from simulation_control_widget.vtk_widget.marker.marker import Marker
from global_variables import GlobalVariables
from MBD_system.MBD_system_items import GroupItem


class Spring(SpringItem):
    """
    classdocs
    """
    __id = itertools.count()

    def __init__(self, _name, body_id_i, body_id_j, u_iP_CAD=np.zeros(2, dtype=float), u_jP_CAD=np.zeros(2, dtype=float), u_iP_LCS=np.zeros(2, dtype=float), u_jP_LCS=np.zeros(2, dtype=float), k=0., c=0., L_i=None, L_j=None, joint_id=None, properties_dict={}, parent=None):
        """
        Create a spring class
        :param _name:       name of a spring as string
        :param spring_type: type of a sprig (linear, rotational)
        :param body_id_i:       id of body i
        :param body_id_j:       id of body j
        :param u_iP_CAD:        in CAD LCS of a body i
        :param u_jP_CAD:        in CAD LCS of a body j
        :param k:               spring stiffness as float (physical units N/length or N/angle)
        :param c:               damping coefficient as float (physical units Ns/length or Ns/angle)
        :param l_0:             undeformed length/angle of a spring as float
        """
        super(Spring, self).__init__(_name, parent)
        #   parent
        self._parent = parent

        #    number
        self.spring_id = self.__id.next()

        #   status
        self.active = True

        #    name as string
        self._name = _name

        #   reference to joint at which the spring is located
        self.joint_id = joint_id

        #   body type list
        self.body_type_list = ["rigid", "rigid"]

        #   direction (compression, tension, None)
        self.direction = None

        #   position of spring in GCS
        self.z_dim = 0.
        #   get q vector of current MBD system at joint object initialization
        if hasattr(self._parent._parent, "evaluate_q0"):
            self.q0 = self._parent._parent.evaluate_q0()
        else:
            self.q0 = None

        #   physical properties
        self.k = k
        self.c = c
        self.F0 = 0.

        #   direction of the spring (unit vector) at simulation time
        self.I_r_ij_0 = None#np.array([1, 0])
        self.r_ij_P = None
        self.I_r_ij_P = None

        #   length of the spring at simulation time
        self.l = None
        self.dl = None

        #   undeformed length
        self.l_0 = None

        #   deformation of spring length
        self.l0l = 0.

        #   spring force at simulation time
        self.F_s = 0.

        #   spring potential energy
        self.energy = 0.

        #   list of properties at simulation time
        self.Q_q_list = []
        self.Q_dq_list = []

        #   flexible body data
        #   flexible body data if a flexible body is in a spring
        self.node_id_i = None
        self.node_id_j = None
        self.node_id_list = [self.node_id_i, self.node_id_j]

        self.element_id_i = None
        self.element_id_j = None
        self.element_id_list = [self.element_id_i, self.element_id_j]

        self.element_ksi_i = None
        self.element_ksi_j = None
        self.element_ksi_list = [self.element_ksi_i, self.element_ksi_j]

        self.L_i = L_i
        self.L_j = L_j
        self.L_list = [self.L_i, self.L_j]

        #   markers list
        self.markers = []

        #   list of forces
        self._Fn_list = []

        #   visualizations properties
        self.color = np.array([1, 0, 0], dtype="float32")
        self.vtk_actor = None
        self.vtk_mapper = None
        self.line = None
        self.points = None
        self.line_grid = None

        # swap body_id that if body is connected to ground that ground is always the last item in list
        if body_id_i == "ground":
            self.body_id_i = body_id_j
            self.body_id_j = body_id_i
            self.u_iP_CAD = u_jP_CAD
            self.u_jP_CAD = u_iP_CAD

        else:
            self.body_id_i = body_id_i
            self.body_id_j = body_id_j
            self.u_iP_CAD = u_iP_CAD
            self.u_jP_CAD = u_jP_CAD

        #   body id list
        self.body_id_list = [self.body_id_i, self.body_id_j]
        self.body_list = []

        #   list of generalized spring force object at each body
        self.Q_e_list = []

        #   list of point vectors to joint constraint in CAD lCS of a body
        self.u_P_CAD_list = [self.u_iP_CAD, self.u_jP_CAD]

        #   predefined empty list to store point vectors of joint constraint in LCS (center of gravity)
        #   of each body
        self.u_iP_LCS = u_iP_LCS
        self.u_jP_LCS = u_jP_LCS
        self.u_P_LCS_list = [self.u_iP_LCS, self.u_jP_LCS]
        #   evaluate vectors in LCS of a body if defined in CAD CS of a body
        self.u_P_LCS_list = self._transform_uP_cad2lcs()
        [self.u_iP_LCS, self.u_jP_LCS] = self.u_P_LCS_list

        #   in GCS
        self.r_P_i = np.zeros(2, dtype="float32")
        self.r_P_j = np.zeros(2, dtype="float32")
        self.r_P_list = [self.r_P_i, self.r_P_j]

        #    set properties dictionary as object property
        self.properties = properties_dict
        self.additional_params_calculated = False

        #   predefined array of signs for body i and body j
        self.signs = np.array([-1, +1])

        #   predefined empty list of spring force matrix to store one force matrix for each body connected with a spring
        self.spring_force_matrix_list = []

        #   init solution containers
        self._solution_containers()

        #   add additional properties
        # self.add_attributes_from_dict(self.properties)
        #
        # #   list of markers
        # self.markers = self._create_markers()
        #
        # #   create pair of forces list
        # self._Fn_list = self._create_Fn_forces()

    def set_vtk_data(self, q=None):
        """

        :param q:
        :return:
        """
        #   points
        self.points = vtk.vtkPoints()
        self.points.SetNumberOfPoints(2)
        for i, node in enumerate(self.r_P_list):
            if len(node) == 2:
                node = np.append(node, 0.)
            self.points.InsertPoint(i, node)

        #   set vtk line object
        self.line = vtk.vtkLine()
        self.line.GetPointIds().SetId(0, 0)
        self.line.GetPointIds().SetId(1, 1)

        self.line_grid = vtk.vtkUnstructuredGrid()
        self.line_grid.Allocate(1, 1)
        self.line_grid.InsertNextCell(self.line.GetCellType(), self.line.GetPointIds())
        self.line_grid.SetPoints(self.points)

        self.vtk_mapper = vtk.vtkDataSetMapper()
        self.vtk_mapper.SetInputData(self.line_grid)

        self.vtk_actor = vtk.vtkActor()
        self.vtk_actor.SetMapper(self.vtk_mapper)
        self.vtk_actor.GetProperty().SetColor(self.color)

    def update_vtk_data(self, q):
        """

        :param q:
        :return:
        """
        #   replace nodes
        for i, rP in enumerate(self.r_P_list):
            rP = np.append(rP, 0.)
            #   set point
            self.points.SetPoint(i, rP)

        #   data is changed
        self.points.Modified()
        self.vtk_mapper.Update()

        for marker, F, rP in zip(self.markers, self._Fn_list, self.r_P_list):
            marker.update_vtk_data(q, rP=rP)

    def _solution_containers(self):
        """

        :return:
        """
        self._step_num_solution_container = []
        self._t_solution_container = []

        #   position of spring ends on each body in GCS
        self._r_P_solution_container = []

        #   spring length
        self._l_solution_container = []

        #   spring force
        self._F_solution_container = []

    def _add_attributes_from_parent(self):
        """
        Add attributes from parent item
        :return:
        """
        d = {}
        if self._parent._typeInfo in ["joint"]:
            d = self._parent.__dict__

        if isinstance(self.joint_id, int):
            d = self._parent._parent.joints[self.joint_id].__dict__

        for key, val in d.iteritems():
            if hasattr(self, key) and key[0] != "_":
                setattr(self, key, val)

                if key == "node_id_i":
                    self.body_type_list[0] = "flexible"

                if key == "node_id_j":
                    self.body_type_list[1] = "flexible"

    def _create_Q_e(self):
        """
        Function creates Q_e list of generalized external forces
        :return:
        """
        Q_e_list = []

        for i, body_id in enumerate(self.body_id_list):
            if isinstance(body_id, int):
                #   size of vector Q_e
                n_Q_e = GlobalVariables.q_i_dim[body_id]
                #   rigid body
                if n_Q_e == 3:
                    Q_e = Force_Q_e_vector(n=n_Q_e)

                #   flexible body
                else:
                    name = "spring_force_on_flexible_body_of_" + self._parent._name

                    #   if parent is object - group of springs, than use this object properties
                    if isinstance(self._parent, GroupItem):
                        Q_e = Force(body_id, force_name=name, node_id=self.node_id_list[i], element_id=self.element_id_list[i], element_ksi=self.element_ksi_list[i], parent=self)

                    else:
                        Q_e = Force(body_id, force_name=name, node_id=self._parent.node_id_list[i], element_id=self._parent.element_id_list[i], element_ksi=self._parent.element_ksi_list[i], parent=self)

            else:
                Q_e = Force_Q_e_vector()

            Q_e_list.append(Q_e)

        return Q_e_list

    def evaluate_rijP(self, q):
        """
        Length of spring at time when vector of absolute coordinates is equal to q
        """
        for i, (body_id, uP) in enumerate(zip(self.body_id_list, self.u_P_LCS_list)):
            self.r_P_list[i] = u_P_lcs2gcs(uP, q, body_id)

        r_ij_P = self.r_P_list[0] - self.r_P_list[1]

        return r_ij_P

    def _transform_uP_cad2lcs(self):
        """
        Function transforms point vectors from cad to lcs of a body to be used during numerical computation
        :return:
        """
        u_P_LCS_list = []
        for body_id, u_P, u_P_LCS in zip(self.body_id_list, self.u_P_CAD_list, self.u_P_LCS_list):
            if body_id == "ground" or body_id == -1 or u_P is None:
                # u_P_LCS = u_P_cad2cm_lcs(body_id, self._parent._parent.ground, _u_P=u_P)
                u_P_LCS = u_P
                # _body = self._parent._parent.ground
            else:
                #   create pointer to body
                u_P_LCS = np.zeros(2, dtype=float)
                if self._parent is not None and (u_P != np.zeros(2, dtype="float")).any():
                    _body = self._parent._parent.bodies[body_id]
                    #   calculate point vector in body LCS (center of gravity)
                    u_P_LCS = u_P_cad2cm_lcs(body_id, _body, _u_P=u_P)

            u_P_LCS_list.append(u_P_LCS)

        return u_P_LCS_list

    def testing(self):
        """

        :return:
        """
        for marker in self.markers:
            print marker
            pprint(vars(marker))

    def _create_markers(self):
        """
        Function create markers
        :return:
        """
        markers = []
        self.r_P_list = []
        for body_id, u_P in zip(self.body_id_list, self.u_P_LCS_list):
            if body_id == "ground" or body_id == -1:
                #   pointer to body object
                body = self._parent._parent.ground
                r_P = u_P

            else:
                r_P = u_P_lcs2gcs(u_P, self.q0, body_id)

                #   pointer to body object
                if self._parent is not None:
                    body = self._parent._parent.bodies[body_id]
                else:
                    body = None

            self.r_P_list.append(r_P)

            if body is not None:
                #   create marker object
                r_node = np.array(np.append(r_P, self.z_dim),  dtype="float32")
                u_node = np.array(np.append(u_P, self.z_dim),  dtype="float32")
                marker = Marker(r_node, uP=u_node, body_id=body, parent=body)
                #   append marker object
                #   to body list of markers
                body.markers.append(marker)
                #   to markers list
                markers.append(marker)

                self.body_list.append(body)

        return markers

    def _create_Fn_forces(self):
        """

        :return:
        """
        Fn_list = []
        for body_id, _u_P in zip(self.body_id_list, self.u_P_CAD_list):
            if body_id == "ground" or body_id == -1:
                u_P_LCS = u_P_cad2cm_lcs(body_id, self._parent._parent.ground, _u_P=_u_P)
                body = self._parent._parent.ground
            else:
                #   create pointer to body
                if self._parent is not None:
                    body = self._parent._parent.bodies[body_id]

                    #   calculate point vector in body LCS (center of gravity)
                    u_P_LCS = u_P_cad2cm_lcs(body_id, body, _u_P=_u_P)

                else:
                    body = None
                    u_P_LCS = np.zeros(2, dtype=float)

            #   create force object
            if body is not None:
                force = Force(body_id, force_name=self._name + "_on_body_" + str(body_id), u_iP_f=u_P_LCS, parent=self)

                #   add pair of contact forces to forces list of MBD system
                # self._parent._parent.forces.append(force)

                #   add pair of contact forces to forces list of spring
                Fn_list.append(force)

                if body_id != "ground":
                    self._parent._parent.bodies[body_id].forces.append(force)
                else:
                    self._parent._parent.ground.forces.append(force)

        return Fn_list

    def _update_F(self, Q_e_list):
        """
        Function updates list of forces
        :return:
        """
        for body_id, F, Q_e in zip(self.body_id_list, self._Fn_list, Q_e_list):
            F.update_Q_e(Q_e=Q_e.Q_e)

    def _track_data(self, step, h, t):
        """

        :param step:
        :param h:
        :param t:
        :return:
        """
        self._step_num_solution_container.append(step)
        self._t_solution_container.append(t)

        #   rPi, rPj - position of spring ends on each body in GCS
        self._r_P_solution_container.append(self.r_P_list)

        #   spring length
        self._l_solution_container.append(self.l)

        #   spring force
        self._F_solution_container.append(self.F_s)

        #   track force
        for force in self._Fn_list:
            force._track_data(step, h, t)

    def evaluate_Q_e(self, q):
        """
        Function evaluates a generalized external force vector
        Defined in subclass
        :return:
        """
        return []

    def evaluate_potential_energy(self, q):
        """
        Defined in subclass
        :param q:
        :return:
        """
        Ep = 0.5 * self.k * np.linalg.norm(self.r_P_list[0] - self.r_P_list[1], ord=2)**2
        return Ep

    def evaluate_F(self, q):
        """

        :param q:
        :return:
        """
        return 0.

    def _evaluate_F(self, l, dl):
        """

        :return:
        """
        return self.k * (l - self.l_0) + self.c * dl + self.F0

    def _evaluate_r_ij_P(self, q):
        """

        :param q:
        :return:
        """
        for i, (body_id, uP) in enumerate(zip(self.body_id_list, self.u_P_LCS_list)):
            self.r_P_list[i] = q2R_i(q, body_id) + Ai_ui_P_vector(uP, q2theta_i(q, body_id))

        [self.r_P_i, self.r_P_j] = self.r_P_list

        #   distance vector
        r_ij_P = self.r_P_i - self.r_P_j

        return r_ij_P

    def _evaluate_r_ij_P_q_list(self, q):
        """

        :param q:
        :return:
        """
        spring_force_matrix_list = []
        for body_id, _u_P in zip(self.body_id_list, self.u_P_LCS_list):
            if body_id == "ground":
                spring_force_matrix_body = Force_Q_e_matrix(body_id, _u_P, 0)

            else:
                _theta = q2theta_i(q, body_id)
                spring_force_matrix_body = Force_Q_e_matrix(body_id, _u_P, _theta)

            #    append spring force matrix object to list
            spring_force_matrix_list.append(spring_force_matrix_body)

        return spring_force_matrix_list

    def reset(self, q):
        """
        Function defined in subclass
        :return:
        """
        for i, (body_id, uP, node_id) in enumerate(zip(self.body_id_list, self.u_P_LCS_list, self.node_id_list)):
            self.r_P_list[i] = u_P_lcs2gcs(uP, q, body_id, node_id=node_id)

    def evaluate_Q_q(self, q):
        """
        Function defined in subclass
        :return:
        """

    def evaluate_Q_dq(self, q):
        """
        Function defined in subclass
        :return:
        """

    def print_F(self):
        """
        Function defined in subclass
        :return:
        """
        print "F ="
        print self._evaluate_F(self.l, self.dl)

