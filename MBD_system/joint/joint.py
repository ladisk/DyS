"""
Created on 18. mar. 2014

@author: lskrinjar
"""
import itertools
from pprint import pprint
import numpy as np


from MBD_system.MBD_system_items import JointItem
from MBD_system.force.force import Force
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.r_ij_P import r_ij_P
from MBD_system.solution_data.solution_data_joint import SolutionDataJoint
from MBD_system.u_P_cad2cm_lcs import u_P_cad2cm_lcs
from MBD_system.transform_cs import cm_lcs2gcs
from MBD_system.spring.spring_rotational import SpringRotational
from MBD_system.spring.spring_translational import SpringTranslational


class Joint(JointItem):
    """
    classdocs
    """
    __id = itertools.count(0)

    #   supported joint types
    _supported_types = ["fixed",
                        "prismatic",
                        "revolute",
                        "revolute joint flexible-flexible",
                        "revolute joint rigid-flexible",
                        "rigid joint flexible-flexible",
                        "rigid joint point mass-flexible",
                        "rigid joint point mass-flexible with rotation",
                        "rigid joint rigid-flexible",
                        "fixed support",
                        "hinged support",
                        "roller support",
                        "slope discontinuity",
                        "fixed joint point mass-rigid"]

    def __init__(self, joint_type, body_id_i, body_id_j, u_iP_CAD=np.array([0, 0]), u_jP_CAD=np.array([0, 0]), u_iQ_CAD=np.array([0, 0]), properties_dict={}, parent=None):
        """
        Create a joint object
        :param joint_type:      type of 2D joint (revolute, prismatic, fixed) as string
        :param body_id_i:       id of body i
        :param body_id_j:       id of body j
        :param u_iP_CAD:        in CAD LCS of a body i
        :param u_jP_CAD:        in CAD LCS of a body j
        :param u_iQ_CAD:        in CAD LCS of a body i (only for prismatic joint)
        :param parent:
        """
        super(Joint, self).__init__(joint_type, parent)
        #    number
        self.joint_id = self.__id.next()  # len(joints_list) + 1
        self._parent = parent

        if self._parent is not None:
            self.joints_list = self._parent._parent.joints
        else:
            self.joints_list = []

        #   joint type
        self.joint_type = joint_type

        #   dictionary of joint properties
        self._dict = properties_dict

        #   boolean for constant C_q matrix of joint
        self.constant = False

        #   visualization (opengl) properties
        self.z_dim = 0.
        #   get q vector of current MBD system at joint object initialization
        self.q0 = self._parent._parent.evaluate_q0()

        #   spring at joint coordinates
        self.spring = None


        if self.joint_type not in self._supported_types:
            raise ValueError, "Joint type not correct!"

        #   body ids
        self.body_id_i = body_id_i
        self.body_id_j = body_id_j
        #   check ids
        if self.joint_type != "slope discontinuity":
            self._check_ids(self.body_id_i, self.body_id_j, uiP=u_iP_CAD, ujP=u_jP_CAD)

        #   body id list
        self.body_id_list = [self.body_id_i, self.body_id_j]
        self.body_list = []

        #   body type list
        self.body_type_list = ["rigid", "rigid"]

        #   flexible body data if a flexible body is in joint
        self.node_id_i = None
        self.node_id_j = None
        self.node_id_list = [self.node_id_i, self.node_id_j]

        self.element_id_i = None
        self.element_id_j = None
        self.element_id_list = [self.element_id_i, self.element_id_j]

        self.element_ksi_i = None
        self.element_ksi_j = None
        self.element_ksi_list = [self.element_ksi_i, self.element_ksi_j]

        #   list of point vectors to joint constraint in CAD lCS of a body
        self.u_P_CAD_list = [u_iP_CAD, u_jP_CAD]

        #   predefined empty list to store point vectors of joint constraint in LCS (center of gravity)
        #   of each body
        self.u_P_LCS_list = []
        self.r_P_GCS_list = [None, None]

        self.u_iQ_CAD = u_iQ_CAD

        #   predefined variable of vector C(q, t) constraint equation
        self.C0 = None
        self.theta0 = None

        #   time dependent attributes and values
        self.C = None

        # if self.body_id_i == "ground" or self.body_id_j == "ground":
        #     _body_id = copy(self.body_id_list)
        #     _body_id.remove("ground")
        #     _body_id = _body_id[0]
        #
        #     indx = self.body_id_list.index(_body_id)
        #     C1 = q2R_i(self.q0, _body_id) + Ai_ui_P_vector(self.u_P_CAD_list[indx], q2theta_i(self.q0, _body_id))
        #
        #     C2 = q2theta_i(self.q0, self.body_id_list[indx])
        #
        #     self.C0 = np.append(C1, C2)
        # else:
        #     self.C0 = np.zeros(3)

        #   predefined empty list
        self.C_q_list = []
        self.Q_d_list = []

        #   list of joint points
        #   1   uPi - body i
        #   2   uQi - body i
        #   3   uPj - body j
        self.u_QP_LCS_list = []
        
        #    vector of lagrange multipliers at simulation time
        self.L = None
        #    list of vectors of constraint reaction forces for each body in joint
        self.Q_c_list = [None, None]

        #   pair of contact force list
        self.contact_bodies_added_to_list = False
        self.list_of_contact_force_objects_constructed = False
        self.force_list = []

        #   calculate u_P vector of every body in body LCS
        for i, (body_id, _u_P) in enumerate(zip(self.body_id_list, self.u_P_CAD_list)):
            if body_id == "ground" or body_id == -1 or body_id is None:
                u_P_LCS = _u_P
                rP = None
                _body = self._parent._parent.ground
            else:
                #   pointer to body object
                _body = self._parent._parent.bodies[body_id]
                if _body.body_type == "rigid body":
                    #   calculate point vector in body LCS (center of gravity)
                    u_P_LCS = u_P_cad2cm_lcs(body_id, _body, _u_P, 0)#_body.theta[2]

                    R_i = q2R_i(self.q0, body_id)
                    theta_i = q2theta_i(self.q0, body_id)
                    rP = cm_lcs2gcs(u_P_LCS, R=R_i, theta=theta_i)
                else:
                    u_P_LCS = None
                    rP = None

            self.u_P_LCS_list.append(u_P_LCS)
            self.r_P_GCS_list[i] = rP

            #   list of body pointers
            self.body_list.append(_body)

        [self.u_iP_LCS, self.u_jP_LCS] = self.u_P_LCS_list

        #   solution options:
        #   Discard (default)
        #   Overwrite (existing file)
        #   Save to new (next available) file
        self._solution_save_options = "discard"
        #   solution file type, optional: .dat, .xlsx, .csv
        self._solution_filetype = ".xlsx"
        
        #   markers
        self.markers = []

        #   create forces
        self._create_Fn_forces()
        
        #   init solution container
        self._solution_containers()
        
        #    add additional properties from dictionary
        self.add_attributes_from_dict(properties_dict)

        #   visualization properties
        self.vtk_actor = None

    def _update_lists(self):
        """

        :return:
        """
        self.body_id_list = [self.body_id_i, self.body_id_j]

        self.node_id_list = [self.node_id_i, self.node_id_j]

        self.element_id_list = [self.element_id_i, self.element_id_j]

        self.element_ksi_list = [self.element_ksi_i, self.element_ksi_j]
    
    def _solution_containers(self):
        """
        
        """
        self.solution_data = SolutionDataJoint(name=self._name, MBD_item=self)

        self._step_num_solution_container = []
        self._t_solution_container = []
        
        self._L_solution_container = []
        
        self._Q_c_solution_container = []
    
    def _track_data(self, step, h, t):
        """

        :param step:
        :param h:
        :param t:
        :return:
        """
        self._step_num_solution_container.append(step)
        self._t_solution_container.append(t)
        
        #    lagrange multipliers vector
        self._L_solution_container.append(self.L)
        
        #    joint constraint force vector
        self._Q_c_solution_container.append(np.hstack(self.Q_c_list))
        
        #    track force objects
        for force in self._Fn_list:
            force._track_data(step, h, t)

    def write_to_file(self):
        """

        :return:
        """
        #    get joint object data
        data = np.hstack((np.array([self._step_num_solution_container]).T,
                          np.array([self._t_solution_container]).T,
                          np.vstack(self._L_solution_container),
                          np.vstack(self._Q_c_solution_container)))
        
        #   create solution data object
        # solution_data = SolutionDataJoint(name=self._name, parent=None)

        #   add solution data to object
        self.solution_data.add_data(data)
        
        #    write data to file
        self.solution_data.write_to_file()

    def _check_ids(self, i, j, uiP=None, ujP=None):
        """
        Function checks if body id i and body id j are different
        :return:
        """
        if i != j:
            pass
        else:
            raise ValueError, "Body i id and body j id are equal!"

        # swap body_id that if body is connected to ground that ground is always the last item in list
        if i == "ground":
            self.body_id_i = j
            self.body_id_j = i
            self.u_iP_CAD = ujP
            self.u_jP_CAD = uiP

        else:
            self.body_id_i = i
            self.body_id_j = j
            self.u_iP_CAD = uiP
            self.u_jP_CAD = ujP

    def _create_Fn_forces(self):
        """
        Function creates list of contact forces as objects
        :return:
        """
        self.Fn = 0
        self._Fn_list = []
        if self._parent is not None:
            for body_id in self.body_id_list:
                if body_id != "ground" and body_id is not None:
                    _Fn = Force(body_id, force_name=self._name + "_Fn_on_body_" + str(body_id), parent=self)
                    #   add pair of contact forces to forces list of MBD system
                    self._parent._parent.forces.append(_Fn)
                    #   add pair of contact forces to forces list of contact
                    self._Fn_list.append(_Fn)
                    self._parent._parent.bodies[body_id].forces.append(_Fn)
            self._Fn_list[0]._visible = False
            self.list_of_contact_force_objects_constructed = True

    def _add_spring(self, _dict):
        """

        :return:
        """
        name = _dict["spring_type"] + "_spring_at_" + self._name
        if _dict["spring_type"] == "rotational":
            spring = SpringRotational(name, self.body_id_i, self.body_id_j, u_iP_LCS=self.u_iP_LCS, k_t=_dict["k_t"], c_t=_dict["c_t"], properties_dict=self._dict, parent=self)

        elif _dict["spring_type"] == "translational":
            spring = SpringTranslational()

        else:
            spring = None

        #   append spring to spring list of MBD system object
        if spring is not None:
            if hasattr(self._parent._parent, "springs"):
                self._parent._parent.springs.append(spring)

        return spring

    def transform_cad2cm_lcs(self, points):
        """
        Function tranforms body points from body CAD CS to body center of mass LCS
        :return:
        """

    def _create_markers(self):
        """
        Function creates markers, defined in subclass
        :return:
        """
        markers = []
        return markers

    def evaluate_C(self, q, t=0.):
        """
        Function defined in subclass
        :param q:
        :param t:
        :return:
        """

    def evaluate_C0(self, q, t=0.):
        """
        Function defined in subclass
        :param q:
        :param t:
        :return:
        """

    def evaluate_rijP(self, q):
        """

        :param q:
        :return:
        """
        return r_ij_P(q2R_i(q, self.body_id_i), q2theta_i(q, self.body_id_i), self.u_iP_LCS, q2R_i(q, self.body_id_j), q2theta_i(q, self.body_id_j), self.u_jP_LCS)

    def evaluate_C_t(self, q):
        """

        :param q:
        :return:
        """
        C_t = np.zeros(3 - self.joint_DOF)
        return C_t

    def evaluate_d(self, q):
        """

        :return:
        """
        rijP = self.evaluate_rijP(q)
        d = np.linalg.norm(rijP)
        return d
    
    def evaluate_Q_c(self, L):
        """
        Function evaluates 
        :param L:    vector of lagrange multipliers of a joint
        """
        self.L = L
        for i, C_q in enumerate(self.C_q_list):
            self.Q_c_list[i] = -np.dot(C_q.matrix.T, L)

    def evaluate_Q_d(self, q):
        """

        :param q:
        :return:
        """
        return None

    def evaluate_C_q(self, q):
        """

        :param q:
        :return:
        """
        return None

    def reset(self, q0):
        """

        :param q0:
        :return:
        """
        self.q0 = q0

        for marker in self.markers:
            marker.theta = marker.theta0

        self._reset(self.q0)

    def _reset(self, q):
        """
        Defined in subclass
        :return:
        """

    def set_vtk_data(self):
        """

        :param q:
        :return:
        """

    def update_vtk_data(self, q):
        """

        :param q:
        :return:
        """

    def evaluate_C_q_L_q(self, L):
        """

        :return:
        """



if __name__ == "__main__":
    #     a = Joint(joint_type = "revolute", body_id_i = 1, body_id_j = 2, u_iP_CAD = np.array([1, 2]), u_jP_CAD = np.array([33, 44]))
    a = Joint(joint_type="prismatic", body_id_i=1, body_id_j=2, u_iP_CAD=np.array([1, 2]), u_jP_CAD=np.array([33, 44]),
              u_iQ=np.array([-2, -2]))
    print a
    print a.create_joint_C_q_matrix()
# b = Joint(joint_type = "fixed", body_id_i = 1, body_id_j = 2, u_iP_CAD = np.array([1, 2]), u_jP_CAD = np.array([33, 44]))
#     pprint(vars(b))
#     b.create_joint_C_q_matrix()
