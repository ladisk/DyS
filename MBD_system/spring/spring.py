'''
Created on 21. feb. 2014

@author: lskrinjar
'''
import time
import itertools

import numpy as np

try:
    from ..MBD_system import *
    from ..cad2cm_lcs import cad2cm_lcs
    from ..r_ij_P import r_ij_P
    from ..force_matrix import Force_Q_e_matrix
    from ..I_A_uP_i import I_A_uP_i_matrix
    from ..check_and_reorder_ids_if_ground import check_and_reorder_ids
    from ..q2R_i import q2R_i
    from ..q2theta_i import q2theta_i
    from ..q2dR_i import q2dR_i
    from ..q2dtheta_i import q2dtheta_i
    from ..dr_ij_P_dq import dr_ij_P_dq
    from ..force_vector import Force_Q_e_vector
    from ..u_P_cad2cm_lcs import u_P_cad2cm_lcs
except:
    None

from MBD_system.force.force import Force
from MBD_system.MBD_system_items import SpringItem
from simulation_control_widget.opengl_widget.marker.marker import Marker
from MBD_system.transform_cs import gcs2lcs_z_axis


class Spring(SpringItem):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, _name, spring_type, body_id_i, body_id_j, u_iP_CAD, u_jP_CAD, k=0, c=0, l_0=None, parent = None):
        """
        Constructor of a spring class
        :param :
            _name - string
            spring_type - translational, rotational
            k - spring stiffness
            c - damping coefficient
            body_i - 
            body_j - 
            u_iP - in CAD LCS of a body
            u_jP - in CAD LCS of a body
            l_0 - undeformed length of a spring
            body_id - number of body on which the force is applied
            F_x - x component of force
            F_y . y component of force
            u_iP_f - vector of acting force in CAD LCS of a body
            
        """
        super(Spring, self).__init__(_name, parent)
        #   parent
        self._parent = parent

        #    number
        self.spring_id = self.__id.next()


        #    name as string
        self._name = _name
        
        #    spring type
        self.spring_type = spring_type

        #   position of spring in GCS
        self.z_dim = 0

        #   physical properties
        self.k = k
        self.c = c

        if self.spring_type == "translational":
            #    spring linear - params

            self.l_0 = l_0
        elif self.spring_type == "torsional":
            #    spring torsional - params
            # self.k_t = k_t
            # self.c_t = c_t
            pass

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

        self.body_id_list = [self.body_id_i, self.body_id_j]

        #   list of point vectors to joint constraint in CAD lCS of a body
        self.u_P_CAD_list = [self.u_iP_CAD, self.u_jP_CAD]

        #   predefined empty list to store point vectors of joint constraint in LCS (center of gravity)
        #   of each body
        self.u_P_list = []

        self.signs = np.array([-1, +1])

        #   predefined empty list of spring force matrix to store one force matrix for each body connected with a spring
        self.spring_force_matrix_list = []
        #   pair of contact force list
        self.force_list = []
        self.markers = []

        for body_id, _u_P in zip(self.body_id_list, self.u_P_CAD_list):
            if body_id == "ground" or body_id == -1:
                u_P_LCS = u_P_cad2cm_lcs(body_id, self._parent._parent.ground, _u_P=_u_P)
                _body = self._parent._parent.ground
            else:
                #   create pointer to body
                _body = self._parent._parent.bodies[body_id]

                #   calculate point vector in body LCS (center of gravity)
                u_P_LCS = u_P_cad2cm_lcs(body_id, _body, _u_P=_u_P)

                #   create force object
                _force = Force(body_id, force_name=self._name + "_on_body_" + str(body_id))
                #   add pair of contact forces to forces list of MBD system
                self._parent._parent.forces.append(_force)

                #   add pair of contact forces to forces list of spring
                self.force_list.append(_force)
                self._parent._parent.bodies[body_id].forces.append(_force)

            #   create markers for u_P points on every body in kinematic constraint
            z_dim_lcs = gcs2lcs_z_axis(_body.R[2], self.z_dim)
            _node = np.array(np.append(u_P_LCS, z_dim_lcs), dtype='float32')
            _marker = Marker(_node, visible=True, parent=_body)

            if body_id == "ground" or body_id == -1:
                self._parent._parent.ground.markers.append(_marker)
            else:
                self._parent._parent.bodies[body_id].markers.append(_marker)

            self.markers.append(_marker)
            self.u_P_list.append(u_P_LCS)

        [self.u_iP, self.u_jP] = self.u_P_list

        self.N_b_ = len(self._parent._parent.bodies)
        self.additional_params_calulated = True

        
        
    def create_spring_force_Q_e_vector(self, bodies=[], q=[]):
        """
        Function calculates force vector from force object properties and vector q (only angle theta is used from q vector)
        in:
            q - body object
        out:
            vector of force acting on a body (Fx, Fy, Mz)
        """
        if self.spring_type == "translational":
            force = self.spring_translational(q)
            return force
        elif self.spring_type == "torsional":
            force = self.spring_rotational(q)
            return force
        else:
            raise IOError, "Spring type not correct"
        
        
    def spring_translational(self, q):
        """

        :param bodies_:
        :param q:
        :return:
        """
        #   predefined list for spring force matrix
        self.spring_force_matrix_list = []
        # theta_ = np.zeros(2)
        # R_ = np.zeros(4)s
        # dq_ = np.zeros(6)
        for body_id, _u_P in zip(self.body_id_list, self.u_P_list):

            if body_id == "ground":
                spring_force_matrix_body = Force_Q_e_matrix(body_id)
                # R_body = np.zeros(2)
                # dq_body = np.zeros(3)
            else:
                _theta = q2theta_i(q, body_id)
                spring_force_matrix_body = Force_Q_e_matrix(body_id, _u_P, _theta)
                R_body = q2R_i(q, body_id)
                dq_body = np.hstack((q2dR_i(q, body_id, self.N_b_), q2dtheta_i(q, body_id)))
            

            #    append spring force matrix object to list
            # print "spring_force_matrix_body ="
            # print spring_force_matrix_body.matrix_
            self.spring_force_matrix_list.append(spring_force_matrix_body)
            # #    append body velocities to create array of velocities of both bodies connected with spring
            # dq_[3 * self.body_id_list.index(body_id):3 * self.body_id_list.index(body_id) + 3] = dq_body
            # R_[2 * self.body_id_list.index(body_id):2 * self.body_id_list.index(body_id) + 2] = R_body
            # theta_[self.body_id_list.index(body_id)] = dq_body[-1]
            

        # [R_body_i_, R_body_j_] = np.split(R_, 2)

        # [theta_i_, theta_j_] = theta_

        R_body_i_ = q2R_i(q, self.body_id_i)
        theta_i_ = q2theta_i(q, self.body_id_i)

        R_body_j_ = q2R_i(q, self.body_id_j)
        theta_j_ = q2theta_i(q, self.body_id_j)


        dq_ = np.array([np.append(R_body_i_, theta_i_), np.append(R_body_j_, theta_j_)]).flatten()
        # print "dq_ =", dq_


        # print "self.l_0 =", self.l_0
        #    if undeformed length of spring is not specified, the undeformed length is calculated from initial position
        if self.l_0 is None:
            r_ij_P_0 = r_ij_P(R_i=R_body_i_, theta_i=theta_i_, u_iP=self.u_iP, R_j=R_body_j_, theta_j=theta_j_, u_jP=self.u_jP)
            self.l_0 = np.linalg.norm(r_ij_P_0, ord=2)
        
        #    vector between points of installed spring element
        r_ij_P_ = r_ij_P(R_i=R_body_i_, theta_i=theta_i_, u_iP=self.u_iP, R_j=R_body_j_, theta_j=theta_j_, u_jP=self.u_jP)

        #    length of vector - spring length
        l = np.linalg.norm(r_ij_P_, ord=2)
        # print "L(t) =", l
        #    unit vector in direction of vector r_ij_P_
        I_r_ij = r_ij_P_ / l
        # print "I_r_ij =", I_r_ij
        #    dr_ij_P_dq
        dr_ij_P_dq_ = dr_ij_P_dq(body_id_i=self.body_id_i, theta_i=theta_i_, u_iP=self.u_iP, body_id_j=self.body_id_j, theta_j=theta_j_, u_jP=self.u_jP)
 
        #    velocity of deformation of spring length
        dl = np.dot(I_r_ij, np.dot(dr_ij_P_dq_, dq_))

        
        #    force value (amplitude) of spring element
        # print "k =", self.k
        # print "dx =", l - self.l_0
        # print self.k * (l - self.l_0)
        f_s = self.k * (l - self.l_0) + self.c * dl
        # print "f_s =", f_s
        #    generalized force on each body that are connected with a spring
        
        self.spring_force_vector_list = []
        for body_id, force_matrix in zip(self.body_id_list, self.spring_force_matrix_list):
            # print "--------------------------------"
            if body_id == "ground":
                spring_force_vector_body = Force_Q_e_vector()
            else:
                # print "force_matrix.matrix_ ="
                # print force_matrix.matrix_
                spring_force_vector_body = Force_Q_e_vector(f_s, force_matrix.matrix_, I_r_ij)
            #    append spring force vector object to list
            # print "id =", spring_force_vector_body.id
            self.spring_force_vector_list.append(spring_force_vector_body)
            

        [Q_i, Q_j] = self.spring_force_vector_list

        # print "Q_i.vector =", Q_i.vector
        # print "Q_j.vector =", Q_j.vector
        # time.sleep(100)
        return Q_i.vector, Q_j.vector
    
    
    def spring_rotational(self, bodies_=[], q=[]):
        if not self.additional_params_calulated:
            
            
            self.N_b_ = len(bodies_)
            #    relative angular displacement between two bodies
            self.theta_0 = q2theta_i(q, self.body_id_i) - q2theta_i(q, self.body_id_j)
            
            self.additional_params_calulated = True

        
        theta_ = np.zeros(2)
        dtheta_ = np.zeros(2)
        for body_id in self.body_id_list:
            if body_id == "ground":
                theta_body = 0
                dtheta_body = 0
            else:
                theta_body = q2theta_i(q, body_id)
                dtheta_body = q2dtheta_i(q, body_id)
                
            #    append theta of a body to array of both body angular velocities
            theta_[self.body_id_list.index(body_id)] = theta_body
            dtheta_[self.body_id_list.index(body_id)] = dtheta_body
        
        theta = theta_[0] - theta_[1]
        dtheta = dtheta_[0] - dtheta_[1]
        

        
        
        torsional_force = self.k_t * (theta - self.theta_0) + self.c_t * dtheta
        
        self.spring_force_vector_list = []
        spring_force_vector_body = Force_Q_e_vector()
        for body_id in self.body_id_list:
            spring_force_vector_body = Force_Q_e_vector()
            
            if body_id != "ground":
                spring_force_vector_body.vector[-1] = self.signs[self.body_id_list.index(int(body_id))] * torsional_force
            #    append spring force vector object to list
            self.spring_force_vector_list.append(spring_force_vector_body)
            
        
        [Q_i, Q_j] = self.spring_force_vector_list
        return Q_i.vector, Q_j.vector
    
    
