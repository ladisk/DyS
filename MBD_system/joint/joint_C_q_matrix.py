"""
Created on 17. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""


from pprint import pprint
import itertools
import numpy as np


from global_variables import GlobalVariables


class Joint_C_q_matrix(object):
    """
    classdocs
    """
    __id = itertools.cycle([0, 1])
    
    
    def __init__(self, joint_type, body_connected_to_ground=False, body_id=None, parent=None):
        """
        Constructor of joint C_q matrix class for each joint the constructor is called twice
        Args:
        joint_type_ - joint type
        
        Returns:
        none
        
        Raises:
        none
        """
        #   parent pointer
        self._parent = parent

        #    joint_C_q_matrix counter:
        #    0 - for body i
        #    1 - for body j
        if body_connected_to_ground:
            self.id = 0
        else:
            self.id = self.__id.next()

        if joint_type == "fixed":
            self.matrix = np.eye(3)
            if self.id == 0:
                self.matrix = self.matrix
            elif self.id == 1:
                self.matrix = - self.matrix
            else:
                raise ValueError, "Matrix for fixed joint not constructed."

        elif joint_type == "revolute":
            if self.id == 0:
                self.matrix = np.hstack((np.eye(2), np.zeros([2, 1])))
            elif self.id == 1:
                self.matrix = -1 * np.hstack((np.eye(2), np.zeros([2, 1])))
            else:
                raise ValueError, "Matrix for revolute joint not constructed."

        elif joint_type == "prismatic":
            self.matrix = np.zeros([2, 3])

        elif joint_type == "fixed support":
            if type(body_id) is int:
                e_n = self._parent._parent._parent.bodies[body_id].mesh.n_NC

                n_NC = self._parent._parent._parent.bodies[body_id].mesh.n_NC
                if e_n == 8:
                    self.matrix = np.zeros([4, e_n])
                if e_n == 12:
                    self.matrix = np.zeros([6, e_n])
            else:
                self.matrix = None

        elif joint_type == "hinged support":
            if type(body_id) is int:
                e_n = self._parent._parent._parent.bodies[body_id].e_n

                q_i_dim = self._parent._parent._parent.bodies[body_id].q_i_dim
                if e_n == 8:
                    self.matrix = np.zeros([4, q_i_dim])
                if e_n == 12:
                    self.matrix = np.zeros([6, q_i_dim])
            else:
                self.matrix = None

        elif joint_type == "fixed joint point mass-rigid":
            if type(body_id) is int:
                q_i_dim = GlobalVariables.q_i_dim[body_id]
                self.matrix = np.zeros([2, q_i_dim])

            else:
                self.matrix = None
        else:
            raise AttributeError, "Joint type not correct!"
            

if __name__ == "__main__":
    from pprint import pprint
    a = Joint_C_q_matrix(joint_type_="revolute")
    pprint(vars(a))
    b = Joint_C_q_matrix(joint_type_="fixed")
    pprint(vars(b))
    
