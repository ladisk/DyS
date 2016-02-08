"""
Created on 17. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""

import itertools
import numpy as np


class Joint_C_q_matrix(object):
    """
    classdocs
    """
    __id = itertools.cycle([0, 1])
    
    
    def __init__(self, joint_type_, body_connected_to_ground = False):
        """
        Constructor of joint C_q matrix class for each joint the constructor is called twice
        Args:
        joint_type_ - joint type
        
        Returns:
        none
        
        Raises:
        none
        """
        #    joint_C_q_matrix counter:
        #    0 - for body i
        #    1 - for body j
        if body_connected_to_ground:
            self.id = 0
        else:
            self.id = self.__id.next()
        
        
        if joint_type_ == "fixed":
            self.matrix = np.eye(3)
            if self.id == 0:
                self.matrix = self.matrix
            elif self.id == 1:
                self.matrix = -1 * self.matrix
            else:
                raise ValueError, "Matrix for fixed joint not constructed."

        elif joint_type_ == "revolute":
            if self.id == 0:
                self.matrix = np.hstack((np.eye(2), np.zeros([2, 1])))
            elif self.id == 1:
                self.matrix = -1 * np.hstack((np.eye(2), np.zeros([2, 1])))
            else:
                raise ValueError, "Matrix for revolute joint not constructed."

        elif joint_type_ == "prismatic":
            self.matrix = np.zeros([2, 3])

        else:
            raise AttributeError, "Joint type not correct!"
            

if __name__ == "__main__":
    from pprint import pprint
    a = Joint_C_q_matrix(joint_type_="revolute")
    pprint(vars(a))
    b = Joint_C_q_matrix(joint_type_="fixed")
    pprint(vars(b))
    
