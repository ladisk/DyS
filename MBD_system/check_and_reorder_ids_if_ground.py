'''
Created on 19. maj 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''


def check_and_reorder_ids(body_i, body_j, u_iP, u_jP):
    """
    Function checks if body i or j or ground and creates and reorders list of body ids (i, j), that ground body is always body j
    Args:
    body_i - body i id
    body_j - body j id
    Returns:
    id_list - list of body ids
    """
    
    if body_i == "ground":
        
        #    swap body ids
        body_i_ = body_j
        body_j_ = body_i
        #    swap u_P vectors
        u_iP_ = u_jP
        u_jP_ = u_iP
    else:
        body_i_ = body_i
        body_j_ = body_j
        
        u_iP_ = u_iP
        u_jP_ = u_jP
        
    id_list_ = [body_i_, body_j_]
    
    return body_i_, body_j_, id_list_, u_iP_, u_jP_

if __name__ == "__main__":
    import numpy as np
    print check_and_reorder_ids(body_i = "ground", body_j = 85, u_iP = np.array([1, 2]), u_jP = np.array([3, 4]))
