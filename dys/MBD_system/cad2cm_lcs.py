'''
Created on 7. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''


def cad2cm_lcs(vector_in_CAD_CS, CM_CAD_LCS):
    """
    Function transforms vector from cad-LCS to center of mass CM-LCS of a body
    Args: 
        vector_in_CAD_LCS - vector in cad-CS
        CM_CAD_LCS - CM vector in cad-CS
    Returns:
        vector_in_CM_LCS - vector in CM-CS
    Raises:
        none
    """
    vector_in_CM_LCS = vector_in_CAD_CS - CM_CAD_LCS
    return vector_in_CM_LCS