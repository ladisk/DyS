'''
Created on 18. feb. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''

try:
    from ....Ai_ui_P import Ai_ui_P_vector
except:
    None

def offset(vertices, CS):
    """
    offsets the vertices of body shape that the body CM (center of mass) has the coordinates = (0, 0, 0)
    coordinates of CM are written in file body_name_data.txt
    in:
        vertices - matrix of points of body
        CS - array of centre of body mass for
    out:
        vertices - matrix of points of body shifted that the body CM has coordinates = (0, 0, 0)
    """

    vertices_shifted = vertices - CS

    return vertices_shifted