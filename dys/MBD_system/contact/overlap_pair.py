__author__ = 'lskrinjar'


import itertools
from operator import attrgetter
from MBD_system.contact.distance.distance import Distance


class Overlap_pair(object):
    '''

    '''
    __id = itertools.count(1)

    def __init__(self, _AABB_i=None, _AABB_j=None, parent=None):
        '''

        '''
        self._parent = parent
        self._AABB_i = _AABB_i
        self._AABB_j = _AABB_j

        #    attribute - list for distance objects
        self.distance_list = []
        self.distance_min = None

        #    add each AABB object to list
        self._AABB_list = [_AABB_i, _AABB_j]

        #    default value
        self.overlap = True  # because the overlap object is created only if two AABB overlap
        self.contact = False


    def check_for_contact_2D(self, q):
        """

        """
        #    update contact nodes and normals in contact object
        self.__update_contact_nodes_normals_2D(q)

        #    predefine empty distance list
        self.distance_list = []

        #    calculate
        self.__distance_node2edge_AABBi_AABBj(self._AABB_i, self._AABB_j)
        self.__distance_node2edge_AABBi_AABBj(self._AABB_j, self._AABB_i)

        self.distance_min = min(self.distance_list, key=attrgetter('_distance'))


    def __update_contact_nodes_normals_2D(self, q):
        """

        """
        #    update nodes GCS position in every AABB inside pair object
        for _AABB in self._AABB_list:
            _AABB.update_nodes_and_normals_GCS_in_AABB_2D(q)


    def __distance_node2edge_AABBi_AABBj(self, AABB_i, AABB_j):
        """
        Function calculates a distance from node to line for all nodes and edges in contact AABB overlap pair
        """
        #   pointers of edge and node body as object attributes
        self.edge_body = AABB_j._parent_body
        self.node_body = AABB_i._parent_body


        for node in AABB_i.nodes_GCS_in_AABB:

            for i, (node_j, normal_) in enumerate(zip(AABB_j.nodes_GCS_in_AABB, AABB_j.normals_GCS_in_AABB)):

                if i < len(AABB_j.nodes_GCS_in_AABB) - 1:
                    _edge_node_1 = node_j
                    _edge_node_2 = AABB_j.nodes_GCS_in_AABB[i + 1]
                    _edge_normal = normal_

                    _dist = Distance(node, _edge_node_1, _edge_node_2, normal_, self.node_body, self.edge_body, self._parent)

                    self.distance_list.append(_dist)


    def build_2D_contact_geometry(self, bodies=[]):
        """

        """
        for body_id_ in self.body_id_list:
            #    assign body object from body list to new variable name
            body_ = bodies[body_id_]
            body_.contact_geometry = ContactGeometry(parent=body_)

            #    add only geometry (triangles and nodes) that meet conditions
            for triangle in body_.geom.geom_data.triangles:
                #    remove faces parallel to contact plane
                if triangle.normal[2] == 1 or triangle.normal[2] == -1:
                    None
                elif not triangle.intersects_plane(self.position_of_plane_in_z_direction):
                    None
                else:
                    for node in triangle.nodes:
                        node.to_2D()
                        body_.contact_geometry.addNode(node)
                    body_.contact_geometry.addTriangle(triangle)

            #    create contact geometry profile as variable (attribute) contact_nodes
            body_.contact_geometry.create_contact_profile()