__author__ = 'lskrinjar'
import time
from pprint import pprint
import itertools
from operator import attrgetter
import numpy as np
import matplotlib.pyplot as plt


from MBD_system.contact.distance.distance import Distance
from MBD_system.contact.distance.distance_line_node import DistanceLineNode
from MBD_system.body.geometry.color_negative import color_negative


class OverlapPair(object):
    """
    classdocs
    """
    __id = itertools.count(-1)

    def __init__(self, _AABB_i=None, _AABB_j=None, parent=None):
        """

        """
        #    parent
        self._parent = parent
        #    id
        self.AABB_id = self.__id.next()
        
        self.AABB_i = _AABB_i
        self.AABB_j = _AABB_j

        #    attribute - list for distance objects
        self.distance_list = []
        self.distance_list_min = []
        self.distance_min = None

        #    add each AABB object to list
        self._AABB_list = [_AABB_i, _AABB_j]

        #    default value
        self.overlap = True  # because the overlap object is created only if two AABB overlap
        self.contact = False
        self.contact_geometry_created = False

        #   visualizations properties
        self._visible = False
        #    color
        self.color = np.random.rand(3)

        for i, AABB in enumerate(self._AABB_list):
            AABB._update_color(self.color)

        # self.max_penetration_depth = self._parent.distance_TOL * 100.
        self.max_penetration_depth = 1E-3

    def check_for_contact_2D(self, q):
        """
        Function checks for contact between two AABB objects that overlap
        """
        #    update contact nodes and normals in contact object
        self._contact_geometry_GCS(q)

        if not self.contact_geometry_created:
            self.distance_list = []
            #    calculate distance of node to edge for all possible combination of nodes and edges in pair of AABBs
            self.distance_list += self._distance_node2edge_AABBi_AABBj()

            self.contact_geometry_created = True

        return self.distance_list
    
    def update_contact_2D(self, q):
        """
        Update
        """

    def get_contact_nodes(self):
        """
        Function finds contact nodes that are present in overlap contact pair object
        """
        _contact_nodes_list = []
        for dist_obj in self.distance_list:
            if dist_obj._inside:
                _contact_nodes_list.append(dist_obj)

        return _contact_nodes_list

    def _contact_geometry_GCS(self, q):
        """

        """
        #    update nodes GCS position in every AABB inside pair object
        for _AABB in self._AABB_list:
            _AABB._contact_geometry_GCS(q)

    def _get_AABB_by_body_id(self, body_id):
        """

        :param body_id:
        :return:
        """
        # Checked, working correctly - OK!
        if self.AABB_i.body_id == body_id:
            AABB = self.AABB_i
        elif self.AABB_j.body_id == body_id:
            AABB = self.AABB_j
        else:
            raise ValueError, "AABB with body_id not found!"
        return AABB
    
    def _contact_geometry_LCS(self, body_id, rP):
        """
        Function gets point coordinate in LCS based on GCS coordinates as input parameter
        :param body_id:
        :param rP:
        :return uP:
        """
        AABB = self._get_AABB_by_body_id(body_id)
        
        #    get index of point
        _indx = np.where(np.all(AABB.nodes_in_AABB_GCS==rP,axis=1) == True)[0][0]
        
        #    point in LCS
        uP = AABB.nodes_in_AABB_LCS[_indx, :]

        return uP

    def _contact_normal_LCS(self, body_id, n):
        """

        :param body_id:
        :param n:
        :return:
        """
        # Checked, working correctly - OK!
        AABB = self._get_AABB_by_body_id(body_id)

        #    get index of point
        _indx = np.where(np.all(AABB.normals_in_AABB_GCS==n,axis=1) == True)[0][0]

        #    point in LCS
        n_LCS = AABB.normals_in_AABB_LCS[_indx, :]
        return n_LCS

    def _contact_tangent_LCS(self, body_id, t):
        """

        :param body_id:
        :param n:
        :return:
        """
        # Checked, working correctly - OK!
        AABB = self._get_AABB_by_body_id(body_id)

        #    get index of point
        _indx = np.where(np.all(AABB.tangents_in_AABB_GCS==t,axis=1) == True)[0][0]

        #    point in LCS
        t_LCS = AABB.tangents_in_AABB_GCS[_indx, :]
        return t_LCS

    def _subdistance_node2edge_AABBi_AABBj(self, AABB_i, AABB_j):
        """

        :return:
        """
        distance_list = []

        #   body ids
        #   i - body id of free node
        #   j - body id of edge nodes
        _body_id_i = AABB_i.body_id
        _body_id_j = AABB_j.body_id
        # print "len(AABB_i.nodes_in_AABB_GCS) =", len(AABB_i.nodes_in_AABB_GCS)
        for i_free_node in xrange(0, len(AABB_i.nodes_in_AABB_GCS)):
            node = AABB_i.nodes_in_AABB_GCS[i_free_node]
            # print "i_free_node =", i_free_node
            #   check every node in AABB except last, as the last node is included in next AABB
            # if i_free_node < (len(AABB_i.nodes_in_AABB_GCS) - 1):
            for i, (node_j, normal_, tangent_, angle_) in enumerate(zip(AABB_j.nodes_in_AABB_GCS, AABB_j.normals_in_AABB_GCS, AABB_j.tangents_in_AABB_GCS, AABB_j.angles_in_AABB_LCS)):
                #   for every node except for edge between first and last node
                if i < len(AABB_j.nodes_in_AABB_GCS) - 1:
                    _edge_node_1 = node_j
                    _edge_node_2 = AABB_j.nodes_in_AABB_GCS[i + 1]
                    _theta_jP = angle_
                    _theta_jR = AABB_j.angles_in_AABB_LCS[i + 1]
                    _dist = DistanceLineNode(node,
                                             _edge_node_1,
                                             normal=normal_,
                                             r_jR=_edge_node_2,
                                             tangent=tangent_,
                                             body_id_i=_body_id_i,
                                             body_id_j=_body_id_j,
                                             theta_jP=_theta_jP,
                                             theta_jR=_theta_jR,
                                             parent=self)

                #   edge data for edge between first and last node
                else:
                    if AABB_j.profile_type=="closed":
                        _edge_node_1 = AABB_j.nodes_in_AABB_GCS[-1]
                        _edge_node_2 = AABB_j.nodes_in_AABB_GCS[0]
                        _theta_jP = AABB_j.angles_in_AABB_LCS[-1]
                        _theta_jR = AABB_j.angles_in_AABB_LCS[0]
                        _dist = DistanceLineNode(node,
                                                 _edge_node_1,
                                                 normal=normal_,
                                                 r_jR=_edge_node_2,
                                                 tangent=tangent_,
                                                 body_id_i=_body_id_i,
                                                 body_id_j=_body_id_j,
                                                 theta_jP=_theta_jP,
                                                 theta_jR=_theta_jR, parent=self)

                #   append distance object to list
                distance_list.append(_dist)

        return distance_list

    def _distance_node2edge_AABBi_AABBj(self):
        """
        Function calculates a distance from node to edge for all nodes and edges in contact AABB overlap pair
        """
        #   pointers of edge and node body as object attributes
        # self.edge_body = self.AABB_j._parent_body
        # self.node_body = self.AABB_i._parent_body

        #   create intertools cycle object to
        _AABB_list = itertools.cycle(self._AABB_list)
        #   call next() only once before loop
        _AABB_list.next()
        #   loop through pair of overlaped AABBs
        #   in the first step
        #   AABBi = AABBi and AABBj = AABBj
        #   in the second (last) step
        #   AABBi = AABBj and AABBj = AABBi
        #   so method in for loop can only be called once but at every step list of distance objects are returned and saved together to one list
        distance_list = []

        if self._parent.body_id_edge is None and self._parent.body_id_node is None:
            for AABB in self._AABB_list:
                AABB_i = AABB
                AABB_j = _AABB_list.next()

                _dist = self._subdistance_node2edge_AABBi_AABBj(AABB_i, AABB_j)
                #   append distance list to object attribute list
                distance_list += _dist

        else:
            self.body_id_list.index("bar")
            AABB_i = self._AABB_list[self.body_id_list.index(self._parent.body_id_edge)]
            AABB_j = self._AABB_list[self.body_id_list.index(self._parent.body_id_node)]
            distance_list = self._subdistance_node2edge_AABBi_AABBj(AABB_i, AABB_j)

        return distance_list

    def build_2D_contact_geometry(self, bodies):
        """

        """
        for body_id_ in self.body_id_list:
            #    assign body object from body list to new variable name
            body_ = bodies[body_id_]
            # body_.contact_geometry = ContactGeometry(parent=body_)

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