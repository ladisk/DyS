__author__ = 'lskrinjar'

from pprint import pprint
import itertools
from operator import attrgetter
from MBD_system.contact.distance.distance import Distance

class OverlapPair(object):
    """

    """
    __id = itertools.count(1)

    def __init__(self, _AABB_i=None, _AABB_j=None, parent=None):
        """

        """
        self._parent = parent
        self.AABB_i = _AABB_i
        self.AABB_j = _AABB_j

        #    attribute - list for distance objects
        self.distance_list = []
        self.distance_min = None

        #    add each AABB object to list
        self._AABB_list = [_AABB_i, _AABB_j]

        #    default value
        self.overlap = True  # because the overlap object is created only if two AABB overlap
        self.contact = False
        self.contact_geometry_created = False

    def check_for_contact_2D(self, q):
        """
        Function checks for contact between two AABB objects that overlap
        """
        #    update contact nodes and normals in contact object
        self._update_contact_nodes_normals_2D(q)

        #    predefine empty distance list
        self.distance_list = []

        if not self.contact_geometry_created:
            #    calculate distance of node to edge
            self._distance_node2edge_AABBi_AABBj(q)
            self._distance_node2edge_AABBi_AABBj(q)

            self.contact_geometry_created = True

            #    find minimum distance
            self.distance_min = min(self.distance_list, key=attrgetter('_distance'))
        else:
            pass
    
    def update_contact_2D(self, q):
        """

        """

    def get_contact_nodes(self):
        """
        Function finds contact nodes that are present in overlap contact pair object
        """
        _contact_nodes_list = []
        for dist_obj in self.distance_list:
            if dist_obj._inside:
                print "----------------------"
                print "name (edge body) =", dist_obj.edge_body._name
                pprint(vars(dist_obj))
                _contact_nodes_list.append(dist_obj)
        return _contact_nodes_list

    def _update_contact_nodes_normals_2D(self, q):
        """

        """
        #    update nodes GCS position in every AABB inside pair object
        for _AABB in self._AABB_list:
            _AABB.update_nodes_and_normals_GCS_in_AABB_2D(q)

    def _distance_node2edge_AABBi_AABBj(self, q):
        """
        Function calculates a distance from node to edge for all nodes and edges in contact AABB overlap pair
        """
        #   pointers of edge and node body as object attributes
        self.edge_body = self.AABB_j._parent_body
        self.node_body = self.AABB_i._parent_body

        for node in self.AABB_i.nodes_LCS_in_AABB:
            for i, (node_j, normal_) in enumerate(zip(self.AABB_j.nodes_LCS_in_AABB, self.AABB_j.normals_LCS_in_AABB)):
                if i < len(self.AABB_j.nodes_LCS_in_AABB) - 1:
                    _edge_node_1 = node_j
                    _edge_node_2 = self.AABB_j.nodes_LCS_in_AABB[i + 1]
                    _edge_normal = normal_

                    #   create distance object
                    _dist = Distance(node, _edge_node_1, n2=_edge_node_2, normal=normal_, q=q, node_body=self.node_body, edge_body=self.edge_body, parent=self._parent)

                    #   append distance object to list
                    self.distance_list.append(_dist)
                    
    def _distance_node2edge_AABBi_AABBj_OLD(self, q):
        """
        Function calculates a distance from node to edge for all nodes and edges in contact AABB overlap pair
        """
        #   pointers of edge and node body as object attributes
        self.edge_body_id = self.AABB_j._parent_body.body_id
        self.node_body_id = self.AABB_i._parent_body.body_id

        for node in self.AABB_i.nodes_GCS_in_AABB:
            for i, (node_j, normal_) in enumerate(zip(AABB_j.nodes_GCS_in_AABB, AABB_j.normals_GCS_in_AABB)):
                if i < len(AABB_j.nodes_GCS_in_AABB) - 1:
                    _edge_node_1 = node_j
                    _edge_node_2 = self.AABB_j.nodes_GCS_in_AABB[i + 1]
                    _edge_normal = normal_

                    _dist = Distance(node, _edge_node_1, n2=_edge_node_2, normal=normal_, node_body=self.node_body_id, edge_body=self.edge_body_id, parent=self._parent)

                    self.distance_list.append(_dist)

    def build_2D_contact_geometry(self, bodies=[]):
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