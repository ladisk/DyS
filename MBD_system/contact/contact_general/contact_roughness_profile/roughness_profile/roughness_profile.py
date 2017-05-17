"""

created by: lskrinjar, aturel
date of creation: 27/03/2016
time of creation: 17:34
"""
import os
import numpy as np
from pprint import pprint
import itertools
import time
from matplotlib import pyplot as plt
import logging
import logging.handlers


from MBD_system.MBD_system_items import TreeItem
from MBD_system.t2n import t2n
from MBD_system.body.geometry.geometry_2D_profile import ContactGeometry2DProfile, ContactGeometry2DClosedProfile, ContactGeometry2DOpenedProfile
from MBD_system.extract_from_dictionary_by_string_in_key import extract_from_dictionary_by_string_in_key
from MBD_system.string2array import string2array


class RoughnessProfile(object):
    """
    classdocs
    """
    __id_list = itertools.cycle(["i", "j"])


    def __init__(self, body_id, uPi=np.zeros(2), uRi=np.zeros(2), dx=1.E-2, theta=0., Ra=0., Rz=0., thetamax=15., filename=None, scale=1., _dict={}, body=None, parent=None):
        """

        [1] Garcia, N.; Stoll, E.: "Monte Carlo Calculation of Electromagnetic-Wave Scattering from Random Rough Surfaces", Physical Review Letters, Volume 52, Issue 20, pp. 1798-1801 (1984).
        [2] FFTW library - free collection of fast C routines for computing discrete Fast Fourier Transforms. Developed at MIT by Matteo Frigo and Steven G. Johnson.
        [3] http://www.mysimlabs.com/surface_generation.html
        :return:
        """
        #   parent
        self._parent = parent

        #   pointer to body
        self.body = body

        #   filename
        self.filename = filename

        #   sub string
        self._substring = "roughness_profile_"+self.__id_list.next()+"."

        #   body properties
        self.body_id = body_id
        self.uPi = uPi
        self.uRi = uRi
        #   theta in degrees
        self.thetamax = thetamax
        self.theta = theta

        #   cistance between two points
        self.dx = dx
        self.dy = dx

        #   length of surface
        self.L = np.linalg.norm(self.uRi-self.uPi, ord=2)

        #   visualization properties
        self.color = np.array([0.2, 0.2, 0.2], dtype="float32")
        self.scale = scale

        #   roughness properties
        self.Ra = Ra
        self.Rz = Rz

        #   status of profile
        self.profile_corrected = False

        #   profile data
        #   nodex - list of nodes in 2D (x, y)
        self.nodes = np.empty([0, 0], dtype="float32")
        #   normals
        self.normals = np.empty([0, 0], dtype="float32")
        #   tangents
        self.tangents = np.empty([0, 0], dtype="float32")
        #   number of nodes
        self.n = len(self.nodes)

        #   profile type (shape)
        self.profile_type = "closed"

        self._profile_shape = None
        self._profile_shapes = ["circle"]
        
        self.properties_dict = extract_from_dictionary_by_string_in_key(_dict, self._substring)

        for key in self.properties_dict:
            setattr(self, key, self.properties_dict[key])

        #   load profile from file is file exists
        self.geometry = None
        if self.filename is not None:
            if os.path.isfile(self.filename):
                self.load_contact_profile(self.filename)
        
                #    contact profile geometry
                if self.profile_type == "opened":
                    self.geometry = ContactGeometry2DOpenedProfile(filename=self.filename, profile=self, body=self.body)

                if self.profile_type == "closed":
                    self.geometry = ContactGeometry2DClosedProfile(filename=self.filename, profile=self, body=self.body)
            
#         self.geometry = ContactGeometry2DProfile(filename=self.filename, profile=self, body=self.body)

                self._parent.body_list[self._parent.body_id_list.index(self.body_id)].geometry_list.append(self.geometry)

    def save_contact_profile(self, nodes=None, normals=None, tangents=None):
        """

        """
        if nodes is None:
            nodes = self.nodes
        if normals is None:
            normals = self.normals
        if tangents is None:
            tangents = self.tangents

        data = np.hstack((nodes, normals, tangents))
        np.savetxt(self.body._name+".prfl", data, "%6.4f", "\t")

    def load_contact_profile(self, filename):
        """
        Load geometry nodes from text file
        """
        # print "Contact_profile", filename, 'loaded for body', self.body._name, 'with id', self.body_id,'.'
        #    load data from file
        data = []
        arc_ = False
        with open(filename, 'r') as file_:
            for line in file_:
                if "#" in line:
                    line = line[:line.index("#")].strip()

                if len(line.strip()) == 0 or line.startswith("#"):
                    pass
                else:
                    #   if string has commas, then convert to array
                    if "," in line and not "arc" in line:
                        line = string2array(line)
                        if data == []:
                            data = np.append(data, line)
                        else:
                            data = np.vstack((data, line))

                    if "arc" in line:
                        arc_ = True
                        if line.startswith("arc R ="):
                            R = float(line[line.find("=")+1:].strip())
                        elif line.startswith("arc x0 ="):
                            x0 = float(line[line.find("=")+1:].strip())
                        elif line.startswith("arc y0 ="):
                            y0 = float(line[line.find("=")+1:].strip())
                        elif line.startswith("arc n ="):
                            n = float(line[line.find("=")+1:].strip())

        if arc_ is True:
            theta = np.arccos(x0 / R)
            arcpoints = np.linspace((np.pi + theta), (np.pi - theta), n)
            arcx = np.cos(arcpoints) * R + x0
            arcy = np.sin(arcpoints) * R + y0
            arc = np.vstack((arcx, arcy)).T

            data = np.concatenate((data, arc), axis=0)

        # data = np.array(np.loadtxt(filename, delimiter=','))
        # if self.body_id == 0:
        #     print data

        #    check size and if normals are included in file
        n, cols = np.shape(data)
        #   cols - number of columns
        #   options (data shape):
        #   2 cols (x, y)
        #   4 cols (x, y, nx, ny)
        #   6 cols (x, y, nx, nz, tx, ty)

        if cols == 2:
            self.nodes = data * self.scale
            self.tangents = self._evaluate_tangents(self.nodes)
            self.normals = self._evaluate_normals(self.tangents)

            logging.getLogger("DyS_logger").info("Profile geometry from file %s read successfully!"%filename)

        if cols == 4:
            self.nodes = data[:, 0:3]
            print "todo"

        if cols == 6:
            print "todo"

    def _evaluate_tangents(self, nodes):
        """

        :param nodes:
        :return:
        """
        self.N = len(nodes)
        tangents = np.zeros([self.N-1, 2])
        for i in xrange(0, self.N - 1):
            tangents[i, :] = (nodes[i+1, :] - nodes[i, :]) / (np.linalg.norm((nodes[i+1, :] - nodes[i, :]), ord=2))

        return tangents

    def _evaluate_normals(self, tangents):
        """

        :param tangents:
        :return:
        """
        normals = np.zeros(np.shape(tangents))
        for i in xrange(0, self.N - 1):
            normals[i, :] = t2n(tangents[i, :])

        return normals

    def create_profile(self, profile_shape=None, R0=None, n=10):
        """

        :return:
        """
        #   shape of profile
        self._profile_shape = profile_shape

        #   radius of circle
        self.R0 = R0

        #   number of nodes
        self.n = int(self.L / self.dx)

        #   x coordinates of profile
        self.x = np.linspace(-self.L / 2., +self.L / 2., self.n)

        #   initial random height
        self.h = self.Ra * np.sqrt(2.0)

        #   create random height of profile points
        self.y0 = self.h * np.random.randn(1, self.n)

        #   correlation of surface using convolution, inverse Fourier transform and normalizing prefactors
        self.y = self._generate_y(self.y0)

        self.xy = np.vstack((self.x, self.y)).T

        i = 0
        thetamax = np.deg2rad(self.thetamax)
        while not self.profile_corrected:
            i += 1
            FLAG = 0
            for j in xrange(1, self.n - 1):
                bc = self.xy[j, :] - self.xy[j-1, :]
                ba = self.xy[j+1, :] - self.xy[j, :]

                #   dot product
                dotbac = np.dot(bc, ba)

                #   theta angle
                theta = np.arccos(dotbac / (np.linalg.norm(bc, ord=2) * np.linalg.norm(ba, ord=2)))

                #   check criteria for angle
                if (np.pi - theta) < thetamax:
                    FLAG += 1
                    self.dy = self.dx * i

            if FLAG == 0:
                self.profile_corrected = True
                break
            else:
                self.xy[:, 1] = self.y = self._generate_y(self.dy)

        if self._profile_shape == "circle":
            self._create_circle()

    def _generate_y(self, dy):
        """

        :return:
        """
        #   reevaluate gaussian filter
        self.filter = np.exp(-self.x ** 2 / (dy ** 2 / 2.))

        y = np.sqrt(2. / np.sqrt(np.pi)) * np.sqrt(self.L / self.n / self.dx) * np.fft.ifft(np.fft.fft(self.y0) * np.fft.fft(self.filter)).T

        return y[:,0].real

    def _create_circle(self):
        """

        :return:
        """
        fi = np.linspace(0, 2*np.pi, self.n)
        self.nodes = np.zeros([self.n, 2])

        for i in xrange(0, self.n):
            self.nodes[i, :] = self.R0 * np.array([np.cos(fi[i]), np.sin(fi[i])])

    def write_contact_profile_to_file(self, filename=None):
        """

        :param filename:
        :return:
        """
        self.filename = filename

        self._name, self._filetype = os.path.splitext(self.filename)

        if self._filetype == ".txt":
            self._write_to_txt_file()

    def _write_to_txt_file(self):
        """

        :return:
        """
        np.savetxt(self.filename, self.nodes, delimiter=',')


if __name__ == "__main__":
    # filename = "SR_test_R.txt"
    # profile = RoughnessProfile(filename=filename)
    # pprint(vars(profile))
    profile = RoughnessProfile(0, uRi=np.array([4., 0]), Ra=0.8)
    n = 7
    profile.create_profile()


    # profile.write_contact_profile_to_file("circle_n="+str(n)+".txt")

    plt.grid()
    plt.plot(profile.x, profile.y, '-x')
    plt.show()
