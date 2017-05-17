# coding=utf-8
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt


from distance_PSCJ import DistancePSCJ
from contact_point_PSCJ import ContactPointPSCJ

if __name__ == "__main__":
    rPi = np.array([-0.000865303862, -0.000306777138])
    rRi = np.array([3.109033787401e-05, 7.829662785387e-06])
    rPj = np.array([1.999999949476e-05, 4.999999873689e-05])

    #   frame of slot in GCS
    frame = np.array([[-0.001130235883, 0.000448081123],
                    [-0.000233841668, 0.000762687929],
                    [ 0.000296022404, -0.000747028579],
                    [-0.00060037181, -0.001061635386]])

    R0_j = 0.75E-3
    h0_iP = 1.6E-3
    R0_i = h0_iP / 2.
    c = R0_i - R0_j
    L = 2E-3
    d = DistancePSCJ(rPi, rRi, rPj, R0_j=R0_j, h0_iP=h0_iP, c=c, L=L)
    d.contact_geometry_GCS()

    cp = ContactPointPSCJ(d)
    pprint(vars(cp))

    fig = plt.figure(num=1, figsize=(6, 5), dpi=72, facecolor='w', edgecolor='g')
    ax = plt.subplot(111, aspect="equal")

    #   pin - body j
    plt.plot(rPj[0], rPj[1], color='r', marker="o")
    circle_j = plt.Circle((rPj[0], rPj[1]), R0_j, color='r', fill=False)
    ax.add_artist(circle_j)

    #   slot - body i
    plt.plot([rPi[0], rRi[0]], [rPi[1], rRi[1]], color='green', marker="o")

    #   frame
    plt.plot(frame[:, 0], frame[:, 1])
    plt.show()