# coding=utf-8
import numpy as np
from mesh.cross_section.cross_section import Izz_rectangle

def deflection_y(F, L, E, I):
    y = (F * L**3) / (3. * E * I)
    return y


if __name__ == "__main__":
    E = 2.07E+11 #N/m2

    w = 0.1    #m
    h = 0.5    #m
    L = 2.      #m

    F_small = (5E+5) * (h**3)  # N
    F_large = (5E+8) * (h**3)  # N

    Izz = Izz_rectangle(w, h)
    print "Izz =", Izz
    print "F_small =", F_small
    print "y_small deformation =", deflection_y(F_small, L, E, Izz)
    print "y_large deformation =", deflection_y(F_large, L, E, Izz)