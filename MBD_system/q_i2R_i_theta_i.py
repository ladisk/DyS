# coding=utf-8


def q_i2R_i_theta_i(q_i):
    """
    Function gets R and theta of a body from vector q (vector of generalized coordinates)
    :param q_i:
    :return:
    """
    R_i = q_i[0:2]
    theta_i = q_i[2]
    return R_i, theta_i


