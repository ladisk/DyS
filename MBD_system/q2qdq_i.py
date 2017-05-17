"""

created by: lskrinjar
date of creation: 23/04/2016
time of creation: 12:26
"""

import numpy as np


from q2q_body import q2q_body
from q2dq_body import q2dq_body


def q2qdq_i(q, id):
    """

    :param q:
    :param id:
    :return:
    """
    q_i = q2q_body(q, id)
    dq_i = q2dq_body(q, id)

    return q_i, dq_i


if __name__ == "__main__":
    q = np.random.rand(12)
    print q
    id = 0
    print q2qdq_i(q, id)
