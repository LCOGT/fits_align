import numpy as np

def cdist_np(XA, XB):


    XA = np.asarray(XA, order='c')
    XB = np.asarray(XB, order='c')

    s = XA.shape
    sB = XB.shape

    if len(s) != 2:
        raise ValueError('XA must be a 2-dimensional array.')
    if len(sB) != 2:
        raise ValueError('XB must be a 2-dimensional array.')
    if s[1] != sB[1]:
        raise ValueError('XA and XB must have the same number of columns '
                         '(i.e. feature dimension.)')

    mA = s[0]
    mB = sB[0]
    n = s[1]
    dm = np.empty((mA, mB), dtype=np.double)

    dm = cdist_euclid(XA, XB)
    return dm


def cdist_euclid(XA, XB):
    distances = []
    for xa in XA:
        dists = np.linalg.norm(xa - XB, axis=1)
        distances.append(dists)
    dist = np.array(distances)
    return dist
