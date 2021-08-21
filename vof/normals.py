import numpy as np

def mycs(c):
    """
        Compute normals using the Mixed Young's centeres method.
        @param c: stencil3x3 of the color function
    """

    NOT_ZERO = 1.e-30
    m1, m2, t0, t1, t2 = 0.0, 0.0, 0.0, 0.0, 0.0

    m = np.zeros([4,3])
    cn = 0

    m1 = c[0,1,0] + c[0,1,2] + c[0,0,1] + c[0,2,1] + c[0,1,1]
    m2 = c[2,1,0] + c[2,1,2] + c[2,0,1] + c[2,2,1] + c[2,1,1]

    if (m1 > m2):
        m[0,0] = 1
    else:
        m[0,0] = -1

    m1 = c[0,0,1] + c[2,0,1] + c[1,0,1]
    m2 = c[0,2,1] + c[2,2,1] + c[1,2,1]
    m[0,1] = 0.5*(m1-m2)

    m1 = c[0,1,0] + c[2,1,0] + c[1,1,0]
    m2 = c[0,1,2] + c[2,1,2] + c[1,1,2]
    m[0,2] = 0.5*(m1-m2);

    m1 = c[0,0,1] + c[0,2,1] + c[0,1,1]
    m2 = c[2,0,1] + c[2,2,1] + c[2,1,1]
    m[1,0] = 0.5*(m1-m2)

    m1 = c[1,0,0] + c[1,0,2] + c[2,0,1] + c[0,0,1] + c[1,0,1]
    m2 = c[1,2,0] + c[1,2,2] + c[2,2,1] + c[0,2,1] + c[1,2,1]

    if (m1 > m2):
        m[1,1] = 1
    else:
        m[1,1] = -1

    m1 = c[1,0,0] + c[1,1,0] + c[1,2,0]
    m2 = c[1,0,2] + c[1,1,2] + c[1,2,2]
    m[1,2] = 0.5*(m1-m2)

    m1 = c[0,1,0] + c[0,1,2] + c[0,1,1]
    m2 = c[2,1,0] + c[2,1,2] + c[2,1,1]
    m[2,0] = 0.5*(m1-m2)

    m1 = c[1,0,0] + c[1,0,2] + c[1,0,1]
    m2 = c[1,2,0] + c[1,2,2] + c[1,2,1]
    m[2,1] = 0.5*(m1-m2)

    m1 = c[0,1,0] + c[2,1,0] + c[1,0,0] + c[1,2,0] + c[1,1,0]
    m2 = c[0,1,2] + c[2,1,2] + c[1,0,2] + c[1,2,2] + c[1,1,2]

    if (m1 > m2):
        m[2,2] = 1
    else:
        m[2,2] = -1

    t0 = np.abs(m[0,0]) + np.abs(m[0,1]) + np.abs(m[0,2])
    m[0,0] = m[0,0]/t0
    m[0,1] = m[0,1]/t0
    m[0,2] = m[0,2]/t0

    t0 = np.abs(m[1,0]) + np.abs(m[1,1]) + np.abs(m[1,2])
    m[1,0] = m[1,0]/t0
    m[1,1] = m[1,1]/t0
    m[1,2] = m[1,2]/t0

    t0 = np.abs(m[2,0]) + np.abs(m[2,1]) + np.abs(m[2,2])
    m[2,0] = m[2,0]/t0
    m[2,1] = m[2,1]/t0
    m[2,2] = m[2,2]/t0

    t0 = np.abs(m[0,0])
    t1 = np.abs(m[1,1])
    t2 = np.abs(m[2,2])

    cn = 0

    if (t1 > t0):
        t0 = t1
        cn = 1

    if (t2 > t0):
        cn = 2

    #Eigen::VectorXd mvect = m.row(3);
    mvect = m[3,:] # check

    fd32(c, mvect)

    for i in range(3):
        m[3,i] = mvect[i]

    #t0 = std::abs(m(3,0)) + std::abs(m(3,1)) + std::abs(m(3,2)) + NOT_ZERO;
    t0 = np.abs(m[3,0]) + np.abs(m[3,1]) + np.abs(m[3,2]) + NOT_ZERO
    m[3,0] = m[3,0]/t0;
    m[3,1] = m[3,1]/t0;
    m[3,2] = m[3,2]/t0;

    t0 = np.abs(m[3,0]);
    t1 = np.abs(m[3,1]);
    t2 = np.abs(m[3,2]);
    if (t1 > t0):
        t0 = t1
    if (t2 > t0):
        t0 = t2
    if (np.abs(m[cn,cn]) > t0):
        cn = 3

    mxyz = np.zeros(3)

    mxyz[0] = m[cn,0]
    mxyz[1] = m[cn,1]
    mxyz[2] = m[cn,2]

    """
    mxyz_dummy = mxyz
    mxyz[0] = mxyz_dummy[1]
    mxyz[1] = -mxyz_dummy[0]
    mxyz[2] = mxyz_dummy[2]
    """

    return mxyz



def fd32(c,mm):

    m1, m2 = 0.0, 0.0

    m1 = c[0,0,0] + c[0,2,0] + c[0,0,2] + c[0,2,2] + 2.0*(c[0,0,1] + c[0,2,1] + c[0,1,0] + c[0,1,2]) + 4.0*c[0,1,1]
    m2 = c[2,0,0] + c[2,2,0] + c[2,0,2] + c[2,2,2] + 2.0*(c[2,0,1] + c[2,2,1] + c[2,1,0] + c[2,1,2]) + 4.0*c[2,1,1]
    mm[0] = m1-m2

    m1 = c[0,0,0] + c[0,0,2] + c[2,0,0] + c[2,0,2] + 2.0*(c[0,0,1] + c[2,0,1] + c[1,0,0] + c[1,0,2]) + 4.0*c[1,0,1]
    m2 = c[0,2,0] + c[0,2,2] + c[2,2,0] + c[2,2,2] + 2.0*(c[0,2,1] + c[2,2,1] + c[1,2,0] + c[1,2,2]) + 4.0*c[1,2,1]
    mm[1] = m1-m2

    m1 = c[0,0,0] + c[0,2,0] + c[2,0,0] + c[2,2,0] + 2.0*(c[0,1,0] + c[2,1,0] + c[1,0,0] + c[1,2,0]) + 4.0*c[1,1,0]
    m2 = c[0,0,2] + c[0,2,2] + c[2,0,2] + c[2,2,2] + 2.0*(c[0,1,2] + c[2,1,2] + c[1,0,2] + c[1,2,2]) + 4.0*c[1,1,2]
    mm[2] = m1-m2

