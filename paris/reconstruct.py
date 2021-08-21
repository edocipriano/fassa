import numpy as np

def al3d(nr, cc):
    """
        compute in the unit cube the plane constant alpha satisfying
        [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
        where the normal is pointing outward from the reference phase
        INPUT: normal coefficients in nr(3) and volume fraction cc.
    """
    athird = 1./3.
    eps0 = 1.e-50

    np1 = abs(nr[0])    # need positive coefficients
    np2 = abs(nr[1])
    np3 = abs(nr[2])
    m1 = min(np1, np2)  # order positive coefficients
    m3 = max(np1, np2)
    if np3 < m1:
        m2 = m1
        m1 = np3
    elif np3 >= m3:
        m2 = m3
        m3 = np3
    else:
        m2 = np3

    # NOTE avoid runTimeWarning in 2D
    m3 += 1.e-12

    denom = max(6.**m2*m3, eps0)
    cch = min(cc, 1.-cc)    # limit to: 0 < cch < 1/2
    c01 = m1*m1*m1/denom    # get cch ranges
    c02  = c01 + 0.5*(m2-m1)/m3
    m12 = m1 + m2
    if m12 <= m3:
        c03 = 0.5*m12/m3
    else:
        numer = m3*m3*(3.0*m12-m3) + m1*m1*(m1-3.0*m3) + m2*m2*(m2-3.0*m3)
        c03 = numer/denom

    # 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12))
    if cch <= c01:
        al3d = (denom*cch)**athird                                  # case (1)
    elif cch <= c02:
        al3d = 0.5*(m1 + np.sqrt(m1*m1 + 8.*m2*m3*(cch-c01)))       # case (2)
    elif cch <= c03:
        p = 2.*m1*m2
        q = 1.5*m1*m2*(m12 - 2.*m3*cch)
        pst = np.sqrt(p)
        arc = athird*np.arccos(q/(p*pst))
        csarc = np.cos(arc)
        al3d = pst*(np.sqrt(3.*(1.-csarc*csarc)) - csarc) + m12     # case (3)
    elif m12 <= m3:
        al3d = m3*cch + 0.5*m12                                     # case (4a)
    else:
        p = m12*m3 + m1*m2 - 0.250
        q = 1.50*m1*m2*m3*(0.50-cch)
        pst = np.sqrt(p)
        arc = athird*np.arccos(q/(p*pst))
        csarc = np.cos(arc)
        al3d = pst*(np.sqrt(3.*(1.-csarc*csarc)) - csarc) + 0.5     # case (4b)

    if cc > 0.5:
        al3d = 1. - al3d

    # compute alpha for the given coefficients
    al3d = al3d + min(0., nr[0]) + min(0., nr[1]) + min(0., nr[2])

    return al3d


def fl3d(nr, alpha, x0, dx):
    """
        compute in the right hexahedron starting at (x01,x02,x03) and of sides
        (dx1,dx2,dx3) the volume cut by the plane
        [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
        INPUT: normal coefficients in nr(3), plane constant alpha, starting
               point x0(3), sides dx(3)
    """
    eps0 = 1.e-50

    # move origin to x0
    al = alpha - nr[0]*nr[0] - nr[1]*x0[1] - nr[2]*x0[2]
    # reflect the figure when negative coefficients
    al = al + max(0., -nr[0]*dx[0]) + max(0.,-nr[1]*dx[1]) + max(0., -nr[2]*dx[2])
    np1 = abs(nr[0]) # need positive coefficients
    np2 = abs(nr[1])
    np3 = abs(nr[2])
    almax = np1*dx[0] + np2*dx[1] + np3*dx[2]
    # NOTE fpe
    almax += 1.e-12
    al = max(0., min(1., al/almax)) # get new al within safe limits
    alh = min(al, 1.-al)            # limit to: 0 < alh < 1/2

    # normalized equation: m1*y1 + m2*y2 + m3*y3 = alh, 0 <= m1 <= m2 <= m3
    # the problem is then solved again in the unit cube
    np1 = np1/almax
    np2 = np2/almax
    np3 = np3/almax
    m1 = min(np1*dx[0], np2*dx[1])  # order coefficients
    m3 = max(np1*dx[0], np2*dx[1])
    top = np3*dx[2]
    if top < m1:
        m2 = m1
        m1 = top
    elif top >= m3:
        m2 = m3
        m3 = top
    else:
        m2 = top

    m12 = m1 + m2
    mm = min(m12, m3)
    denom = max(6.*m1*m2*m3, eps0)

    # 1: al<=m1; 2: m1<=al<=m2; 3: m2<=al<=mm; 4: mm<=al<=1/2 (a:m12<=m3; b:m3<m12))
    if alh <= m1:
        frac = alh*alh*alh/denom                              # case (1)
    elif alh <= m2:
        frac = 0.5*alh*(alh-m1)/(m2*m3) +  m1*m1*m1/denom     # case (2)
    elif alh <= mm:
        top = alh*alh*(3.*m12-alh) + m1*m1*(m1-3.*alh)        # case (3)
        frac = (top + m2*m2*(m2-3.*alh))/denom
    elif m12 <= m3:
        frac = (alh - 0.5*m12)/m3                             # case (4a)
    else:
        top = alh*alh*(3.-2.*alh) + m1*m1*(m1-3.*alh)
        frac = (top + m2*m2*(m2-3.*alh) + m3*m3*(m3-3.*alh))/denom # case (4b)

    top = dx[0]*dx[1]*dx[2]
    if al <= 0.5:
        fl3d = frac*top
    else:
        fl3d = (1.-frac)*top

    return fl3d

