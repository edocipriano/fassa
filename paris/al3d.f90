!===============================================================================
! DESCRIPTION OF FUNCTION AL3D:
! compute in the unit cube the plane constant alpha satisfying
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! where the normal is pointing outward from the reference phase
! INPUT: normal coefficients in nr(3) and volume fraction cc 
!-------------------------------------------------------------------------------
FUNCTION AL3D(nr,cc)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),cc
  REAL(8) :: AL3D
  REAL(8) :: cch,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc
  REAL(8), PARAMETER :: athird=1.d0/3.d0,eps0=1.d-50 
  INTRINSIC DABS,DMIN1,DMAX1,DSQRT,DACOS,DCOS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
  m1 = DMIN1(np1,np2)                              ! order positive coefficients
  m3 = DMAX1(np1,np2)
  if (np3 < m1) then
     m2 = m1
     m1 = np3
  else if (np3 >= m3) then
     m2 = m3
     m3 = np3
   else
     m2 = np3
  endif

  denom = DMAX1(6.d0*m1*m2*m3,eps0)                           
  cch = DMIN1(cc,1.d0-cc)                              ! limit to: 0 < cch < 1/2
  c01 = m1*m1*m1/denom                                          ! get cch ranges
  c02  = c01 + 0.5d0*(m2-m1)/m3
  m12 = m1 + m2
  if (m12 <= m3) then
     c03 = 0.5d0*m12/m3
  else
     numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
     c03 = numer/denom
  endif

! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
  if (cch <= c01) then
     AL3D = (denom*cch)**athird                                       ! case (1)
  else if (cch <= c02) then 
     AL3D = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch-c01)))          ! case (2)
  else if (cch <= c03) then
     p = 2.d0*m1*m2
     q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     AL3D = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12        ! case (3)
  else if (m12 <= m3) then
     AL3D = m3*cch + 0.5d0*m12                                       ! case (4a)
  else                                  
     p = m12*m3 + m1*m2 - 0.25d0                                     
     q = 1.5d0*m1*m2*m3*(0.5d0-cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     AL3D = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0     ! case (4b)
  endif

  if (cc > 0.5d0)  AL3D = 1.d0 - AL3D
  
! compute alpha for the given coefficients  
  AL3D = AL3D + DMIN1(0.d0,nr(1)) + DMIN1(0.d0,nr(2)) + DMIN1(0.d0,nr(3))
  
END FUNCTION AL3D
