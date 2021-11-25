!===============================================================================
! DESCRIPTION OF FUNCTION AREA3D:
! compute in the unit cube the area cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT: normal coefficients in nr(3) and volume fraction cc 
! NOTE : the cut area A is invariant with respects to reflections, ordering  
!        of the coefficients and midpoint (i.e., A(cc) = A(1-cc))
!-------------------------------------------------------------------------------
FUNCTION AREA3D(nr,cc)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),cc
  REAL(8) :: AREA3D
  REAL(8) :: cch,ccr,al,c00,c01,c02,c03,np1,np2,np3
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

  denom = DMAX1(2.d0*m1*m2*m3,eps0)                           
  cch = DMIN1(cc,1.d0-cc)                              ! limit to: 0 < cch < 1/2
  ccr = DMAX1(cch,eps0)                                     ! avoid cch = m1 = 0  
  c01 = m1*m1*m1/(3.d0*denom)                                   ! get cch ranges
  c02  = c01 + 0.5d0*(m2-m1)/m3
  m12 = m1 + m2
  if (m12 <= m3) then
     c03 = 0.5d0*m12/m3
  else
     numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
     c03 = numer/(3.d0*denom)
  endif
  c00 = DSQRT(m1*m1 + m2*m2 + m3*m3)

! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
  if (ccr <= c01) then
     al = (3.d0*denom*cch)**athird                                         
     AREA3D = c00*al*al/denom                                         ! case (1)
  else if (ccr <= c02) then 
    al = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch - c01)))           
    AREA3D = c00*(2.d0*al-m1)/(2.d0*m2*m3)                            ! case (2)
  else if (ccr <= c03) then
     p = 2.d0*m1*m2                                                   
     q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12
     AREA3D = c00*(-al*al + m1*(2.d0*al-m1) + m2*(2.d0*al-m2))/denom  ! case (3)
  else if (m12 <= m3) then                                  
     al = m3*cch + 0.5d0*m12                                         
     AREA3D = c00/m3                                                 ! case (4a)
  else                                  
     p = m12*m3 + m1*m2 - 0.25d0                                     
     q = 1.5d0*m1*m2*m3*(0.5d0-cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0 
     AREA3D = c00*(2.d0*al*(1.d0-al) - c00*c00)/denom                ! case (4b)
  endif

END FUNCTION AREA3D
