SUBROUTINE CHAPFUN0(RE, MU0, ZKM, TAU, NB, NL, CHF0)
!===============================================================================
! PURPOSE:
!   To compute the Chapman function for direct solar beam making an angle
!   acos(mu0) with surface normal at every given height (nadir observation).
!
! INPUT:
!   RE    D(1)    Radius of the Earth, Re ~= 6356.8 km
!   MU0   D(1)    cos(SZA), same at every ZKM; mu0 > 0.0
!   ZKM   D(NB)   Array of heights, TOA to BOA: ZKM(I) > ZKM(I+1)
!   TAU   D(NB)   Optical depth profile: DTAU(I) = TAU(I+1)-TAU(I) > 0.0
!   NB    I(1)    Number of boundaries
!   NL    I(1)    Number of layers, NB-1
!
! OUTPUT:
!   CHF0   D(NB)   Values of the Chapman function at each boundary, CHF0(1)=0.0
!
! TREE:
!   -
!
! COMMENTS:
!   See Fig.B1 in [1] for geometry.
!   When RE->Inf => X2->X1 => X2-X1 becomes INACCURATE (see X2 & X1 below).
!
! REFERENCES:
!   1. Dahlback A and Stamnes K, 1991: Planet. Space Sci, 39(5), p.671
!   2. Adams CN and Kattawar GW, 1978: Icarus, 35, p.139
!===============================================================================
!
! DECLARATION OF VARIABLES
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NB, NL
    REAL*8, INTENT(IN) :: RE, MU0
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NB), INTENT(IN) :: ZKM, TAU
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NB), INTENT(OUT) :: CHF0
!
! LOCAL SCALARS
    INTEGER &
        IB, IL
    REAL*8 &
        SIN2, RS, X1, X2
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NL) :: &
        EXT
    REAL*8, DIMENSION(NB) :: &
        R2, R
!===============================================================================
!
    SIN2 = 1.0D0 - MU0*MU0
    R = RE + ZKM
    R2 = R*R
    DO IL = 1, NL; EXT(IL) = (TAU(IL+1) - TAU(IL))/(ZKM(IL) - ZKM(IL+1)); END DO
!
!   Accumulate the Chapman function at each boundary
    CHF0(1) = 0.0D0
    DO IB = 2, NB
        CHF0(IB) = 0.0D0
        RS = R2(IB)*SIN2
        X2 = DSQRT(R2(1) - RS)
        DO IL = 1, IB-1
            X1 = DSQRT(R2(IL+1) - RS)
            CHF0(IB) = CHF0(IB) + EXT(IL)*(X2 - X1) ! X2 > X1
            X2 = X1
        END DO ! IL = 1, IB-1
    END DO ! IB = 1, NB
!
END SUBROUTINE CHAPFUN0
!===============================================================================
! 27May17 - SS for SZA > 0 & VZA = 0 [2: Table I]
!
! 15Jan16 - RE = 6356.8, ZKM = [15 10 6 3 1 0], TKM = [0.0 0.1 0.3 0.6 1.0 1.5]
!           SZA = 0 (zer) => CHF0 == TKM - ok for both RE & RE*100 & RE*0.1
!           SZA = 1 (one) => CHF0 ~= TKM (very close but slightly higher) - ok
!           SZA = 60,RE = 6356.8  => CHF0 > TKM - ok
!           SZA = 60, RE = RE*1000 => CHF0 ~= 2.0*TKM
!
!           Note, RE->Inf => X2->X1 => X2-X1 becomes INACCURATE.
!
!           NL=30: Rayleigh profile from IPRT. Reasonable values of CHF0 for
!           SZA = 0.0, 0.1, 1.0, 10, 80.
!===============================================================================