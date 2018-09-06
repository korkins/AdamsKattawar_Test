! RJD SPURR: see email !!!
!
PROGRAM TEST_SGLSCAT_TOA_PS
!
    IMPLICIT NONE
!
    INTEGER, PARAMETER :: &
        NL = 750, &
        NB = NL+1
    REAL*8, PARAMETER :: &
        RE = 6371.0D0
!
    INTEGER IL, ICASE
    REAL*8 :: MU0, SMU0, SZA, VZA, AZA, CAZ, TAU0, EXT, MUMU0, MM0, SMU, MU, X, R11, I1P, I1S, DH, HMAX, HMIN, &
              SZAP, MU0P, VZAP, MUP, CHF, U, SR1, D, I1P_BM, I1S_BM
!
    REAL*8, DIMENSION(NB) :: HKM, TAUE, CHF0, R, R2
    REAL*8, DIMENSION(NL) :: TAU
!
    ICASE = 1
!
    TAU0 = 0.25D0
    HMAX = 100.0D0
    HMIN = 0.0D0
!
    SELECT CASE (ICASE)
        CASE(1)
            SZA  = 84.26D0
            VZA  = 0.0D0
            AZA  = 0.0D0
            I1P_BM = 0.0161D0
            I1S_BM = 0.0193D0
        CASE(2)
            SZA  = 0.0D0
            VZA  = 70.0D0
            AZA  = 0.0D0
            I1P_BM = 0.0975D0
            I1S_BM = 0.102D0
        CASE(3)
            SZA  = 84.26D0
            VZA  = 70.0D0
            AZA  = 0.0D0
            I1P_BM = 0.0738D0
            I1S_BM = 0.0962D0
        CASE DEFAULT
            SZA  = 84.26D0
            VZA  = 70.0D0
            AZA  = 180.00D0
            I1P_BM = 0.0790D0
            I1S_BM = 0.0872D0
    END SELECT
!
    DH = HMAX/NL
    HKM(1) = HMAX
    DO IL = 2, NB
        HKM(IL) = HKM(IL-1) - DH
    END DO
    R = RE + HKM
    R2 = R*R
!
    EXT = TAU0/(HKM(1) - HKM(NB))
    DO IL = 1, NL; TAU(IL) = EXT*( HKM(IL) - HKM(IL+1) ); END DO
!
!   I1P:
    CAZ  = COSD( AZA )
    MU0  = COSD( SZA )
    SMU0 = SIND( SZA )
!
    MU = -COSD( VZA ) ! VZA' = PI - VZA
    SMU = SIND( VZA )
    MUMU0  = MU*MU0
    X = MUMU0 + SMU*SMU0*CAZ
    R11 = 0.75D0*(1.0D0 + X*X)
    MM0 = MU0/(MU0 - MU) ! > 0; MU < 0
    I1P = MM0*R11*( 1.0D0 - EXP( TAU(1)/MU/MM0 ) )  ! Initialize on TOA                 
    TAUE(1) = 0.0D0
!
    DO IL = 2, NL
        TAUE(IL) = TAUE(IL-1) + TAU(IL-1)
        I1P = I1P + MM0*R11*( 1.0D0 - EXP( TAU(IL)/MU/MM0 ) )*EXP( TAUE(IL)/MU )*EXP( -TAUE(IL)/MU0 )
    END DO
    TAUE(NB) = TAU0
!
    WRITE(*, *) 'SK - myself;  BM - benchmark Adams & Kattawar, Icarus(1978)'
    WRITE(*, *) 'P - plain-parallel;  S - pseudosphere'
    WRITE(*, *) 'Single scattering only'
    WRITE(*, *)
    WRITE(*, 10) 'I1P_SK = ', 0.25D0*I1P
    WRITE(*, 10) 'I1P_BM = ', I1P_BM
    WRITE(*, *)
!
    VZAP = VZA ! By definition of geometry
    MUP = -COSD( VZAP )
    SZAP = SZA ! By definition 
    MU0P = COSD( SZAP )
    MM0 = MU0P/(MU0P - MUP)
!
    I1S = MM0*R11*(1.0D0 - EXP( TAU(1)/MUP/MM0 )) ! Initialize on TOA
!
    SR1 = SIND( VZA )*R(1) ! sin(VZA)*(h_max+Re)
    DO IL = 2, NL
        VZAP = ASIND( SR1/R(IL) ) ! VZAP > VZA
        MUP = -COSD( VZAP )
        U = COSD( VZAP - VZA )
        MU0P = SQRT( 1.0D0 - U*U )*SMU0*CAZ + U*MU0
        MM0 = MU0P/(MU0P - MUP)
        ! D = R(IL)*SIND( VZAP - VZA )/SIND( VZA )        ! R(IL) = RE + HKM(IL); R*EXT = units of Tau
        D = SQRT( R2(IL) + R2(1) - 2.0D0*R(IL)*R(1)*U )   ! avoids 0 = sin(VZA=0)
        CHF = EXT*D
        CALL CHAPFUN0(RE, MU0P, HKM(1:IL), TAUE(1:IL), IL, IL-1, CHF0(1:IL))
        I1S = I1S + MM0*R11*(1.0D0 - EXP( TAU(IL)/MUP/MM0) )*DEXP( -(CHF + CHF0(IL)) )
    END DO
!
    WRITE(*, 20) 'NB = ', NB, 'DH(KM) = ', DH
    WRITE(*, 10) 'I1S_SK = ', 0.25D0*I1S
    WRITE(*, 10) 'I1S_BM = ', I1S_BM
!
!
10  FORMAT(1X, A, 3X, 1F8.4)
20  FORMAT(1X, A, 1I3, 3X, A, 1F8.4)
30  FORMAT(1X, A, 2X, 5F8.4)
    WRITE(*, *)
    WRITE(*, *) 'DONE! Press Enter'
    READ(*, *)
!
END PROGRAM