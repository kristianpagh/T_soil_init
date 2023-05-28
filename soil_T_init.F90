SUBROUTINE SOIL_T_INIT(TMEAN,TDELTA,DAY,PTG)

! ======================================================================
! Purpose:
! --------
! *SOIL_T_INIT* Computes an initial soil vertical temperature profile as 
! a function of average TMEAN, the Winter to Summer temperature span 
! TDELTA and the Julian day. This is done for the 14 ISBA DIFF layers 
! defined by Descharme et al. (2013); doi:10.1002/jgrd.50631
!
! This subroutine is valid for latitudes northward of the Tropic of 
! Cancer (approx. 23.5 deg. N). For latitudes southward of the Tropic of 
! Capricorn (approx. 23.5 deg. S) the cosine term in the main equation 
! should be positive. For tropical latitudes a different equation should 
! be used.
!
! Author:     Kristian Pagh Nielsen,                          2023-05-25
! -------

! ======================================================================
! 1. Declarations
! ======================================================================

IMPLICIT NONE

INTEGER, PARAMETER              :: JPIM = SELECTED_INT_KIND(9)
INTEGER, PARAMETER              :: JPRM = SELECTED_REAL_KIND(6,37)

INTEGER(KIND=JPIM), PARAMETER   :: NLMAX      = 20
INTEGER(KIND=JPIM), PARAMETER   :: NLAYER     = 14
! Termal diffusivity (could be an input variable in stead):
REAL(KIND=JPRM), PARAMETER      :: KAPPA      = 0.0432 ! m²/day
REAL(KIND=JPRM), PARAMETER      :: OMEGA      = 0.01720242383 ! 2*pi/365.25 1/day
REAL(KIND=JPRM), PARAMETER      :: DZ(NLAYER) = (/ 0.01, 0.04, 0.1, & 
 & 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0 /)
REAL(KIND=JPRM), PARAMETER      :: ZERO       = 0.0
REAL(KIND=JPRM), PARAMETER      :: ONE        = 1.0
REAL(KIND=JPRM), PARAMETER      :: TWO        = 2.0
REAL(KIND=JPRM), PARAMETER      :: THREE      = 3.0
REAL(KIND=JPRM), PARAMETER      :: FOUR       = 4.0
REAL(KIND=JPRM), PARAMETER      :: FIVE       = 5.0
REAL(KIND=JPRM), PARAMETER      :: SIX        = 6.0
REAL(KIND=JPRM), PARAMETER      :: SEVEN      = 7.0
REAL(KIND=JPRM), PARAMETER      :: TEN        = 10.0
REAL(KIND=JPRM), PARAMETER      :: HUNDRED    = 100.0

! ----------------------------------------------------------------------
! 1.1 Input variables
! ----------------------------------------------------------------------

REAL(KIND=JPRM), INTENT(IN)     :: TMEAN
REAL(KIND=JPRM), INTENT(IN)     :: TDELTA
INTEGER(KIND=JPIM), INTENT(IN)  :: DAY

! ----------------------------------------------------------------------
! 1.2 Output variables
! ----------------------------------------------------------------------

REAL(KIND=JPRM), INTENT(OUT)    :: PTG(NLAYER)

!-----------------------------------------------------------------------
! 1.3   LOCAL VARIABLES:
!-----------------------------------------------------------------------

INTEGER(KIND=JPIM)              :: JL
REAL(KIND=JPRM)                 :: FACTOR
REAL(KIND=JPRM)                 :: RDAY
REAL(KIND=JPRM)                 :: ZTG(NLMAX)
REAL(KIND=JPRM)                 :: ZZ

! ======================================================================
! 2. Produce vertical profiles for the 14 ISBA DIFF layers
! ======================================================================

PTG=ZERO
ZTG=ZERO
FACTOR = SQRT(OMEGA/(TWO*KAPPA))
RDAY = DAY*ONE
DO JL=0,10
   ZZ = JL/HUNDRED
   ZTG(JL+1) = TMEAN - TDELTA*EXP(-ZZ*FACTOR)* &
    &                  COS(OMEGA*RDAY-ZZ*FACTOR)
ENDDO
PTG(1) = (ZTG(1)+ZTG(2))/TWO
PTG(2) = SUM(ZTG(2:5))/FOUR
PTG(3) = SUM(ZTG(5:11))/SEVEN
ZTG=ZERO
DO JL=1,20
   ZZ = JL/TEN
   ZTG(JL) = TMEAN - TDELTA*EXP(-ZZ*FACTOR)* &
    &                COS(OMEGA*RDAY-ZZ*FACTOR)
ENDDO
PTG(4) = (ZTG(1)+ZTG(2))/TWO
PTG(5) = SUM(ZTG(2:4))/THREE
PTG(6) = SUM(ZTG(4:6))/THREE
PTG(7) = SUM(ZTG(6:8))/THREE
PTG(8) = SUM(ZTG(8:10))/THREE
PTG(9) = SUM(ZTG(10:15))/SIX
PTG(10) = SUM(ZTG(15:20))/SIX
ZTG=ZERO
DO JL=2,12
   ZZ = JL*ONE
   ZTG(JL-1) = TMEAN - TDELTA*EXP(-ZZ*FACTOR)* &
    &                  COS(OMEGA*RDAY-ZZ*FACTOR)
ENDDO
PTG(11) = (ZTG(1)+ZTG(2))/TWO
PTG(12) = SUM(ZTG(2:4))/THREE
PTG(13) = SUM(ZTG(4:7))/FOUR
PTG(14) = SUM(ZTG(7:11))/FIVE

!For test plots:
!DO JL=1,14
!   WRITE(*,*)DZ(JL),PTG(JL)
!ENDDO
!WRITE(*,*)

END SUBROUTINE SOIL_T_INIT
