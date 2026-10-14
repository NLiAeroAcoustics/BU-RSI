        subroutine eigen (mm, sigmaR, Xmn, Amn, Bmn, ier)

c***********************************************************************
c*
c* Routine:   eigen (ver 1.00)
c*
c* Purpose:   Evaluate eigenvalues
c*
c* Method:    See C.S. Ventres, et. al., "Turbofan Noise Generation"
c*
c* I/O:       Input File = None
c*            Output File = None 
c*
c* Called by:
c*
c*     Argument List:
c*
c* Author:    Unknown
c* Date:      circa 1990
c* Revisions:
c*
c*
c***********************************************************************
c...+....1....+....2....+....3....+....4....+....5....+....6....+....7..

      implicit real*8 (a-h, o-z)
      integer ier(4)
 
 
      if (Xmn .eq. 0.D0) go to 22
       n = mm
        if (n .lt. 0) n = -n
         XmnHUB = Xmn * sigmaR
          s = sigmaR**2

c     Bessel J (@ tip)
       call besj (Xmn, n, BJM, BJPRIM, 'EIGEN ', ier(1))


c     Bessel Y (@ tip)
       call besy (Xmn, n, BYM, BYPRIM, 'EIGEN ', ier(2))


c     Bessel J (@ hub)
       call besj (XmnHUB, n, BJMH, BJMHPR, 'EIGEN ', ier(3))


c     Bessel Y (@ hub)
       call besy (XmnHUB, n, BYMH, BYMHPR, 'EIGEN ', ier(4))

      if (ier(2) .eq. 3 .or. ier(4) .eq. 3) write(3,100) mm

 100  format(1X, //,
     . 25X,'* *  N  O  T  E  * *', ///,
     . 12X, 'Y Bessel function set to minus machine infinity ', //,
     . 12X, ' while calculating the radial mode shape weighting ', //,
     . ' factors Amn and Bmn for m = ', I4, //)

c  End of Bessel function calls

      if (DABS(BJPRIM) .lt. DABS(BYPRIM)) go to 10
       Amn = -BYPRIM / BJPRIM
        Bmn = 1.D0
         go to 20

10      Amn = 1.D0
         Bmn = -BJPRIM / BYPRIM

20    PART1 = (1.D0 -(mm / Xmn)**2) *
     .         (Amn * BJM + Bmn * BYM)**2 / (1.D0 - s)
     
      PART2 = (s - (mm / Xmn)**2) *
     .         (Amn * BJMH + Bmn * BYMH)**2 / (1.D0 - s)
     
      FACTOR = DSQRT(PART1 - PART2)
            
      Amn = Amn / FACTOR
       Bmn = Bmn / FACTOR
       
      return

22    continue
 
      Amn = 1.D0
       Bmn = 0.D0
 
      do 24 j = 1, 4
       ier(j) = 0
24    continue

      return
       end
