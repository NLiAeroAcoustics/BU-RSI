        subroutine rmode (mm, X, Amn, Bmn, PSI, ierror)

c***********************************************************************
c*
c* Routine:   rmode (ver 1.00)
c*
c* Purpose:   Evaluate mode shapes
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
c* Author:    unknown
c* Date:      circa 1990
c* Revisions:
c*
c*
c***********************************************************************
c...+....1....+....2....+....3....+....4....+....5....+....6....+....7..

      implicit real*8 (a-h, o-z)
      integer ierror(2)


      if (X .eq. 0.D0) go to 300
       FACTOR = 1.D0
        m = mm
          
      if (m .lt. 0) go to 100
        
200   continue

      call besj (X, m, BJ, DUM, 'rmode ', ierror(1))
       call besy (X, m, BY, DUM, 'rmode ', ierror(2))

      if (ierror(2) .eq. 3) write(3,400) m
        
 400  format(1X, //,
     . 25X, '* *  N  O  T  E  :  * *',///,
     . 12X, 'Y Bessek function set to minues machine infinity ', //,
     . 12X, 'while calculating the radial mode shape, PSI, for m = ',
     . I4, //)
 
      PSI = FACTOR * (Amn * BJ + Bmn * BY)
       return

100   m = -m
       k = m / 2
        if (DFLOAT(m) / 2.D0 .gt. DFLOAT(k)) FACTOR = -1.D0
        go to 200

300   continue
       PSI = Amn
        ierror(1) = 0
         ierror(2) = 0

      return
       end
