        subroutine anfu (m, X, s, f, ieb)

c***********************************************************************
c*
c* Routine:   anfu (ver 1.00)
c*
c* Purpose:   Evaluate f = J' * Y - J * Y'
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
c*     m - Cicumferential mode order (input)
c*     x - Guess for Xmn (input)
c*     s - Hub-to-tip radius ratio sigma (input)
c*     f - Function value (output)
c*     eb - Array of Bessel funciton errors (output)
c*
c* Author:    Unknown
c* Date:      circa 1980
c* Revisions:
c*
c*
c***********************************************************************
c...+....1....+....2....+....3....+....4....+....5....+....6....+....7..

      implicit real*8 (a-h, o-z)
      integer ieb(4)

      
c  Bessel J 

      call besj (X, m, BJM, ZJ, 'anfu  ', ie1)
       ieb(1) = ie1
        call besj (s*X, m, BJMS, ZJS, 'anfu  ', ie2)
         ieb(2) = ie2

c  Bessel Y

      call besy (X, m, BYM, ZY, 'anfu  ', ie3)
       ieb(3) = ie3
        call besy (s*X, m, BYMS, ZYS, 'anfu  ', ie4)
         ieb(4) = ie4

      f   = ZJS * ZY - ZJ * ZYS

      return
       end
