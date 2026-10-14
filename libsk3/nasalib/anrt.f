      subroutine anrt (m, n, s, ec, Xmn, ieb, iec)

c***********************************************************************
c*
c* Routine:   anrt (ver 1.00)
c*
c* Purpose:
c*
c* Method:    See C.S. Ventres, et. al., "Turbofan Noise Generation"
c*
c* I/O:       Input File = None
c*            Output File = None 
c*
c* Called By: ANRT
c*
c*     Argument List:
c*
c*     m - Cicumferential mode order (input)
c*     n - Radial mode order + 1 (input)
c*     ec - Convergance criteria for xmn
c*     xmn - Non-dimensional root (output)
c*     F - Function value (output)
c*     ieb - Sum of Bessel funciton errors (output)
c*     iec - error flag: iec = 0 (no error); iec = 1 - xmn (no convergance)
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

      Pi = 4.D0 * DATAN(1.D0)
      
      if (m .eq. 0 .and. n .eq. 1) go to 7

      Xmn = m - Pi

       do 6 nn = 1, n
        dX = Pi / 10.D0
         is = (-1)**(nn+1)
          X  = Xmn + Pi

        if (m .eq. 0 .and. nn .eq. 1) X = 0.D0
         j = 1
 
  	  if (X .le. 0.D0) then
	   F = 0.D0
	    go to 8
	  endif
	    
        call anfu (m, X, s, f, ieb)

  8      dX = -isig(F*is) * dX
          TE = F

  1      if (DABS(F) .le. ec) go to 3
          if (j .gt. 1000) go to 2

          X  = X + dX

          call anfu (m, X, s, F, ieb)

           if (F*TE .lt. 0.D0) dX = -dX / 2.D0
            j  = j + 1
             TE = F
              go to 1

  2        iec = 1
            go to 5
  3          iec = 0
  5           Xmn = X

  6        continue
  66      continue

          mu = n - 1
          
      if (ieb(3) .eq. 3 .or. ieb(4) .eq. 3) write (3,100) m, mu
        
 100  format(1X,//,
     . 25X, '* *  N  O  T  E  * *', ///,
     . 12X, 'Y Bessel function set to minus machine infinity while', //,
     . 12X, 'calculating xmn for (m, mu) =', 12X, I4, ',', 1X, I4, //)

      if (iec .eq. 1) write (3,200) m, mu
        
200   format(1X, //,
     . 25X, '* *  W  A  R  N  I  N  G  * *',///,
     . 10X, 'Calculation of xmn did not converge', /, 10X,
     . 'for (m, mu) = ', 1X, I4, ',', 1X, I4, //)
     
       return
       
 7    continue
 
c       write(11, '(1x, 1H(, I2, 1H, ,I2, 1H), E12.4, 2H *)' ) m, nn, Xmn

        do 20 i = 1, 4
         ieb(i) = 0
   20   continue
 
         iec = 0
          Xmn = 0.D0

      return
       end
