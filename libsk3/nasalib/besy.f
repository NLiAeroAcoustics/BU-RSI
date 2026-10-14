      subroutine besy (X, n, BY, BYPRIM, SUBNAM, ier)

c***********************************************************************
c*
c* Routine:   besy (ver 1.00)
c*
c* Purpose:   Evaluate Y & Y'
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
c* Author:    Donald B. Hanson
c* Date:      circa 1990
c* Revisions:
c*
c*
c***********************************************************************
c...+....1....+....2....+....3....+....4....+....5....+....6....+....7..

      implicit real*8 (a-h, o-z)
      real*8 Y(0:1030)
      character*6 SUBNAM


      BY = 0.D0
       BYPRIM = 0.D0
        ier = 0
         if (n .gt. 1029) then
          ier = 4
           write(3,100) SUBNAM, n, X
 
 100  format(1X,//,
     . 20X, '* *  F  A  T  A  L      E  R  R  O  R  * *', ///,
     . 16X, 'Order, n of the Y Bessel function is out of ' //
     . 16X, 'range while being called from subroutine ', A6, //,
     . 16X, 'n =', I4, ' X =', D10.4, //,
     . 16X, 'Program will be terminated.', //)
      
       stop
  
       endif
        if (n .LT. 0) then
         ier = 1
          write(3,200) SUBNAM, n, X
 
 200  format(1X, //,
     . 20X, '* *  F  A  T  A  L      E  R  R  O  R  * *', ///,
     . 16X, 'Order, n of the Y Bessel function is negative', //,
     . 16X, 'while being called from subroutine ', A6, //,
     . 16X, 'n =', I4, ' X =', D10.4, //,
     . 16X, 'Program will be terminated.', //)

      stop

      endif
       if (X .le. 0.D0) then
        ier = 2
         write(3,300) SUBNAM, n, X
 
 300  format(1X, //,
     . 20X, '* *  F  A  T  A  L      E  R  R  O  R  * *',///,
     . 16X, 'Argument, X of the Y Bessel function is negative or ', //
     . 16X, 'zero while being called from subroutine ', A6, //,
     . 16X, 'n =', I4, ' X =', D10.4, //,
     . 16X, 'Program will be terminated.', //)
      
      stop

      endif

c     Branch if X less than or equal to 4

   20 if (X - 4.D0) 40, 40, 30

c     Compute Y0 and Y1 for X greater than 4

   30 T1 = 4.D0 / X
       T2 = T1 * T1
      
      P0 = ((((-0.0000037043D0  * T2 + 0.0000173565D0) * T2 -
     .          0.0000487613D0) * T2 + 0.00017343D0)   * T2 -
     .          0.001753062D0)  * T2 + 0.3989423D0
 
      Q0 = ((((0.0000032312D0  * T2 - 0.0000142078D0) * T2 +
     .         0.0000342468D0) * T2 - 0.0000869791D0) * T2 +
     .         0.0004564324D0) * T2 - 0.01246694D0
 
      P1 = ((((0.0000042414D0  * T2 - 0.0000200920D0) * T2 +
     .         0.0000580759D0) * T2 - 0.000223203D0)  * T2 +
     .         0.002921826D0)  * T2 + 0.3989423D0
 
      Q1 = ((((-0.0000036594D0  * T2 + 0.00001622D0)   * T2 -
     .          0.0000398708D0) * T2 + 0.0001064741D0) * T2 -
     .          0.0006390400D0) * T2 + 0.03740084D0
      
      A = 2.D0 / DSQRT(X)
       B = A * T1
        C = X - 0.7853982D0

      Y0 = A * P0 * DSIN(C) + B * Q0 * DCOS(C)
       Y1 = -A * P1 * DCOS(C) + B * Q1 * DSIN(C)

      go to 90

c       COMPUTE Y0 AND Y1 FOR X LESS THAN OR EQUAL TO 4

   40 XX = X / 2.D0
       X2 = XX * XX
        T = DLOG(XX) + 0.5772157D0

      SUM = 0.D0
       TERM = T
        Y0 = T

      do 70 L = 1, 15
       if (L - 1) 50, 60, 50

   50 SUM = SUM + 1.D0 / DFLOAT(L-1)

   60 FL = L
       TS = T - SUM
        TERM = (TERM * (-X2) / FL**2) * (1.D0 - 1.D0 / (FL*TS))

   70 Y0 = Y0 + TERM
       TERM = XX * (T - 0.5D0)
         SUM = 0.D0
          Y1 = TERM

      do 80 L = 2, 16
       SUM = SUM + 1.D0 / DFLOAT(L-1)
        FL = L
         FL1 = FL - 1.D0
          TS = T - SUM
      TERM = (TERM * (-X2) / (FL1*FL)) * ((TS - 0.5D0 / FL) /
     .        (TS + 0.5D0 / FL1))

   80 Y1 = Y1 + TERM
       PI2 = 0.6366198D00
        Y0 = PI2 * Y0
         Y1 = -PI2 / X + PI2 * Y1

c     CHECK IF Y0 IS DESIRED

   90 continue
      if (n .EQ. 0) then

c     Return Y0 and Y'0

      BY = Y0
       BYPRIM = (-1.D0 * Y1)

      go to 170
       endif

c   Perform recursion to find Yn(X)

  130 continue
      Y(0) = Y0
       Y(1) = Y1

      do 145 k = 1, n
       T = DFLOAT(2 * k) / X
        Y(k + 1) = T * Y(k) - Y(k - 1)
         if (DABS(Y(k + 1)) - 1.7D33) 145, 145, 141

  141     IER = 3
           BY = Y(k)
            BYPRIM = 0.5D0 * (Y(k - 1) - Y(k + 1))
            
      go to 170
  145  continue

  150 BY = Y(n)
       BYPRIM = 0.5D0 * (Y(n - 1) - Y(n + 1))

  170 return
       end
