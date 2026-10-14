      subroutine besj (X, n, BJ, BJPRIM, SUBNAM, ier)

c***********************************************************************
c*
c* Routine:   besj (ver 1.00)
c*
c* Purpose:   Evaluate J & J'
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
      real*8 J(0:1030)                         
      character*6 SUBNAM

      BJ = 0.D0
       BJPRIM = 0.0D0
        nplus = n + 1
         ier = 0

      if (X .eq. 0.D0) return 
       if (X .lt. 1.0D-04 .and. n .gt. 10) return 
        if (X .lt. 0.D0) then
         ier = 2
          write(3,100) SUBNAM, n, X
          
 100  format(1X, //,
     . 20X, '* *  F  A  T  A  L      E  R  R  O  R  * *', ///,
     . 16X, 'Argument, X of the J Bessel function is negative', //,
     . 16X, 'while being called from subroutine ', A6, //,
     . 16X, 'n =', I4, ' X =', D10.4, //,
     . 16X, 'Program will be terminated.', //)

       stop

      elseif (n .lt. 0) then
       ier = 1
        write(3,200) SUBNAM, n, X
        
 200  format(1X, //,
     . 20X, '* *  F  A  T  A  L      E  R  R  O  R  * *', ///,
     . 16X, 'Order, n of the J Bessel function is negative', //,
     . 16X, 'while being called from subroutine ', A6, //,
     . 16X, 'n =', I4, ' X =', D10.4, //,
     . 16X, 'Program will be terminated.', //)
      
       stop

      endif

c   Start backward recursion 30 orders above the J order desired

      mmax = MAX(INT(X), nplus) + 30

      if (MMAX .gt. 1028) then
       ier = 3
        write(3,300) SUBNAM, n, X

300   format(1X, //,
     . 20X, '* *  F  A  T  A  L      E  R  R  O  R  * *', ///,
     . 16X,'Order, n or argument, X of the J Bessel function is ', //,
     . 16X,'out of range while being called from subroutine ', A6, //,
     . 16X, 'n =', I4, ' X =', D10.4, //,
     . 16X, 'Program will be terminated.', //)
     
        stop
      endif

c   Backward recursion starts here
             
      J(mmax + 2) = 0.D0
       J(mmax + 1) = 1.D-30

        do 20 i = mmax, 0, -1
         J(i) = (2.D0 * (i+1) / X) * J(i+1) - J(i+2)

   20   continue
                                                                              
      SUM = J(2)
       do 30 i = 4, mmax, 2
        SUM = SUM + J(i)
   30  continue
                                                                    
      SUM = J(0) + 2.D0 * SUM
       BJ = J(N) / SUM 
                                                    
c   Backward recursion ends here     

c   If Bessel function is near zero and less than machine accuracy
c   then set Bessel function equal to zero
         
      if (DABS(BJ) .le. 1.D-35) BJ = 0.D0

c   Calculate Bessel function derivative. If derivative near zero and less 
c   than machine accuracy set derivative of Bessel function equal to zero 
                                                                     
      if (n .EQ. 0) then
       BJPRIM = (-1.D0 * J(1)) / SUM

      else
       BJPRIM = 0.5D0 * (J(N-1) - J(N+1)) / SUM

      if (DABS(BJPRIM) .le. 1.D-14) BJPRIM = 0.D0

      endif

      return
       end
