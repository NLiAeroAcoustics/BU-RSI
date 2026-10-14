        subroutine kernel (z, K, gamma, MR, h1, h2, error, CKF)

c***********************************************************************
c*
c* Routine:   kernel (ver 1.00)
c*
c* Purpose:   Evaluate cascade response kernel function
c*
c* Method:    See C.S. Ventres, et. al., "Turbofan Noise Generation"
c*
c* I/O:       Input File = None
c*            Output File = None 
c*
c* Called by:
c*
c* Author:    Unknown
c* Date:      circa 1980
c* Revisions: Ed Envia, NASA GRC, 1995, 2007
c*
c***********************************************************************
c...+....1....+....2....+....3....+....4....+....5....+....6....+....7..

        implicit real*8 (a-h, o-z)
        real*8 MR, K
        complex*16 i, CKF, an, anp, anm, Rn, root

        Rn(an) = (beta**2 * h2 / (2.D0 * d**2)) * (an**2 - K**2) * 
     .             CDEXP(i * (an + K * MR) * z) /
     .              ((an - Gn * h1 / d**2) * (an + K / MR))

        data i/(0.D0, 1.D0)/, zmin/1.D-07/

        Pi = 4.D0 * DATAN(1.D0)
         beta = DSQRT(1.D0 - MR**2)

        if(DABS(z) .gt. zmin) go to 10

         CKF = beta * (-1.D0 / z + i * (K / MR) * DLOG(DABS(z))) *
     .          CDEXP(i * K * MR * z) / (2.D0 * Pi)
     
        return

 10     continue

        d = DSQRT(h1**2 + (beta * h2)**2)
         R = d * DSQRT(K**2 + (d * DLOG(error) / (beta * h2 * z))**2)
        
        N1 = (gamma - R) / (2.D0 * Pi)
         N2 = (gamma + R) / (2.D0 * Pi)

        CKF = (0.D0, 0.D0)

        if(z .ge. 0.D0) go to 50

        do 40 n = N1, N2
         Gn = gamma - 2.D0 * Pi * n
          S = (K * d)**2 - Gn**2

         if(S .ge. 0.D0) root = -(1.D0, 0.D0) * DSQRT(S)
          if(S .lt. 0.D0) root = -i * DSQRT(-S)

            anm = (Gn * h1 + beta * h2 * root) / d**2
             CKF = CKF + Rn(anm)

 40     continue

         return

 50     continue

        do 70 n = N1, N2
         Gn = gamma - 2.D0 * Pi * n
          S = (K * d)**2 - Gn**2

         if (S .ge. 0.D0) root = (1.D0, 0.D0) * DSQRT(S)
          if (S .lt. 0.D0) root = i * DSQRT(-S)

           anp = (Gn * h1 + beta * h2 * root) / d**2
            CKF = CKF - Rn(anp)

 70     continue
       
        CKF = CKF - beta**2 * K / (2.D0 * MR) *
     .         DSINH(beta**2 * K * h2 / MR) * 
     .          CDEXP(-i * beta**2 * K * z / MR) / 
     .          (DCOSH(beta**2 * K * h2 / MR) -
     .            DCOS(gamma + K * h1 / MR))

        return
         end
