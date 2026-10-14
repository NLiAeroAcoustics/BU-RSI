c***********************************************************************
c*
c* Routine:   CASC (ver 2.00)
c*
c* Purpose:   This program calculates complex pressure distribution
c*            across the blade/vane chord for incident turbulence.
c*
c* Method:    See C.S. Ventres, et. al., "Turbofan Noinse Generation"
c*
c* I/O:       Input File = None
c*            Output File = None 
c*
c* Called By: STATOR
c*
c* Author:    C.S. Ventres, et. al., Bolt Beranek Newman
c* Date:      circa 1980
c* Revisions: Ed Envia, NYMA, Inc. --  08/25/96
c*            Julian Winkler, BU -- May 2011
c*            (added k3 via Graham's similarity rules)
c*            Nuo Li, BU Mar 2024
c*            (added sweep modeling)
c***********************************************************************
c...+....1....+....2....+....3....+....4....+....5....+....6....+....7..

      subroutine casc_k3_no_r
     &(MR_3D,k1_3D,k3_3D,sweep,sigma_3D,h,alpha,N,F,ier)

cjw   ---------------------------------------------------------
cjw   NOTES:
cjw   ------ 
cjw   (1) MR_3D = U
cjw   (2) k1_3D = kz
cjw   (3) sigma_3D > 0 (make sure that a positive value is passed)
cjw   (5) K_3D < 0, but sigma_3D > 0 --> recover M.L.'s results 
cjw   ----------------------------------------------------------

      implicit real*8 (a-h, o-z)
      
      integer iPVT(N)
      real*8 B(N),K_3D,k1,k1_3D,k2,k2_3D,k3,k3_3D
      real*8 MR_3D,beta,beta3D,comu,comu3D, sweep
      real*8 THETA,h,h1_3D,h2_3D,sigma_3D,sigma,zlocs(N)
      real*8 K, MR, pi, k1s_3D, k3s_3D, MR_S_3D, beta_s_3D
c      complex*16 K, MR
      complex*16 A(N,N), F(N), Zi(N), CKF, i
      complex*16 Cpcorr(N), G, G0, G1, G2, G3, lift


      i = (0.D0, 1.D0)
	pi = dacos(-1.d0)
      error = 1.D-04

cjw   Calculate missing parameters from input:	
      beta3D = DSQRT(1.D0 - MR_3D**2.0)
	  
 
cjw   Similarity parameter
cNL   include negative k3 and sweep effect	 
!	  if (k3_3D .lt. 0.D0) then
	  k1s_3D = k1_3D*DCOS(sweep)-k3_3D*DSIN(sweep)
	  k3s_3D = k3_3D
	  MR_S_3D = MR_3D*DCOS(sweep)
	  beta_s_3D = DSQRT(1.D0 - MR_s_3D**2.0)
!	  else
!	  k1_3D = k1_3D*DCOS(sweep)-k3_3D*DSIN(sweep)
!	  k3_3D = k1_3D*DSIN(sweep)+k3_3D*DCOS(sweep)
!	  end if
      THETA = k1s_3D / k3s_3D * MR_S_3D / beta_s_3D	  
	  
cjw   In the original casc K is computed using omega_b, which
cjw   is passed as -omega_b to casc --> make K_3D negative!
      K_3D = -k1s_3D * MR_S_3D / beta_s_3D**2
      h1_3D = h * DSIN(alpha)
      h2_3D = h * DCOS(alpha)
      comu3D = sigma_3D/(k1s_3D*h2_3D)-h1_3D/h2_3D
      k2_3D = k1s_3D * comu3D


c	write(*,*) 'in here,   MR_3D,  k1_3D, k3_3D, '
c      	write(*,*) 'sigma_3D,h,alpha'
c     	write(*,*) MR_3D,k1s_3D,k3_3D
c     &		sigma_3D,h,alpha

cjw   -------------------------
cjw   Graham's similarity rules
cjw   -------------------------


     
      MR = + MR_3D*DSQRT(1.D0-1.D0/THETA**2)


c	write(*,*)'THETA',THETA 
c        write(*,*) 'MR',DSQRT(1.D0-1.D0/THETA**2)
cjw   Equiv. 2-D compressible gust:
      beta = SQRT(1.D0 - MR**2)
      k1 = k1s_3D * (1.D0 + k3s_3D**2 / k1s_3D**2)
      k3 = 0.0
      h1 = h1_3D
      h2 = h2_3D * beta_s_3D / beta
      comu = beta_s_3D / beta * comu3D
      k2 = k1 * comu
      sigma = k1 * h1 + k2 * h2
cjw   if(sigma.gt.2*pi) sigma = sigma - INT(sigma/(2*pi))*2*pi
cjw   Make K negative again...
      K = - k1 * MR / beta**2
cjw   Note that K < 0 will affect gamma:
      gamma = sigma - K * MR * h1

c      write(*,*) 'K', K
c      write(*,*) 'gamma', gamma

cjw   Calculate the GX with neg. K (as before):
      G0 = i*beta*K/(2.0*pi*MR)
      G1 = beta*K**2/(2.0*pi)*(1.0/MR**2-0.5)
      G2 = i*beta*K**3/(4.0*pi)*(1.0/(2.0 *MR)-1./MR**3)
      G3 = -beta*K**4/(12.0*pi)*(1.0/MR**4-1.0/(2.0*MR**2)-0.125)

      do 5 l = 1, N
        B(l) = 1.D0
    5 continue
        B(N) = 0.5D0         

      do 10 m = 1, N
         z = DCOS((m - 0.5D0) * Pi / N)
	 zlocs(m) = z
cjw   write(40,*) z
         F(m) = - CDEXP(i * k1 * z)
   10 continue

      do 20 m = 1, N
        do 20 l = 1, N
	  S = -ALOG(2.0)
          do 40 ir=1,n
            bir=1.
            if(ir.eq.n) bir=0.5
      S = S-2.0*bir/ir*DCOS((m-0.5)*ir*pi/N)*DCOS(pi*ir*l/N)
40    continue
		
      z = DCOS((m-0.5D0)*Pi/N)-DCOS(1.D0*l*Pi/N)
	  
c ... Cascade kernel function calculation
cjw   The kernel has not been modified.

      call kernel (z, K, gamma, MR, h1, h2, error, CKF)


      G = (G0+G1*z+G2*z**2+G3*z**3)*CDEXP(i*K*MR*Z)
      A(m,l) = B(l)*(CKF + G*(S-LOG(abs(z))))*Pi/N
       
   20 continue

c ... Matrix inversion

       call zgeco (A, N, N, iPVT, RCOND, Zi)

       if (RCOND + 1.D0 .eq. 1.D0) then
c        write (3,100)

 100     FORMAT(1X,//,
     .20X,'* *  F  A  T  A  L      E  R  R  O  R  * *',///,
     .12X,'Inversion subroutine detected a singular matrix during',//,
     .12X,'the pressure distribution calculation.',//)

      stop
       endif


      if (RCOND .le. 1.D-05) then ier = 1
       job = 0
        call zgesl (A, N, N, iPVT, F, job)


	do 300 m = 1,N
cjw     This correction is =1 for 2-D gusts.
	Cpcorr(m) = beta/beta_s_3D*CDEXP(+i*(MR_S_3D*K_3D - MR*K)*zlocs(m))
        F(m)=Cpcorr(m)*F(m)  
 300	continue


	
	return
       end

