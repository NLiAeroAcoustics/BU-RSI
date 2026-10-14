!!
!!  Calculates the spectral integrals given the surface unsteady pressure at all  frequency-wave number 
!!	combos of interest
!!  Using Liepmann weighting currently
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine get_Pmns(iFRQ,m2,n2,NC,NRcalc,NV,NF,Ng,dk3,datainterp,p_casc,Xmn2, Amn2, Bmn2,propmat,  &
		M_ann,specfreq,RDUCT,sos,hub_to_tip,tipc_to_tipr,sweep, iSTDYL,power_nondim)

!!!!
! this computes the contribution of cut on modes (m < 400, n < 50) to the power at the given frequency : specfreq
!!!!

	implicit none
	integer, parameter :: rk = selected_real_kind(p=12,r=307)

! variables common to rsi_aero

        integer  NV, NB, NRdat, NRcalc, NC, NF, iR, jR , Ng, jk, ik3,  &
                m2, ic, ierp, iFRQ, iSTDYL

        real(rk)   T, rho, Mr_tip, M_ann, tipc_to_tipr, hub_to_tip, &
                sweepdeg, RDUCT, BW, conv_tol_xmn, specfreq, sos
        real(rk)  datainterp(NRcalc,11)
        real(rk) R(NRcalc), c_ctip(NRcalc), stag(NRcalc), Maxial(NRcalc), Mcirc(NRcalc) , &
                c_ctip_over2(NRcalc),  Mach_ang(NRcalc) , Mach(NRcalc), &
                beta(NRcalc), sweep(NRcalc)
        complex(rk)  p_casc(NC,NRcalc,NV,Ng+1)
        real(rk) Ts, T_corr, pi,   k1, k3, sigma

! variables specific to rsi_spec
	integer  n2, mg
        real(rk) burms(NRcalc), murms(NRcalc), wurms(NRcalc), &
                        Lx(NRcalc), Lt(NRcalc), Lr(NRcalc)
        real(rk) Xmn2(m2,n2), Amn2(m2,n2), Bmn2(m2,n2), propmat(m2,n2)
	complex(rk)  Pmns(2)

! variable specific to get_Pmns
	integer  iRs(NRcalc), m, n, ii, q1max, q1, jk1, jk2, ierPsi(2), mode_loc, is(NC-1), muse, mode2
	real(rk)  kmns(m2,n2), gmnsup(m2,n2), gmnsdown(m2,n2), dk3, wr(NRcalc), dr, arg
	real(rk)  om_R_oversos, mthmode(m2,n2), power_nondim(m2,n2,2), muup, mudown
	real(rk)  effk1, effk2, effk3, h, h_chalf, h1, h2, innerblade, newinnerblade, allk
	real(rk)  phiweight,  mu, zn(NC-1), zn1(NC-1), g0(NC-1), g1(NC-1), Psi, sinc
	real(rk) mthmodeup(m2,n2), mthmodedown(m2,n2), lweight
	real(rk), allocatable :: smnup(:,:), smndown(:,:), smntemp(:,:)
	complex(rk)  cc,  g2(NC-1), cx(NC-1), cmns, capr, z2loop(NC), k3int(NC,NC)
	complex(rk),  allocatable :: rintup(:,:,:), rintdown(:,:,:)


	pi = dacos(-1.d0)
	cc = dcmplx(0.d0,1.d0)
!!!!!
! put input data at all radial calculation positions into arrays with explanatory names
!!!!!
	R = datainterp(:,1)

        c_ctip = datainterp(:,2)
        c_ctip_over2 = c_ctip/2.d0
        stag = datainterp(:,3)
        Maxial = datainterp(:,4)
        Mcirc = datainterp(:,5)


        Mach_ang = datan(Mcirc/ Maxial)   ! absolute flow angle
        Mach = dsqrt(Maxial**2 + Mcirc**2) * &  !absolute Mach number
                !dcos(sweep) * &                 ! sweep correction
                dcos(iSTDYL*(stag - Mach_ang))  ! turning correction


        burms = datainterp(:,6)
        murms = datainterp(:,7)
        wurms = datainterp(:,8)
        Lx = datainterp(:,9)
        Lt = datainterp(:,10)
        Lr = datainterp(:,11)

	om_R_oversos = 2*pi*specfreq*RDUCT/12/sos    !omega R_tip/ sos
 
	kmns = dsqrt(om_R_oversos**2- ( 1.d0 - M_ann**2)*Xmn2**2)
	gmnsup = (M_ann*om_R_oversos - kmns)/(1.d0 - M_ann**2)
	gmnsdown = (M_ann*om_R_oversos + kmns)/(1.d0 - M_ann**2)


	iRs = (/ (ii, ii = 1,NRcalc) /)
	wr = 1.d0 + (-1.d0)**iRs/3.d0  ! Simpson rule weights for radial integration
	wr(1) = wr(1)/2.d0
	wr(NRcalc) = wr(NRcalc)/2.d0

	
	mthmodeup = 0.d0
	mthmodedown = 0.d0
		
!	q1max = int(4*pi*R(NRcalc)/NV) + 1
	q1max = int(4*pi*R(NRcalc)/NV) + 1  + 2
	allocate (rintup(2*q1max+1,n2,NRcalc),rintdown(2*q1max+1,n2,NRcalc), smnup(2*q1max+1,n2), smndown(2*q1max+1,n2))
	allocate(smntemp(n2,NRcalc))

!!!!!!!
! start mode loop  -- parallelize this part
!!!!!!
!  !$omp parallel  do  &
!  !$omp default(shared)    &
!  !$omp firstprivate( m,  muse, q1max, rintup, rintdown, smnup, smndown)   &
!  !$omp firstprivate( iR, effk1, effk2, q1, h, h_chalf, h1, h2, innerblade)  &
!  !$omp firstprivate( newinnerblade, mode_loc, k3int, jk2, jk1, ik3, effk3, allk)  &
!  !$omp firstprivate( phiweight, arg, Psi, ierPsi, muup, mudown, is, zn, zn1, z2loop)  &
!  !$omp firstprivate( g0, g1, g2, cx, cmns, capr )


!$omp parallel  do  &
!$omp default(private)   & 
!$omp shared(  mthmodeup, mthmodedown)  
! !$omp shared( m2, n2, NRcalc, om_R_oversos, Maxial, NV, R, c_ctip, tipc_to_tipr )&
! !$omp shared( c_ctip_over2, stag, pi, cc, Ng, p_casc,  Amn2, Bmn2, Xmn2, propmat ) &
! !$omp shared( gmnsup, gmnsdown, kmns   ,Mach, burms, Lx, dk3, wr, dr, q1max, NC)&
! !$omp shared(  mthmodeup, mthmodedown)  

	do m = -m2+1, m2-1	
		if (m .lt. 0)  then 
			muse = -m  ! for locating values in the Xmn, Amn arrays
		else
			muse = m
		endif

!!!!!!!
! start radial integration loop 
!!!!!!
	rintup = 0.0
	rintdown = 0.0
	smnup = 0.0
	smndown = 0.0

	do iR = 1, NRcalc
	dr = (R(2) - R(1))
!!!!!!
! k1 is basically fixed , k2 changes with q1, k3 is integration variable
! q1max is limit for the summation in the final equation for power of a given mode
!  k1 and k2 nondim with chord/2
	!  effk1 and effk2 nondim with R_tip
!!!!!!
	effk1 = om_R_oversos/Maxial(iR)    ! omega/U1 * R_tip   assumes flow aligned with l.e.
		do q1 = -q1max , q1max
			effk2 = (q1*NV - m )/R(iR)

				mg = m  - q1*NV   ! this is actual mg for innerblade phase angle 
						 ! must do mod and if negative change signs
						 !  will be an integer. 
						 !   from  m = mg + q*Nv
						 !  matches Posson's def of k2 being -mg/r

				mode_loc = mod(mg,NV)
				if  ( mode_loc .lt. 0 )  mode_loc = mode_loc + NV

		k3int = 0.d0
		do jk2 = 1,NC
		do jk1 = 1,NC
		
			if (sweep(iR) .ne. 0) then
				do ik3 = 1, 2*Ng + 1
				!if (ik3 .lt. Ng+1) then
					effk3 = dk3*(ik3-Ng-1)/tipc_to_tipr/c_ctip_over2(iR)

					allk = sqrt(effk1**2 + effk2**2 + effk3**2)
				!!! Liepmann
					phiweight = (effk1**2 + effk3**2)/( 1 + (allk*Lx(iR))**2)**3

					k3int(jk1,jk2) =  k3int(jk1,jk2) +   &
					  p_casc(jk1,iR,mode_loc+1,ik3)* &
										  dconjg(p_casc(jk2,iR,mode_loc+1,ik3))  *&
												phiweight* sinc(effk3*Lr(iR)/2)
				end do
			else
			
			do ik3 = 1, Ng + 1
				effk3 = dk3*(ik3-1)/tipc_to_tipr/c_ctip_over2(iR)


				allk = sqrt(effk1**2 + effk2**2 + effk3**2)
			!!! Liepmann
				phiweight = (effk1**2 + effk3**2)/( 1 + (allk*Lx(iR))**2)**3


				if (ik3 .ne. 1) then    !accounts for negative values of k3
					k3int(jk1,jk2) =  k3int(jk1,jk2) +   &
                                            p_casc(jk1,iR,mode_loc+1,ik3)* &
				          dconjg(p_casc(jk2,iR,mode_loc+1,ik3))  * &
						phiweight*2* sinc(effk3*Lr(iR)/2)


				else
					k3int(jk1,jk2) =  k3int(jk1,jk2) +   &
					  p_casc(jk1,iR,mode_loc+1,ik3)* &
                                          dconjg(p_casc(jk2,iR,mode_loc+1,ik3))  *&
                                                phiweight* sinc(effk3*Lr(iR)/2)
				endif



			end do				
			end if ! if non-zero sweep
		end do  !jk1
		end do  !jk2

				k3int = k3int*dk3*2/pi**2*Lx(iR)**5*(2*pi)**3  !extra (2pi)^3 from transform of phi
										!! finish Liepmann


		!!! now integrate with z2 and then z1 using Filon's integration rule for oscillatory functions

		Psi = 0.d0
			do n = 1, n2

			if (propmat(muse+1,n) .ne. 0) then 
				arg = Xmn2(muse+1,n)*R(iR)
				call rmode(m, arg, Amn2(muse+1,n), Bmn2(muse+1,n), Psi,ierPsi)



				muup = -(-m*dsin(stag(iR))/R(iR) + gmnsup(muse+1,n)*dcos(stag(iR)))*tipc_to_tipr*c_ctip_over2(iR)
					! so it is normalized by c/2 like z1 and z2
				mudown = -(-m*dsin(stag(iR))/R(iR) + gmnsdown(muse+1,n)*dcos(stag(iR)))*tipc_to_tipr*c_ctip_over2(iR)


!		muup = muup/2
!		mudown = mudown/2
!	if (m .eq.1)	write(*,*) 'mus up', m, muse, n, muup, gmnsup(muse+1,n)
!	if (m .eq. 1) write(*,*) 'mus down', m, R(iR), muse, n, mudown, gmnsdown(muse+1,n)
!!! prev divisions only to match old code
!!!

		is = (/ (ii, ii = 1,NC-1)/)
		zn = (is -0.5d0)*pi/NC
	 	zn1 = (is + 1- 0.5d0)*pi/NC
		z2loop = 0.d0
! do upstream
		do jk1 = 1,NC
		g0 = 0.5d0*muup*(dcos(zn) + dcos(zn1))
		g1 = 0.5d0*muup*(dcos(zn) - dcos(zn1))
		g2 = 0.5d0*(dsin(g1) + cc*(dsin(g1)/g1 - dcos(g1)))/g1	

		cx = (dconjg(g2)*k3int(jk1,1:NC-1)+g2*k3int(jk1,2:NC)) *(zn1-zn)
		z2loop(jk1) =  sum( cx*cdexp(cc*g0))

		end do

		g0 = 0.5d0*(-muup)*(dcos(zn) + dcos(zn1))
		g1 = 0.5d0*(-muup)*(dcos(zn) - dcos(zn1))
		g2 = 0.5d0*(dsin(g1) + cc*(dsin(g1)/g1 - dcos(g1)))/g1	
		cx = (dconjg(g2)*z2loop(1:NC-1)+g2*z2loop(2:NC))*(zn1-zn)	

		cmns = sum(cx*cdexp(cc*g0))
!		capr = Psi*(m*dcos(stag(iR))/R(iR) + gmnsup(muse+1,n)*dsin(stag(iR)))*burms(iR)*c_ctip_over2(iR)*Mach(iR)
! changing for nondim stuff  6/25/14
!		capr = Psi*(m*dcos(stag(iR))/R(iR) + gmnsup(muse+1,n)*dsin(stag(iR)))*burms(iR)*c_ctip_over2(iR)*tipc_to_tipr
!                rintup(q1+q1max+1,n,iR) = cmns*capr**2*Mach(iR)**3* wr(iR)*dr  !put Simpson's weight in
! changind for delta r stuff  7/17/14
		capr = Psi*(m*dcos(stag(iR))/R(iR) + gmnsup(muse+1,n)*dsin(stag(iR)))
		rintup(q1+q1max+1,n,iR) = cmns*capr**2*Mach(iR)*(Mach(iR))**2*burms(iR)**2*Lr(iR)*c_ctip_over2(iR)*tipc_to_tipr* wr(iR)*dr

! do downstream

		z2loop = 0.d0
		do jk1 = 1,NC
		g0 = 0.5d0*mudown*(dcos(zn) + dcos(zn1))
                g1 = 0.5d0*mudown*(dcos(zn) - dcos(zn1))
                g2 = 0.5d0*(dsin(g1) + cc*(dsin(g1)/g1 - dcos(g1)))/g1

                cx = (dconjg(g2)*k3int(jk1,1:NC-1)+g2*k3int(jk1,2:NC))*(zn1-zn)
                z2loop(jk1) =  sum( cx*cdexp(cc*g0))
		end do 

                g0 = 0.5d0*(-mudown)*(dcos(zn) + dcos(zn1))
                g1 = 0.5d0*(-mudown)*(dcos(zn) - dcos(zn1))
                g2 = 0.5d0*(dsin(g1) + cc*(dsin(g1)/g1 - dcos(g1)))/g1
                cx = (dconjg(g2)*z2loop(1:NC-1)+g2*z2loop(2:NC))*(zn1-zn)

                cmns = sum(cx*cdexp(cc*g0))
!                capr = Psi*(m*dcos(stag(iR))/R(iR) + gmnsdown(muse+1,n)*dsin(stag(iR)))*burms(iR)*c_ctip_over2(iR)*Mach(iR)
! changing for nondim stuff
!                capr = Psi*(m*dcos(stag(iR))/R(iR) + gmnsdown(muse+1,n)*dsin(stag(iR)))*burms(iR)*c_ctip_over2(iR)*tipc_to_tipr
!                rintdown(q1+q1max+1,n,iR) = cmns*capr**2*Mach(iR)**3*wr(iR)*dr  !put Simpson's weight in
! changind for delta r stuff  7/17/14
		capr = Psi*(m*dcos(stag(iR))/R(iR) + gmnsdown(muse+1,n)*dsin(stag(iR)))
		rintdown(q1+q1max+1,n,iR) = cmns*capr**2*Mach(iR)*(Mach(iR))**2*burms(iR)**2*Lr(iR)*c_ctip_over2(iR)*tipc_to_tipr * wr(iR)*dr


	
			endif



			end do  ! n loop

		end do   ! q1 loop

	end do  ! iR loop

! finish integration wtih respect to r

	smnup = sum(dreal(rintup), dim = 3)   !should be size 2*q1max+1 , n2

	smndown = sum(dreal(rintdown), dim = 3)

	
	
! store the m'th mode's integrated and  q1 summed value   ---- put negative m together with positive m

	mthmodeup(muse+1,:) = sum(smnup,dim=1) + mthmodeup(muse+1,:)

	mthmodedown(muse+1,:) = sum(smndown,dim=1) + mthmodedown(muse+1,:)



	end do  !  m loop
!$omp end parallel  do


	power_nondim(:,:,1) = mthmodeup*((1.d0-M_ann**2)*NV)**2* om_R_oversos/8.d0 /pi**2 &
		/kmns/ (1-hub_to_tip**2)/(om_R_oversos+M_ann*kmns)**2
	power_nondim(:,:,2) = mthmodedown*((1.d0-M_ann**2)*NV)**2* om_R_oversos/8.d0 /pi**2 &
		/kmns/ (1-hub_to_tip**2)/(om_R_oversos-M_ann*kmns)**2


end

!***********************************************************************************************************************************
!  SINC
!
!  Sine cardinal (sinc) function.
!***********************************************************************************************************************************

      FUNCTION SINC (X) RESULT (Y)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: X
      DOUBLE PRECISION :: Y

      IF (X .EQ. 0.0D0) THEN
         Y = 1.0D0
      ELSE
         Y = SIN(X)/X
      END IF

      RETURN

      END FUNCTION SINC
