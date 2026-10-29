!!
!!
!!  BU - RSI code 
!!  Computation of exhaust broadband noise due to wake turbulence interacting with an exit guide vane
!!
!!  Originally based on work by Nallasamy and Envia, this code is now quite different than the NASA RSI code
!!  The theory that is represented in this code is described in the 2014 JSV paper found in the DOCS directory
!!
!!  This code was created by Sheryl Grace in 2015 and updated by Nuo Li in 2024 with python
!!  wrapper and sweep modeling
!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none

	integer, parameter :: rk = selected_real_kind(p=12,r=307)

	integer  NV, NB, NRdat, NRcalc, NC, NF, iR, jR , Ng, jk, ik3,  &
		m, ic, ierp, iFRQ, iSTDYL, nmax, mmax, ii,jj

	real(rk)   T, rho,  M_ann, tipc_to_tipr, hub_to_tip, &
     		RDUCT, BW,  conv_tol_xmn, sos, dk3
	real(rk), allocatable :: radius(:), chord(:),sweepdeg(:), &
     		stagger(:), Mx(:), Mt(:),  R(:), freq(:)
	real(rk), allocatable  :: data(:,:), datainterp(:,:)
	complex(rk), allocatable ::  p_casc(:,:,:,:), deltap(:)
	real(rk), allocatable  :: c_ctip(:), stag(:), Maxial(:), Mcirc(:) , &
                c_ctip_over2(:), h(:), h_chalf(:), Mach_ang(:) , Mach(:), &
                beta(:), sweep(:), beta_s(:)
	real(rk) Ts, T_corr, pi, dR, weight, theta, k1, k3, sigma
	character(50) :: line
	integer  ieb(4), iec, ier(4), n, mloc
	integer, allocatable :: ansloc(:)
	real(rk) xmn_dat, amn_dat, bmn_dat, Rg, gamma, gc, dB_ref, BWc,  Xmax
	real(rk)  sPWRup, sPWRdown, sPWRdown_dim
	real(rk), allocatable :: sPWRdbup(:), sPWRdbdown(:)
	real(rk), allocatable :: backurms(:), maxurms(:), urmsgausswidth(:), &
		Lxdat(:), Ltdat(:), Lrdat(:), sweep2tip_sweep(:)
	real(rk), allocatable :: Xmn(:,:), Amn(:,:), Bmn(:,:), mattest(:,:)
	real(rk), allocatable :: Xmn2(:,:), Amn2(:,:), Bmn2(:,:), propmat(:,:)
	real(rk), allocatable ::   power_nondim(:,:,:), sPWRmup(:), sPWRmdown(:)
	namelist /inputconst/ NV, NB, NRdat, NRcalc, NC,  NF, Ng  &
                      T, rho,  M_ann, iSTDYL, tipc_to_tipr,    &
                        hub_to_tip, RDUCT, BW, conv_tol_xmn
	 namelist /inputmat/ radius, chord, sweepdeg, stagger, Mx, Mt

	namelist /inputturbmat/ backurms, maxurms, urmsgausswidth,  &
		Lxdat, Ltdat, Lrdat
	integer :: iunit, ierr



      open (unit = 10, file = 'aero95const.input', form = 'formatted',  &
            status = 'unknown')
      open (unit = 11, file = 'aero95mat.input', form = 'formatted',  &
            status = 'unknown')
      open (unit = 14, file = 'turb95mat.input', form = 'formatted',  &
            status = 'unknown')
open (unit = 22, file = 'fs.dat', form = 'formatted',  &
            status = 'unknown')


        read(10,inputconst)     !reads size paramaters NV, NB, NRdat
                                !    NRcalc, NC as well as other params


        allocate(radius(NRdat),chord(NRdat),stagger(NRdat),  &
                Mx(NRdat),Mt(NRdat), sweepdeg(NRdat))

        read(11,inputmat)       !reads in the geom and flow info
		
		
	
    if (NF .eq. 0)  then
		ierr = 0
		do
			if (ierr /= 0) exit
			read(22,*, iostat=ierr) line
			NF = NF + 1
		end do
	endif
	NF = NF - 1
	!write(*,*) 'nf', NF	
	allocate( freq(NF) )
	rewind(22)
	read(22,*) freq
	!write(*,*) 'nf, freq',SIZE(freq), freq
!!!!!!!!!!!
! adjust input if needed
!!!!!!!!!!!
        pi = dacos(-1.d0)      
        !sweep = sweepdeg(:) * pi/180.d0
        Ts = 518.00      !standard day temperature - use to correct things
!       T_corr = dsqrt(Ts/T)
        T_corr = 1.d0
!       M_ann = M_ann*T_corr
!       M_rtip = M_rtip*T_corr

	gc =  3.2174d+01
	Rg = 1.716d+03
	gamma = 1.4d0
	sos = dsqrt(gamma*Rg*T)
	allocate(backurms(NRdat), maxurms(NRdat), urmsgausswidth(NRdat), &
		Lxdat(NRdat), Ltdat(NRdat), Lrdat(NRdat) )

	read(14,inputturbmat)  !reads in the turbulence info

!!!!!!!!!!!
! interpolate specified values onto radial locations needed for calculation
!!!!!!!!!!!


        allocate(data(NRdat,12))
!NL 	Added sweepdeg as function of radius		
	data(:,12) = sweepdeg(:) * pi/180.d0
        data(:,1) = radius(:)
        data(:,2) = chord(:) !* dcos(data(:,12)) !NL sweep coordinate transformation. c_s = c*cos(psi)
        data(:,3) = stagger(:)*pi/180.d0
        data(:,4) = Mx(:)*T_corr
        data(:,5) = Mt(:)*T_corr
	data(:,6) = backurms(:)
	data(:,7) = maxurms(:)
	data(:,8) = urmsgausswidth(:)
	data(:,9) = Lxdat(:)
	data(:,10) = Ltdat(:)
	data(:,11) = Lrdat(:)
	



	        dR = (1.d0 - hub_to_tip)/(NRcalc-1)
        allocate( R(NRcalc) )
        allocate( datainterp(NRcalc,12) )
        R = hub_to_tip + (/ ((iR -1)*dR, iR = 1, NRcalc )  /)

        do iR = 1,NRcalc
                jR = 1
             do while  (radius(jR) .lt. R(iR) .and. jR .lt. NRdat)
                        jR = jR +  1
                end do

	if (jR .eq. 1)  then
                weight = (R(iR) - radius(2))/(radius(1)-radius(2))
        else  if (jR .eq. NRdat .and. R(iR) .gt. radius(jR)) then
        	weight = (R(iR) - radius(jR-1))/(radius(jR)-radius(jR-1))
        else
        	weight = (R(iR) - radius(jR-1))/(radius(jR)-radius(jR-1))
        endif
	if (weight .ne. 1.d0) then
                if (jR .eq. 1) then
                datainterp(iR,:) = weight*(data(1,:)-data(2,:)) + data(2,:)
		else if (jR .eq. NRdat .and. R(iR) .gt. radius(jR)) then
		datainterp(iR,:) = weight*(data(jR,:) - data(jR-1,:)) + data(jR-1,:)
		else
                datainterp(iR,:) = (1.d0-weight)*data(jR-1,:)+weight*data(jR,:)
		endif
        else
                datainterp(iR,:) = data(jR,:)
        end if
        end do

        deallocate(data)

	
        allocate( c_ctip(NRcalc), stag(NRcalc), Maxial(NRcalc),Mcirc(NRcalc),sweep2tip_sweep(NRcalc) )
        stag = datainterp(:,3)
        Maxial = datainterp(:,4)
        Mcirc = datainterp(:,5)
		sweep = datainterp(:,12)
		!sweep2tip_sweep = sweep(:) / sweep(NRcalc)
		c_ctip = datainterp(:,2) 
        allocate(c_ctip_over2(NRcalc), h(NRcalc), h_chalf(NRcalc) )
        c_ctip_over2 = c_ctip/2.d0
        h = 2.d0 * pi * R/NV  / c_ctip /tipc_to_tipr  ! gap to chord
        h_chalf = 2.d0 * pi * R/NV  / c_ctip_over2 /tipc_to_tipr  ! gap to semichord
		
        allocate(Mach_ang(NRcalc), Mach(NRcalc) ,beta(NRcalc), beta_s(NRcalc) )
        Mach_ang = datan(Mcirc/ Maxial)   ! absolute flow angle
        Mach = dsqrt(Maxial**2 + Mcirc**2) * &  !absolute Mach number
                !dcos(sweep) * &                 ! sweep correction
                dcos(iSTDYL*(stag - Mach_ang))  ! turning correction
        beta = sqrt(1 - Mach**2)
		beta_s = sqrt(1 - (Mach*dcos(sweep))**2)

!!!!!!!!!!!
!  read in frequencies of interest
!!!!!!!!!!!
 

!!!!!!!!!!!!
! find the duct eigenvalues and functions  do for circ and axial modes much greater than will need
!!!!!!!!!!!!
	mmax = 801
	nmax = 100
	allocate(  Xmn(mmax,nmax), Amn(mmax,nmax), Bmn(mmax,nmax) )

!!!
! only compute if less than max for cut off for highest frequency
!!!
		Xmax = maxval(freq)/sos*2*pi*RDUCT/12/sqrt(1-M_ann**2)

	Xmn = 1000000    ! initialized like this so when looking for cut-on modes can find them
	Amn = 0
	Bmn = 0
	m = 0
	n = 3
	do while ( m .lt. mmax .and. n .ne. 2)	
		n = 1
		xmn_dat = 0
		do while (xmn_dat .lt. Xmax .and. n .le. nmax)
			call anrt(m,n,hub_to_tip,conv_tol_xmn,xmn_dat,ieb,iec)
			call eigen(m,hub_to_tip,xmn_dat,amn_dat,bmn_dat,ier)
				if (xmn_dat .lt. Xmax) then
			Xmn(m+1,n) = xmn_dat
			Amn(m+1,n) = amn_dat
			Bmn(m+1,n) = bmn_dat
	write(*,*) 'm,n,xmn', m,n,xmn_dat
				endif
		n = n + 1
		end do
	m = m + 1
	end do

!!!!!!!!!!
! set the reference dB value (nondim by rho c^2 R_tip^3) and bandwidth correction
!!!!!!!!!	
	!1.3558d+12 == slug/lbm/(10^(-12)W) or slug per lbm over reference power in watt
	dB_ref = 10.d0 * dlog10(1.3558d+12 * rho/gc * (gamma*Rg*T)*(RDUCT/12.d0)**3)
	BWc = 10.d0 * dlog10(BW)


!!!!!!!!!!!!!
! for each frequency, add all axial and circular mode power contributions for cut=on modes
!!!!!!!!!!
        allocate( deltaP(NC) , p_casc(NC, NRcalc, NV ,  2*Ng+1 ) )
	allocate (mattest(mmax,nmax), ansloc(mmax))
	allocate(sPWRdbup(NF), sPWRdbdown(NF))

	dk3 = 15.29*RDUCT*tipc_to_tipr/Ng

	do iFRQ = 1, NF
		
		p_casc = 0.d0   ! clear p_casc for each frequency 

		do  iR = 1, NRcalc

                do m = 0,NV-1
		k1 = 2.d0*pi*freq(iFRQ)*c_ctip_over2(iR)*tipc_to_tipr*RDUCT/12/Mach(iR)/sos
                                  !gives omega c/2 U

                        sigma = 2.d0*pi*m/NV
!NL 						k3 goes from -ng to 0 to + ng
					if (sweep(iR) .ne. 0) then
                        do ik3 = 1, 2*Ng + 1
						! k3new = k1*dsin(sweep(iR))+k3*dcos(sweep(iR)), k1new = k1*dcos(sweep(iR))-k3*dsin(sweep(iR)) after rotation
                                ierp = 0
                                k3 = dk3*(ik3-Ng-1)	
								!write(*,*) 'ik3, k3',ik3, k3
                                if (k3 .eq. 0) then
								
								theta = 10.d0
								!else if (k3 .lt. 0) then
								!theta = -k1/k3*Mach(iR)/beta(iR)
								else
                                theta = (k1*dcos(sweep(iR))-k3*dsin(sweep(iR)))/k3*Mach(iR)*dcos(sweep(iR))/beta_s(iR)
								!theta = k1/(k1*dsin(sweep(iR))-k3*dcos(sweep(iR)))*Mach(iR)/beta(iR)
                                endif
                                                ! set all subcritical gust responses
                                                !  to zero												
                                if (theta**2 .lt. 1.d0) then
                                        deltaP = (0.d0,0.d0)
                                else
								!write(*,*) 'ir, k1, k3',iR, k1, k3
								!write(*,*) 'ik3, k3, theta', ik3, (k1*dsin(sweep(iR))-k3*dcos(sweep(iR))), theta	
								!write(*,*) 'k1, k3, theta', k1*dcos(sweep(iR))-k3*dsin(sweep(iR)), (k1*dsin(sweep(iR))+k3*dcos(sweep(iR))), theta
                                        call casc_k3_no_r(Mach(iR),  &
                                          k1, k3, sweep(iR), sigma, h_chalf(iR), &
                                          stag(iR), NC, deltaP, ierp)
                                                        ! passing sigma or -sigma
                                                        ! may make difference!!! orig had -sigma
                                endif
								!write(*,*) 'k1, k3, theta', k1*dcos(sweep(iR))+k3*dsin(sweep(iR)), (k1*dsin(sweep(iR))+k3*dcos(sweep(iR))), theta
									
				   p_casc(:,iR,m+1,ik3) = deltaP(:)
						end do   !ik3 Neg.
					else !If sweep = 0, compute only half k3
                        do ik3 = 0, Ng

                                ierp = 0
                                k3 = dk3*ik3
                                if (k3 .gt. 0) then
                                theta = k1/k3*Mach(iR)/beta(iR)
                                                ! set all subcritical gust responses
                                                !  to zero	
                                else
                                theta  = 10.d0
                                endif
                                if (theta**2 .lt. 1.d0) then
                                        deltaP = (0.d0,0.d0)
                                else
                                        call casc_k3_no_r(Mach(iR),  &
                                          k1, k3, sweep(iR), sigma, h_chalf(iR), &
                                          stag(iR), NC, deltaP, ierp)
								
                                endif
				   p_casc(:,iR,m+1,ik3+1) = deltaP(:)
                        end do   !ik3					
					endif !if sweep loop	

                if (ierp .ne. 0) write(6, 1100) ierp


                end do    ! m

        end do   ! iR



		sPWRup = 0.d0    !upstream value
		sPWRdown = 0.d0    ! downstream value 
		mattest = 0.d0
		ansloc = 0.d0
		sPWRdown_dim = 0.d0
		Xmax = freq(iFRQ)*2*pi*RDUCT/12/sos/sqrt(1-M_ann**2)
		mattest = -Xmn + Xmax	
		ansloc = minloc(mattest,dim=2, mask = mattest  .gt. 0)  ! will have the nth locations of cutoff

		if (ansloc(mmax) .eq. nmax) then 
			mloc = mmax
		else 
		mloc = minloc(ansloc, dim=1, mask  = ansloc .eq. 0)  - 1! will have the mth location of cutoff
			if (mloc .eq. 0) mloc = 1
		endif


	allocate (Xmn2(mloc,ansloc(1)), Amn2(mloc, ansloc(1)), Bmn2(mloc,ansloc(1)), propmat(mloc, ansloc(1)))
	Xmn2 = Xmn(1:mloc, 1:ansloc(1))
	Amn2 = Amn(1:mloc, 1:ansloc(1))
	Bmn2 = Bmn(1:mloc, 1:ansloc(1))
	propmat = 0
	do m = 1,mloc
		Xmn2(m,ansloc(m)+1:ansloc(1)) = 0.d0
		propmat(m,1:ansloc(m)) = 1
		Amn2(m, ansloc(m) + 1: ansloc(1)) = 0.d0
                Bmn2(m, ansloc(m) + 1: ansloc(1)) = 0.d0
	end do

	allocate(power_nondim(mloc,ansloc(1),2), sPWRmup(mloc), sPWRmdown(mloc))	
	power_nondim = 0.d0
		call get_Pmns(iFRQ,mloc,ansloc(1),NC,NRcalc,NV,NF,Ng,dk3,datainterp,p_casc,Xmn2, Amn2, Bmn2,propmat, &
				M_ann,freq(iFRQ),RDUCT,sos,hub_to_tip,tipc_to_tipr,sweep,iSTDYL,power_nondim)

	deallocate(propmat, Xmn2, Amn2, Bmn2)
        sPWRmup = sum(power_nondim(:,:,1),dim=2)
        sPWRmdown = sum(power_nondim(:,:,2),dim=2)
	sPWRup = sum(sPWRmup)
	sPWRdown = sum(sPWRmdown)
	sPWRdown_dim = sPWRdown* rho * (gamma*Rg*T)*(RDUCT/12.d0)**3/gc*1.3558
        sPWRdbup(iFRQ) = 10.d0*dlog10(sPWRup)+dB_ref+BWc
        sPWRdbdown(iFRQ) = 10.d0*dlog10(sPWRdown)+dB_ref+BWc
	
	write(*,*) 'before', sPWRup, sPWRdown
	write(*,*) freq(iFRQ), sPWRdbup(iFRQ), sPWRdbdown(iFRQ)

	open (unit = 33, file = 'fort.33', form = 'formatted', status = 'unknown', position = 'append')
	write(33,*) freq(iFRQ), sPWRdbup(iFRQ), sPWRdbdown(iFRQ)
	close(33)
	deallocate(power_nondim, sPWRmup, sPWRmdown)
	end do

1100 format(' ** Warning **: Error code from inversion routine = ', I3)

end
