!* Uses parts of the mainspec.f and mainaero.f codes to get the Xmn
!* needed to compute the edge frequencies in the duct
!*
!***********************************************************************
!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8

      program main

	implicit none
	integer  NV, NB, NRdat, NRcalc, NC, NF,  Ng,  &
		m, n, iFRQ, iSTDYL, ieb(4), iec

	integer, parameter :: rk = selected_real_kind(p=12,r=307)

	real(rk)   T, rho,  M_ann, tipc_to_tipr, hub_to_tip, &
     		sweepdeg, RDUCT, BW, conv_tol_xmn, xmn_dat
	real(rk), allocatable :: radius(:), chord(:), &
     		stagger(:), Mx(:), Mt(:),  R(:)

	 namelist /inputconst/ NV, NB, NRdat, NRcalc, NC NF, Ng  &
                      T, rho,  M_ann, iSTDYL, tipc_to_tipr,    &
     			hub_to_tip, sweepdeg, RDUCT, BW, conv_tol_xmn


      open (unit = 10, file = 'aero95const.input', form = 'formatted',  &
            status = 'unknown')


	read(10,inputconst)     !reads size paramaters NV, NB, NRdat
 				!    NRcalc, NC as well as other params

      do 20 iFRQ = 1, NF 	! begin frequency sweep



       do 30 m = 0, 801     	! count up in m untill m limit is reached
        do 40 n = 1, 100		! count up in n untill cutoff is reached

         call anrt (m, n, hub_to_tip, conv_tol_xmn, xmn_dat, ieb, iec)
		write(102,*) m,n,xmn_dat

   40  continue
   30  continue
   20  continue

      stop
       end
