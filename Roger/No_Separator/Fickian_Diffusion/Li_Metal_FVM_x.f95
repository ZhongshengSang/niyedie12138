!  1-D, transient, rectangular, electrode problem solved by the method of Newman.
! Simple Diffusion
! Nicholas Brady 4/20/2017

!   The user must specify
!     1. parameters for the reaction kinetics, operating conditions, and electrode structure
!     2. an initial guess
!     3. the coefficients aik, bik, dik, gi, pik, eik, fi
!       (see J.S. Newman, Electrochemical Systems, p. 540-)

!________________________________________________
! program description
!       Finite Volume Formulation
!       Diffusion through a rectangular medium
!       First Order Reaction
!       Fixed Concentrations at the boundaries
!
! ********** Governing Equations **********
!       electrolyte concentration (1): D*d2cdx2 = dcdt

!
! ********** Boundary Conditions **********
!       at x=0,     dcdx=flux,   where flux is proportional to the current
!       at x=xmax,  c=1      
!       at t=0,     c=1      
!________________________________________________
!******************************INPUT MODULE******************************************

module user_input
SAVE    
     parameter(N=1,NJ=12,Numbertimesteps=3.6d3*1000,tmax=3.6d3*1000)
     ! 3.6d3 = 1 hour
     
          ! *********************** ELECTRODE PARAMETERS ***********************
     double precision :: xmax                                     ! Distance between separator and anode

          ! *********************** GENERAL  PARAMETERS ***********************
     double precision :: Rigc = 8.314,Temp = 298, Fconst = 96485  ! Ideal gas constant [J/(mol*K)], Temperature [K], Faraday's Constant [C/mol]
     double precision :: rxn_k=5.62d-9, alpha_a = 0.5, alpha_c = 0.5             ! cathodic and anodic charge transfer coefficients
     real             :: PI=3.141592654                           ! Geometric constant

          ! *********************** ELECTROLYTE  ***********************
     double precision :: eps = 0.26             ! electrode porosity, diffusion coefficient
     double precision :: c0, cbulk=0.001, c0min=1.0d-22           ! electrolyte concentration, bulk electrolyte conc, and minimum conc
                    ! CATION
     double precision :: z_cat = 1.0, diff_cat = 1.00d-6
                    ! ANION
     double precision :: z_an = -1.0, diff_an = 1.00d-6/2.0
          ! *********************** ELECTRODE MATERIAL PHYSICAL PARAMETERS ***********************
     double precision :: sigma=4.269d-2, cimax=1.788337d-1        ! solid-state electronic conductivity
     double precision :: vmin=1.0                                 ! minimum voltage allowed (numerically)
     double precision :: molar_mass_Fe3O4 = 231.533, density_Fe3O4 = 5.15 ! [g/mol],  [g/cm3]

          ! *********************** INITIAL CONDITIONS ***********************
     double precision :: Phi_1_init=3.3, c0_init=0.001, cs_init = 1.0d-5   ! initial conditions: voltage, electrolyte conc, solid-state conc.

          !****************** CHANGES BY EXPERIMENTAL CONDITIONS ***********************
     double precision :: Discharge_Time, Discharge_Relax_Time, Charge_Time, Charge_Relax_Time           !
     double precision :: i_applied

          !****************** COULOMB COUNTING ******************
     double precision :: mAhg = 0.0, LixFe3O4

          !****************** VARIABLE: STATE OF EXPERIMENT ***********************
     CHARACTER(LEN=1) :: state
     CHARACTER(LEN=4) :: electron_equiv
     integer(kind=8)  :: Discharge


        ! N: number of dependent variables
        ! NJ: number of node points (for 1-D, 101 is typical)
        ! tmax: duration of the simulation [=] seconds
        ! c_on: time when current is applied [=] seconds

        ! Rigc: ideal gas constant [=] J/mol/K
        ! Temp: temeperature [=] K
        ! c_applied: applied current density [=] A/cm^2
        ! c0min_fail: lower bound on unreacted concentration in pores [=] mol/cm^3
        ! cimax_fail: upper bound on reacted concentration [=] mol/cm^3
        
        ! eps: porosity
        ! diff0: diffusivity coefficent [=] cm^2/s
        ! sigma: electronic conductivity [=] S/cm
        ! spec_a: specificy surface area [=] cm^2/cm^3
        ! xmax: radius of the particle [=] cm
        ! cbulk: concentration of lithium ions in the electroltyte [=] mol/cm^3
        
        ! rxn_k: reaction rate [=] 1/(mol*cm)**(0.5)/s
        ! alpha_a/c: charge transfer coefficient
        ! cimax: maximum concentration of reacted lithium in Fe3O4 crystal [=] mol/cm^3
!_______________________

               
end module user_input

!************************OPEN CIRCUIT POTENTIAL MODULE**********************************
 module mod_OSP
  use user_input
  implicit double precision(a-h,o-z)
  SAVE

contains 
 FUNCTION OSP(c00)
  double precision :: OSP, cii, xish, c00
  double precision :: U_ref = 3.3

      xish = cii/cimax

      OSP = U_ref + Rigc*Temp/Fconst * DLOG(c00/cbulk)

      !! Nernstian Potential Lithium Metal vs RHE


RETURN
END FUNCTION OSP

! functions to find the derivatives of the open-circuit voltage

!  FUNCTION DOSPDC(yyy,xxx)
!   double precision :: DOSPDC, xxx,yyy
!   DOSPDC=(OSP(yyy,xxx*1.01)-OSP(yyy,xxx*0.99))/(0.02*xxx)

! RETURN
! END FUNCTION  DOSPDC

!  FUNCTION DOSPDC0(yyy,xxx)
!   double precision :: DOSPDC0, xxx,yyy
!   DOSPDC0=(OSP(yyy*1.01,xxx)-OSP(yyy*0.99,xxx))/(0.02*yyy)

! RETURN
! END FUNCTION  DOSPDC0

end module mod_OSP

!************************ ELECTROCHEMICAL REACTION MODULE **********************************

 module MOD_echem_rxn
  use user_input
  use mod_OSP
  SAVE

  contains
    FUNCTION EX_I(c00,cii) ! exchange current
      double precision :: EX_I, c00, cii
    ! c0 is the electrolyte concentration
    ! ci is the solid-state concentration

    ! rxn_k is the reaction rate [=] 1/(mol*cm)**(0.5)/s
      

      EX_I = Fconst*rxn_k*(c00**alpha_a)*((cimax-cii)**alpha_a)*(cii**alpha_c)

    RETURN
    END FUNCTION EX_I


    ! Electrochemical Reaction
    ! ex_i * [exp(alpha_a*F*eta/RT) - exp(-alpha_c*F*eta/RT)]
    ! where eta = U

    FUNCTION echem_rxn(c00, cii, Phi1, Phi2) ! (c00,cii,eta) , (c00,cii,phi_1, phi_2)
      double precision :: echem_rxn, eta, U0, Phi1, Phi2
      double precision :: c00, cii

    U0 = OSP(c00)
!     eta = -0.05/2
    eta = Phi1 - Phi2 - U0
!     eta = Phi_1 + Phi_2 - OCP

    ex_curr = EX_I(c00,cii)

    echem_rxn = ex_curr * ( DEXP(alpha_a*Fconst*eta/(Rigc*Temp)) - DEXP(-alpha_c*Fconst*eta/(Rigc*Temp)) )



    RETURN
    END FUNCTION echem_rxn

 end module MOD_echem_rxn


!******************************VARIABLE MODULE******************************************
! do not edit this module

module variables
use user_input
SAVE
 
! Declare array/matrix sizes depending on the user inputted values.
     double precision, dimension(N,N+1,NJ) :: E
     double precision, dimension(N,2*N+1) :: D
     double precision, dimension(N,N) :: A,B,X,Y,smA,smB,smD,smE,smP
     double precision, dimension(N,NJ) :: cprev, delC
     double precision, dimension(N) :: dcdx,d2cdx2,G,smG,smF,ID
     double precision, dimension(NJ) :: xx,delx ! Investigate if this line is necessary: delx, xx
     double precision :: delT,time                     ! Investigate if this line is necessary: delT, time
     integer :: NP1,t1,t2,clock_rate,clock_max
     !****************** FVM VARIABLES ***********************
     double precision, dimension(N,N) :: fE, fW, dE, dW, rj
     double precision, dimension(N)   :: cE, cW, dcdxE, dcdxW
     double precision :: alphaE, alphaW, betaE, betaW

end module variables

!**********************SUBROUTINES FOR WRITING DATA******************************
!____________________________________________________________________________
!   This section of code stores all of the subroutines used to write data files.
!   It is up to the user to edit the "call" in the main program to get the right data.
!____________________________________________________________________________

subroutine write_all_voltage(it)
      use user_input
      use variables
      use mod_OSP
      use MOD_echem_rxn
      implicit double precision(a-h,o-z)
      double precision :: i_rxn_0
      double precision :: Phi_1, Phi_2
      double precision :: to_electrons ! multiply by this to convert mol/cm3 to electrons (mol/mol)

      ! Write data at all positions
      open(56, file = 'Time_Voltage_Position.txt', status = 'unknown')

      ! Write just edge information
      open(57, file = 'Time_Voltage.txt', status = 'unknown')

      t_write=time/float(3600)

! Just the  information from node NJ
        c0 = cprev(1,1)
        
        U0 = OSP(c0)
        
!         i_rxn_0 = echem_rxn(c0,cs,Phi_1,Phi_2)

        LixFe3O4 = mAhg*molar_mass_Fe3O4*3.6/Fconst ! 1 mAh = 3.6 C

        to_electrons = 1.0/density_Fe3O4*molar_mass_Fe3O4

! **** TERMINAL SCREEN INFORMATION ****
        if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
          write(*,17) 'State', 'Time', 'Voltage', 'Anode c0', 'Edge cs'
          ! need to write units 

        end if
        write(*,15) state, t_write, U0, c0*1000, cprev(1,2)*1000, cprev(1,3)*1000, cprev(1,NJ)*1000

! **** TEXT FILE INFORMATION ****
        if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
          write(56,17) 'State', 'Time', 'Voltage', 'Position', 'mAhg', 'LixFe3O4', 'Elect_Conc', 'Solid_Conc', 'Reaction'

        end if

        do j = 1,NJ
            c0 = cprev(1,j)

            U0 = OSP(c0)
!             i_rxn_0 = echem_rxn(c0,cs,Phi_1,Phi_2)
          
          write(56,15) state, t_write, Phi_1, xx(j), mAhg, LixFe3O4, c0, cs, i_rxn_0
        
        end do

        if (it.EQ.1) then
          write(57,17) 'State', 'Time', 'Voltage', 'mAhg', 'LixFe3O4', 'Elect_Conc', 'Solid_Conc', 'Reaction'
          write(57,17) 'CDR',   'hours', 'Volts',  'mAhg', 'LixFe3O4', 'mol/cm3'  , 'LixFe3O4'   , ','
        end if  

        write(57,15) state,     t_write, Phi_1,   mAhg,   LixFe3O4,   c0,   cs*to_electrons,   i_rxn_0
!         write(*,15)  state,     t_write, Phi_1,   mAhg,   LixFe3O4,   c0,   cs,   i_rxn_0

! Format for terminal
! numbers
 15    format(1(A5,1X), 3(F12.5,1X), 20(ES15.5,1X))
! headers
 17    format(1(A5,1X), 3(A12,1X),   20(A15,1X))
 

! Format for text file
 16    format(1(A7,1X), 10(ES15.5,1X))
 18    format()

!  16    format(1(A5,1X),  2(A12,1X),  2(A7,1X), 12(A12,1X))
!  17    format(1(F5.2,1X),2(ES12.5,1X),2(F7.5,1X),12(ES12.5,1X))

end subroutine write_all_voltage
      
!****************************** MAIN PROGRAM ******************************************  

program unsteady  
use user_input  
use variables
use mod_OSP
implicit double precision(a-h,o-z)

! *************************************************************************************************************
! ****************************************** EXPERIMENTAL CONDITIONS ******************************************
! *************************************************************************************************************
! Under this line, conditions that change for each experimental run are defined.
! Using conditional statements the simulated parameters are adjusted for different experiments
! Usually these conditional statements are triggered through the correspond shell script

! Experimental specific parameters include
! xmax            - distance from anode to separator
! c_applied       - the applied current in A/g (per gram of active material) ~ related to flux

! Discharge_ON    - time that current is being passed during discharge
! Discharge_OFF   - time that current is OFF after discharge
! Charge_ON       - time that current is passed during charge; usually this time is irrelevant because charge is done until a cutoff voltage is reached
! Charge_OFF      - time that current is OFF after charge
! Const_Volt_Time - Amount of time for the constant voltage hold (at the end of charge)

! Sometimes, the initial conditions will change as well for different experiments or simulations
! Phi_1_init, c0_init, cs_init, mAhg


xmax  = 1.0d-0 ! cm
state = 'D'
i_applied = -6.36801e-5

      
call system_clock(t1,clock_rate,clock_max)
call initial_condition(h)

write_time = 0
write_time_step = 0.5 * 3600.0 ! write data every 0.5 hours

do 1 it=1,Numbertimesteps

    U0 = OSP(cprev(1,1))
!     print *, cprev(1,1), U0
!   call write_all_voltage(it)

! !       write one of the early data points
      if (it.EQ.1) then 
          call write_all_voltage(it)
          write_time = write_time + write_time_step  
      
      elseif (time.GE.write_time) then 
          call write_all_voltage(it)
          write_time = write_time + write_time_step
      
      ! Write a data point near the very end of the program
      elseif (it.GE.(Numbertimesteps)) then 
          call write_all_voltage(it)

      end if


15    format(9(ES20.5,1X))
16    format(9(A20,1X))
20    format(1(F5.2,1X),2(ES12.5,1X),2(F7.5,1X),12(ES12.5,1X))


  call bound_val(h)

! dynamic time-stepping so we can quickly move through the relaxation period.  
! Significantly reduces simulation time during long relaxation periods
      if (state.EQ.'R') then 
        delT = delT * 1.0001 ! * 1.00005

      else
        delT=tmax/float(Numbertimesteps)

      end if

      time=time+delT

 1    continue
call system_clock(t2,clock_rate,clock_max)
write ( *, * ) 'Elapsed real time =', real(t2-t1)/real(clock_rate )
end program unsteady  

!*********************************BOUND VAL**************************************      
!____________________________________________________________________________
!   This subroutine assumes solution of a linear boundary-value problem at 
!   a given time step.  There is thus no need for iteration.
!     c=cprev+delC, we solve for delC using forward time central difference method.
!     It is assumed that the time step is small enough that the linearization is exact.
!     Thus, iterations are not necessary.  This assumption is of course tested by 
!     computer experiments on the effect of time-step size.
!____________________________________________________________________________

subroutine bound_val(h)
      use user_input
      use variables
      implicit double precision(a-h,o-z)
      
! initialize all of the finite volume (matrix) coefficients to zero

      do 150 j=1,nj
      do 123 ic=1,n
      do 123 kc=1,n
        dE(ic,kc) = 0.0
        dW(ic,kc) = 0.0
        fE(ic,kc) = 0.0
        fW(ic,kc) = 0.0
        rj(ic,kc) = 0.0
123      continue

call fillmat(h,j)
! Fillmat is determining small coefficents based on the user input. Must go to
! fillmat and change coefficent definitions based on governing differential equations.

call ABDGXY(h,j)
! ABDGXY equates the large coefficents based on the small coefficents. These equations
! can be found in Newman appendix C.
 
call BAND(J)
! BAND(J) computes delC and calls MATINV to solve the problem using gaussian elimination.

! for all the dependent variables, the values are update as c = cprev + delC

150   continue
        do 12 k=1,n
         do 12 j=1,nj
             cprev(k,j) = cprev(k,j) + delC(k,j)
12       continue    
      return
end subroutine bound_val
    


!******************************INITIAL GUESS****************************************

subroutine initial_condition(h)
      use user_input
      use variables
      implicit double precision(a-h,o-z)

! Define delta t       
      delT=tmax/float(Numbertimesteps)
      h=xmax/float(nj-2)
      do j=1,NJ

        if (j.EQ.1) then
          xx(j) = 0.0
        else if (j.EQ.NJ) then
          xx(NJ) = xmax
        else
          xx(j)=h*float(j-1) - h/2.0
        end if

      end do

      do 11 j=2,NJ-1
         delx(j)=h
11     continue
      ! it is common in the finite volume code to set the control volumes at the boundaries to zero
         delx(1)=0.0d0
         delx(NJ)=0.0d0
          
      do 1 j=1,NJ
        cprev(1,j) = c0_init

1     continue
      return
end subroutine initial_condition


!***********************************FILLMAT****************************************

subroutine fillmat(h,j)
      use user_input
      use variables
      use mod_OSP
      use MOD_echem_rxn
      implicit double precision(a-h,o-z)
      double precision :: i_rxn_0
      double precision :: cs, Phi_1, U
      double precision :: step_c0, step_cs, step_phi_1, step_phi_2
      double precision :: N_W, N_E
      double precision :: u_an, u_cat ! anion and cation mobilities
      double precision :: rxn_coef, diff_eff, max_i_applied

      c0 = cprev(1,j)     ! electrolyte concentration
      
      U = OSP(c0)

      u_cat    = diff_cat/(Rigc*Temp)
      u_an     = diff_an /(Rigc*Temp)

      diff_eff = (z_cat*u_cat*diff_an - z_an*u_an*diff_cat)/(z_cat*u_cat - z_an*u_an)

      max_i_applied = cbulk*diff_eff/xmax

      if ( (j.EQ.1).AND.(time.EQ.0) ) then
        print*, max_i_applied*Fconst, xmax**2/diff_eff/3600.0
      end if
      

! Butler-Volmer Equation and its partial derivatives
! ________ ELECTROCHEMICAL REACTION RATE [A/cm2]


      ! It is important to benchmark the numerical derivatives and have a good way to determine
      ! the appropriate step size that will produce the most accurate results 
      ! In addition this method of multiplying a certain fraction does not work well when the variable
      ! is close to zero
      ! Use a conditional statement when close to zero to use additive steps
      ! if ( ABS(variable).LE.1.0d-10 ) then
      !     step  = 1.0d-10
      !     df_dv = (f(v+step) - f(v-step))/(2*step)
      ! end if

      ! for the concentrations, they are in log() terms and log(-#) is undefined
      ! for these variables it might be more useful to instead do something like:
      ! if ( variable.LE.1.0d-10 ) then
      !     step  = 1.0d-10
      !     df_dv = (f(2*step) - f(step))/(step)
      ! end if

! it is often necessary to have "write" statements inside of fillmat when debugging the code      

 19    format(2(A5,1X),  2(A12,1X), 12(A12,1X))
 20    format(1(I5,1X), 1(F5.2,1X),2(ES12.5,1X),12(ES12.5,1X))

!______________________Boundary-condition at x = 0
! At the boundaries it is usually useful to set the volume to zero (delX(j) = 0)
!     dE(i,k)*dcdxE + rj(i,k)*c = g(i)

      if (j.eq.1) then 

        alphaE=delx(j)/(delx(j+1)+delx(j))
        betaE=2.0/(delx(j)+delx(j+1))

        do ic=1,N
          cE(ic)    = alphaE*cprev(ic,j+1) + (1.d0 - alphaE)*cprev(ic,j)
          dcdxE(ic) = betaE * (cprev(ic,j+1) - cprev(ic,j))
        end do

! *************** Electrolyte Concentration ****************
! Both Fixed Concentration and No Flux work

!  No Flux Condition


! ### ROGER ### !
        dE(1,1) = 1.0
        smG(1)  = +i_applied/Fconst/diff_eff + dE(1,1)*dcdxE(1)

        return
      end if


!______________________Boundary-condition at x=xmax
! At the boundaries it is usually useful to set the volume to zero (delX(j) = 0)
! The general boundary condition equation looks like: dW(i,k)*dcdxW + rj(i,k)*delC = g(i)

      if (j.eq.NJ) then

        alphaW=delx(j-1)/(delx(j-1)+delx(j))
        betaW=2.0/(delx(j-1)+delx(j))

        do ic=1,N
          cW(ic)    = alphaW*cprev(ic,j) + (1.d0 - alphaW)*cprev(ic,j-1)
          dcdxW(ic) = betaW * (cprev(ic,j) - cprev(ic,j-1))
        end do

! *************** Electrolyte Concentration ****************
!  Fixed Concentration Condition
!     C = cbulk at x = xmax

        rj(1,1) = 1.0
        smG(1)  = cbulk*1.0 - rj(1,1)*cprev(1,j)

        return
      end if


! !______________________Governing Equations
! Finite Difference Method (Control Volume Formulation)
! A material-balance for species i at node k

!       dc(i,k)/dt * delX(k) = Flux_west(i) - Flux_east(i) + R(i)*delX(k)

! where R(i) is the production rate per unit volume, Flux_west(i) is the flux of species i
! at the "west" face of node k, and Flux_east(i) is the flux of species i at the "east" face 
! of node k.

! In general the flux of species i can be affected by all of the variable and their gradients
      
!       Flux_west(i) = d_west(j,k) * dcdx_west(j) + f_west*c_west(j)
!       Flux_east(i) = d_east(j,k) * dcdx_east(j) + f_east*c_east(j)

! The general formulat to apply to each node point can be written as:

!      g(i) = rj(i,k)*c + fW(i,k)*cW + dW(i,k)*dcdxW - fE(i,k)*cE - dE(i,k)*dcdxE

! The optimal interpolation formula to use to evaluate the variable c and their derivatives
! depend on the local Peclet number. The general form is can be seen below
! For a low Peclet number, a linear interpolation is most appropriate
! Linear interpolation, alpha = 0.5, beta = 1/delX

! *** NOTE *** For large Peclet numbers this central-difference approximation causes an oscillatory
! behavior in the concentration. This can be eliminated by an upwind scheme, in which betaW and betaE
! remain the same, but 
! alphaE = 0, alphaW = 0 when flow is from west to east
! alphaE = 1, alphaW = 1 when flow is from east to west

        alphaW=delx(j-1)/(delx(j-1)+delx(j))
        alphaE=delx(j)/(delx(j+1)+delx(j))
        betaW=2.0/(delx(j-1)+delx(j))
        betaE=2.0/(delx(j)+delx(j+1))

        do ic=1,N
          cW(ic)    = alphaW*cprev(ic,j) + (1.d0 - alphaW)*cprev(ic,j-1)
          cE(ic)    = alphaE*cprev(ic,j+1) + (1.d0 - alphaE)*cprev(ic,j)
          dcdxW(ic) = betaW * (cprev(ic,j) - cprev(ic,j-1))
          dcdxE(ic) = betaE * (cprev(ic,j+1) - cprev(ic,j))
        end do

        

        ! -D*dcdxW - (-D*dcdxE) = delX * delC/delt
        ! *** LINEARIZED *** !
          ! (-D*dcdxW - D*dcdxW_prev) - (-D*dcdxE - D*dcdxE_prev) = delX * delC/delt
        !!*******************!
        dW(1,1) = -diff_eff
        dE(1,1) = -diff_eff
        fW(1,1) = 0.d0
        fE(1,1) = 0.d0

        rj(1,1) = +delx(j) * (1.0/delT)

        smG(1)  =  0.0 &
        &         + (fW(1,1)*cW(1) + dW(1,1)*dcdxW(1)) &
        &         - (fE(1,1)*cE(1) + dE(1,1)*dcdxE(1))

        return



end subroutine fillmat




! Below this point in code are the fundamental subroutines. Do not edit anything below.
!************************************************************************************
!************************************************************************************
!************************************************************************************


!************************************ABDGXY******************************************

subroutine ABDGXY(h,j)
      use user_input
      use variables
      implicit double precision(a-h,o-z)

      if(j.eq.1) then
          do 1 ii=1,n
          do 10 kk=1,n
             X(ii,kk)=0.d0
             B(ii,kk)=rj(ii,kk) - (1.d0 - alphaE)*fE(ii,kk) + betaE*dE(ii,kk)
             D(ii,kk)= -alphaE*fE(ii,kk) - betaE*dE(ii,kk)

10         continue
             G(ii)=smG(ii)
1        continue
          return
      end if
      if (j.eq.NJ) then 
          do 2 ii=1,n
          do 20 kk=1,n
             Y(ii,kk)=0.d0
             A(ii,kk)=(1.d0 - alphaW)*fW(ii,kk) - betaW*dW(ii,kk)
             B(ii,kk)=rj(ii,kk) + betaW*dW(ii,kk) + alphaW*fW(ii,kk)
20         continue
             G(ii)=smG(ii)
2         continue
          return
      end if
      do 3 ii=1,n
      do 30 kk=1,n
             A(ii,kk)=(1.d0 - alphaW)*fW(ii,kk) - betaW*dW(ii,kk)
             B(ii,kk)=rj(ii,kk) + betaW*dW(ii,kk) + alphaW*fW(ii,kk) &
             &        - (1.d0 - alphaE)*fE(ii,kk) + betaE*dE(ii,kk)
             D(ii,kk)= -alphaE*fE(ii,kk) - betaE*dE(ii,kk)
30     continue
             G(ii)=smG(ii)
3     continue
      return
end subroutine ABDGXY

!***********************************MATINV*****************************************

SUBROUTINE MATINV(N,M,DETERM)     
use variables, only: A,B,delC,D,ID
implicit double precision (A-H,O-Z)    

      DETERM=1.0
      DO 1 I=1,N
1     ID(I)=0
      DO 18 NN=1,N
      BMAX=1.1
      DO 6 I=1,N
      IF (ID(I).NE.0) GOTO 6
      BNEXT=0.0
      BTRY=0.0
      DO 5 J=1,N
      IF (ID(J).NE.0) GOTO 5
      IF (DABS(B(I,J)).LE.BNEXT) GOTO 5
      BNEXT=DABS(B(I,J))
      IF (BNEXT.LE.BTRY) GOTO 5
      BNEXT=BTRY
      BTRY=DABS(B(I,J))
      JC=J
5     CONTINUE 
      IF (BNEXT.GE.BMAX*BTRY) GOTO 6
      BMAX=BNEXT/BTRY
      IROW=I
      JCOL=JC 
6     CONTINUE
      IF (ID(JC).EQ.0) GOTO 8
      DETERM=0.0
      RETURN
8     ID(JCOL)=1
      IF (JCOL.EQ.IROW) GOTO 12
9     DO 10 J=1,N
      SAVE=B(IROW,J)
      B(IROW,J)=B(JCOL,J)
10    B(JCOL,J)=SAVE
      DO 11 K=1,M
      SAVE=D(IROW,K)
      D(IROW,K)=D(JCOL,K)
11    D(JCOL,K)=SAVE
12    F=1.0/B(JCOL,JCOL)
      DO 13 J=1,N
13    B(JCOL,J)=B(JCOL,J)*F
      DO 14 K=1,M
14    D(JCOL,K)=D(JCOL,K)*F
      DO 18 I=1,N
      IF (I.EQ.JCOL) GOTO 18
      F=B(I,JCOL)
      DO 16 J=1,N
16    B(I,J)=B(I,J)-F*B(JCOL,J)
      DO 17 K=1,M
17    D(I,K)=D(I,K)-F*D(JCOL,K)
18    CONTINUE
      RETURN
      END

!*************************************BAND******************************************

SUBROUTINE BAND(J)
use variables, only: A,B,delC,D,G,X,Y,NP1,E
use user_input, only: N,NJ
implicit double precision (A-H,O-Z) 


101   FORMAT(15H DETERM=0 AT J=,I4)
      IF (J-2) 1,6,8
1     NP1=N+1
      DO 2 I=1,N
      D(I,2*N+1)=G(I)
      DO 2 L=1,N
      LPN=L+N
2     D(I,LPN)=X(I,L)
      CALL MATINV(N,2*N+1,DETERM)
      IF (DETERM) 4,3,4
3     PRINT 101,J
4     DO 5 K=1,N
      E(K,NP1,1)=D(K,2*N+1)
      DO 5 L=1,N
      E(K,L,1)=-D(K,L)
      LPN=L+N
5     X(K,L)=-D(K,LPN)
      RETURN
6     DO 7 I=1,N
      DO 7 K=1,N
      DO 7 L=1,N
7     D(I,K)=D(I,K)+A(I,L)*X(L,K)
8     IF (J-NJ) 11,9,9
9     DO 10 I=1,N
      DO 10 L=1,N
      G(I)=G(I)-Y(I,L)*E(L,NP1,J-2)
      DO 10 M=1,N
10    A(I,L)=A(I,L) + Y(I,M)*E(M,L,J-2)
11    DO 12 I=1,N
      D(I,NP1)=-G(I)
      DO 12 L=1,N
      D(I,NP1)=D(I,NP1)+A(I,L)*E(L,NP1,J-1)
      DO 12 K=1,N
12    B(I,K)=B(I,K) + A(I,L)*E(L,K,J-1)
      CALL MATINV(N,NP1,DETERM)
      IF (DETERM) 14,13,14
13    PRINT 101,J
14    DO 15 K=1,N
      DO 15 M=1,NP1
15    E(K,M,J)=-D(K,M)
      IF (J-NJ) 20,16,16
16    DO 17 K=1,N
17    delC(K,J)=E(K,NP1,J)
      DO 18 JJ=2,NJ
      M=NJ-JJ+1
      DO 18 K=1,N
      delc(K,M)=E(K,NP1,M)
      DO 18 L=1,N
18    delC(K,M)=delC(K,M) +E(K,L,M)*delC(L,M+1)  
      DO 19 L=1,N
      DO 19 K=1,N
19    delC(K,1)=delC(K,1)+X(K,L)*delC(L,3)
20    RETURN
      END
