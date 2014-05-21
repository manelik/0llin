
!
! code 0llin
!
! 0+1 code that solves einstein (friedmann) equations in a simplified 3+1 scenario
! when the spacetime is homogeneous and isotropic
!

program main
! Main routine of program 0llin (that is Ollin written with a zero replacing the "O")

  implicit none

! Variables declaration

! geometric quantities (conformal exponent, trace of extrinsic curvature, lapse
  real(8) :: phi,trK,alpha
  real(8) :: phi_k1,trK_k1,alpha_k1
! Constraint
  real(8) :: ham
! Matter
  real(8) :: rho,trS
  real(8) :: rho_k1,trS_k1

! Cosmological constant
  real(8) :: cosmo_lambda

! Dust array
  real(8) :: dust_rho0,dust_rho0_k1
  logical :: has_dust
!  integer :: Nfluid

! Space curvature
  real(8) :: k_curvature

! Time and related
  real(8) :: t, dt
  integer :: Nt, ft

! An auxiliar real that can be reused many times in the code

  real(8) :: aux

! Sign of the expansion
  integer :: exp_sign


! Output directory

  character*30 :: out_dir

! counters
  
  integer :: i,j,k


! Numbers
  real(8) :: zero,one,two,half,third,smallpi

  zero = 0.d0
  one  = 1.d0
  two  = 2.d0
  half = .5d0
  third= one/3.d0

  smallpi =  acos(-1.d0)

! Greeting

  write(*,*)"******************************************"
  write(*,*)"*   _______  .__  .__  .__               *"
  write(*,*)"*   \   _  \ |  | |  | |__| ____         *"
  write(*,*)"*   /  /_\  \|  | |  | |  |/    \        *"
  write(*,*)"*   \  \_/   \  |_|  |_|  |   |  \       *"
  write(*,*)"*    \_____  /____/____/__|___|  /       *"
  write(*,*)"*          \/                  \/        *"
  write(*,*)"*                                        *"
  write(*,*)"*   A simple 0+1 relativistic code       *"
  write(*,*)"*     First version: Yesterday           *"
  write(*,*)"*     Last version:  Today               *"
  write(*,*)"*                                        *"
  write(*,*)"*     Author:       Jose Manuel          *"
  write(*,*)"*                                        *"
  write(*,*)"******************************************"



! Parse parameters

  write(*,*) "Time steps:"
  read(*,*) Nt
  write(*,*) "(initial) time step"
  read(*,*) dt
  write(*,*) "Screen output frequency"
  read(*,*) ft
  write(*,*) "Space-time curvature (0,1,-1)"
  read(*,*) k_curvature
  write(*,*) "Cosmological constant"
  read(*,*) cosmo_lambda
  write(*,*) "Initial scale_factor"
  read(*,*) aux
  phi = log(aux)

  write(*,*) "Is dust included (.false./.true.)"
  read(*,*) has_dust

  if(has_dust) then
     write(*,*) "Initial dust density"
     read(*,*) dust_rho0
  end if

  write(*,*) "Sign of expansion (+1 expanding, -1 contracting)"
  read(*,*) exp_sign 

  write(*,*) "Output directory name"
  read(*,*) out_dir
 
! Sanity check
  if(exp_sign/=1.and.exp_sign/=-1) then
     print*, "Wrong sign for expansion. can be just +1/-1"
     print*,
     print*, "You were advised"
     stop
  end if

! ******************
!    Initialize
! ******************

! Initialize time
  t= 0.d0
  alpha = 1.d0
  alpha_k1 = 1.d0
! Add matter contributions
! starting with cosmological constant  


  rho = cosmo_lambda/(8.d0*smallpi)

  trS = - 3.d0*cosmo_lambda/(8.d0*smallpi)

! dust
  rho = rho + dust_rho0
!  dust is pressureless, it does not add to trS
  

!  rho_k1 =  cosmo_lambda/(8.d0*smallpi)
!  trS_k1 = -3.d0*cosmo_lambda/(8.d0*smallpi)

! NOTE: alpha and matter sources are kept constant for the time being
! and because of this I calculate the midpoint values outside main loop just once
! For other matter choices they must be updated


! The extrinsic curvature must be set from the hamiltonian constraint
! once initial rho has been calculated

  aux = 24.d0*smallpi*rho - 9.d0*k_curvature*exp(-4.d0*phi) 
  trK = -exp_sign*sqrt(aux)

! Hamiltonian constraint, initially vanishes identically

  ham = two*third*trK**2 + 6.d0*k_curvature*exp(-4.d0*phi)-16.d0*smallpi*rho

! ***********************
!    Open outputfiles
! ***********************

  call system('mkdir -p '//trim(out_dir))

! Save parameters
! It is convenient to save a copy of the parameters that produced the simulation
  open(unit=88,file=trim(out_dir)//"/"//trim(out_dir)//".par",status="unknown")

  write(88,*) Nt, "#Time steps:"
  write(88,*) dt, "#(initial) time step"
  write(88,*) ft, "#Screen output frequency"
  write(88,*) k_curvature, "#Space-time curvature (0,1,-1)"
  write(88,*) cosmo_lambda,"#Cosmological constant"
  write(88,*) exp(phi), "#Initial scale_factor"
  write(88,*) has_dust, "#dust_included?"
  if(has_dust)   write(88,*) dust_rho0, "#dust_initial density"
  write(88,*) exp_sign, "#Sign of expansion (+1 expanding, -1 contracting)"
  write(88,*) out_dir, "#Output directory name"

  close(88)

! Open outfiles
  open(unit=88,file=trim(out_dir)//"/phi.tl",status="unknown")
  open(unit=89,file=trim(out_dir)//"/trK.tl",status="unknown")
  open(unit=90,file=trim(out_dir)//"/ham.tl",status="unknown")
  open(unit=91,file=trim(out_dir)//"/scale_factor.tl",status="unknown")
  open(unit=92,file=trim(out_dir)//"/Hubble_factor.tl",status="unknown")

  if(has_dust) open(unit=93,file=trim(out_dir)//"/dust_rho0.tl",status="unknown")

! Write starting values

  write(88,"(2ES16.8)") t, phi
  write(89,"(2ES16.8)") t, trK
  write(90,"(2ES16.8)") t, ham
  write(91,"(2ES16.8)") t, exp(2.d0*phi)
  write(92,"(2ES16.8)") t, -alpha*third*trK
  if(has_dust)  write(93,"(2ES16.8)") t, dust_rho0

! Screen output
  write(*,*)"+-------------------------------+"
  write(*,*)"|  Time        |  Hamiltonian   |"
  write(*,*)"+-------------------------------+"
  write(*,"(2ES16.8)") t, ham



! ***********************
!    Integrate
! ***********************

  do i=1,Nt
     !Runge-Kutta2
     ! First evaluation, midpoint estimation
     phi_k1 = phi + half*dt*(-half*third*alpha*trK) 
     trK_k1 = trK + half*dt*(alpha*(third*trK**2+ 4.d0*smallpi*(rho+trS)))

     if(has_dust) dust_rho0_k1 = dust_rho0 +half*dt*( dust_rho0*alpha*trK)

     ! Lapse for the moment is constant
     ! alpha_k1 = whatever...

     ! Matter sources should be also updated at midpoint
     ! cosmological constant is constant and boring
     rho_k1 =  cosmo_lambda/(8.d0*smallpi)
     trS_k1 = -3.d0*cosmo_lambda/(8.d0*smallpi)
     ! Dust
     if(has_dust) rho_k1 = rho_k1 + dust_rho0_k1

     ! Final evaluation: Full step but with sources calculated at midpoints
     phi = phi + dt*(-half*third*alpha_k1*trK_k1)
     trK = trK + dt*(alpha*(third*trK_k1**2 + 4.d0*smallpi*(rho_k1+trS_k1)))

     if(has_dust) dust_rho0 = dust_rho0 +dt*( dust_rho0_k1*alpha_k1*trK_k1)

     ! alpha = whatever
     ! rho = 
     ! trS = 

     t = t + dt

     ! update matter
     ! Cosmological constant
     rho =  cosmo_lambda/(8.d0*smallpi)
     trS = -3.d0*cosmo_lambda/(8.d0*smallpi)
     ! Dust
     if(has_dust) rho = rho + dust_rho0

     ! update hamiltonian constraint
     ham = two*third*trK**2 + 6.d0*k_curvature*exp(-4.d0*phi)-16.d0*smallpi*rho

     write(88,"(2ES16.8)") t, phi
     write(89,"(2ES16.8)") t, trK
     write(90,"(2ES16.8)") t, ham
     write(91,"(2ES16.8)") t, exp(2.d0*phi)
     write(92,"(2ES16.8)") t, -alpha*third*trK
     if(has_dust) write(93,"(2ES16.8)") t, dust_rho0

! Output to screen
     if(MOD(i,ft)==0) then
        write(*,"(2ES16.8)") t, ham
     end if
  end do

  close(88)
  close(89)
  close(90)
  close(91)
  close(92)
  if(has_dust) close(93)

  write(*,*)"+-------------------------------+"
  print*,
  print*, 
  print*, "Finished! Have a nice day!"
  print*,
  print*, 

end program main
