


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

! Space curvature
  real(8) :: k_curvature

! Time and related
  real(8) :: t, dt
  integer :: Nt

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

! Parse parameters

  write(*,*) "Time steps:"
  read(*,*) Nt
  write(*,*) "(initial) time step"
  read(*,*) dt
  write(*,*) "Space-time curvature (0,1,-1)"
  read(*,*) k_curvature
  write(*,*) "Cosmological constant"
  read(*,*) cosmo_lambda
  write(*,*) "Initial scale_factor"
  read(*,*) aux
  phi = log(aux)
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

  trS = - cosmo_lambda/(8.d0*smallpi)

  rho_k1 =  cosmo_lambda/(8.d0*smallpi)
  trS_K1 = -cosmo_lambda/(8.d0*smallpi)

! NOTE: alpha and matter sources are kept constant for the time being
! and because of this I calculate the midpoint values outside main loop just once
! For other matter choices they must be updated


! The extrinsic curvature must be set from the hamiltonian constraint
! once initial rho has been calculated

  aux = 24.d0*smallpi*rho - 9.d0*k_curvature*exp(-4.d0*phi) 
  trK = -exp_sign*sqrt(aux)

! ***********************
!    Open outputfiles
! ***********************

  call system('mkdir -p '//trim(out_dir))

  open(unit=88,file=trim(out_dir)//"/phi.tl",status="unknown")
  open(unit=89,file=trim(out_dir)//"/trK.tl",status="unknown")

! Write starting values

  write(88,"(2ES16.8)") t, phi
  write(88,"(2ES16.8)") t, trK

! ***********************
!    Integrate
! ***********************

  do i=1,Nt
     !Runge-Kutta2
     ! First evaluation, midpoint estimation
     phi_k1 = phi + half*dt*(-half*third*alpha*trK) 
     trK_k1 = trK + half*dt*(4.d0*smallpi*(rho+trS))

     ! Lapse for the moment is constant
     ! alpha_k1 = whatever...

     ! Matter sources should be also updated at midpoint
     ! cosmological constant is constant and boring
     !rho_k1 =  cosmo_lambda/(8.d0*smallpi)
     !trS_K1 = -cosmo_lambda/(8.d0*smallpi)

     ! Final evaluation: Full step but with sources calculated at midpoints
     phi = phi + dt*(-half*third*alpha_k1*trK_k1)
     trK_k1 = trK + half*dt*(4.d0*smallpi*(rho_k1+trS_k1))

     ! alpha = whatever
     ! rho = 
     ! trS = 

     t = t + dt

     write(88,"(2ES16.8)") t, phi
     write(89,"(2ES16.8)") t, trK


  end do

  close(88)
  close(89)

  print*, 
  print*, "Finished! Have a nice day!"

end program main
