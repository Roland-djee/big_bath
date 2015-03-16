!------------------------------------------------------------------------------
! interactions module
!------------------------------------------------------------------------------
!
! MODULE: interactions
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> This module contains the subroutines needed to compute the interactions 
!> (hyperfine, dipolar) between the spin systems considered.
!
! REVISION HISTORY:
! 11-03-2015 - Initial Version
!------------------------------------------------------------------------------

module interactions
  use types 
  use constants
  use write
  implicit none

contains

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Selects the pairs for nuclear decoherence.
  !> @brief
  !> Selects among the total number of pairs the flip-flopping ones 
  !> likely to interact with the nucleus.
  !
  ! REVISION HISTORY:
  ! TODO_11_03_2015 - Break the subroutine into smaller ones - TODO_select_nuc
  !
  !> @param[in]  C12, HYPER
  !> @param[out] --      
  !> @return     C12, DC, DJ
  !--------------------------------------------------------------------------- 

  subroutine select_nuc
    implicit none
    double precision :: CA1(nb_imp - 1)
    ! Local variables
    integer(kind=8) :: loc, k, l, m, n
    integer(kind=8) :: Np
    double precision, allocatable :: C12_tmp(:)
    double precision, allocatable :: DC_tmp(:)
    double precision, allocatable :: DJ_tmp(:)
    
    logical, allocatable :: mask(:)
  
    Np = (nb_imp - 1) * (nb_imp - 2) / 2
    allocate(C12_tmp(Np), DC_tmp(Np), DJ_tmp(Np))
    allocate(DC_n(Np))
    
    ! locate the nearest impurity to the input J value
    ! loc  = minloc(abs(HYPER - Jnuc), 1)
    ! Jnuc = HYPER(loc)
    !> locate the nucleus 
    loc = 1
    Jnuc = HYPER(loc)
    
    !> select corresponding C12 values between the nucleus and
    !> all other impurities
    m = 0
    n = 0
    do k=1,nb_imp - 1
       do l=k + 1,nb_imp
          m = m + 1
          if ((k .ne. loc) .and. (l == loc)) then
             n = n + 1
             CA1(n) = C12(m)
             C12(m) = 0.d0
          else if (k == loc) then
             n = n + 1
             CA1(n) = C12(m)
             C12(m) = 0.d0
          end if
       end do
    end do
  
    !> creating a temporary array of C12 differences between the nucleus and
    !> all impurity pairs
    m = 0
    do k=1,nb_imp - 2
       do l=k + 1,nb_imp - 1
          m = m + 1
          DC_tmp(m) = CA1(k)   - CA1(l)
       end do
    end do

    !> creating an array of J differences for the impurity pairs
    m = 0
    do k=2,nb_imp - 1
       do l=k + 1,nb_imp
          m = m + 1
          DJ_tmp(m) = HYPER(k) - HYPER(l)
          DC_n(m)   = C12_n(k) - C12_n(l)
       end do
    end do

    !> creating a temporary C12 array excluding previous CA1 values
    C12_tmp = pack (C12, C12 .ne. 0.d0)
    
    deallocate (C12)
    allocate (C12(Np), DJ(Np), DC(Np))
    allocate (mask(Np))
    
    !> reset C12, J arrays and nb_pairs for IFF
    C12      = C12_tmp
    DJ       = DJ_tmp
    DC       = DC_tmp
    nb_pairs = Np 

    write(*,10),'Number of nuclear pairs :',nb_pairs
10  format( (a, i) )

    !> selecting the pairs above a defined C12 threshold if so
    if (threshold) then
       mask = .false.
       where (abs(C12) .gt. threshold_val * 2.d0 * pi) mask = .true.  
       
       nb_pairs = count(mask)

       write(*,20),'Number of selected pairs:',nb_pairs
20     format( (a, i) )

       deallocate(C12_tmp, DC_tmp, DJ_tmp)
       allocate(C12_tmp(nb_pairs), DC_tmp(nb_pairs), DJ_tmp(nb_pairs))
    
       C12_tmp = pack(C12, mask)
       DC_tmp  = pack(DC, mask)
       DJ_tmp  = pack(DJ, mask)

       deallocate(C12,DC,DJ)
       allocate(C12(nb_pairs),DC(nb_pairs),DJ(nb_pairs))

       C12 = C12_tmp
       DC  = DC_tmp
       DJ  = DJ_tmp
    end if

    deallocate(C12_tmp, DJ_tmp, DC_tmp)
  
  end subroutine select_nuc

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Ask for the computation of couplings.
  !> @brief
  !> Directs to the computing subroutines for each interaction type.
  !
  ! REVISION HISTORY:
  ! TODO_11_03_2015 - - TODO_couplings
  !
  !> @param[in]  --
  !> @param[out] --      
  !> @return     --
  !--------------------------------------------------------------------------- 

  subroutine couplings
    implicit none
    
    !> define the number of impurity pairs
    nb_pairs = nb_imp * (nb_imp - 1) / 2

    write(*,10),'Total number of pairs   :',nb_pairs
10  format( (a, i) )
   
    !> calculate hyperfine coupling values between the electron
    !> and each impurity
    call hyperfine
     
    !> calculate dipolar coupling values between members of all pairs
    call dipolar

    !> calculate dipolar coupling values between the nuclei
    !> and each impurity
    if (e_n_coupling) call dipolar_n
      
  end subroutine couplings

  subroutine dipolar_n
    implicit none    
    ! Local variables and arrays
    double precision :: Orient_norm
    integer(kind=8) :: proj(nb_imp)
    integer(kind=8) :: X_dist(nb_imp)
    integer(kind=8) :: Y_dist(nb_imp)
    integer(kind=8) :: Z_dist(nb_imp)
    double precision :: dist(nb_imp)
    double precision :: cos_t12(nb_imp)

    ! allocations
    allocate (C12_n(nb_imp))

    ! Orientation norm
    Orient_norm = dsqrt(dble(Orientation%x**2) + &
                        dble(Orientation%y**2) + &
                        dble(Orientation%z**2))
    
    ! coordinate difference between every impurity spin pairs
    ! avoiding double counting
 
    X_dist = Imp_x
    Y_dist = Imp_y
    Z_dist = Imp_z
   
    ! distance between every spin partners of impurity pairs
    dist    = dsqrt(dble(X_dist**2) + dble(Y_dist**2) + dble(Z_dist**2))
    ! Projection
    proj    = Orientation%x*X_dist + Orientation%y*Y_dist + Orientation%z*Z_dist
    ! cos(theta_12)
    cos_t12 = proj / (Orient_norm * dist)
    
    ! array of C12 values for each impurity pair
    C12_n = pref2 * (1.d0 - 3.d0 * cos_t12**2) / ((scaling * dist)**3)
  
  end subroutine dipolar_n

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Calculates the dipolar couplings.
  !> @brief
  !> Calculate the dipolar value between the members of pairs formed by all 
  !> impurities.
  !
  ! REVISION HISTORY:
  ! TODO_12_03_2015 -  - TODO_F
  !
  !> @param[in]  Imp_x,Imp_y,Imp_z
  !> @param[out] --      
  !> @return     C12
  !---------------------------------------------------------------------------

  subroutine dipolar
    implicit none    
    ! Local variables and arrays
    integer(kind=8) :: i, j, k
    double precision :: Orient_norm
    integer(kind=8) :: proj(nb_pairs)
    integer(kind=8) :: X_dist(nb_pairs)
    integer(kind=8) :: Y_dist(nb_pairs)
    integer(kind=8) :: Z_dist(nb_pairs)
    double precision :: dist(nb_pairs)
    double precision :: cos_t12(nb_pairs)

    !> allocation
    allocate (C12(nb_pairs))

    !> B field norm
    Orient_norm = dsqrt(dble(Orientation%x**2) + &
                        dble(Orientation%y**2) + &
                        dble(Orientation%z**2))
    
    !> coordinate difference between every impurity spin pairs
    !> avoiding double counting
    k = 0
    do i=1,nb_imp - 1
       do j=i + 1,nb_imp
          k = k + 1
          X_dist(k) = Imp_x(j) - Imp_x(i)
          Y_dist(k) = Imp_y(j) - Imp_y(i)
          Z_dist(k) = Imp_z(j) - Imp_z(i)
       end do
    end do
  
    !> distance between every spin partners of impurity pairs
    dist    = dsqrt(dble(X_dist**2) + dble(Y_dist**2) + dble(Z_dist**2))
    !> Projection
    proj    = Orientation%x*X_dist + Orientation%y*Y_dist + Orientation%z*Z_dist
    !> cos(theta_12)
    cos_t12 = proj / (Orient_norm * dist)
    
    !> array of C12 values for each impurity pair
    C12 = pref * (1.d0 - 3.d0 * cos_t12**2) / ((scaling * dist)**3)
  
  end subroutine dipolar

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Calculates the hyperfine values.
  !> @brief
  !> Calculate the hyperfine value between the electron and all impurities.
  !
  ! REVISION HISTORY:
  ! TODO_12_03_2015 -  - TODO_F
  !
  !> @param[in]  Imp_x,Imp_y,Imp_z
  !> @param[out] --      
  !> @return     HYPER
  !--------------------------------------------------------------------------- 

  subroutine hyperfine
    implicit none
    ! Local variables and arrays
    double precision :: X_imp(nb_imp), Y_imp(nb_imp), Z_imp(nb_imp)
    double precision :: term1(nb_imp), term2(nb_imp), term3(nb_imp)
 
    allocate (HYPER(nb_imp))

    !> scale impurities coordinates in angstroms
    X_imp = scaling * Imp_x
    Y_imp = scaling * Imp_y
    Z_imp = scaling * Imp_z
   
    !> compute the terms of the KLW
    term1 = F(X_imp, Y_imp, Z_imp) * dcos(k0 * X_imp)
    term2 = F(Y_imp, Z_imp, X_imp) * dcos(k0 * Y_imp)
    term3 = F(Z_imp, X_imp, Y_imp) * dcos(k0 * Z_imp)

    !> array of hyperfine values
    HYPER = term1 + term2 + term3
    HYPER = p * HYPER**2
               
  end subroutine hyperfine

  !---------------------------------------------------------------------------  
  !> @author 
  !> Dr. Roland Guichard University College London
  !
  ! DESCRIPTION: 
  !> Calculates the F part of the hyperfine Kohn-Luttinger wavefunction.
  !> @brief
  !> Calculate the F part of the hyperfine value between the electron and 
  !> all impurities.
  !
  ! REVISION HISTORY:
  ! TODO_12_03_2015 -  - TODO_F
  !
  !> @param[in]  x,y,z
  !> @param[out] --      
  !> @return     F
  !--------------------------------------------------------------------------- 
 
  function F(x, y, z)
    implicit none
    double precision :: F(nb_imp)
    double precision, intent(in) :: x(nb_imp), y(nb_imp), z(nb_imp)
    ! Local array
    double precision :: arg(nb_imp)

    arg = dsqrt((x / (n * b))**2 + (y**2 + z**2) / ((n * a)**2))
    F   = dexp(-arg) / dsqrt(pi * (n * a)**2 * n * b)

  end function F

end module interactions
