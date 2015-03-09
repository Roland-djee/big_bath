module generate
  use types
  use constants
  implicit none

contains

subroutine generate_impurities
  use types
  implicit none
  ! local variables
  integer(kind=8) :: i
  integer(kind=8), allocatable :: imp_temp(:), index(:)
  double precision :: rand2(1)
  double precision, allocatable :: rand1(:)
  logical, allocatable :: mask(:)
  
  ! define the number of impurities
  nb_imp = nint(lattice%imp*nb_sites/1.D2)
  
  ! allocate the impurities , random 
  allocate (imp(nb_imp))
  allocate (rand1(nb_imp))
  allocate (mask(nb_imp))

  ! define randomly nb_imp sites as impurities among the nb_sites
  ! excluding 0
  call random_seed()
  call random_number(rand1)
  imp = nint( rand1 * (nb_sites - 1) + 1)
    
  ! if so, remove the origin and replace by another random site
  call random_number(rand2)
  where (imp == origin) 
     imp = nint ( rand2 * (nb_sites - 1) + 1)
  end where

  ! check and remove duplicates if any
  mask = .true.
  do i=nb_imp,2,-1
     mask(i) = .not.(any(imp(:i-1) == imp(i)))
  end do
  allocate (index, source=pack([(i,i=1,nb_imp)], mask))
  allocate (imp_temp, source=imp(index))
  
  ! reset parameters
  nb_imp = size(imp_temp)
  print*,'nb imp', nb_imp
  deallocate(imp)
  allocate(imp(nb_imp))
  imp = imp_temp
  
  ! allocate coordinates arrays
  allocate (Imp_x(nb_imp),Imp_y(nb_imp),Imp_z(nb_imp))
  
  ! impurities coordinates
  Imp_x = Sites_x(imp)
  Imp_y = Sites_y(imp)
  Imp_z = Sites_z(imp)

  ! Arbitrary set the first impurity as the nuclear spin qubit for J=1Mrad

  Imp_x(1) = nuclear%x
  Imp_y(1) = nuclear%y
  Imp_z(1) = nuclear%z
  
  deallocate (rand1)
  deallocate (mask)
  deallocate (imp_temp)
  deallocate (index)
  deallocate (Sites_x,Sites_y,Sites_z)
  
end subroutine generate_impurities

subroutine generate_crystal
  implicit none
  ! local variables
  integer(kind=8) :: i, j, k, l, m, N0, N1
  integer(kind=8) :: volume
  integer, allocatable :: tmp_x(:), tmp_y(:), tmp_z(:)
  
  integer(kind=8) :: x, y, z
  double precision :: r, ph, t
  
  ! define the lattice volume and the number of sites
  volume    = (2*lattice%range)**3
  if (lattice%range == 0) volume = 1
  nb_sites = nb_basisvec*volume
  
  ! allocate the sites coordinates array
  !allocate (lattice_sites(nb_sites))
  allocate (Sites_x(nb_sites))
  allocate (Sites_y(nb_sites))
  allocate (Sites_z(nb_sites))
  
  ! generate the crystal sites coordinates
  m = 0
  N0 = lattice%range
  N1 = lattice%range
  if (lattice%range == 0) then
     N0 = 0
     N1 = 1
  end if
  do i=-N0,N1-1
     do j=-N0,N1-1
        do k=-N0,N1-1
           !   print*,'i,j,k',i,j,k
           do l=1,nb_basisvec
              x = basis_vectors(l)%x + i*lattice%modulo
              y = basis_vectors(l)%y + j*lattice%modulo
              z = basis_vectors(l)%z + k*lattice%modulo
              
              r  = scaling * sqrt(dble(x)**2 + dble(y)**2 + dble(z)**2)
              ph = datan(abs(dble(y)/dble(x)))
              if (dble(y) .gt. 0.d0 .and. dble(x) .le. 0.d0) then
                 ph = pi - ph
              else if (dble(y) .le. 0.d0 .and. dble(x) .le. 0.d0) then
                 ph = pi + ph
              else if (dble(y) .le. 0.d0 .and. dble(x) .gt. 0.d0) then
                 ph = 2.d0*pi- ph
              end if
              t  = acos(dble(z) / r)
              t  = pi/2.d0 - t
              
              !print*,'r,ph,t',r,ph,t
              
              if (sphere)then
                 if ((r .le. R1) .and. (r .gt. R0) &
                      & .and. (ph .ge. 0.d0) .and. (ph .le. phi) &
                      & .and. ((pi/2.d0 - t) .lt. theta)) then
                    m = m + 1
                    Sites_x(m) = x
                    Sites_y(m) = y
                    Sites_z(m) = z
                    !lattice_sites(m)%x = basis_vectors(l)%x + i*lattice%modulo
                    !lattice_sites(m)%y = basis_vectors(l)%y + j*lattice%modulo
                    !lattice_sites(m)%z = basis_vectors(l)%z + k*lattice%modulo
                    ! Select the origin indice
                    if (Sites_x(m) == 0 .and. &
                         Sites_y(m) == 0 .and. &
                         Sites_z(m) == 0) origin = m
                    
                 end if
              else
                 m = m + 1
                    Sites_x(m) = x
                    Sites_y(m) = y
                    Sites_z(m) = z
                    !lattice_sites(m)%x = basis_vectors(l)%x + i*lattice%modulo
                    !lattice_sites(m)%y = basis_vectors(l)%y + j*lattice%modulo
                    !lattice_sites(m)%z = basis_vectors(l)%z + k*lattice%modulo
                    ! Select the origin indice
                    if (Sites_x(m) == 0 .and. &
                         Sites_y(m) == 0 .and. &
                         Sites_z(m) == 0) origin = m

              end if
           end do
           !10 continue
        end do
     end do
  end do
  
  nb_sites = m 
  
  print*,'nb sites',nb_sites
  
  allocate(tmp_x(m), tmp_y(m), tmp_z(m))
  
  tmp_x = Sites_x(1:m)
  tmp_y = Sites_y(1:m)
  tmp_z = Sites_z(1:m)
  
  deallocate(Sites_x,Sites_y,Sites_z)
  allocate(Sites_x(m),Sites_y(m),Sites_z(m))
  
  Sites_x = tmp_x
  Sites_y = tmp_y
  Sites_z = tmp_z

  deallocate(tmp_x, tmp_y, tmp_z)
  deallocate(basis_vectors)
  
end subroutine generate_crystal

end module generate
