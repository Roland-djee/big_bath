module test_generate
  use types
  use generate
  implicit none

contains

  subroutine unit_test_generate_impurities
    logical :: test

    call generate_impurities

  end subroutine unit_test_generate_impurities

  subroutine unit_test_generate_crystal
    integer :: i,j,k,l,m
    logical :: test

    lattice%range  = 5
    lattice%modulo = 4 
    nb_basisvec    = 3

    basis_vectors(1)%x = 1 
    basis_vectors(1)%y = 2
    basis_vectors(1)%z = 3

    basis_vectors(2)%x = 3
    basis_vectors(2)%y = 2
    basis_vectors(2)%z = 1

    basis_vectors(3)%x = 2
    basis_vectors(3)%y = 1
    basis_vectors(3)%z = 3

    sphere = .false.

    call generate_crystal

    test = .true.
    
    allocate(basis_vectors(3))
    basis_vectors(1)%x = 1 
    basis_vectors(1)%y = 2
    basis_vectors(1)%z = 3

    basis_vectors(2)%x = 3
    basis_vectors(2)%y = 2
    basis_vectors(2)%z = 1

    basis_vectors(3)%x = 2
    basis_vectors(3)%y = 1
    basis_vectors(3)%z = 3

    if (nb_sites == 3000) then
       write(*,*)'nb_sites...ok'
    else
       test = .false.
    end if

    m = 0
    do i=-5,4
       do j=-5,4
          do k=-5,4
             do l=1,3
                m = m + 1
                if (Sites_x(m) == basis_vectors(l)%x + i*4) then
                else
                   test = .false.
                end if
                if (Sites_y(m) == basis_vectors(l)%y + j*4) then
                else
                   test = .false.
                end if
                if (Sites_z(m) == basis_vectors(l)%z + k*4) then
                else
                   test = .false.
                end if
             end do
          end do
       end do
    end do

    if (test) then
       write(*,*)'Sites...ok'
       write(*,*)'Test generate successful'
    else
       write(*,*)'Test generate failed'
    end if
    
  end subroutine unit_test_generate_crystal

end module test_generate
