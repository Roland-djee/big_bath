module test_read
  use types
  use read
  implicit none

contains

  subroutine unit_test_read
    implicit none
    logical :: test

    call read_input

    test = .true.

    if (Orientation%x == 0) then
       write(*,*)'Orientation%x...ok'
    else
       test = .false.
    end if

    if (Orientation%y == 0) then
       write(*,*)'Orientation%y...ok'
    else
       test = .false.
    end if

    if (Orientation%z == 1) then
       write(*,*)'Orientation%z...ok'
    else
       test = .false.
    end if

    if (Systype == 'P') then
       write(*,*)'Systype...ok'
    else
       test = .false.
    end if

    if (nb_basisvec == 8) then
       write(*,*)'nb_basisvec...ok'
    else
       test = .false.
    end if

    if (basis_vectors(1)%x == 0 .and. &
        basis_vectors(1)%y == 0 .and. &
        basis_vectors(1)%z == 0) then
       write(*,*)'basis_vectors(1)...ok'
    else
       test = .false.
    end if

    if (basis_vectors(2)%x == 0 .and. &
        basis_vectors(2)%y == 2 .and. &
        basis_vectors(2)%z == 2) then
       write(*,*)'basis_vectors(2)...ok'
    else
       test = .false.
    end if

    if (basis_vectors(3)%x == 2 .and. &
        basis_vectors(3)%y == 0 .and. &
        basis_vectors(3)%z == 2) then
       write(*,*)'basis_vectors(3)...ok'
    else
       test = .false.
    end if

    if (basis_vectors(4)%x == 2 .and. &
        basis_vectors(4)%y == 2 .and. &
        basis_vectors(4)%z == 0) then
       write(*,*)'basis_vectors(4)...ok'
    else
       test = .false.
    end if

    if (basis_vectors(5)%x == 3 .and. &
        basis_vectors(5)%y == 3 .and. &
        basis_vectors(5)%z == 3) then
       write(*,*)'basis_vectors(5)...ok'
    else
       test = .false.
    end if

    if (basis_vectors(6)%x == 3 .and. &
        basis_vectors(6)%y == 1 .and. &
        basis_vectors(6)%z == 1) then
       write(*,*)'basis_vectors(6)...ok'
    else
       test = .false.
    end if

    if (basis_vectors(7)%x == 1 .and. &
        basis_vectors(7)%y == 3 .and. &
        basis_vectors(7)%z == 1) then
       write(*,*)'basis_vectors(7)...ok'
    else
       test = .false.
    end if

    if (basis_vectors(8)%x == 1 .and. &
        basis_vectors(8)%y == 1 .and. &
        basis_vectors(8)%z == 3) then
       write(*,*)'basis_vectors(8)...ok'
    else
       test = .false.
    end if

    if (lattice%modulo == 4) then
       write(*,*)'lattice%modulo...ok'
    else
       test = .false.
    end if

    if (lattice%imp == 4.67d0) then
       write(*,*)'lattice%imp...ok'
    else
       test = .false.
    end if

    if (lattice%impdist == 'Uniform') then
       write(*,*)'lattice%imp...ok'
    else
       test = .false.
    end if

    if (lattice%range == 10) then
       write(*,*)'lattice%range...ok'
    else
       test = .false.
    end if

    if (lattice%param == 5.43d0) then
       write(*,*)'lattice%param...ok'
    else
       test = .false.
    end if

    if (sphere == .true.) then
       write(*,*)'lattice%param...ok'
    else
       test = .false.
    end if

    if (R0 == 10.D0) then
       write(*,*)'R0...ok'
    else
       test = .false.
    end if

    if (R1 == 50.D0) then
       write(*,*)'R1...ok'
    else
       test = .false.
    end if

    if (phi == 180.D0 * convert) then
       write(*,*)'phi...ok'
    else
       test = .false.
    end if

    if (theta == 90.D0 * convert) then
       write(*,*)'theta...ok'
    else
       test = .false.
    end if

    if (nuclear%x == 0 .and. &
        nuclear%y == 0 .and. &
        nuclear%z == 0) then
       write(*,*)'nuclear...ok'
    else
       test = .false.
    end if

    if (Tmax == 10.D0) then
       write(*,*)'Tmax...ok'
    else
       test = .false.
    end if

    if (nb_pts_t == 1000) then
       write(*,*)'nb_pts_t...ok'
    else
       test = .false.
    end if

    if (power == 8) then
       write(*,*)'power...ok'
    else
       test = .false.
    end if

    if (threshold == .true.) then
       write(*,*)'threshold...ok'
    else
       test = .false.
    end if

    if (threshold_val == 1.d-2) then
       write(*,*)'threshold_val...ok'
    else
       test = .false.
    end if

    if (e_n_coupling == .true.) then
       write(*,*)'e_n_coupling...ok'
    else
       test = .false.
    end if

    if (test) then
       write(*,*)'Test read_input successful'
    else
       write(*,*)'Test read_input failed'
    end if


  end subroutine unit_test_read

end module test_read
