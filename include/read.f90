module read
  use types
  use constants
  implicit none

contains

subroutine read_inputs_generate
  implicit none
  ! local variables
  integer(kind=8) :: i, System=10

  character (len=*), parameter :: fmt1 = "(t50, i, i, i)"
  character (len=*), parameter :: fmt2 = "(t50, i)"
  character (len=*), parameter :: fmt3 = "(t50, es)"
  character (len=*), parameter :: fmt4 = "(t50, a)"
  character (len=*), parameter :: fmt5 = "(t50, l)"
  
  open(unit=System, file="../input/System.inp")
  read(System, fmt1) Orientation%x, Orientation%y, Orientation%z
  read(System, fmt4) Systype
  read(System, fmt2) nb_basisvec
  allocate(basis_vectors(nb_basisvec))
  do i=1,nb_basisvec
     read(System, fmt1)basis_vectors(i)%x, basis_vectors(i)%y, &
          basis_vectors(i)%z
  end do
  read(System, fmt2) lattice%modulo
  read(System, fmt3) lattice%imp
  read(System, fmt4) lattice%impdist
  read(System, fmt2) lattice%range
  read(System, fmt3) lattice%param
  read(System, fmt5) sphere
  read(System, fmt3) R0
  read(System, fmt3) R1
  read(System, fmt3) phi
  read(System, fmt3) theta
  read(System, fmt1) nuclear%x, nuclear%y, nuclear%z
  read(System, fmt3) Tmax
  read(System, fmt2) nb_pts_t
  read(System, fmt2) power
  close(System)

  phi   = phi * convert
  theta = theta * convert

end subroutine read_inputs_generate

end module read