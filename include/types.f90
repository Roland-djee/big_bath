module types
  implicit none
  
  ! IO units, inputs 
  character(len=20) :: Systype
  integer(kind=8) :: nb_basisvec, nb_sites, nb_imp, nb_pairs
  integer(kind=8) :: origin
  integer(kind=8) :: en_dipolar
  integer(kind=8) :: power

  integer(kind=8) :: nb_pts_t, CP_seq
  double precision :: Tmax
  double precision :: expval_up, expval_down
  double precision :: polar_up, polar_down
  double precision :: Jnuc
  double precision :: R0, R1, theta, phi 

  logical :: crys_sites
  logical :: crys_imp
  logical :: coup_val
  logical :: deco_val
  logical :: sphere

  ! Arrays
  integer(kind=8), allocatable :: imp(:), Sites_x(:), Sites_y(:), Sites_z(:)
  integer(kind=8), allocatable :: Imp_x(:), Imp_y(:), Imp_z(:)
  double precision, allocatable :: HYPER(:)
  double precision, allocatable :: C12(:), C12_n(:)
  double precision, allocatable :: DJ(:), DC(:), DC_n(:)

  ! Vector type
  type vectors
     integer(kind=8) :: x, y, z
  end type vectors
  type (vectors) :: nuclear
  type (vectors) :: Orientation
  type (vectors), allocatable :: basis_vectors(:)

  ! Lattice type
  type lat
     integer(kind=8)   :: modulo, range   
     double precision  :: imp, param
     character(len=20) :: impdist
  end type lat
  type (lat) :: lattice

end module types
