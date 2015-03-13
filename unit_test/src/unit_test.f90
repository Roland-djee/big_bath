!------------------------------------------------------------------------------
! unit_test code for big_bath
!------------------------------------------------------------------------------
!
! MAIN: unit_test
!
!> @author
!> Dr. Roland Guichard University College London
!
! DESCRIPTION: 
!> This part of the unit_test code is the main to unit test big_bath. 
!
! REVISION HISTORY:
! 09-03-2015 - Initial Version
! TODO_09-03-2015 - Complete main - TODO_main
! TODO_09-03-2015 - Deallocate arrays - TODO_main
!------------------------------------------------------------------------------

program unit_test
  use test_read
  use test_generate
  use test_interactions
  
  implicit none

  call unit_test_read

  call unit_test_generate_crystal

  call unit_test_generate_impurities

  deallocate(Imp_x,Imp_y,Imp_z)

  call unit_test_F

  call unit_test_hyperfine

end program unit_test
