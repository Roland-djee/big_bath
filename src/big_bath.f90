program big_bath
  use types
  use constants
  use read
  use generate
  use write
  use interactions
  implicit none

  call read_inputs_generate

  call generate_crystal
  
  call crystal_outputs

  call generate_impurities

  call impurities_outputs

  call couplings
  
  call couplings_outputs

  call select_nuc
 
end program big_bath
