! subsidiary program to do fourier space helmholtz decomp stuff
program spr_fh

  use common
  use input
  use setup
  use k_scatter
  implicit none
  include 'rev.inc'

  call read_input
  call initial_setup

  call write_output

  stop

end program spr_fh
