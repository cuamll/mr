

subroutine update_canonical()
  use common
  implicit none
  integer :: incdec, comp

  ! do we want something like this???????
  e(x,y,z,1) ! x component
  e(x,y,z,2) ! y component
  e(x,y,z,3) ! z component


  comp = rand3()
  incdec = pm1()

  ! now we can do some cool clever shit

  ! need to figure out necessary checks to keep it canonical

  e(x,y,z,comp) + incdec * delta

  ! and for rotational updates:
  e(x,y,z,comp) + delta
  e(x,y,z,(mod(comp + 1, 3) + 1))

  ! HAHA! SEE WHAT I HAVE DONE HERE LOL

end subroutine update_canonical

! int casts in fortran do greatest number. i.e. ceiling
function rand3()
  integer, intent(out) :: result
  result = int(rand() * 3)
end function rand3

function pm1()
  integer, intent(out) :: result
  result = 2 * int(2 * rand() - 1) - 1
end function pm1
