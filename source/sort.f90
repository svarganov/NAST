module sort

contains
!---------------------------------------------------
subroutine insert_sort(myarray,fun)
implicit none

double precision, intent(inout) :: myarray(:)
double precision, optional      :: fun
integer                         :: i,j
double precision                :: key

do j = 2, size(myarray)
 key = myarray(j)
 i = j - 1
 do while ( (i > 0) .and. (myarray(i) > key) )
    myarray(i+1) = myarray(i)
    i = i - 1
 end do
 myarray(i+1)=key
end do

end subroutine insert_sort
!---------------------------------------------------
subroutine shake_sort(myarray)
implicit none

double precision, intent(inout) :: myarray(:)
integer                         :: left, right, i,j
double precision                :: tmp

left = 2
right = size(myarray)

do while (left <= right)
  i = right
  j = left
  do while (i>=left)
     if (myarray(i-1) > myarray(i)) then
        tmp = myarray(i)
        myarray(i) = myarray(i-1)
        myarray(i-1) = tmp
     end if
     i = i-1
  end do
  left = left+1

  do while (j<=right)
    if (myarray(j-1) > myarray(j)) then
       tmp = myarray(j)
       myarray(j) = myarray(j-1)
       myarray(j-1) = tmp
    end if
    j =j+1
  end do
  right = right-1
end do

end subroutine shake_sort
!--------------------------------------------------
end module
