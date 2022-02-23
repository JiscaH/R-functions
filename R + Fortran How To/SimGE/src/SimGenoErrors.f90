! Array storage in Fortran is column-major:
! the first index varies fastest

subroutine mkerrors(nind, nsnp, genofr, eprobfr)
implicit none

integer, intent(IN) :: nind, nsnp
integer, intent(INOUT) :: genofr(nind*nsnp)
double precision, intent(IN) :: eprobfr(nsnp*3*3)
integer :: Genos(nSnp, nInd), l, i, x, j, h
double precision :: EProb(3,3,nSnp), p(3)
real :: r(nSnp, nInd)

Genos = -9
j = 0
do l=1,nSnp
  do i=1, nInd
    j = j+1
    if (GenoFR(j)/=-9) then
      Genos(l,i) = GenoFR(j) +1
    endif
  enddo
enddo

j=0
do i=1,3
  do h=1,3
    do l=1,nSnp
      j = j+1
      EProb(h,i,l) = EProbFR(j)
    enddo
  enddo
enddo

call init_random_seed()
call RANDOM_NUMBER(r)

j=0
do l=1, nSnp
  do i=1,nInd
    if (Genos(l,i)==-9)  cycle
    p = Eprob(Genos(l,i), :, l) / sum(Eprob(Genos(l,i), :, l))
    if (dble(r(l,i)) < p(1)) then
      x = 0
    else if (dble(r(l,i)) < (p(1) + p(2))) then
      x = 1
    else
      x = 2
    endif
    
    j = j+1
    GenoFR(j) = x
  enddo
enddo

end subroutine MkErrors


subroutine init_random_seed()
! Copyright (c) 2012, pplab (shenyu828@gmail.com)
! random_seed based on system time
    implicit none
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)
end subroutine init_random_seed
