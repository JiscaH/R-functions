subroutine distsF(E1, N1, L1, E2, N2, L2, DV)
implicit none

integer, intent(IN) :: L1, L2
double precision, intent(IN) :: E1(L1), N1(L1), E2(L2), N2(L2)
double precision, intent(INOUT) :: DV(L1 * L2)
integer :: i, j

do i=1,L1
    do j=1, L2
        DV((j-1)*L1 + i) = SQRT((E1(i) - E2(j))**2 + (N1(i) - N2(j))**2) 
    enddo
enddo

end subroutine distsF