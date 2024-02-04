subroutine degrees_und(a, n, d)
    ! degree centrality for undirected network
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n), intent(out) :: d
    real(8), dimension(n, n) :: a_bi
    external :: binarize

    call binarize(a, n, a_bi)
    d = sum(a_bi, 1)
end subroutine degrees_und


subroutine strengths_und(a, n, d)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n), intent(out) :: d
    d = sum(a, 1)
end subroutine strengths_und


subroutine density_und(a, n, d)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: a
    real(8), intent(out) :: d
    integer :: i, j

    d = 0
    do i = 1, n
        do j = 1, n
            if (a(i, j) > 0.) then
                d = d + 1
            end if
        end do
    end do

    d = d / (dble(n) * (dble(n) - 1) )
end subroutine density_und
