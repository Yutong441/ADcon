module params
    implicit none
    real(8), parameter :: inf = 1e10
end module params

subroutine binarize(a, n, d)
    ! binarize an array thresholded at 0
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n, n), intent(out) :: d
    integer :: i, j
    do j = 1, n
        do i = 1, n
            if (a(i, j) /= 0.) then
                d(i, j) = 1
            else
                d(i, j) = 0
            end if
        end do
    end do
end subroutine binarize


subroutine eig_sym(a, n, lambda, eig_vec)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n), intent(out) :: lambda
    real(8), dimension(n, n), intent(out) :: eig_vec
    real(8), allocatable :: work(:)
    integer :: lwork, info
    external :: dsyev  ! from lapack

    eig_vec = a
    lwork = 3*n - 1
    allocate(work(lwork))
    call dsyev("v", "u", n, eig_vec, n, lambda, work, lwork, info)
    deallocate(work)
end subroutine eig_sym


subroutine delete_sym_mat(a, n, k, nrow_d, d2)
    ! delete the rows and columns of `a` if `k` > 0
    ! a: n x n array
    implicit none
    integer, intent(in) :: n, nrow_d
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n), intent(in) :: k
    integer :: i

    integer(8), dimension(nrow_d) :: d1
    real(8), dimension(nrow_d, nrow_d), intent(out) :: d2

    d1 = pack([(i, i=1, n)], k > 0)
    d2 = a(d1, d1)
end subroutine delete_sym_mat


subroutine invert(a, n, d)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n, n), intent(out) :: d
    integer :: i, j

    do j = 1, n
        do i = 1, n
            if (a(i, j) /= 0) then
                d(i, j) = 1/a(i, j)
            else
                d(i, j) = 0
            end if
        end do
    end do
end subroutine invert
