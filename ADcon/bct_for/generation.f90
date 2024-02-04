subroutine random_choice(k, k_out)
    ! k: maximum integer
    implicit none
    integer, intent(in) :: k
    integer, intent(out) :: k_out
    real(8) :: u
    intrinsic :: random_number

    call random_number(u)
    k_out = floor(u*k) + 1
end subroutine random_choice


subroutine randmio_und(g, n, itr, r)
    ! equivalent to bct.randmio_und but up to 100 times faster
    implicit none
    integer, intent(in) :: n, itr
    real(8), dimension(n, n), intent(in) :: g
    real(8), dimension(n, n), intent(out) :: r
    integer :: k, pos_k, max_attempts, itr_all, e1, e2, it, ia, a, b, c, d
    integer, dimension(:), allocatable :: i, j
    real(8) :: u
    external :: random_choice

    ! the following block is equivalent to np.where(np.tril(a))
    k = count(g /= 0)/2
    if (allocated(i)) deallocate(i)
    if (allocated(j)) deallocate(j)
    allocate(i(k))
    allocate(j(k))
    pos_k = 1
    do it = 1, n
        do ia = 1, n
            if (it >= ia .and. g(it, ia) /= 0) then
                i(pos_k) = it
                j(pos_k) = ia
                pos_k = pos_k + 1
            end if
        end do
    end do

    ! maximum number of rewiring attempts per iteration
    max_attempts = n * k / (n * (n - 1))
    itr_all = itr*k
    r(:, :) = g(:, :)

    iteration: do it = 1, itr_all
        attempt: do ia = 1, max_attempts  ! while not rewired
            random: do
                call random_choice(k, e1)
                call random_choice(k, e2)
                do
                    if (e1 == e2) then
                        call random_choice(k, e2)
                    else
                        exit
                    end if
                end do

                a = i(e1)
                b = j(e1)
                c = i(e2)
                d = j(e2)

                if (a /= c .and. a /= d .and. b /= c .and. b /= d) then
                    exit random ! all 4 vertices must be different
                end if
            end do random

            call random_number(u)
            if (u > .5) then
                i(e2) = d
                j(e2) = c  ! flip edge c-d with 50% probability
                c = i(e2)
                d = j(e2)  ! to explore all potential rewirings
            end if

            ! rewiring condition
            if (.not. (r(a, d) /= 0 .or. r(c, b) /= 0)) then
                r(a, d) = r(a, b)
                r(a, b) = 0
                r(d, a) = r(b, a)
                r(b, a) = 0
                r(c, b) = r(c, d)
                r(c, d) = 0
                r(b, c) = r(d, c)
                r(d, c) = 0

                j(e1) = d
                j(e2) = b  ! reassign edge indices
                exit attempt
            end if
        end do attempt
    end do iteration
    deallocate(i)
    deallocate(j)
end subroutine randmio_und
