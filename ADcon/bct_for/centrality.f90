subroutine eigenvector_centrality_und(a, n, d)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n, 1), intent(out) :: d
    real(8), dimension(n) :: lambda
    real(8), dimension(n, n) :: eig_vec
    integer, dimension(1) :: i
    external :: eig_sym

    call eig_sym(a, n, lambda, eig_vec)
    i = maxloc(lambda)
    d = abs(eig_vec(:, i))
end subroutine eigenvector_centrality_und


subroutine clustering_coef_wu(a, n, d)
    use :: params
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n), intent(out) :: d

    real(8), dimension(n, n) :: a_bi
    real(8), dimension(n, n) :: ws
    real(8), dimension(n) :: k
    integer :: i
    external :: binarize

    call binarize(a, n, a_bi)
    k = dble(sum(a_bi, 2))
    ws = a**(1.0/3.0)
    ws = matmul(ws, matmul(ws, ws))
    do i = 1, n
        d(i) = ws(i, i)
        if (d(i) == 0) then
            k(i) = inf
        end if
    end do
    d = d / (k * (k - 1))
end subroutine clustering_coef_wu


subroutine rich_club_wu(a, n, k, rw)
    implicit none
    integer, intent(in) :: n, k
    real(8), dimension(n, n), intent(in) :: a
    real(8), intent(out) :: rw

    real(8), dimension(n) :: deg, deg_k
    real(8), dimension(n*n) :: wrank
    integer :: i, nrow_d, er, info
    real(8), dimension(:, :), allocatable :: cutcij

    external :: degrees_und, binarize, delete_sym_mat, find_k, dlasrt

    call degrees_und(a, n, deg)
    ! sort the weights of the network, with the strongest connection first
    wrank = reshape(a, (/n*n/))
    call dlasrt('D', n*n, wrank, info)

    do i = 1, n
        deg_k(i) = deg(i) - k + 1 ! small nodes
    end do
    nrow_d = count(deg_k > 0)
    if (nrow_d <= 0) then
        rw = 1e10
    else
        ! remove small nodes with node degree < k
        allocate(cutcij(nrow_d, nrow_d))
        call delete_sym_mat(a, n, deg_k, nrow_d, cutcij)
        ! total number of connections in subset E>r
        er = count(cutcij /= 0)
        ! E>r number of connections with max weight in network
        ! total weight of connections in subset E>r / weighted rich-club coefficient
        rw = sum(cutcij) / sum(wrank(1:er))
    end if
end subroutine rich_club_wu


subroutine rich_club_all(a, n, rc)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n), intent(out) :: rc
    real(8), dimension(n) :: deg
    integer :: info, max_deg, i
    external :: degrees_und, rich_club_wu, dlasrt

    call degrees_und(a, n, deg)
    max_deg = maxval(deg)
    do i = 1, max_deg
        call rich_club_wu(a, n, i, rc(i))
    end do
end subroutine rich_club_all


subroutine find_k(a, n, top_perc, k)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: top_perc
    real(8), dimension(n, n), intent(in) :: a
    integer, intent(out) :: k
    real(8), dimension(n) :: deg
    integer :: info, top_degree
    external :: dlasrt, degrees_und

    call degrees_und(a, n, deg)
    call dlasrt('d', n, deg, info)
    top_degree = int(top_perc*n)
    k = int(deg(top_degree))
end subroutine find_k
