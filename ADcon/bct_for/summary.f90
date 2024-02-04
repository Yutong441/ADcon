subroutine graph_metrics(m, n, k, globa, noda)
    implicit none
    integer, intent(in) :: n, k
    real(8), dimension(n, n), intent(in) :: m
    real(8), dimension(12), intent(out) :: globa
    real(8), dimension(5, n), intent(out) :: noda

    real(8), dimension(n) :: cent_eig, cent_str, cent_deg, clust, cent_bet, ecc
    real(8), dimension(n, n) :: d
    real(8) :: rc, rc_indiv, dense, lambd, ge, radius, diam
    integer :: k_indiv
    real(8) :: top_perc = 0.12
    external :: eigenvector_centrality_und, strengths_und, degrees_und
    external :: find_k, clustering_coef_wu, rich_club_wu, density_und
    external :: distance_wei

    call eigenvector_centrality_und(m, n, cent_eig)
    call strengths_und(m, n, cent_str)
    call degrees_und(m, n, cent_deg)
    call clustering_coef_wu(m, n, clust)
    call density_und(m, n, dense)
    call rich_club_wu(m, n, k, rc)
    call find_k(m, n, top_perc, k_indiv)
    call rich_club_wu(m, n, k_indiv, rc_indiv)
    call distance_wei(m, n, d, cent_bet)
    call charpath(d, n, lambd, ge, ecc, radius, diam)

    globa(1) = sum(cent_deg)/dble(n)
    globa(2) = sum(cent_eig)/dble(n)
    globa(3) = sum(cent_str)/dble(n)
    globa(4) = sum(clust)/dble(n)
    globa(5) = dense
    globa(6) = rc
    globa(7) = rc_indiv
    globa(8) = sum(cent_bet)/dble(n)
    globa(9) = ge
    globa(10) = lambd
    globa(11) = radius
    globa(12) = diam

    noda(1, :) = cent_deg
    noda(2, :) = cent_eig
    noda(3, :) = cent_str
    noda(4, :) = clust
    noda(5, :) = cent_bet
end subroutine graph_metrics


subroutine normalize_metrics(m, n, k, itr, globa, norm_noda)
    ! Args:
    ! `m`: adjacency matrix (of weight)
    ! `n`: first dimension of the matrix
    ! `itr`: generate how many random networks with the same degree distribution 
    ! to compare with `m`
    ! Return: 
    ! `globa`: 11 x 2 matrix: first column is unnormalized graph metric statistics
    ! second column is the normalized ones
    ! `norm_noda`: 5 x n matrix
    use omp_lib
    implicit none
    integer, intent(in) :: n, k, itr
    real(8), dimension(n, n), intent(in) :: m
    real(8), dimension(12, 2), intent(out) :: globa
    real(8), dimension(5, n), intent(out) :: norm_noda

    real(8), dimension(12) :: norm_globa, rand_globa, denom_globa, denom_partial
    real(8), dimension(5, n) :: rand_noda
    real(8), dimension(n, n) :: m_norm, mr
    integer :: i

    external :: graph_metrics, randmio_und

    ! normalize the matrix
    m_norm = m
    call graph_metrics(m_norm, n, k, norm_globa, norm_noda)

    !$OMP PARALLEL PRIVATE(mr, rand_globa, denom_partial, rand_noda) SHARED(m_norm, denom_globa)
    ! generate random matrices
    denom_globa(:) = 0
    denom_partial(:) = 0

    !$OMP DO
    do i = 1, itr
        call randmio_und(m_norm, n, 1, mr)
        call graph_metrics(mr, n, k, rand_globa, rand_noda)
        denom_partial = denom_partial + rand_globa
    end do

    !$OMP CRITICAL
    denom_globa = denom_globa + denom_partial
    !$OMP END CRITICAL

    !$OMP END PARALLEL

    denom_globa = denom_globa / itr
    globa(:, 1) = norm_globa(:)
    globa(:, 2) = norm_globa(:) / denom_globa(:)
end subroutine normalize_metrics


subroutine rich_club_norm(a, n, itr, rc)
    use omp_lib
    implicit none
    integer, intent(in) :: n, itr
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n), intent(out) :: rc

    real(8), dimension(n) :: rc_ori, rc_partial, rc_globa, rc_rand
    real(8), dimension(n, n) :: a_rand
    integer :: i

    external :: rich_club_all, randmio_und

    ! normalize the matrix
    call rich_club_all(a, n, rc_ori)

    !$OMP PARALLEL PRIVATE(a_rand, rc_rand, rc_partial) SHARED(a, rc_globa)
    ! generate random matrices
    rc_globa(:) = 0
    rc_partial(:) = 0

    !$OMP DO
    do i = 1, itr
        call randmio_und(a, n, 1, a_rand)
        call rich_club_all(a_rand, n, rc_rand)
        rc_partial = rc_partial + rc_rand
    end do

    !$OMP CRITICAL
    rc_globa = rc_globa + rc_partial
    !$OMP END CRITICAL

    !$OMP END PARALLEL

    rc_globa = rc_globa / itr
    rc = rc_ori/rc_globa
end subroutine rich_club_norm
