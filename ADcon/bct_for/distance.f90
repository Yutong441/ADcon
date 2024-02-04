subroutine init_distance(n, d)
    ! Dijkstra's algorithm starts with infinite distance
    use :: params
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(out) :: d
    integer :: i

    d(:, :) = inf
    do i = 1, n
        d(i, i) = 0
    end do

end subroutine init_distance


subroutine distance_wei(a, n, d, bc)
    ! return the shortest distance between all pairs of vertices + betweenness centrality
    use :: params
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n, n), intent(out) :: d
    real(8), dimension(n), intent(out) :: bc

    ! for distance
    real(8), dimension(n, n) :: a1, a_inv
    integer(2), dimension(n) :: status_vec
    integer :: u, i, j, len_v
    integer, dimension(:), allocatable :: v  ! intermediate nodes 
    real(8), dimension(n) :: w  ! final nodes
    real(8) :: new_d, min_d

    ! for betweenness
    integer, dimension(n) :: np, q  ! np = number of shortest path
    real(8), dimension(n) :: dp
    integer, dimension(n, n) :: p
    integer :: q_ind, z

    external :: init_distance, invert

    call invert(a, n, a_inv)
    call init_distance(n, d)
    bc(:) = 0
    all_node: do u = 1, n
        ! distance permanence (true is temporary)
        a1 = a_inv
        status_vec(:) = 1

        ! initialize v (intermediate node list)
        len_v = 1
        if (allocated(v)) then
            deallocate(v)
        end if
        allocate(v(len_v))
        v(1) = u

        ! initialize arrays for betweenness
        np(:) = 0
        np(u) = 1
        p(:, :) = 0
        q(:) = 1  ! indices
        q_ind = n

        whileloop: do
            intermediate: do j = 1, len_v
                if (q_ind > 0) then
                    q(q_ind) = v(j)
                end if
                q_ind = q_ind - 1

                status_vec(v(j)) = 0  ! distance u->V is now permanent
                a1(:, u) = 0  ! no in-edges as already shortest
                w = reshape(a1(v(j), :), (/n/))

                do i = 1, n
                    ! among the terminal nodes
                    if (w(i) /= 0) then
                        new_d = d(u, v(j)) + a1(v(j), i)
                        if (d(u, i) > new_d) then
                            d(u, i) = new_d
                            np(i) = np(v(j))
                            p(i, :) = 0
                            p(i, v(j)) = 1
                        else if (d(u, i) == new_d) then
                            np(i) = np(i) + np(v(j)) ! NP(u->w) sum of old and new
                            p(i, v(j)) = 1  ! v is also predecessor
                        end if
                    end if
                end do
            end do intermediate

            if (sum(status_vec) == 0) then
                exit whileloop  ! all nodes reached
            end if

            ! same outcome as np.min(D[u, status_vec])
            min_d = minval(pack(d(u, :), status_vec == 1))
            if (min_d == inf) then  ! some nodes cannot be reached
                exit whileloop
            end if

            ! same outcome as np.where(d[u, :] == min_d)
            len_v = count(d(u, :) == min_d)
            if (len_v == 0) then
                exit whileloop
            else
                deallocate(v)
                allocate(v(len_v))
                v = pack([(i, i=1, n)], d(u, :) == min_d)
            end if
        end do whileloop

        dp(:) = 0.
        do j = 1, n -1
            z = q(j)
            bc(z) = bc(z) + dp(z)
            do i = 1, n
                if (p(z, i) /= 0) then
                    dp(i) = dp(i) + (1 + dp(z)) * np(i) / np(z)
                end if
            end do
        end do
    end do all_node
    bc(:) = bc(:) / dble((n - 1)*(n - 2))
end subroutine distance_wei


subroutine charpath(a, n, lambda, efficiency, ecc, radius, diameter)
    use :: params
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: a
    real(8), dimension(n, n) :: a_mask
    real(8), intent(out) :: lambda, efficiency, radius, diameter
    real(8), dimension(n), intent(out) :: ecc
    integer :: i, j, num

    lambda = 0.
    efficiency = 0.
    num = 0
    do i = 1, n
        do j = 1, n
            if (a(i, j) /= inf .and. i /= j) then
                lambda = lambda + a(i, j)
                efficiency = efficiency + 1/a(i, j)
                num = num + 1
                a_mask(i, j) = a(i, j)
            else
                a_mask(i, j) = 0
            end if
        end do
    end do

    lambda = lambda/num
    efficiency = efficiency/num
    ecc = maxval(a_mask, 2)
    radius = minval(ecc)
    diameter = maxval(ecc)
end subroutine charpath
