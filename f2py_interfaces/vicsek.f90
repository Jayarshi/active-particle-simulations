subroutine evolve(x, theta, n, num_part, rad_int, noise, v, dt, boxL)
    implicit none
    
    integer, intent(in) :: n, num_part
    real(8), intent(inout) :: x(num_part, 2)
    real(8), intent(inout) :: theta(num_part, 1)
    real(8), intent(in) :: rad_int, noise, v, dt, boxL
    
    !f2py intent(in,out) x
    !f2py intent(in,out) theta
    !f2py intent(in) n
    !f2py intent(hide), depend(x) :: num_part=shape(x,0)
    !f2py intent(in) rad_int
    !f2py intent(in) noise
    !f2py intent(in) v
    !f2py intent(in) dt
    !f2py intent(in) boxL
    
    integer :: step, i, j, nc, cx, cy, ncx, ncy, dxcell, dycell
    real(8) :: cell_size, r2, dx, dy, sx, sy, s
    integer, allocatable :: head(:,:), linked(:)
    real(8), allocatable :: avgt(:)
    
    ! Calculate actual grid size and restrict dynamically
    nc = int(boxL / rad_int) ! Nearby cells
    if (nc < 1) nc = 1
    cell_size = boxL / dble(nc)
    r2 = rad_int * rad_int
    
    ! Pre-allocate temporary memory OUTSIDE the step loop to maximize speed
    allocate(head(0:nc-1, 0:nc-1))
    allocate(linked(num_part))
    allocate(avgt(num_part))
    
    do step = 1, n
        ! --- 1. Cell Linking ---
        head = -1
        linked = -1
        do i = 1, num_part
            cx = int(x(i, 1) / cell_size)
            cy = int(x(i, 2) / cell_size)
            
            ! Ensure floating point boundary safety
            if (cx < 0) cx = 0
            if (cx >= nc) cx = nc - 1
            if (cy < 0) cy = 0
            if (cy >= nc) cy = nc - 1
            
            linked(i) = head(cx, cy)
            head(cx, cy) = i
        end do
        
        ! --- 2. Calculate New Average Thetas ---
        do i = 1, num_part
            sx = 0.0d0
            sy = 0.0d0
            cx = int(x(i, 1) / cell_size)
            cy = int(x(i, 2) / cell_size)
            
            if (cx < 0) cx = 0
            if (cx >= nc) cx = nc - 1
            if (cy < 0) cy = 0
            if (cy >= nc) cy = nc - 1
            
            ! Iterate 9 neighboring cells with periodic conditions
            do dxcell = -1, 1
                do dycell = -1, 1
                    ncx = modulo(cx + dxcell, nc)
                    ncy = modulo(cy + dycell, nc)
                    
                    j = head(ncx, ncy)
                    do while (j /= -1)
                        ! Check pairwise distances exactly as the Python equivalent does
                        dx = x(i, 1) - x(j, 1)
                        dy = x(i, 2) - x(j, 2)
                        
                        if (dx*dx + dy*dy <= r2) then
                            sx = sx + cos(theta(j, 1))
                            sy = sy + sin(theta(j, 1))
                        end if
                        j = linked(j)
                    end do
                end do
            end do
            avgt(i) = atan2(sy, sx)
        end do
        
        ! --- 3. Update Positions & Orientations ---
        do i = 1, num_part
            call random_number(s)
            s = s - 0.5d0
            theta(i, 1) = avgt(i) + noise * s * dt
            
            x(i, 1) = x(i, 1) + v * cos(theta(i, 1)) * dt
            x(i, 2) = x(i, 2) + v * sin(theta(i, 1)) * dt
            
            ! Periodic Boundary Conditions
            x(i, 1) = modulo(x(i, 1), boxL)
            x(i, 2) = modulo(x(i, 2), boxL)
        end do
    end do
    
    deallocate(head)
    deallocate(linked)
    deallocate(avgt)

end subroutine evolve
