subroutine abp_store(dr,mu,dt,ens,n,ntimes,xstore,ystore,cstore,sstore,times)
implicit none

integer ens,n,ntimes,intrvl
integer i,j,k
real dr,mu,dt,sdth
real :: xstore(ens,ntimes), ystore(ens,ntimes)
real :: cstore(ens,ntimes), sstore(ens,ntimes)
real :: times(ntimes)
real u1,u2,q,r,x,y,theta
real, parameter :: pi=3.14159265358979323846264338327950288419716939937510
real, parameter :: eps=1e-12

!f2py intent(in) 2D,mu,dt,ens,n,ntimes
!f2py intent(in,out) xstore,ystore,cstore,sstore,times
!f2py depend(ens,ntimes) xstore,ystore,cstore,sstore
!f2py depend(ntimes) times

sdth = sqrt(2.0*dr*dt)
intrvl = max(1,n/ntimes)

do j = 1,ens
    x = 0.0
    y = 0.0
    call random_number(q)
    theta = 2.0*pi*q
    k = 1
    do i = 1,n
        x = x + (-mu*x + cos(theta))*dt
        y = y + (-mu*y + sin(theta))*dt
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1,eps)
        r = sqrt(-2.0*log(u1))*cos(2.0*pi*u2)
        theta = theta + r*sdth
        if (mod(i,intrvl) == 0) then
            if (k <= ntimes) then
                xstore(j,k) = x
                ystore(j,k) = y
                cstore(j,k) = cos(theta)
                sstore(j,k) = sin(theta)
                if (j == 1) times(k) = i*dt
                k = k + 1
            endif
        endif
    enddo
enddo

end subroutine abp_store