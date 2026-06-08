program RTP
implicit none

integer :: i,n,j,ens
real :: v,v0,tmax,p,u1,u2,q,r,a,g,d2
real, parameter :: dt=1e-2, pi=4*atan(1.0), o=1e-12
real, dimension(:), allocatable :: x,t,avgx,avgx2

print*, "Input initial position a, speed v0, gamma, d2, tmax"
read*, a,v0,g,d2,tmax

n = int(tmax/dt)
ens = 10000

allocate(x(1:n+1))
allocate(avgx(1:n+1))
allocate(avgx2(1:n+1))
allocate(t(1:n+1))

avgx  = 0.0
avgx2 = 0.0
do i = 1,n+1
    t(i) = (i-1)*dt
enddo
do j = 1,ens
    x(1) = a
    call random_number(q)
    if (q <= 0.5) then
        v = v0
    else
        v = -v0
    endif
    avgx(1)  = avgx(1)  + x(1)
    avgx2(1) = avgx2(1) + x(1)*x(1)
    do i = 1,n
        call random_number(u1)
        call random_number(u2)
        u1 = max(o,u1)
        r = sqrt(-2.0*log(u1))*cos(2.0*pi*u2)
        x(i+1) = x(i) + v*dt + r*sqrt(d2*dt)
        call random_number(p)
        if (p <= g*dt) v = -v
        avgx(i+1)  = avgx(i+1)  + x(i+1)
        avgx2(i+1) = avgx2(i+1) + x(i+1)*x(i+1)
    enddo
enddo

open(unit=10,file="RTP_data.dat")
do i = 1,n+1
    write(10,*) t(i), avgx(i)/ens, avgx2(i)/ens
enddo
close(10)

end program RTP