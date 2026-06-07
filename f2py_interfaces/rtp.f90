subroutine rtp_store(g,d2,dt,ens,n,ntimes,xstore,times)
implicit none
integer ens,n,ntimes,intrvl
integer i,j,k
real g,d2,dt,gt,sd2t
real :: xstore(ens,ntimes)
real :: times(ntimes)
real u1,u2,p,q,r,v,x
real, parameter :: pi=3.14159265358979323846264338327950288419716939937510

!f2py intent(in) g,d2,dt,ens,n,ntimes
!f2py intent(in,out) xstore,times
!f2py depend(ntimes) times
!f2py depend(ens,ntimes) xstore

gt = g*dt
sd2t = sqrt(d2*dt)
intrvl = int(n/ntimes)

do j = 1, ens
    x = 0.0
    call random_number(q)
    if (q <= 0.5) then
        v = 1.0
    else
        v = -1.0
    endif
    k = 1
    do i = 1, n
        if (i <= n) then
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1,1e-12)
            r = sqrt(-2.0*log(u1))*cos(2.0*pi*u2)
            x = x + v*dt + r*sd2t
            call random_number(p)
            if (p <= gt) v = -v
            if (mod(i,intrvl) == 0) then
                if (k <= ntimes) then
                    xstore(j,k) = x
                    if (j == 1) times(k) = i*dt
                    k = k + 1
                endif
            endif
        endif
    enddo
enddo

end subroutine
