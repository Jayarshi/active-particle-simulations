program ABP
implicit none
!Parameter declaration
integer :: i,n,j,ens
real :: u1,u2,r,m,d2,sd2t,q
real, parameter :: dt=0.01,pi=4*atan(1.0),tmax=300, o=1e-12
real, dimension(:),allocatable :: x,avgx,avgx2,y,avgy,avgy2,theta,t
print*, "Input the values of initial D_r and Mu"
read*, d2, m
ens=10000
sd2t=sqrt(2*d2*dt)
!Allocate array size
n=int(tmax/dt)
allocate(x(1:n+1))
allocate(avgx(1:n+1))
allocate(avgx2(1:n+1))
allocate(y(1:n+1))
allocate(avgy(1:n+1))
allocate(avgy2(1:n+1))
allocate(t(1:n+1))
allocate(theta(1:n+1))
!Update
avgx=0.0
avgx2=0.0
avgy=0.0
avgy2=0.0
do j=1,ens
    t(1)=0
    x(1)=0.0
    y(1)=0.0
    call random_number(q)
    theta(1)=2*pi*q
    do i=1,n
        call random_number(u1)
        call random_number(u2)
        u1 = max(o,u1)
        r=sqrt(-2*log(u1))*cos(2*pi*u2) !Random Gaussian number
        avgx(i)=avgx(i)+x(i) !Different x(i) stored for different ensemble
        avgx2(i)=avgx2(i)+x(i)*x(i)
        avgy(i)=avgy(i)+y(i) !Different y(i) stored for different ensemble
        avgy2(i)=avgy2(i)+y(i)*y(i)
        theta(i+1)=theta(i)+r*sd2t
        x(i+1)=x(i)+(-m*x(i)+cos(theta(i)))*dt
        y(i+1)=y(i)+(-m*y(i)+sin(theta(i)))*dt
        t(i+1)=t(i)+dt
    enddo
enddo
!Generating data points
open(unit=10,file="ABP_data.dat")
do i=1,n+1
    if (mod(i,60)==0) then
        write(10,*) t(i), avgx(i)/ens, avgy(i)/ens, &
                    avgx2(i)/ens, avgy2(i)/ens, &
                    (avgx2(i)+avgy2(i))/ens
    endif
enddo
close(10)
end program ABP