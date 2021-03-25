implicit none
double precision::t,tf,ht,hx,xi,xf,hx3,hx2,d,x(-2:10002),u(-2:10002),k0(-2:10002),k1(-2:10000),k2(-2:10000),&
k3(-2:10002),u0(-2:10002),u1(-2:10002),u2(-2:10002),u3(-2:10002),sp
real::finish,start
integer::nt,nx,ti,tfi
integer::i,j,a1,p

!write(*,*) "deltasquare,tf,ht,hx"
!read(*,*) d,tf,ht,hx


hx=0.1
nx=10/hx

do i=-2,nx+2
x(i)=i*hx
end do




d=1
tf=1000
ht=0.0001
a1=1
p=1
sp=16


!hx=0.2 used

!specify initial condition
do i=-2,nx+2
!u(i)=3*d**(1/3)*sp/(dcosh(0.5*sqrt(sp)*(d**(-1/3)*(x(i)-x(3*nx/5)))))**2!&
!+3*d**(1/3)*50/(cosh(0.5*sqrt(50.0)*(d**(-1/3)*(x(i)-x(2*nx/5)))))**2
!write(*,*) u(i),x(i)
!u(i)=-2*exp(-10*(x(i)-x(nx/2))**2)
u(i)=sin(2*3.1419*i/nx)
!u(i)=10*dfloat(i-nx/2)*exp(-0.1*dfloat(i-nx/2)**2)


end do


!d3ydx3=u(i+2)-2*u(i+1)+2*u(i-1)-u(i-2)
!dydx=(u(i+1)-u(i-1))/(hx2)

open(1,file="sima.dat")

t=0
do while(t.le.tf)




do j=0,nx-1
	
	u0(j)=u(j)
	k0(j)=-u(j)*(u(j+1)-u(j-1))/(2*hx)-d*(u(j+2)-2*u(j+1)+2*u(j-1)-u(j-2))/(2*hx**3)
	u1(j)=u0(j)+ht*k0(j)/2
	
end do

k0(-1)=k0(nx-1)
k0(-2)=k0(nx-2)
k0(nx)=k0(0)
k0(nx+1)=k0(1)
k0(nx+2)=k0(2)
u1(-1)=u1(nx-1)
u1(-2)=u1(nx-2)
u1(nx+1)=u1(1)
u1(nx+2)=u1(2)
u1(nx)=u1(0)

do j=0,nx-1

	
	k1(j)=-u1(j)*(u1(j+1)-u1(j-1))/(2*hx)-d*(u1(j+2)-2*u1(j+1)+2*u1(j-1)-u1(j-2))/(2*hx**3)
	u2(j)=u0(j)+ht*k1(j)/2

end do
k1(-1)=k1(nx-1)
k1(-2)=k1(nx-2)
k1(nx+1)=k1(1)
k1(nx+2)=k1(2)
k1(nx)=k1(0)
u2(-1)=u2(nx-1)
u2(-2)=u2(nx-2)
u2(nx+1)=u2(1)
u2(nx+2)=u2(2)
u2(nx)=u2(0)

do j=0,nx-1

	k2(j)=-u2(j)*(u2(j+1)-u2(j-1))/(2*hx)-d*(u2(j+2)-2*u2(j+1)+2*u2(j-1)-u2(j-2))/(2*hx**3)
	u3(j)=u0(j)+ht*k1(j)

end do
k2(-1)=k2(nx-1)
k2(-2)=k2(nx-2)
k2(nx+1)=k2(1)
k2(nx+2)=k2(2)
k2(nx)=k2(0)
u3(-1)=u3(nx-1)
u3(-2)=u3(nx-2)
u3(nx+1)=u3(1)
u3(nx+2)=u3(2)
u3(nx)=u3(0)


do j=0,nx-1
	
	k3(j)=-u3(j)*(u3(j+1)-u3(j-1))/(2*hx)-d*(u3(j+2)-2*u3(j+1)+2*u3(j-1)-u3(j-2))/(2*hx**3)
end do
k3(-1)=k3(nx-1)
k3(-2)=k3(nx-2)
k3(nx+1)=k3(1)
k3(nx+2)=k3(2)
k3(nx)=k3(0)

if(mod(a1,100000).eq.0) then
write(*,*) t
p=p+1
end if


	do j=-2,nx+2
	u(j)=u(j)+ht*(k0(j)+2*k1(j)+2*k2(j)+k3(j))/6
		if(mod(a1,2000).eq.0 .and. mod(j,1).eq.0 .and. j.ge.0 .and. j.le.nx) then
		write(1,*) a1/2000,x(j),u(j)
		end if	
	end do

t=t+ht
a1=a1+1


end do
close(1)


end
