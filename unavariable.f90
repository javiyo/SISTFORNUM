program main 
use sistemas
implicit none
integer :: i,j,np
real(8),allocatable :: x(:),m(:,:),u(:),r(:)
real(8) :: a, b,h,alpha,beta
real(8) :: val(2),besselmatrix(2,2),coeficientes(2)
	besselmatrix(1,:)=[bessel_j0(1.0d0),bessel_y0(1.0d0)]
	besselmatrix(2,:)=[bessel_j0(2.0d0),bessel_y0(2.0d0)]
	val=[0,1]
	np=100
	alpha=0
	beta=1
	b=2.0
	a=1.0
	h=(b-a)/(np-1)
	allocate(x(np),m(np,np),u(np),r(np))
	m=0
	m(1,1)=1
	m(np,np)=1


	do i=1,np
	x(i)=a+(real(i)-1.0)*h
	end do


	do i=2,np-1
	m(i,i-1)=(x(i)*x(i))-x(i)*(h/2.0)			!Bucle que aproxima los puntos i
	m(i,i)=-2.0*(x(i)*x(i))+((x(i)*x(i))*(h*h))
	m(i,i+1)=(x(i)*x(i))+x(i)*(h/2.0)
	r(i)=0	
	end do
	
	
	r(1)=alpha
	r(np)=beta
	do i=1,np
	print *, 'Valores de x', x(i)
	end do
	call gausspivote(besselmatrix,val,coeficientes)
	print *, 'bessel en 1',bessel_j0(1.0d0),bessel_y0(1.0d0)
	print *, 'bessel en 2',bessel_j0(2.0d0),bessel_y0(2.0d0)
	print *, 'coef, bessel'
	print *, coeficientes
	call gausspivote(m,r,u)
	do i=1,np
	print *, 'Valores de x', x(i), 'Valores de u', u(i),'valor de bessel',coeficientes(1)*bessel_j0(x(i))+&
	coeficientes(2)*bessel_y0(x(i))
	end do
end program main
