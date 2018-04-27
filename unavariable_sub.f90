program main 
use sistemas
use contorno
implicit none
integer :: i,j,np
real(8),allocatable :: x(:),m(:,:),u(:),r(:)
real(8) :: a, b,h,alpha,beta
real(8) :: val(2),besselmatrix(2,2),coeficientes(2)
	besselmatrix(1,:)=[bessel_j0(1.0d0),bessel_y0(1.0d0)]    !Nos crea la matriz para hallar los coeficientes que resuelven
	besselmatrix(2,:)=[bessel_j0(2.0d0),bessel_y0(2.0d0)]	 !el sistema de forma precisa con las funciones de bessel
	val=[0,1]						 !vector resultado dado por las condiciones de contorno
	np=100
	alpha=0.0d0
	beta=1.0d0
	b=2.0d0
	a=1.0d0
	h=(b-a)/(np-1)
	allocate(x(np),m(np,np),u(np),r(np))
	
	call invoca_matvectdif(m,r,x,np,h,a,alpha,beta)		!Llama a la subrutina que crea matriz, vector resultado y valores de x
	call gausspivote(besselmatrix,val,coeficientes)		!resuelve el sistema que hemos creado

	print *, 'bessel en 1',bessel_j0(1.0d0),bessel_y0(1.0d0)
	print *, 'bessel en 2',bessel_j0(2.0d0),bessel_y0(2.0d0)
	print *, 'coef, bessel'
	print *, coeficientes
	call gausspivote(m,r,u)
	do i=1,np
	print *, 'Valores de x', x(i), 'Valores de u', u(i),'valor de bessel',coeficientes(1)*bessel_j0(x(i))+&
	coeficientes(2)*bessel_y0(x(i))
	end do


	open(unit=13, file='solucion_unavariable.dat')
	write(13,*) '#Numero puntos:', np
	write(13,*) '#X','                U'
	do i=1,np
	write(13,*) x(i),u(i)
	end do
	close(11)

end program main
