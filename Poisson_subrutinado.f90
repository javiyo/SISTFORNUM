program main 
use sistemas
use contorno
implicit none
integer :: i,j,k,cont,np,np2
real(8),allocatable :: m(:,:),u(:),r(:),uapro(:,:),error(:)
real(8) :: a, b,h,alpha,beta,errormedio
real(8),parameter :: pi=dacos(-1.0d0),tolconj=1.0d-6
integer,parameter :: iteraconj=100
	np=20
	np2=np**2 !np elevado al cuadrado, ya que np significa el numero de puntos en una dimension (util para la matriz)
	alpha=0
	beta=1
	b=pi
	a=0
	h=(b-a)/(np-1)
	allocate(m(np2,np2),u(np2),r(np2),error(np2),uapro(np,np)) !Alocatar es importante

	call invoca_matvect(m,r,h,np)				!Llama a la subrutina que saca la matriz y vector solucion
	call rverdadero(uapro,h,a,np)				!Llama a la subrutina que nos devuelve la aproximación
	call grad_conj(m,r,u,tolconj,iteraconj)			!Resuelve el sistema por el gradiente conjugado
	
	cont=0
	open(unit=12, file='poisson_data.dat')			!Genera el fichero que vamos a representar en GNUPLOT
	write(12,*) '#Numero puntos:', np
	write(12,*) '#Numero de iteraciones gradiente conjugado', iteraconj
	write(12,*) '#Tolerancia gradiente conjugado', tolconj
	do j=1,np
		do i=1,np
		cont = cont +1
		write(12,*) cont,u(cont),uapro(j,i) 		!Escribe el numero de punto por el que vamos, aproximacion y real
		end do
	end do
	close(12)

	cont=0
	do j=1,np						!Este bucle es de comprobación por pantalla, nos hace el error medio
		do i=1,np
		cont =cont+1
		error(cont)=abs(u(cont)-uapro(j,i))
		errormedio=error(cont)+errormedio
		print *,cont,'Valores de u', u(cont), 'val aprox',uapro(j,i),'error',abs(u(cont)-uapro(j,i))
		end do
	end do
	errormedio=errormedio/np2
	print *, errormedio


	!open(unit=11, file='matriz.txt')                       !/en desuso/, comprobación de que hace la matriz correcta
	!do i=1,np2
	!write(11,"(25I2)")int(m(i,:))
	!end do
	!close(11)
	
end program main
