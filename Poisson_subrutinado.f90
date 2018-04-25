program main 
use sistemas
use contorno
implicit none
integer :: i,j,k,cont,np,np2
real(8),allocatable :: m(:,:),u(:),r(:),uapro(:,:)
real(8) :: a, b,h,alpha,beta
real(8),parameter :: pi=dacos(-1.0d0),tolconj=1.0d-6
integer,parameter :: iteraconj=100
	np=50
	np2=np**2 !np elevado al cuadrado, ya que np significa el numero de puntos en una dimension
	alpha=0
	beta=1
	b=pi
	a=0
	h=(b-a)/(np-1)
	allocate(m(np2,np2),u(np2),r(np2),uapro(np,np))

	call invoca_matvect(m,r,h,np)
	call rverdadero(uapro,h,a,np)
	call grad_conj(m,r,u,tolconj,iteraconj)

	do i=1,np2
	cont =cont+1
	print *,  'Valores de u', u(i),'/////', cont
	end do	
	
	open(unit=12, file='resultado_exacto(aprox).txt')
	write(12,*) 'Numero puntos:', np
	do j=1,np
		do i=1,np
		write(12,*) uapro(j,i)
		end do
	end do
	close(11)

	open(unit=13, file='solucion_numerica.txt')
	write(13,*) 'Numero puntos:', np
	write(13,*) 'Numero de iteraciones gradiente conjugado', iteraconj
	write(13,*) 'Tolerancia gradiente conjugado', tolconj
	do i=1,np2
	write(13,*) u(i)
	end do
	close(11)

	cont=0
	do j=1,np
		do i=1,np
		cont =cont+1
		print *, 'val aprox',uapro(j,i),'/////', cont
		end do
	end do

	!open(unit=11, file='matriz.txt')
	!do i=1,np2
	!write(11,"(25I2)")int(m(i,:))
	!end do
	!close(11)
	
end program main
