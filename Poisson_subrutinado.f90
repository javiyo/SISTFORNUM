program main 
use sistemas
use contorno
implicit none
integer :: i,j,k,cont,np,np2
real(8),allocatable :: m(:,:),u(:),r(:),uapro(:,:)
real(8) :: a, b,h,alpha,beta
real(8),parameter :: pi=dacos(-1.0d0)
	np=20
	np2=np**2 !np elevado al cuadrado, ya que np significa el numero de puntos en una dimension
	alpha=0
	beta=1
	b=pi
	a=0
	h=(b-a)/(np-1)
	allocate(m(np2,np2),u(np2),r(np2),uapro(np,np))

	call invoca_matvect(m,r,h,np)
	call rverdadero(uapro,h,a,np)
	call grad_conj(m,r,u,1.0d-6,100)

	do i=1,np2
	cont =cont+1
	print *,  'Valores de u', u(i),'/////', cont
	end do	
	cont=0
	do j=1,np
		do i=1,np
		cont =cont+1
		print *, 'val aprox',uapro(j,i),'/////', cont
		end do
	end do

	open(unit=11, file='matriz.txt')
	do i=1,np2
	write(11,"(25I2)")int(m(i,:))
	end do
	close(11)
	
end program main
