module contorno
implicit none
contains


subroutine invoca_matvect(m,r,h,np)
	integer, intent(in) :: np
	real(8), intent(in) :: h
	real(8), intent(inout) :: m(:,:),r(:)
	integer :: i,np2
	m=0.0d0
	np2=np**2
	m=0
	do i=1,np
	m(i,i)=1.0d0
	m(np2-i+1,np2-i+1)=1.0d0
	r(i)=0.0d0
	r(np2-i+1)=0.0d0
	end do
	
	do i=np+1,np2,np
	m(i-1,i-1)=1.0d0
	m(i,i)=1.0d0
	r(i-1)=0.0d0
	r(i)=0.0d0
	end do

	do i=np,np2-np
	if (mod(i,np)/=0.AND.mod(i-1,np)/=0) then
		m(i,i-np)=-1.0d0
		m(i,i-1)=-1.0d0
		m(i,i)=4.0d0
		m(i,i+1)=-1.0d0
		m(i,i+np)=-1.0d0
		r(i)=h**2
	end if
	end do	
	
end subroutine invoca_matvect

subroutine rverdadero(uapro,h,a,np)
integer, intent(in) :: np
real(8), intent(in) :: h,a
real(8), intent(inout) :: uapro(:,:)
real(8),dimension(np) :: x,y
integer :: xcont, ycont, i, j
real(8),parameter :: pi=dacos(-1.0d0)
	do i=1,np
	x(i)=a+(real(i)-1.0d0)*h
	y(i)=a+(real(i)-1.0d0)*h
	end do
		
	uapro=0.0d0
	do xcont=1,np
	do ycont=1,np
		do j=1,100,2
			do i=1,100,2
	uapro(xcont,ycont)=uapro(xcont,ycont)+&
	&(16.0d0/((real(j)**2+real(i)**2)*real(i)*real(j)*pi**2))*dsin(real(j)*x(xcont))*dsin(real(i)*y(ycont))
			end do
		end do
	end do
	end do

end subroutine rverdadero





end module contorno
