module contorno
implicit none
contains


subroutine invoca_matvect(m,r,h,np) !Invoca la matriz y el vector resultado del problema de poisson
	integer, intent(in) :: np
	real(8), intent(in) :: h
	real(8), intent(inout) :: m(:,:),r(:)
	integer :: i,np2
	!Declaracion de variables
	m=0.0d0				!Puesta de la matriz principal a 0
	np2=np**2
	do i=1,np			!Este bucle coloca en la matriz los resultados que conocemos correspondientes al contorno
	m(i,i)=1.0d0			!Colocaría los valores correspondientes a la condicion que nos dan de la primera fila y ultima
	m(np2-i+1,np2-i+1)=1.0d0	!Escribe por encima y por abajo los 1 en la diagonal, para los que u es trivial
	r(i)=0.0d0
	r(np2-i+1)=0.0d0		!De igual manera este bucle nos da la condicion de que en los extremos el resultado es 0
	end do				
	
	do i=np+1,np2,np 		!Este bucle de igual manera nos coloca los otros vectores que conocemos, correspondiente
	m(i-1,i-1)=1.0d0		!a los extremos que no estan situados en la primera fila ni en la ultima, los intermedios
	m(i,i)=1.0d0
	r(i-1)=0.0d0
	r(i)=0.0d0
	end do

	do i=np,np2-np
	if (mod(i,np)/=0.AND.mod(i-1,np)/=0) then	!Este bucle condicionado coloca los valores no triviales que debemos resolver
		m(i,i-np)=-1.0d0			!La condicion es que no esté en un extremo, de los que ya conocemos
		m(i,i-1)=-1.0d0
		m(i,i)=4.0d0
		m(i,i+1)=-1.0d0				!Se colocan filas de estos valores que vienen de despejar al resolver la
		m(i,i+np)=-1.0d0			!ecuacion empleando la derivacion numérica
		r(i)=h**2				!El vector resultado sale h² de multiplicar todo por h²
	end if
	end do	
	
end subroutine invoca_matvect

subroutine rverdadero(uapro,h,a,np) 			!Esta subrutina es una subrutina de validacion, se trata simplemente de
integer, intent(in) :: np				!una subrutina que nos da una solucion muy cercana a la exacta, ya que
real(8), intent(in) :: h,a				!la exacta nos la han dado en forma de suma infinita.
real(8), intent(inout) :: uapro(:,:)
real(8),dimension(np) :: x,y
integer :: xcont, ycont, i, j
real(8),parameter :: pi=dacos(-1.0d0)
	do i=1,np
	x(i)=a+(real(i)-1.0d0)*h			!Sacamos el vector x e y que tienen valores equiespaciados
	y(i)=a+(real(i)-1.0d0)*h		
	end do
		
	uapro=0.0d0
	do xcont=1,np
	do ycont=1,np
		do j=1,100,2
			do i=1,100,2			!Realizamos la suma de la funcion 100 veces ya que el 100 es un numero  
	uapro(xcont,ycont)=uapro(xcont,ycont)+&		!bastante grande y guay, lo hacemos para cada x e y (xcont,ycont)
	&(16.0d0/((real(j)**2+real(i)**2)*real(i)*real(j)*pi**2))*dsin(real(j)*x(xcont))*dsin(real(i)*y(ycont))
			end do
		end do
	end do
	end do

end subroutine rverdadero

subroutine invoca_matvectdif(m,r,x,np,h,a,alpha,beta)		!Crea la matriz y el vector que usaremos para solucionar la 
	integer, intent(in) :: np				!diferencial ordinaria
	real(8), intent(in) :: h,alpha,beta,a
	real(8), intent(inout) :: m(:,:),r(:),x(:)
	integer :: i
	m=0.0d0							
	m(1,1)=1.0d0
	m(np,np)=1.0d0						!Introducimos en la matriz los 1 correspondientes a nuestro contorno
	
	do i=1,np
	x(i)=a+(real(i)-1.0d0)*h 				!Generamos el vector x equiespaciado para resolver el problema
	end do

	do i=2,np-1						
	m(i,i-1)=(x(i)*x(i))-x(i)*(h/2.0d0)			!Establecemos las ecuaciones correspondientes a haber hecho la 
	m(i,i)=-2.0d0*(x(i)*x(i))+((x(i)*x(i))*(h*h))		!derivación numérica, que contiene las relaciones correspondientes
	m(i,i+1)=(x(i)*x(i))+x(i)*(h/2.0d0)
	r(i)=0							!Como nos dicen que el vector resultado es 0 pues lo imponemos
	end do
	r(1)=alpha
	r(np)=beta						!Alpha y beta son los dos valores que nos dicen
end subroutine invoca_matvectdif




end module contorno
