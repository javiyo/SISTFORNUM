module sistemas
implicit none
contains


        subroutine grad_conj(A,b,x,tol,nmax)
            implicit none
            real(8), intent(in) :: A(:,:), b(:),tol
	    integer, intent(in) :: nmax
            real(8), intent(inout) :: x(:)
            real(8), allocatable :: r(:), v(:), r1(:), v1(:), t(:,:)
            real(8) :: tao, s
            integer :: n,i !iteraciones
            n=size(x)
            allocate(r(n))
            allocate(r1(n))
            allocate(v(n))
            allocate(v1(n))
            allocate(t(n,1))
            r1=matmul(A,x)-b !residuo inicial
            v=-r1            !direccion inicial de descenso
            do i=1,nmax
                t= reshape(v,[1,n]) !vector traspuesto a v para el calculo del paso(tao)
                tao = real(norma(r1)**2)/real(norma(matmul(matmul(t,A),v))) ! calculo del paso
                x= x + tao*v !calculo del nuevo itinerante para cada ciclo
                r=r1 !cambio de variable ya que al calcular el parametro s, se necesitan tanto r1 como su anterior; r
                r1 = r + tao*matmul(A,v) !calculo del nuevo residuo
                if (norma(r1)<tol*norma(b)) EXIT ! condicion puesta para asegurar exactitud en los resultados
                s=(norma(r1)**2)/(norma(r)**2) ! calculo del parametro s
                v=(-1)*r1 + (s*v) !calculo de la nueva direccion de descenso
            end do
        end subroutine grad_conj
        function norma (v)
                    implicit none
                    integer:: n,i
                    real(8):: v(:), norma
                    n=size(v)
                    norma=0.0
                    do i=1,n
                        norma=(v(i)**2)+norma
                    end do
                    norma=sqrt(norma)
end function norma

end module sistemas
