module bairstow_module
    contains
        subroutine bairstow(ai, grau, r0, s0, tol, x1, x2, bi_out, convergiu)
            implicit none
            integer, intent(in) :: grau
            real, intent(inout) :: r0, s0
            real, intent(in) :: tol
            real, intent(in) :: ai(0:grau)
            real, intent(out) :: x1, x2
            real, intent(out) :: bi_out(0:grau)
            logical, intent(out) :: convergiu
            real :: delr, dels, delta, det
            real, allocatable :: bi(:), ci(:)
            integer :: j, iter, max_iter
            
            allocate(bi(0:grau))
            allocate(ci(0:grau))
            
            delr = 1.0
            dels = 1.0
            iter = 0
            max_iter = 100
            convergiu = .false.
            
            do while (((abs(delr) >= tol) .or. (abs(dels) >= tol)) .and. (iter < max_iter))
                bi(0) = ai(0)
                ci(0) = bi(0)
                
                if (grau >= 1) then
                    bi(1) = ai(1) + r0 * bi(0)
                    ci(1) = bi(1) + r0 * ci(0)
                end if
                
                do j = 2, grau
                    bi(j) = ai(j) + r0 * bi(j-1) + s0 * bi(j-2)
                    ci(j) = bi(j) + r0 * ci(j-1) + s0 * ci(j-2)
                end do
                
                if (grau < 2) exit
                
                det = ci(grau-1)**2 - ci(grau-2)*ci(grau)
                if (abs(det) < 1e-12) then
                    write(*,*) "Determinante muito pequeno - parando iterações"
                    exit
                end if
                
                delr = (-bi(grau-1)*ci(grau-1) + bi(grau-2)*ci(grau)) / det
                dels = (-bi(grau)*ci(grau-2) + bi(grau-1)*ci(grau-1)) / det
                
                r0 = r0 + delr
                s0 = s0 + dels
                iter = iter + 1
            end do
            
            if (iter < max_iter) convergiu = .true.
            
            delta = r0**2 + 4.0 * s0
            if (delta >= 0.0) then
                x1 = (r0 + sqrt(delta)) / 2.0
                x2 = (r0 - sqrt(delta)) / 2.0
            else
                write(*,*) "Raízes complexas encontradas."
                x1 = 0.0
                x2 = 0.0
            end if
            
            bi_out = bi
            
            deallocate(bi)
            deallocate(ci)
        end subroutine bairstow
end module bairstow_module

!#####################################################################
program code4
    use bairstow_module
    implicit none
    
    integer :: n, i, grau, k, num_raizes
    real :: r0, s0, tol, rinit, sinit
    real, allocatable :: ai(:), bi(:), xr(:)
    logical :: convergiu
    
    write(*,*) "Grau do polinômio:"
    read(*,*) n
    
    allocate(ai(0:n))
    allocate(bi(0:n))
    allocate(xr(n+2))
    
    do i = 0, n
        if (i == 0) then
            write(*,*) "Coef. do termo x^0 (constante):"
        else
            write(*,*) "Coef. do termo x^", i, ":"
        end if
        read(*,*) ai(i)
    end do
    
    write(*,*) "Chute inicial r0:"
    read(*,*) rinit
    write(*,*) "Chute inicial s0:"
    read(*,*) sinit
    
    tol = 1.0e-6
    grau = n
    k = 1
    
    do while (grau > 2)
        r0 = rinit
        s0 = sinit
        call bairstow(ai(0:grau), grau, r0, s0, tol, xr(k), xr(k+1), bi(0:grau), convergiu)
        
        if (.not. convergiu) then
            write(*,*) "Método não convergiu para grau", grau
            exit
        end if
        
        grau = grau - 2
        ai(0:grau-2) = bi(2:grau)
        k = k + 2
    end do
    
    if (grau == 2) then
        call raiz_quad(ai(2), ai(1), ai(0), xr(k), xr(k+1))
        num_raizes = k + 1
    else if (grau == 1) then
        if (abs(ai(1)) > 1e-6) then
            xr(k) = -ai(0)/ai(1)
        else
            write(*,*) "Divisão por zero em polinômio linear"
            xr(k) = 0.0
        end if
        num_raizes = k
    else
        num_raizes = k - 1
    end if
    
    write(*,*) "Raízes encontradas:"
    do i = 1, num_raizes
        write(*,'(A,I2,A,F12.6)') "Raiz ", i, " = ", xr(i)
    end do
    
    deallocate(ai)
    deallocate(bi)
    deallocate(xr)
    
contains
    
    subroutine raiz_quad(a, b, c, x1, x2)
        real, intent(in) :: a, b, c
        real, intent(out) :: x1, x2
        real :: delta
        
        delta = b**2 - 4.0*a*c
        if (delta >= 0.0) then
            x1 = (-b + sqrt(delta))/(2.0*a)
            x2 = (-b - sqrt(delta))/(2.0*a)
        else
            write(*,*) "Raízes complexas no polinômio quadrático"
            x1 = 0.0
            x2 = 0.0
        end if
    end subroutine raiz_quad
    
end program code4