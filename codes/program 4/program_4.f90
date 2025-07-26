program bairstow

implicit none
integer :: grau, ind, grau_original
real, allocatable :: ai(:), bi(:)
real :: r, s, tol
integer :: k
real, allocatable :: xr(:)

! Variáveis para varredura
real :: rmin, rmax, smin, smax, dr, ds, rr, ss
logical :: convergiu
integer :: i, n_iter
real, allocatable :: temp_ai(:), temp_bi(:)

write(*,*) "Digite o grau do polinômio (>=2):"
read(*,*) grau

if (grau < 2 ) then
    write(*,*) "Grau deve ser >= 2."
    stop
end if

grau_original = grau
allocate(ai(0:grau))
allocate(xr(0:grau - 1))  ! No máximo, grau raízes reais ou complexas

tol = 1.0e-6

write(*,*) "Entre com os coeficientes do polinômio:"
do ind = 0, grau
    if (ind == 0) then
        write(*,*) "Coeficiente livre a0:"
    else
        write(*,*) "Coeficiente a", ind, ":"
    end if
    read(*,*) ai(ind)
end do

write(*,*) "Entre com o valor inicial de r0:"
read(*,*) r
write(*,*) "Entre com o valor inicial de s0:"
read(*,*) s

k = 0
do while (grau > 2)
    allocate(bi(0:grau))
    call bairstow_sub(ai, tol, grau, r, s, xr(k), xr(k+1), bi, n_iter)
    grau = grau - 2
    ai(0:grau) = bi(2:grau+2)  ! Redução do polinômio
    deallocate(bi)
    k = k + 2
end do

! Tratamento do restante do polinômio
if (grau == 2) then
    xr(k)   = (-ai(1) + sqrt(ai(1)**2 - 4*ai(2)*ai(0))) / (2*ai(2))
    xr(k+1) = (-ai(1) - sqrt(ai(1)**2 - 4*ai(2)*ai(0))) / (2*ai(2))
    k = k + 2
else if (grau == 1) then
    xr(k) = -ai(0)/ai(1)
    k = k + 1
end if

! Impressão das raízes encontradas
write(*,*) "Raízes encontradas:"
do ind = 0, k - 1
    write(*,'(A,I2,A,F12.6)') "x(", ind, ") = ", xr(ind)
end do

! -------------------
! Parte adicional: Varredura para mapa
! -------------------
write(*,*) "\n--- Varredura para gerar mapa de convergência ---"
write(*,*) "Digite rmin, rmax, dr:" 
read(*,*) rmin, rmax, dr
write(*,*) "Digite smin, smax, ds:" 
read(*,*) smin, smax, ds

open(unit=20, file="mapa.dat", status="replace")

do ss = smin, smax, ds
    do rr = rmin, rmax, dr
        grau = grau_original
        allocate(temp_ai(0:grau))
        allocate(temp_bi(0:grau))
        temp_ai = ai
        r = rr
        s = ss
        convergiu = .true.

        do while (grau > 2 .and. convergiu)
            call bairstow_sub(temp_ai, tol, grau, r, s, xr(0), xr(1), temp_bi, n_iter)
            if (n_iter >= 5000) then
                convergiu = .false.
                exit
            end if
            grau = grau - 2
            temp_ai(0:grau) = temp_bi(2:grau+2)
        end do

        if (convergiu) then
            write(20,'(F7.3,1X,F7.3,1X,I5)') rr, ss, n_iter
        else
            write(20,'(F7.3,1X,F7.3,1X,I5)') rr, ss, -1
        end if

        deallocate(temp_ai)
        deallocate(temp_bi)
    end do
    write(20,*)   ! <<< linha em branco para separar blocos
end do

close(20)
write(*,*) "Arquivo 'mapa.dat' gerado com sucesso."


end program bairstow


!------------------------------------------------------------------

subroutine bairstow_sub(ai, tol, grau, r, s, x1, x2, bo, iter_out)
    implicit none

    integer, intent(in) :: grau
    real, intent(inout) :: r, s
    real, intent(in) :: tol
    real, intent(in) :: ai(0:grau)
    real, intent(out) :: x1, x2
    real, intent(out) :: bo(0:grau)
    integer, intent(out) :: iter_out
    integer :: i, iter, max_iter
    real, allocatable :: bi(:), ci(:)
    real :: deter, delr, dels, delta

    allocate(bi(0:grau))
    allocate(ci(0:grau))

    delr = 1.0
    dels = 1.0
    iter = 0
    max_iter = 5000

    do while (((abs(delr) >= tol) .or. (abs(dels) >= tol)) .and. (iter < max_iter))
        bi(grau) = ai(grau)
        bi(grau-1) = ai(grau-1) + r*bi(grau)

        ci(grau) = 0.0
        ci(grau-1) = bi(grau)

        do i = grau - 2, 0, -1
            bi(i) = ai(i) + r*bi(i+1) + s*bi(i+2)
            ci(i) = bi(i+1) + r*ci(i+1) + s*ci(i+2)
        end do

        deter = ci(1)**2 - ci(0)*ci(2)

        if (abs(deter) < 1.0e-12) then
            exit
        end if

        delr = (-bi(1)*ci(1) + bi(0)*ci(2)) / deter
        dels = (-bi(0)*ci(1) + bi(1)*ci(0)) / deter

        r = r + delr
        s = s + dels
        iter = iter + 1
    end do

    delta = r**2 + 4*s
    if (delta >= 0.0) then
        x1 = (r + sqrt(delta)) / 2
        x2 = (r - sqrt(delta)) / 2
    else
        x1 = r / 2
        x2 = sqrt(abs(delta)) / 2  ! Parte imaginária
    end if

    bo = bi  ! Retorna novo polinômio reduzido

    iter_out = iter

    deallocate(bi)
    deallocate(ci)
end subroutine bairstow_sub
