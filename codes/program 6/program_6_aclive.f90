program aclive_maximo
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: k, max_iter, unit
    real(dp) :: x, y, alpha, fx, fy, grad_norm,x_con,y_con
    real(dp), parameter :: tol = 1.0e-8_dp
    character(len=30) :: arq_saida

    ! Entrando com valores iniciais
    max_iter = 100 ! Número máximo de iterações (mais que suficiente para convergir mas tem uma condição de parada abaixo)
    x = 0.0_dp
    y = 0.0_dp

    ! Abrindo arquivo para salvar trajetória
    arq_saida = "trajetoria.txt"
    open(newunit=unit, file=arq_saida, status="replace", action="write")

    do k = 1, max_iter
        ! Gradientes de f em (x, y)
        fx = 4.0_dp + 2.0_dp*x - 8.0_dp*x**3 + 2.0_dp*y
        fy = 2.0_dp + 2.0_dp*x - 6.0_dp*y

        ! Verifica norma do gradiente (condição de parada pelo gradiente, se for pequeno o suficiente, considera que convergiu)
        grad_norm = sqrt(fx**2 + fy**2)
        ! write(*,*) 'Iteração:', k, 'Gradiente:', grad_norm
        if (grad_norm < tol) exit

        ! Busca o melhor alpha ( o alpha aqui é o h, passo do aclive máximo )
        call busca_alpha(x, y, fx, fy, alpha)

        ! Imprime na tela
        print '(A,I2,A,F10.5,A,F10.5,A,F10.5)', 'Iteração ', k,': x = ', x, ' y = ', y, ' f(x,y) = ', f(x,y)

        ! Salva no arquivo
        write(unit,'(F12.6,1X,F12.6,1X,F12.6)') x, y, f(x,y)

         ! Atualiza x e y(por ultimo pra captar o ponto inicial   )
        x_con=x
        y_con=y
        x = x + alpha * fx
        y = y + alpha * fy
        if (abs(x - x_con) < tol .and. abs(y - y_con) < tol) exit ! condição de parada adicional ( verifica que convergiu )
    end do

    close(unit)
    print *, 'Trajetória salva em ', trim(arq_saida)

contains

    function f(x, y) result(val)
        real(dp), intent(in) :: x, y
        real(dp) :: val
        val = 4.0_dp*x + 2.0_dp*y + x**2 - 2.0_dp*x**4 + 2.0_dp*x*y - 3.0_dp*y**2
    end function

    subroutine busca_alpha(x, y, fx, fy, alpha_otimo) ! método bruto, uso de um passo pequeno pra achar o alpha ótimo
        real(dp), intent(in) :: x, y, fx, fy
        real(dp), intent(out) :: alpha_otimo
        real(dp) :: alpha, fval, fmax
        real(dp), parameter :: a_min = 0.0_dp, a_max = 1.0_dp, passo = 0.001_dp

        fmax = -1.0e30_dp
        alpha_otimo = 0.0_dp

        do alpha = a_min, a_max, passo
            fval = f(x + alpha * fx, y + alpha * fy)
            if (fval > fmax) then
                fmax = fval
                ! write(*,*) 'Alpha:', alpha, 'fval:', fval
                alpha_otimo = alpha
            end if
        end do
    end subroutine

end program aclive_maximo
