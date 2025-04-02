program sedimentation
    implicit none
    integer, parameter :: n = 1000  ! Número de passos de tempo
    real(8), parameter :: t_max = 250.0d0  ! Tempo final adimensional
    real(8), parameter :: dt = t_max / n  ! Discretizando passo de tempo
    real(8), parameter :: St_values(3) = (/ 1.25d0, 1.005d0, 10.0d0 /)  ! Diferentes valores de Stokes
    real(8) :: v_star,v_i, t_star, v_analytical, St,epsilon_v,v(n) !demais variáveis
    integer :: i, j !iterar

    !Loop sobre diferentes valores de St
    do j = 1, size(St_values)
        St = St_values(j)
        v_star = 0.0d0  !Condição inicial: v*(0) = 0, Re=0
        t_star = 0.0d0
        v_i = 1.0d0
        epsilon_v=3/80

        ! Abrir arquivo para salvar resultados
        open(unit=10, file='results_St' // trim(adjustl(num2str(St))) // '.dat', status='replace')
        write(10,*) "t_star, v_star_numeric,v(i), v_star_analytical"

        !loop de integração numérica usando RK4
        do i = 1, n
            v_analytical = 1.0d0 - exp(-t_star / St)  !Solução exata
            v(i) = v_i * (1.0d0 + (1.0d0 + 2.0d0 * epsilon_v) / &
                ((1.0d0 - epsilon_v) * EXP((1.0d0 + 2.0d0 * epsilon_v) * t_star / St) - epsilon_v))
            write(10,*) t_star, v_star,2-v(i), v_analytical  !Salvar no arquivo
            
            ! chama Runge-Kutta 4ª ordem
            call rk4(v_star, t_star, dt, St)

            t_star = t_star + dt  ! Avança no tempo
        end do
        close(10)
    end do
contains

    !fyunção para a equação diferencial: dv/dt = (1 - v) / St
    !função para a equação diferencial : dv/dt = (1- v- (3/8)*Re_s*v**2)/St
    function dvdt(v, St) result(dv)
        real(8), intent(in) :: v, St
        real(8) :: dv
        real(8) :: Re_s
        Re_s = 0.00000000001
        !dv = (1.0d0-v)/St
        dv = (1- v- (3/8)*Re_s*v**2)/St
    end function dvdt

    ! Rotina de Runge-Kutta de quarta ordem
    subroutine rk4(v, t, dt, St)
        real(8), intent(inout) :: v  ! velocidade adimensional
        real(8), intent(in) :: t, dt, St
        real(8) :: k1, k2, k3, k4

        k1 = dt * dvdt(v, St)
        k2 = dt * dvdt(v + 0.5d0 * k1, St)
        k3 = dt * dvdt(v + 0.5d0 * k2, St)
        k4 = dt * dvdt(v + k3, St)

        v = v + (k1 + 2.0d0*k2 + 2.0d0*k3 + k4) / 6.0d0
    end subroutine rk4

    ! Função auxiliar para converter número para string
    function num2str(x) result(str)
        real(8), intent(in) :: x
        character(20) :: str
        write(str, '(F6.2)') x
    end function num2str

end program sedimentation
