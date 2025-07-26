program calor_1d_geracao

  implicit none
  !  integer, parameter :: n = 150            ! Número de nós
  integer, parameter :: n = 50            ! Número de nós
  integer :: i, j, nt
  real(8), parameter :: L = 0.01d0        ! Comprimento [m]
  real(8), parameter :: tfinal = 60.0d0   ! tempo total [s]
  real(8), parameter :: alpha = 5d-6      ! Difusividade térmica [m²/s]

  real(8), parameter :: k = 30.0d0        ! condutividade [W/m·K]
  real(8), parameter :: q = 1d7          ! Geração interna [W/m³]
  real(8), parameter :: h = 1100.0d0      ! coef. de convecção [W/m²·K]
  real(8), parameter :: Tinf = 523.15d0   ! Temperatura do fluido [K]
  real(8), parameter :: T0 = 298.15d0     ! temperatura inicial [K]
  real(8), parameter :: rhocp = k / alpha ! rho * cp

  real(8) :: dx, dt, Fo, Bi, Aa 
  real(8), dimension(n) :: T, Tnew
  real(8), allocatable :: a(:), b(:), c(:), d(:)
  real(8), allocatable :: Tmap(:,:) ! Matriz para armazenar T(x,t) e plotar depois

  allocate(a(n), b(n), c(n), d(n)) 
  ! acima aloca vetores para o método de Thomas, aloca pra o tamanho da matriz tridiagonal dos coeficientes que tem n elementos (quantidade de nós)
  
  dx = L / (n - 1) ! Espaçamento entre nós
  dt = 0.5d0 * dx**2 / alpha ! Passo de tempo baseado na relação de fourier
  nt = int(tfinal / dt) + 1 ! Número de passos de tempo

  Fo = alpha * dt / dx**2 ! Número de Fourier
  Bi = h * dx / k ! Número de Biot
  Aa = q * dt / rhocp ! Aa aqui é o termo de geração dividido por rho*cp

  allocate(Tmap(n, nt))

  ! Condiçãp inicial
  T = T0
  Tmap(:, 1) = T

   ! Aqui resolve pelo método de Thomas
  do j = 2, nt
    ! primeiro nó, como especificado na aula 9
     a(1) = 0.0d0
     b(1) = 1.0d0 + 2.0d0 * Fo
     c(1) = -2.0d0 * Fo
     d(1) = T(1) + Aa
      ! para os nós internos
     do i = 2, n-1
        a(i) = -Fo
        b(i) = 1.0d0 + 2.0d0 * Fo
        c(i) = -Fo
        d(i) = T(i) + Aa
     end do

    ! para o último nó, como especificado na aula 9
     a(n) = -2.0d0 * Fo
     b(n) = 1.0d0 + 2.0d0 * Fo + 2.0d0 * Fo * Bi
     c(n) = 0.0d0
     d(n) = T(n) + Fo * (q * dx**2 / k) + 2.0d0 * Fo * Bi * Tinf
     ! write(*,*) 'd(n) = ', d(n)

     call thomas(a, b, c, d, Tnew, n)
      ! write(*,*) 'Tnew = ', Tnew
     T = Tnew
     Tmap(:, j) = T
  end do

  ! Salva resultado final
  open(unit=10, file="temperatura_final.dat", status="unknown")
  do i = 1, n
     write(10,*) (i - 1) * dx, T(i)
  end do
  close(10)

  ! Salva matriz T(x,t), uma linha para cada passo de tempo
  ! eu to plotando no tempo o mapa de todas as temperaturas em cada nó, ou seja, ao invés de plotar T(x) para cada t, eu ploto T(x,t)
  ! assim fica mais fácil visualizar a evolução temporal
  open(unit=20, file="mapa_temperatura.dat", status="unknown")
  do j = 1, nt
     write(20,'(50F12.5)') (Tmap(i,j), i = 1, n)
  end do
  close(20)

  print *, "Simulação finalizada."
  print *, "Arquivo 'temperatura_final.dat' contém o perfil final."
  print *, "Arquivo 'mapa_temperatura.dat' contém T(x,t) para visualização temporal."

end program calor_1d_geracao

subroutine thomas(a, b, c, d, x, n) ! o método de thomas para resolver sistemas tridiagonais
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: a(n), b(n), c(n)
  real(8), intent(inout) :: d(n)
  real(8), intent(out) :: x(n)
  real(8), dimension(n) :: c_, d_
  integer :: i

  c_(1) = c(1) / b(1)
  d_(1) = d(1) / b(1)

  do i = 2, n
     c_(i) = c(i) / (b(i) - a(i) * c_(i-1))
     d_(i) = (d(i) - a(i) * d_(i-1)) / (b(i) - a(i) * c_(i-1))
  end do

  x(n) = d_(n)
  do i = n-1, 1, -1
     x(i) = d_(i) - c_(i) * x(i+1)
  end do
end subroutine thomas
