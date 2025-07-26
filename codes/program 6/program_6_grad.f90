program gradiente_conjugado
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: i, max_iter
  real(dp) :: x, y, fx, fy, dx, dy, alpha, g2_old, g2_new, beta
  real(dp) :: tol, norm_grad

  max_iter = 100
  tol = 1.0e-6_dp

  ! Condição inicial
  x = 0.0_dp
  y = 0.0_dp

  call gradiente(x, y, fx, fy)
  dx = -fx
  dy = -fy
  g2_old = fx*fx + fy*fy

  do i = 1, max_iter
     alpha = busca_alpha(x, y, dx, dy)

     ! Atualiza x e y
     x = x + alpha * dx
     y = y + alpha * dy

     ! Calcula novo gradiente
     call gradiente(x, y, fx, fy)
     g2_new = fx*fx + fy*fy
     norm_grad = sqrt(g2_new)

     print '(A,I3,A,F10.5,A,F10.5,A,F10.5)', 'Iteração ', i, ': x = ', x, ' y = ', y, ' f(x,y) = ', f(x,y)

     if (norm_grad < tol) exit

     ! Atualiza direção com Fletcher-Reeves
     beta = g2_new / g2_old
     dx = -fx + beta * dx
     dy = -fy + beta * dy
     g2_old = g2_new
  end do

contains

  function f(x, y) result(val)
    real(dp), intent(in) :: x, y
    real(dp) :: val
    val = 4.0_dp*x + 2.0_dp*y + x**2 - 2.0_dp*x**4 + 2.0_dp*x*y - 3.0_dp*y**2
  end function f

  subroutine gradiente(x, y, fx, fy)
    real(dp), intent(in) :: x, y
    real(dp), intent(out) :: fx, fy
    fx = 4.0_dp + 2.0_dp*x - 8.0_dp*x**3 + 2.0_dp*y
    fy = 2.0_dp + 2.0_dp*x - 6.0_dp*y
  end subroutine gradiente

  function busca_alpha(x, y, dx, dy) result(alpha)
    real(dp), intent(in) :: x, y, dx, dy
    real(dp) :: alpha, dalpha, f0, df, alpha_old
    integer :: j, max_inner
    real(dp) :: eps

    alpha = 0.0_dp
    max_inner = 50
    eps = 1.0e-8_dp

    do j = 1, max_inner
       f0 = f(x + alpha*dx, y + alpha*dy)
       df = derivada_f_alpha(x, y, dx, dy, alpha)
       dalpha = derivada_segunda_f_alpha(x, y, dx, dy, alpha)

       if (abs(df) < eps) exit

       alpha_old = alpha
       alpha = alpha - df / dalpha
       if (abs(alpha - alpha_old) < eps) exit
    end do
  end function busca_alpha

  function derivada_f_alpha(x, y, dx, dy, alpha) result(df)
    real(dp), intent(in) :: x, y, dx, dy, alpha
    real(dp) :: df, xp, yp, fxp, fyp

    xp = x + alpha * dx
    yp = y + alpha * dy
    call gradiente(xp, yp, fxp, fyp)
    df = fxp * dx + fyp * dy
  end function derivada_f_alpha

  function derivada_segunda_f_alpha(x, y, dx, dy, alpha) result(d2f)
    real(dp), intent(in) :: x, y, dx, dy, alpha
    real(dp) :: d2f, xp, yp, dxx, dxy, dyy

    xp = x + alpha * dx
    yp = y + alpha * dy

    ! Segunda derivadas parciais
    dxx = 2.0_dp - 24.0_dp*xp**2
    dxy = 2.0_dp
    dyy = -6.0_dp

    ! d²f/dα² = d^T * H * d
    d2f = dx**2 * dxx + 2.0_dp * dx * dy * dxy + dy**2 * dyy
  end function derivada_segunda_f_alpha

end program gradiente_conjugado
