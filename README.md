# Repositório de Métodos numéricos
### André Brandão

The present repository aims to serve as a place to hold all the fortran codes i write during the discipline of Numerical Methods in mechanical sciences.
Below it's explained how it's organized with the explanation of each code presented in the [code folder](./codes).

Along with the codes, in the [Analysis](./Analysis) folder is the pdf files with some analysis made with every code developed.

## Code 1 - 4th order Runge kutta method

A sphere sedimenting in a fluid is studied here. Considering the drag force the second law of newton applied to this problem stays: 

$$
    m_p \frac{dv_z}{dt}=-6\pi\eta av_z-\frac{9}{4}\pi\rho_ra^2v_z^2+\frac{4}{3}\pi a^3\Delta\rho g
$$

or in its adimensional version:

$$
    St\frac{dv_z^{\ast}}{dt^{\ast}} = -v_z^{\ast}-\frac{3}{8}Re_s v_z^{2^\ast}+1
$$

Its version without the drag version as shown below:

$$
    St\frac{dv_z^{\ast}}{dt^{\ast}} = 1-v_z^{\ast}
$$

This code aims to calculate the velocity of the sphere in this process of sedimentation using the 4th order runge kutta method.

The selected lines in the code, showed below, shows the idea of this method:

```fortran
    do i = 1,n ! iterando no passo de tempo
        call rk4(v_star, t_star, dt, St)
        t_star = t_star + dt  ! Avança no tempo
    function dvdt(v, St) result(dv)
        real(8), intent(in) :: v, St
        real(8) :: dv
        real(8) :: Re_s
        Re_s = 0.001
        ! funcoes a serem consideradas
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
```

An example of comparisson of the results using this method and the analitical result is shown below:

![temporal development of the adimensional velocity to St=1.005](figures/St1005.PNG)
  
*Figura 1 – temporal development of the adimensional velocity to St=1.005.*
## Code 2 - Methods of the bissection and the false position

The numerical methods presented here aim to find the zero of functions. Both methods have some strong qualities and some defects. they're based in the idea that if you enter with a interval, and use its convergence method, this initial interval narrow around the root of the function studied.

In the bissection method the idea is that the interval is always narrowed by the half of the previous interval. then it's checked in which new interval the root is, then the process repeats untill you reach the tolerence you want.
Considering $x_l$ the lower value of the interval and $x_u$ the upper value, then the new value will be:
$$
    x_m = \frac{x_l+x_u}{2}
$$
then it's checked if $f(x_m)*f(x_l)\le 0$, if yes then $x_u$ is updated to $x_m$, if not then $x_u$ is updated to $x_m$, and the process repeats.
