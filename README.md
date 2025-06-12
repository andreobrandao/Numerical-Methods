# Repositório de Métodos Numéricos - André Brandão

<p align="center">
  <img src="figures/logo-unb.png" width="600"/>
</p>



The present repository aims to serve as a place to hold all the fortran codes i write during the discipline of Numerical Methods in mechanical sciences.
Below it's explained how it's organized with the explanation of each code presented in the [code folder](./codes).

Along with the codes, in the [Analysis](./Analysis) folder is the pdf files with some analysis made with some of the codes developed.

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

### How to use

To run the program, run "ifx program_1.f90" and then run the executable file with "./a.out".

## Code 2 - Methods of the bissection and the false position

The numerical methods presented here aim to find the zero of functions. Both methods have some strong qualities and some defects. they're based in the idea that if you enter with a interval, and use its convergence method, this initial interval narrow around the root of the function studied.

In the bissection method the idea is that the interval is always narrowed by the half of the previous interval. then it's checked in which new interval the root is, then the process repeats untill you reach the tolerence you want.
Considering $x_l$ the lower value of the interval and $x_u$ the upper value, then the new value will be:

$$
    x_m = \frac{x_l + x_u}{2}
$$

then it's checked if $f(x_m)*f(x_l)\le 0$, if yes then $x_u$ is updated to $x_m$, if not then $x_u$ is updated to $x_m$, and the process repeats.

The false position method takes in consideration the proximity of the function analyzed in the values of the interval, and narrow this interval considering this proximity.
The new value of the interval is defined using the next equation:

$$
    x_m = x_u - \frac{f(x_u)\cdot (x_l - x_u)}{f(x_l) - f(x_u)}
$$

then just like the previous method it's checked if $f(x_m)*f(x_l)\le 0$, if yes then $x_u$ is updated to $x_m$, if not then $x_u$ is updated to $x_m$, and the process repeats.

### How to use

To run the program, run "ifx program_2.f90" and then run the executable file with "./a.out".

## Code 3 - Methods of the Secant and Müller
The third program written in fortran makes a comparisson between the secant method and the Müller method to find the roots of a given polynomial. This comparisson is given in terms of the relative error in % to every iteration untill the defined tolerance is achieved.
The secant method follows the idea of the newton-raphson method, with the difference that it uses two different points from the function or available tabeled numbers and its f(x) values to approximate the root of the studied equation using a line. The following equation shows this approximation:

$$
  x_{i+1}=x_i - \frac{f(x_i)\cdot(x_{i-1} - x_i)}{f(x_{i-1}) - f(x_i)}
$$

The Müller method makes the approximation to the root of the given function defining another curve that coincides with the studied function in three points (initial guesses). The root is approximated by the following relation:

$$
  x_{i+1} = x_i - \frac{2c}{b \pm \sqrt{b^2 - 4ac}}
$$

where:

$$
  a = \frac{\delta_1 - \delta_0}{h_1 - h_o},\ b = ah_1 + \delta_1,\ c = f(x_2) 
$$

and:

$$
  h_0 = x_1 - x_0,\ h_1 = x_2-x_1,\ \delta_0 = \frac{f(x_1) - f(x_0)}{x_1 - x_0},\ \delta_1 = \frac{f(x_2) - f(x_1)}{x_2 - x_1}
$$

The following figure shows the difference between the two methods:

![Comparisson between the Secant and Müller method](figures/secant_muller.png)
  
*Figura 2 – Comparisson between the Secant and Müller method.*

### How to use

To run the program, run "ifx program_3.f90" and then run the executable file with "./a.out". The file "grafico.plt" has the gnuplot commands to plot the graph "Relative error (%) x iteration". So you must have the gnuplot installed in your pc.

The usual output graph of this code is shown below:

![Relavtive error x iteration](figures/RE_code3.png)
  
*Figura 3 – Relavtive error x iteration.*

## Code 4 - Method of Bairstow

The code 4 implements Bairstow's method to find the roots of any polynomial function. The method is based on dividing the polynomial by a second-degree polynomial of the form:

$$
f(x)=x^2 -r\cdot x - s
$$

The goal is to adjust the values of r and s until the remainder of the division is zero. This adjustment is performed using the Newton-Raphson method.

Once the division yields no remainder, the quadratic factor is considered a root (or pair of roots, in the case of complex conjugates), and the remaining quotient becomes the new polynomial to analyze. This process repeats until all roots of the original polynomial are found.

Depending on the initial guesses for r and s, the method may fail to converge to a solution. To explore how different initial values affect convergence, this program also generates a fractal map that visualizes the convergence behavior of Bairstow's method over a selected range of r and s values.

A fractal, in this context, is a colorful graphical map where each point represents an initial guess (r,s), and the color indicates the number of iterations required to converge. Points that do not converge within a set number of iterations are typically colored black or left blank.

These fractals help to visualize the sensitivity of the method to the initial values and reveal intricate structures in the convergence landscape.

Examples of such fractals are shown below:

<p align="center">
  <img src="figures/mapa_1.png" width="45%" />
  <img src="figures/mapa_2.png" width="45%" />
</p>
*Figura 4 – Example of fractal maps to polynomials of the 3th and 5th degree*

## How to use

To run the program, run "ifx program_4.f90" and then run the executable file with "./a.out". The file "mapa.gnu" has the gnuplot commands to plot the fractal map. So after the file "mapa.dat" is generated run "gnuplot -persist mapa.gnu" to create the map.

When you run the program it will be asked to enter with the coefficients of the polynomial, the initial guess of r and s, and then a range of r and s to create the fractal. Inside the file "mapa.gnu", there are two lines written set xrange [-100:100] and set yrange [-100:100], what means that it's set to create the map in this range. If you choose another range in the fortran code then adjust these lines in the gnuplot file so it will run as it should.
