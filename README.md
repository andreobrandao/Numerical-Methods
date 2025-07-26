
# Numerical Methods Repository - André Brandão

<p align="center">
  <img src="figures/logo-unb.png" width="600" />
</p>

This repository serves as a collection of all the Fortran codes I wrote during the Numerical Methods course in Mechanical Sciences.  
Below is an overview of the folder structure and explanations for each code located in the [codes folder](./codes).

---

## Code 1 - 4th Order Runge-Kutta Method

This code studies a sphere sedimenting in a fluid. Considering the drag force, Newton’s second law applied to this problem reads:

$$
m_p \frac{dv_z}{dt} = -6\pi \eta a v_z - \frac{9}{4} \pi \rho_r a^2 v_z^2 + \frac{4}{3} \pi a^3 \Delta \rho g
$$

Or in its dimensionless form:

$$
St \frac{dv_z^{\ast}}{dt^{\ast}} = -v_z^{\ast} - \frac{3}{8} Re_s (v_z^{\ast})^{2} + 1
$$

The drag-free version is:

$$
St \frac{dv_z^{\ast}}{dt^{\ast}} = 1 - v_z^{\ast}
$$

This code calculates the velocity of the sphere during sedimentation using the 4th order Runge-Kutta method.

The core part of the implementation is shown below:

```fortran
do i = 1, n  ! Iterate over time steps
    call rk4(v_star, t_star, dt, St)
    t_star = t_star + dt  ! Advance in time
end do

function dvdt(v, St) result(dv)
    real(8), intent(in) :: v, St
    real(8) :: dv
    real(8) :: Re_s
    Re_s = 0.001
    ! Functions to consider
    ! dv = (1.0d0 - v) / St
    dv = (1 - v - (3.0/8.0)*Re_s*v**2) / St
end function dvdt

subroutine rk4(v, t, dt, St)
    real(8), intent(inout) :: v  ! dimensionless velocity
    real(8), intent(in) :: t, dt, St
    real(8) :: k1, k2, k3, k4

    k1 = dt * dvdt(v, St)
    k2 = dt * dvdt(v + 0.5d0 * k1, St)
    k3 = dt * dvdt(v + 0.5d0 * k2, St)
    k4 = dt * dvdt(v + k3, St)

    v = v + (k1 + 2.0d0*k2 + 2.0d0*k3 + k4) / 6.0d0
end subroutine rk4
```

An example comparison between this numerical method and the analytical result is shown below:

![Temporal development of the dimensionless velocity for St=1.005](figures/St1005.PNG)  
*Figure 1 – Temporal development of the dimensionless velocity for \( St=1.005 \).*

### How to use

Compile the program with:

```bash
ifx program_1.f90
```

and run the executable:

```bash
./a.out
```

---

## Code 2 - Bisection and False Position Methods

These numerical methods find the root (zero) of functions. Both have advantages and drawbacks. They are based on narrowing an initial interval that contains the root.

In the **bisection method**, the interval is halved each iteration. Given lower bound $(x_l\)$ and upper bound $(x_u\)$, the midpoint is:

$$
x_m = \frac{x_l + x_u}{2}
$$

If $(f(x_m) \times f(x_l) \leq 0\)$, then the root lies between $(x_l\)$ and $\(x_m)$, so update $(x_u = x_m)$; otherwise, update $(x_l = x_m\)$. Repeat until the desired tolerance is reached.

The **false position method** takes into account the function values at the interval ends to better approximate the root:

$$
x_m = x_u - \frac{f(x_u) \cdot (x_l - x_u)}{f(x_l) - f(x_u)}
$$

Then the interval is updated similarly based on the sign of $(f(x_m) \times f(x_l))$.

### How to use

Compile:

```bash
ifx program_2.f90
```

Run:

```bash
./a.out
```

---

## Code 3 - Secant and Müller Methods

This program compares the Secant and Müller methods for root finding of a given polynomial. The comparison is based on the relative error (%) at each iteration until the desired tolerance is met.

The **Secant method** approximates the root using two previous points:

$$
x_{i+1} = x_i - \frac{f(x_i) \cdot (x_{i-1} - x_i)}{f(x_{i-1}) - f(x_i)}
$$

The **Müller method** fits a parabola through three points and approximates the root by:

$$
x_{i+1} = x_i - \frac{2c}{b \pm \sqrt{b^2 - 4ac}}
$$

where:

$$
a = \frac{\delta_1 - \delta_0}{h_1 - h_0}, \quad b = a h_1 + \delta_1, \quad c = f(x_2)
$$

and

$$
h_0 = x_1 - x_0, \quad h_1 = x_2 - x_1, \quad \delta_0 = \frac{f(x_1) - f(x_0)}{x_1 - x_0}, \quad \delta_1 = \frac{f(x_2) - f(x_1)}{x_2 - x_1}
$$

Comparison plot:

![Comparison between the Secant and Müller methods](figures/secant_muller.png)  
*Figure 2 – Comparison between the Secant and Müller methods.*

### How to use

Compile:

```bash
ifx program_3.f90
```

Run:

```bash
./a.out
```

The file `grafico.plt` contains gnuplot commands to plot "Relative error (%) vs iteration". You need gnuplot installed.

Output example:

![Relative error vs iteration](figures/RE_code3.png)  
*Figure 3 – Relative error vs iteration.*

---

## Code 4 - Bairstow's Method

This program implements Bairstow’s method to find roots of any polynomial by dividing it by a quadratic polynomial:

$$
f(x) = x^2 - r x - s
$$

The method iteratively adjusts $\(r\)$ and $\(s\)$ using Newton-Raphson until the remainder is zero.

Each quadratic factor corresponds to one or two roots (including complex conjugates). The polynomial is deflated by dividing out this factor and the process repeats until all roots are found.

Since convergence depends on initial guesses for $\(r\)$ and $\(s\)$, this program can generate fractal maps visualizing convergence over ranges of $\(r\)$ and $\(s\)$.

A fractal map is a colorful image where each pixel corresponds to an initial guess $\((r,s)\)$, colored by the number of iterations to converge. Non-convergent points are black.

Examples:

<p align="center">
  <img src="figures/mapa_1.png" width="45%" />
  <img src="figures/mapa_2.png" width="45%" />
</p>

*Figure 4 – Example fractal maps for 3rd and 5th degree polynomials.*

### How to use

Compile:

```bash
ifx program_4.f90
```

Run:

```bash
./a.out
```

The program will ask for polynomial coefficients, initial guesses for $\(r\)$ and $\(s\)$, and a range for fractal generation.

The gnuplot script `mapa.gnu` plots the fractal map. After generating `mapa.dat`, run:

```bash
gnuplot -persist mapa.gnu
```

Make sure to adjust `set xrange` and `set yrange` in `mapa.gnu` if you change the fractal range.

---

## Code 5 - Linear Systems Solutions

Two programs solve two different problems.

The first code [program5.f90](./codes/program%205/program5.f90) solves a 5x5 system that arises from studying concentrations in interconnected reactors. Based on input concentrations and flow rates, it calculates the output concentrations.

This system is solved using LU decomposition. The coefficient matrix $\( \mathbf{A} \)$ is decomposed into:

$$
\mathbf{A} = \mathbf{L} \mathbf{U}
$$

where $\(\mathbf{L}\)$ is lower triangular and $\(\mathbf{U}\)$ is upper triangular. Then the following systems are solved in order:

$$
\mathbf{L} \cdot \mathbf{d} = \mathbf{b}
$$

$$
\mathbf{U} \cdot \mathbf{x} = \mathbf{d}
$$

### How to use

Compile:

```bash
ifx program_5.f90
```

Run:

```bash
./a.out
```

The second code is associated with the determination of the transient one-dimensional temperature distribution in a solid fuel associated with a nuclear reactor, subject to a heat conduction process with internal generation, with a surroundind flow of heated air.

To get the time distribution of the one-dimensional temperature it's used the finite difference method to the x direction. Considering the symmetry of the problem only half of the geometry is simulated, the other half follows the same idea.

This code plots the temperature from the center of the reactor to its edge over time, or in other words, creates a file that contains the node and the temperature over time. Plus it plots a textfile withe final temperature. There is a python file that creates the temporal profile of temperature development along the one-dimensional grid.

An example of this temporal profile plotted by the python code is shown below:

![Relative error vs iteration](figures/program5.png)  
*Figure 5 – temporal profile*
### How to use

Compile:

```bash
ifx program_5_2.f90
```

Run:

```bash
./a.out
```
and to run the graph

```bash
python3 mapa_program5.py
```
## Code 6 - multidimensional optimization

Three programs solve the same problem of multidimensional optimization here.

The idea of the programs here is to find a critical point of f(x,y) using three different methods: **Random search**, **Maximum slope method (aclive máximo)** and **Conjugate gradient method**.

$$
f(x,y) = 4x + 2y + x^2 - 2x^4 + 2xy - 3y^2
$$

the program [program6_random.f90](./codes/program%206/program6_random.f90) that solves using the random search method uses an equal and random search space for the following number of points [50, 200, 350, 500, 700, 1000], and compares the error of the value found for these two cases when compared to the previous case with fewer number of search points. It plots in the terminal window the tables of comparison.

the program [program6_aclive.f90](./codes/program%206/program6_aclive.f90) that solves using the maximum slope method starts from the starting point (0,0). It plots a table with the values needed to plot the contour lines and the search path for the optimal point in the file **trajetoria.txt**, as shown in the image below.

![Relative error vs iteration](figures/program6.png)  
*Figure 6 – contour lines and the search path*

the program [program6_grad.f90](./codes/program%206/program6_grad.f90) that solves using the conjugate gradient method also starts from the initial point (0,0). And plots the values found in each iteration in the terminal window itself.

### How to use

Compile:

```bash
ifx program_6_random.f90
ifx program_6_aclive.f90
ifx program_6_grad.f90
```

Run to each one:

```bash
./a.out
```

and to run the contour lines and the search path

```bash
python3 mapa_program6.py
```

## Code 7 - Numerical Integration

Here include two programs, one to run manually and one to create a graph T_med x n.

The idea here is to implement numerical integration techniques to approximate the definite integral of a given function f(x,y) within a specific interval $0/le x/le 8$ and $0/le y/le 6$. The methods available include the **Trapezoidal Rule**, **Simpson’s 1/3 Rule** and **Simpson’s 3/8 Rule**.

There is a two-dimensional domain of dimensions 8 x 6, and the temperature at each point in this domain is given by the expression:

$$
f(x,y) = 2xy + 2x - x^2 - 2y^2 + 72
$$

The idea is to discretize this space into n x n nodes and associate the temperature with each node using the expression above, and then determine the average temperature given by the following expression:

$$
\overline{T} = \frac{1}{8 \cdot 6} \int_0^8 \int_0^6 T(x, y) dy dx
$$

The program [program7.f90](./codes/program%207/program7.f90) works manually, where you enter the grid size nxn, and it will compute the average temperature for all three cases.

Program [program7_varre.f90](./codes/program%207/program7_varre.f90) solves for different grid sizes and generates a comparison table between the three methods, which can be viewed by running the Python program in the same directory.

### How to use

Compile:

```bash
ifx program7.f90
ifx program7_varre.f90
```

Run to each one:

```bash
./a.out
```

and to run the contour lines and the search path

```bash
python3 graph7.py
```
