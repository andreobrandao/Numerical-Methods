program code3
implicit none

! variaveis para o metodo da secante
real :: fx1,fx2
real :: x_1,x_2 
real :: xm 
integer, parameter :: m = 30
integer :: iter
real :: tol, erro
real, allocatable :: errorelat1(:)

! variaveis para o metodo de muller
complex :: x0,x1,x2,xn
complex :: fx0m,fx1m,fx2m
complex :: h0,h1,d0,d1,a,b,c
complex :: check1,check2
integer, parameter :: n = 30
integer :: h,iter2
real :: erro2
real, allocatable :: errorelat2(:)

! Alocando espaco para os vetores de erro relativo
allocate(errorelat1(m))
allocate(errorelat2(n))

!Abre o arquivo de saida
open (unit=10, file="saida.dat", status="unknown")

! Titulo
write(*,*) "***************************************************************************"
write(*,*) "*                                                                         *"
write(*,*) "*                                                                         *"
write(*,*) "*                Program_3 - metodo da secante e de Muller                *"
write(*,*) "*                                                                         *"
write(*,*) "*                   Andre Brandao - Metodos numericos                     *"
write(*,*) "*                                                                         *"
write(*,*) "*                                                                         *"
write(*,*) "***************************************************************************"
write (*,*) "            "
write (*,*) "METODO DA SECANTE"
write (*,*) "Entre com o chute inicial x1: "
read (*,*) x_1
write (*,*) "Entre com o chute inicial x2: "
read (*,*) x_2
write (*,*) "            "

write (*,*) "METODO DE MULLER"
write (*,*) "Entre com o chute inicial x0: "
read (*,*) x0
write (*,*) "Entre com o chute inicial x1: "
read (*,*) x1
write (*,*) "Entre com o chute inicial x2: "
read (*,*) x2

write(10,*) "iter       errosecante       erromuller"

tol= 1.0e-6 !tolerancia
iter=1
erro = 1.0
! fx = (x-1)*(x-3)*(x-5)*(x-7)*(x-9)

! Metodo da secante
do while (erro.GE.tol)
    call funcao (x_1,fx1) ! aplica a funcao para o primeiro chute
    call funcao (x_2,fx2) ! aplica a funcao para o segundo chute

    ! Verifica se a divisao fx2-fx1 é proximo de zero
    if (abs(fx2-fx1).LE.1.0e-44) then
        exit
    end if
    ! Atualiza o valor de xm
    xm = x_1 - ((fx1*(x_2 - x_1))/(fx2 - fx1))
    if (abs(xm).GE.1.0e-12) then
        erro = abs((xm - x_1)/xm)
        errorelat1(iter)= abs((xm - x_2)/xm)*100
    else
        erro = abs(xm - x_1)
    end if
    ! Atualiza os dados para a proxima iteracao
    x_1 = x_2
    x_2 = xm
    iter=iter + 1
end do

! Escreve na tela os resultados
write (*,*) "Metodo da secante"
write (*,*) "Raiz encontrada: ", xm
write (*,*) "Numero de iteracoes: ", iter
write (*,*) "Erro: ", erro
write (*,*) "            "


! Metodo de muller
erro2 = 1.0
iter2 = 1

do while (erro2.GE.tol)
    ! Aplica a funcao para os chutes iniciais
    call funcaoc(x0,fx0m)
    call funcaoc(x1,fx1m)
    call funcaoc(x2,fx2m)

    ! Variaveis utilizadas no metodo
    h0 = x1-x0
    h1 = x2-x1
    d0 = (fx1m - fx0m)/(x1-x0)
    d1 = (fx2m - fx1m)/(x2-x1)

    a = (d1 - d0)/(h1 + h0)
    b = a*h1 + d1
    c = fx2m

    ! Verifica se a divisao (x1-x0) ou (x2-x1) é proximo de zero
    if (abs(x1-x0).LE.1.0e-44.OR.abs(x2-x1).LE.1.0e-44) then
        exit
    end if
    ! Verifica qual maximiza o denominador
    check1 = b+sqrt(b**2 - 4*a*c)
    check2 = b-sqrt(b**2 - 4*a*c)
    ! atualiza o valor de xn
    if (abs(check1).GE.abs(check2)) then
        xn = x2 - (2*c / (b + csqrt(b**2 - 4*a*c)))
    else
        xn = x2 - (2*c / (b - csqrt(b**2 - 4*a*c)))
    end if

    ! Determina o erro
    if (abs(xn).GE.1.0e-12) then
        erro2 = abs((xn - x1)/xn)
        errorelat2(iter2) = abs((real(xn) - real(x2))/real(xn))*100
    else
        erro2 = abs(xn - x1)
    end if

    ! Atualiza os dados para a proxima iteracao
    x0=x1
    x1=x2
    x2=xn
    iter2=iter2 + 1
end do

! Escreve os resultados no arquivo de saida
if (iter.GE.iter2) then
    do h = 1, iter
        write(10,'(I2,2ES20.5)') h, errorelat1(h), errorelat2(h)
    end do
else
    do h = 1, iter2
        write(10,'(I2,2ES20.5)') h, errorelat1(h), errorelat2(h)
    end do
end if



! Escreve na tela os resultados
write (*,*) "Metodo de muller"
write (*,*) "Raiz encontrada: ", real(xn)
write (*,*) "Numero de iteracoes: ", iter2
write (*,*) "Erro: ", erro2

!Plota o grafico
call system("gnuplot grafico.plt")
end program code3

! Subrotina que define a funcao
subroutine funcaoc(x,fx)
    implicit none
    complex, intent(in) :: x
    complex, intent(out) :: fx
    fx = x**5 - 25*x**4 + 230*x**3 - 950*x**2 + 1689*x - 945
  end subroutine funcaoc

  subroutine funcao(x,fx)
    implicit none
    real, intent(in) :: x
    real, intent(out) :: fx
    fx = x**5 - 25*x**4 + 230*x**3 - 950*x**2 + 1689*x - 945
  end subroutine funcao