program num_root_methods

! Variaveis
implicit none

real :: x_l,x_u,x_m,x_antigo
integer :: metodo,funcao,iter
real :: fxl, fxu,fxm,erro
character (len=100) :: f1

! Titulo
write(*,*) "***************************************************************************"
write(*,*) "*                                                                         *"
write(*,*) "*                                                                         *"
write(*,*) "*            Program_2 - método da bissecção e da posição falsa           *"
write(*,*) "*                                                                         *"
write(*,*) "*                   André Brandao - Métodos numericos                     *"
write(*,*) "*                                                                         *"
write(*,*) "*                                                                         *"
write(*,*) "***************************************************************************"

!Escolher funcao
write(*,*)"Escolha a funcao:"
write(*,*)"1 - f(x) = -0.5x**2 + 2.5x + 4.5"
write(*,*)"2 - f(x) = 5x**3 - 5x**2 + 6x - 2"
write(*,*)"3 - f(x) = -25 + 82x - 90x**2 + 44x**3 - 8x**4 + 0.7x**5"
write(*,*)"4 - sen(x) = x**3"
write(*,*)"5 - ln(x**4) = 0.7"
read(*,*) funcao

if (funcao.eq.1) then
f1 = "1 - f(x) = -0.5x**2 + 2.5x + 4.5"
end if
if (funcao.eq.2) then
f1 = "2 - f(x) = 5x**3 - 5x**2 + 6x - 2"
end if
if (funcao.eq.3) then
f1 = "3 - f(x) = -25 + 82x - 90x**2 + 44x**3 - 8x**4 + 0.7x**5"
end if
if (funcao.eq.4) then
f1 = "4 - sen(x) = x**3"
end if
if (funcao.eq.5) then
f1 = "5 - ln(x**4) = 0.7"
end if

!Escolha do metodo
write(*,*) "Qual método deseja utilizar?"
write(*,*) "1 - Método da bisseccao"
write(*,*) "2 - Método da falsa posicao"
read(*,*) metodo

!Seleciona os limites
write(*,*) "Determine os limites inferiores (x_l) e superiores (x_u)."
write(*,*)"x_l: "
read(*,*) x_l
write(*,*)"x_u: "
read(*,*) x_u
erro=100
x_m=1

!Metodo da bissecao
if (metodo.EQ.1) then
 iter=1
 open(unit=10, file='saida_bisec.txt', status='unknown')
 write(*,*) "iteracao    x_m    erro"
 write(10,*) "funcao escolhida:",f1
 write(10,*) "iteracao    x_m    erro"
 !Comeca o loop
 do while (erro.GE.10E-06)
  !Usa a funcao escolhida
  x_antigo = x_m
  if(funcao.EQ.1) then
   call funcao1(x_l, fxl)
   call funcao1(x_u, fxu)
  end if
  if(funcao.EQ.2) then
    call funcao2(x_l, fxl)
    call funcao2(x_u, fxu)
  end if
  if(funcao.EQ.3) then
    call funcao3(x_l, fxl)
    call funcao3(x_u, fxu)
  end if
  if(funcao.EQ.4) then
    call funcao4(x_l, fxl)
    call funcao4(x_u, fxu)
  end if
  if(funcao.EQ.5) then
    call funcao5(x_l, fxl)
    call funcao5(x_u, fxu)
  end if
  !checa se tem raiz
  if (fxl*fxu.LE.0) then
   x_m= (x_l+x_u)/2
   !Atualiza o x para x_m para a funcao sendo utilizada
   if(funcao.eq.1) then
    call funcao1(x_m, fxm)
   end if
   if(funcao.eq.2) then
    call funcao2(x_m, fxm)
   end if
   if(funcao.eq.3) then
    call funcao3(x_m, fxm)
   end if 
   if(funcao.eq.4) then
    call funcao4(x_m, fxm)
   end if
   if(funcao.eq.5) then
    call funcao5(x_m, fxm)
   end if  
   !Decide aqui qual limite precisa ser atualizado
   if (fxl*fxm.LE.0) then
    x_u=x_m
    iter=iter+1
   else
    x_l=x_m
    iter=iter+1
   end if
  !escreve na tela o numero da iteracao, o xm e o erro
  erro=abs((x_m-x_antigo)/x_m)
  write(*,'(I3,2X,F10.6,F10.6)') iter, x_m, erro
  write(10,'(I3,2X,F10.6,F10.6)') iter, x_m, erro
  else
   write(*,*) "intervalo nao contem raiz"
   stop
  end if
end do
else
 open(unit=11, file='saida_fals.txt', status='unknown')
 iter=1

 write(*,*) "iteracao    x_m    erro"
 write(11,*) "funcao escolhida:",f1
 write(11,*) "iteracao    x_m    erro"
 !Comeca o loop
 do while (erro.GE.10E-06)
  !Usa a funcao escolhida
  x_antigo = x_m
  if(funcao.EQ.1) then
    call funcao1(x_l, fxl)
    call funcao1(x_u, fxu)
   end if
   if(funcao.EQ.2) then
     call funcao2(x_l, fxl)
     call funcao2(x_u, fxu)
   end if
   if(funcao.EQ.3) then
     call funcao3(x_l, fxl)
     call funcao3(x_u, fxu)
   end if
   if(funcao.EQ.4) then
     call funcao4(x_l, fxl)
     call funcao4(x_u, fxu)
   end if
   if(funcao.EQ.5) then
     call funcao5(x_l, fxl)
     call funcao5(x_u, fxu)
   end if
  !Verifica se tem raiz no intervalo
   
  if (fxl*fxu.LE.0) then
   x_m = x_u - fxu*(x_l-x_u)/(fxl-fxu)
   !atualiza o xm para a funcao escolhida
   if(funcao.eq.1) then
    call funcao1(x_m, fxm)
   end if
   if(funcao.eq.2) then
    call funcao2(x_m, fxm)
   end if
   if(funcao.eq.3) then
    call funcao3(x_m, fxm)
   end if 
   if(funcao.eq.4) then
    call funcao4(x_m, fxm)
   end if
   if(funcao.eq.5) then
    call funcao5(x_m, fxm)
   end if 
   !Aqui decide que extremo vai trocar
   if (fxm*fxl.le.0) then
    x_u = x_m
    iter = iter+1
   else
    x_l = x_m
    iter=iter+1
   end if 
   
  !escreve na tela o numero da iteracao, o xm e o erro
  erro=abs((x_m-x_antigo)/x_m)
  write(*,'(I3,2X,F10.6,F10.6)') iter, x_m, erro
  write(11,'(I3,2X,F10.6,F10.6)') iter, x_m, erro
  else
   write(*,*) "intervalo nao contem raiz"
   stop
  end if   
  
  end do
end if

end program num_root_methods

subroutine funcao1(x,fx)
  real, intent(in) :: x
  real, intent(out) :: fx
  fx = -0.5*(x**2) + (2.5*x) + 4.5
end subroutine funcao1

subroutine funcao2(x,fx)
  implicit none
  real, intent(in) :: x
  real, intent(out) :: fx
  fx = 5*(x**3) - 5*(x**2) + (6*x) - 2
end subroutine funcao2

subroutine funcao3(x,fx)
  implicit none
  real, intent(in) :: x
  real, intent(out) :: fx
  fx = -25 + 82*x - 90*(x**2) + 44*(x**3) - 8*(x**4) + 0.7*(x**5)
end subroutine funcao3

subroutine funcao4(x,fx)
  implicit none
  real, intent(in) :: x
  real, intent(out) :: fx
  fx = sin(x)-x**3
end subroutine funcao4

subroutine funcao5(x,fx)
  implicit none
  real, intent(in) :: x
  real, intent(out) :: fx
  fx = 4*log(x)-0.7
end subroutine funcao5
