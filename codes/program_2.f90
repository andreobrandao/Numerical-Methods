program num_root_methods

!Variaveis
implicit none

real :: x_l,x_u,x_m,x_antigo
integer :: metodo,funcao,iter
real :: fxl, fxu,fxm,erro
character (len=100) :: funcao1

!Titulo
write(*,*) "***************************************************************************"
write(*,*) "*                                                                         *"
write(*,*) "*                                                                         *"
write(*,*) "*            Program_2 - metodo da bisseccao e da posicao falsa           *"
write(*,*) "*                                                                         *"
write(*,*) "*                   Andre Brandao - Metodos numericos                     *"
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
funcao1 = "1 - f(x) = -0.5x**2 + 2.5x + 4.5"
end if
if (funcao.eq.2) then
funcao1 = "2 - f(x) = 5x**3 - 5x**2 + 6x - 2"
end if
if (funcao.eq.3) then
funcao1 = "3 - f(x) = -25 + 82x - 90x**2 + 44x**3 - 8x**4 + 0.7x**5"
end if
if (funcao.eq.4) then
funcao1 = "4 - sen(x) = x**3"
end if
if (funcao.eq.5) then
funcao1 = "5 - ln(x**4) = 0.7"
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
 write(10,*) "funcao escolhida:",funcao1
 write(10,*) "iteracao    x_m    erro"
 !Comeca o loop
 do while (erro.GE.10E-06)
  !Usa a funcao escolhida
  x_antigo = x_m
  if(funcao.EQ.1) then
   fxl = -0.5*(x_l**2) + (2.5*x_l) + 4.5
   fxu = -0.5*(x_u**2) + (2.5*x_u) + 4.5
  end if
  if(funcao.EQ.2) then
   fxl = 5*(x_l**3) - 5*(x_l**2) + (6*x_l) - 2
   fxu = 5*(x_u**3) - 5*(x_u**2) + (6*x_u) - 2
  end if
  if(funcao.EQ.3) then
   fxl = -25 + 82*x_l - 90*(x_l**2) + 44*(x_l**3) - 8*(x_l**4) + 0.7*(x_l**5)
   fxu = -25 + 82*x_u - 90*(x_u**2) + 44*(x_u**3) - 8*(x_u**4) + 0.7*(x_u**5)
  end if
  if(funcao.EQ.4) then
   fxl = sin(x_l)-x_l**3
   fxu = sin(x_u)-x_u**3
  end if
  if(funcao.EQ.5) then
   fxl = 4*log(x_l)-0.7
   fxu = 4*log(x_u)-0.7
  end if	
  !checa se tem raiz
  if (fxl*fxu.LE.0) then
   x_m= (x_l+x_u)/2
   !Atualiza o x para x_m para a funcao sendo utilizada
   if(funcao.eq.1) then
    fxm = -0.5*(x_m**2) + (2.5*x_m) + 4.5
   end if
   if(funcao.eq.2) then
    fxm = 5*(x_m**3) - 5*(x_m**2) + (6*x_m) - 2
   end if
   if(funcao.eq.3) then
    fxm = -25 + 82*x_m - 90*(x_m**2) + 44*(x_m**3) - 8*(x_m**4) + 0.7*(x_m**5)
   end if 
   if(funcao.eq.4) then
    fxm = sin(x_m)-x_m**3
   end if
   if(funcao.eq.5) then
    fxm = 4*log(x_m)-0.7
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
  end if
  end do		
else
 open(unit=11, file='saida_fals.txt', status='unknown')
 iter=1

 write(*,*) "iteracao    x_m    erro"
 write(11,*) "funcao escolhida:",funcao1
 write(11,*) "iteracao    x_m    erro"
 !Comeca o loop
 do while (erro.GE.10E-06)
  !Usa a funcao escolhida
  x_antigo = x_m
  if(funcao.EQ.1) then
   fxl = -0.5*(x_l**2) + (2.5*x_l) + 4.5
   fxu = -0.5*(x_u**2) + (2.5*x_u) + 4.5
  end if
  if(funcao.EQ.2) then
   fxl = 5*(x_l**3) - 5*(x_l**2) + (6*x_l) - 2
   fxu = 5*(x_u**3) - 5*(x_u**2) + (6*x_u) - 2
  end if
  if(funcao.EQ.3) then
   fxl = -25 + 82*x_l - 90*(x_l**2) + 44*(x_l**3) - 8*(x_l**4) + 0.7*(x_l**5)
   fxu = -25 + 82*x_u - 90*(x_u**2) + 44*(x_u**3) - 8*(x_u**4) + 0.7*(x_u**5)
  end if
  if(funcao.EQ.4) then
   fxl = sin(x_l)-x_l**3
   fxu = sin(x_u)-x_u**3
  end if
  if(funcao.EQ.5) then
   fxl = 4*log(x_l)-0.7
   fxu = 4*log(x_u)-0.7
  end if
  !Verifica se tem raiz no intervalo
  if (fxl*fxu.LE.0) then
   x_m = x_u - fxu*(x_l-x_u)/(fxl-fxu)
   !atualiza o xm para a funcao escolhida
   if(funcao.eq.1) then
    fxm = -0.5*(x_m**2) + (2.5*x_m) + 4.5
   end if
   if(funcao.eq.2) then
    fxm = 5*(x_m**3) - 5*(x_m**2) + (6*x_m) - 2
   end if
   if(funcao.eq.3) then
    fxm = -25 + 82*x_m - 90*(x_m**2) + 44*(x_m**3) - 8*(x_m**4) + 0.7*(x_m**5)
   end if 
   if(funcao.eq.4) then
    fxm = sin(x_m)-x_m**3
   end if
   if(funcao.eq.5) then
    fxm = 4*log(x_m)-0.7
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
  end if   
  
  end do
end if

end program num_root_methods
