program comparacao_distribuicao_pontos
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307) !definindo alta precisão 
  integer :: i, j, nx, ny, n_pontos, total_pontos
  real(dp) :: x, y, fx, fmin_reg, fmin_rand, f_comp1, erro_rand,erro_reg, f_comp2
  real(dp) :: dx, dy
  real(dp) :: xmin_reg, ymin_reg, xmin_rand, ymin_rand
  integer, dimension(6) :: pontos = [50, 200, 350, 500, 700, 1000]
  integer :: k
  real(dp), parameter :: xinf = -3.0_dp, xsup = 3.0_dp
  real(dp), parameter :: yinf = -3.0_dp, ysup = 3.0_dp
  call random_seed()

  print *, "---------------------------------------------------------------------------------------"
  print *, " n_pontos | Mín (Regular)|    erro(%)        | Mín (Randômica)    |  erro(%)          "
  print *, "--------------------------------------------------------------------------------------"

  do k = 1, size(pontos) !varredura para cada quantidade pontos
     n_pontos = pontos(k)

     ! --------- Distribuição REGULAR ----------
     ! Aproxima um arranjo nx x ny ≈ n_pontos
     nx = int(sqrt(real(n_pontos))) !em x
     ny = n_pontos / nx !em y
     total_pontos = nx * ny

     dx = (xsup - xinf) / real(nx - 1)
     dy = (ysup - yinf) / real(ny - 1)

     fmin_reg = -1.0e30_dp !pra iniciar com um valorpequeno

     do i = 0, nx-1
        do j = 0, ny-1
           x = xinf + i*dx
           y = yinf + j*dy
           fx = 4*x + 2*y + x**2 - 2*x**4 + 2*x*y - 3*y**2

           if (fx > fmin_reg) then 
            ! armazenando o mínimo encontrado para cada laço
              fmin_reg = fx
              xmin_reg = x
              ymin_reg = y
           end if
        end do
     end do
     if (k==1) then
         f_comp2 = fmin_reg !primeiro valor de referência
         erro_reg = 0.0_dp
      else
         erro_reg = abs(fmin_reg - f_comp2) / abs(f_comp2)
      end if
      f_comp2 = fmin_reg
     ! --------- Distribuição RANDÔMICA ----------
     fmin_rand = -1.0e30_dp
    
     do i = 1, n_pontos
        call random_number(x)
        call random_number(y)
        x = xinf + (xsup - xinf) * x
        y = yinf + (ysup - yinf) * y
        fx = 4*x + 2*y + x**2 - 2*x**4 + 2*x*y - 3*y**2

        if (fx > fmin_rand) then
           fmin_rand = fx
           xmin_rand = x
           ymin_rand = y
        end if
     end do
       ! --------- Comparação dos resultados ----------
     
     if (k==1) then
         erro_rand = 0.0_dp
      else
         erro_rand = abs(fmin_rand - f_comp1) / abs(f_comp1)
      end if
      f_comp1=fmin_rand
     ! --------- Tabela de resultados ----------
     write(*,'(I5, 2X, F12.5, " @ (",F5.2,",",F5.2,")",2X,F6.3, 2X, F12.5, " @ (",F5.2,",",F5.2,")",3X,F6.3)') &
           n_pontos, fmin_reg, xmin_reg, ymin_reg, erro_reg*100.0_dp, fmin_rand, xmin_rand, ymin_rand, erro_rand*100.0_dp

  end do

  print *, "------------------------------------------------------"

end program comparacao_distribuicao_pontos
