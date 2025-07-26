program lu_solver
  implicit none
  integer, parameter :: n = 5
  integer :: i, j, k
  real(8) :: A(n,n), L(n,n), U(n,n), b(n), d(n), c(n), factor

  ! Matriz A
A = reshape([6.0d0,  0.0d0, -1.0d0,  0.0d0,  0.0d0, &
             -3.0d0,  3.0d0,  0.0d0,  0.0d0,  0.0d0, &
              0.0d0, -1.0d0,  9.0d0,  0.0d0,  0.0d0, &
              0.0d0, -1.0d0, -8.0d0, 11.0d0, -2.0d0, &
             -3.0d0, -1.0d0,  0.0d0,  0.0d0,  4.0d0], &
             shape=[n,n], order=[2,1])

  ! Vetor b
  b = [50.0d0, 0.0d0, 160.0d0, 0.0d0, 0.0d0]

  ! Inicializa L como identidade e U como A
  L = 0.0d0
  U = A
  do i = 1, n
     L(i,i) = 1.0d0
  end do

  ! Decomposição LU sem pivoteamento
  do k = 1, n-1
     do i = k+1, n ! varre as linhas abaixo da diagonal
        factor = U(i,k) / U(k,k)
        !write(*,*) 'Fator de escala para linha ', i, ': ', factor
        L(i,k) = factor
        do j = k, n ! varre a coluna aqui
           U(i,j) = U(i,j) - factor * U(k,j)
        end do
     end do
  end do
!do i=1,n
!    do j=1,n
!     print *, 'U(', i, ',', j, ') = ', U(i,j)
!    end do
!end do
  ! Substituição progressiva: L * d = b
  d(1) = b(1)
  do i = 2, n
     d(i) = b(i)
     do j = 1, i-1
        d(i) = d(i) - L(i,j) * d(j)
     end do
  end do
 

  ! Substituição regressiva: U * c = d
  c(n) = d(n) / U(n,n)
  do i = n-1, 1, -1
     c(i) = d(i)
     do j = i+1, n
        c(i) = c(i) - U(i,j) * c(j)
     end do
     c(i) = c(i) / U(i,i)
  end do

  ! Imprime o vetor solução c
  print *, 'Solução do sistema:'
  do i = 1, n
     print *, 'c(',i, ') = ', c(i)
  end do

end program lu_solver
