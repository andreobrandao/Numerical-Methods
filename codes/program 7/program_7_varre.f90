program integ_varredura
  implicit none
  integer :: n, i, j, unit, nf, idx
  integer, parameter :: n_vals(6) = [7, 13, 19, 25, 31, 37]
  real, allocatable :: x(:), y(:), T(:)
  real :: dx, dy
  real :: Ix1, Ixn, Ixint, Ixsoma
  real :: Iy, T_med
  real :: Ix_simpson, Iy_simpson, T_med_simpson
  real :: Ix_simpson38, Iy_simpson38, T_med_simpson38
  logical :: usar_simpson
  character(len=30) :: arq_saida

  ! Arquivo para salvar os resultados
  open(newunit=nf, file="resultados_temp.txt", status="replace", action="write")
  write(nf,'(A)') "n  Trapézio   Simpson13   Simpson38"

  do idx = 1, size(n_vals)
     n = n_vals(idx)

     if (mod(n,2) == 0) then
        usar_simpson = .false.
     else
        usar_simpson = .true.
     end if

     allocate(x(n), y(n), T(n*n))
     dx = 8.0 / (n - 1)
     dy = 6.0 / (n - 1)

     do i = 1, n
        x(i) = (i - 1) * dx
        y(i) = (i - 1) * dy
     end do

     do i = 1, n
        do j = 1, n
           T((i-1)*n + j) = 2.0 * x(i) * y(j) + 2.0 * x(i) - x(i)**2 - 2.0 * y(j)**2 + 72.0
        end do
     end do

     ! ------------------ Trapézio ------------------
     Ix1 = (dx / 2.0) * (T(1) + T(n) + 2.0 * sum(T(2:n-1)))
     Ixn = (dx / 2.0) * (T((n-1)*n + 1) + T((n-1)*n + n) + 2.0 * sum(T((n-1)*n + 2:(n-1)*n + n - 1)))

     Ixsoma = 0.0
     do i = 2, n-1
        Ixint = (dx / 2.0) * (T((i-1)*n + 1) + T((i-1)*n + n) + 2.0 * sum(T((i-1)*n + 2:(i-1)*n + n - 1)))
        Ixsoma = Ixsoma + Ixint
     end do

     Iy = (dy / 2.0) * (Ix1 + Ixn + 2.0 * Ixsoma)
     T_med = Iy / (8.0 * 6.0)

     ! -------------- Simpson 1/3 --------------------
     T_med_simpson = 0.0
     if (usar_simpson) then
        Iy_simpson = 0.0
        do j = 1, n
           Ix_simpson = T(1 + (j-1)*n) + T(n + (j-1)*n)
           do i = 2, n-1
              if (mod(i,2) == 0) then
                 Ix_simpson = Ix_simpson + 4.0 * T(i + (j-1)*n)
              else
                 Ix_simpson = Ix_simpson + 2.0 * T(i + (j-1)*n)
              end if
           end do
           Ix_simpson = dx / 3.0 * Ix_simpson

           if (j == 1 .or. j == n) then
              Iy_simpson = Iy_simpson + Ix_simpson
           else if (mod(j,2) == 0) then
              Iy_simpson = Iy_simpson + 4.0 * Ix_simpson
           else
              Iy_simpson = Iy_simpson + 2.0 * Ix_simpson
           end if
        end do
        Iy_simpson = dy / 3.0 * Iy_simpson
        T_med_simpson = Iy_simpson / (8.0 * 6.0)
     end if

     ! -------------- Simpson 3/8 --------------------
     T_med_simpson38 = 0.0
     Iy_simpson38 = 0.0
     do j = 1, n
        Ix_simpson38 = T(1 + (j-1)*n) + T(n + (j-1)*n)
        do i = 2, n-1
           select case (mod(i-1,3))
           case (0)
              Ix_simpson38 = Ix_simpson38 + 2.0 * T(i + (j-1)*n)
           case default
              Ix_simpson38 = Ix_simpson38 + 3.0 * T(i + (j-1)*n)
           end select
        end do
        Ix_simpson38 = 3.0 * dx / 8.0 * Ix_simpson38

        if (j == 1 .or. j == n) then
           Iy_simpson38 = Iy_simpson38 + Ix_simpson38
        else if (mod(j-1,3) == 0) then
           Iy_simpson38 = Iy_simpson38 + 2.0 * Ix_simpson38
        else
           Iy_simpson38 = Iy_simpson38 + 3.0 * Ix_simpson38
        end if
     end do
     Iy_simpson38 = 3.0 * dy / 8.0 * Iy_simpson38
     T_med_simpson38 = Iy_simpson38 / (8.0 * 6.0)

     ! Salva resultado no arquivo
     write(nf,'(I3,3F12.6)') n, T_med, T_med_simpson, T_med_simpson38

     deallocate(x, y, T)
  end do

  close(nf)

  print *, "Arquivo 'resultados_temp.txt' gerado com sucesso."
end program integ_varredura
