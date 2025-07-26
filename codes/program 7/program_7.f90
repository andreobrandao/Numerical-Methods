program integ
  implicit none
  integer :: n
  real, allocatable :: x(:), y(:), T(:)
  real :: dx, dy
  integer :: i, j, unit
  real :: Ix1, Ixn, Ixint, Ixsoma
  real :: Iy, T_med
  real :: Ix_simpson, Iy_simpson, T_med_simpson
  real :: Ix_simpson38, Iy_simpson38, T_med_simpson38
  logical :: usar_simpson, usar_simpson38
  character(len=30) :: arq_saida

  write(*,*) "Enter with the amount of nodes in each direction of your grid:"
  read(*,*) n

  if (mod(n,2) == 0) then
     usar_simpson = .false.
     print *, "Nota: n é par, então a Regra de Simpson 1/3 não será usada."
  else
     usar_simpson = .true.
     print *, "Nota: n é ímpar, a Regra de Simpson 1/3 será usada também."
  end if

  if (mod(n-1,3) == 0) then
     usar_simpson38 = .true.
     print *, "Nota: n-1 é múltiplo de 3, a Regra de Simpson 3/8 será usada também."
  else
     usar_simpson38 = .false.
     print *, "Nota: n-1 não é múltiplo de 3, a Regra de Simpson 3/8 será ignorada."
  end if

  allocate(x(n), y(n), T(n*n))
  dx = 8.0 / (n - 1)
  dy = 6.0 / (n - 1)

  ! Geração das coordenadas
  do i = 1, n
    x(i) = (i - 1) * dx
    y(i) = (i - 1) * dy
  end do

  ! Geração do campo de temperatura
  do i = 1, n
    do j = 1, n
       T((i-1)*n + j) = 2.0 * x(i) * y(j) + 2.0 * x(i) - x(i)**2 - 2.0 * y(j)**2 + 72.0
    end do
  end do

  ! Salvar o campo de temperatura
  arq_saida = "temp_field_program7.txt"
  open(newunit=unit, file=arq_saida, status="replace", action="write")
  do i = 1, n
    do j = 1, n
      write(unit, '(F10.5,1X,F10.5,1X,F10.5)') x(i), y(j), T((i-1)*n + j)
    end do
  end do
  close(unit)

  ! ---------------------------
  ! Regra do trapézio
  ! ---------------------------
  Ix1 = (dx / 2.0) * (T(1) + T(n) + 2.0 * sum(T(2:n-1)))   ! y = 0
  Ixn = (dx / 2.0) * (T((n-1)*n + 1) + T((n-1)*n + n) + 2.0 * sum(T((n-1)*n + 2:(n-1)*n + n - 1))) ! y = 6

  Ixsoma = 0.0
  do i = 2, n-1
    Ixint = (dx / 2.0) * (T((i-1)*n + 1) + T((i-1)*n + n) + 2.0 * sum(T((i-1)*n + 2:(i-1)*n + n - 1)))
    Ixsoma = Ixsoma + Ixint
  end do

  Iy = (dy / 2.0) * (Ix1 + Ixn + 2.0 * Ixsoma)
  T_med = Iy / (8.0 * 6.0)

  write(*,'(A,F10.5)') "Average temperature (Trapézio): ", T_med

  ! ---------------------------
  ! Regra de Simpson 1/3 (se possível)
  ! ---------------------------
  if (usar_simpson) then
     Iy_simpson = 0.0
     do j = 1, n
        Ix_simpson = T(1 + (j-1)*n) + T(n + (j-1)*n)
        do i = 2, n - 1
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

     write(*,'(A,F10.5)') "Average temperature (Simpson 1/3): ", T_med_simpson
  end if

  ! ---------------------------
  ! Regra de Simpson 3/8 (se possível)
  ! ---------------------------
  if (usar_simpson38) then
     Iy_simpson38 = 0.0
     do j = 1, n
        Ix_simpson38 = T(1 + (j-1)*n) + T(n + (j-1)*n)
        do i = 2, n - 1
           if (mod(i-1,3) == 0) then
              Ix_simpson38 = Ix_simpson38 + 2.0 * T(i + (j-1)*n)
           else
              Ix_simpson38 = Ix_simpson38 + 3.0 * T(i + (j-1)*n)
           end if
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

     write(*,'(A,F10.5)') "Average temperature (Simpson 3/8): ", T_med_simpson38
  end if

end program integ
