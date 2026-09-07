!
!   cbp_ej9_qinf.f90
!------------------------------------------------------------------------------------------------
!   3 - Examen QINF, Ejercicio 9. 
!   Introducción a la Información y Computación Cuánticas. Máster en Física Avanzada UNED.
!
!   Autor: Cayetano Bayona Pacheco, Junio 2026.
!------------------------------------------------------------------------------------------------
!   
!   Programa para resolver la ecuación de Lindblad junto con su aproximación de Monte Carlo
!   para un sistema cuántico de 1 qubit sometido a un Hamiltoniano H = Gamma * sigma_z, 
!   donde sigma_z es la matriz de Pauli Z.
!   Este sistema está sometido a un decaimineto |+> -> |-> con una constante de decaimiento gamma.
!
!   ---> ES MUY IMPORTANTE CREAR LAS CARPETAS /out_lindblad Y /out_monte_carlo ANTES DE EJECUTAR
!   ---> EL PROGRAMA PARA QUE FORTRAN EXPORTE LOS DATOS SIN DAR ERRORES.
!
!------------------------------------------------------------------------------------------------
!   Compilación y ejecución (opciones):
!       gfortran cbp_ej9_qinf.f90 -llapack -lblas -fopenmp -o cbp_ej9_qinf.exe
!       gfortran -O3 -march=native -fopenmp -fno-stack-arrays cbp_ej9_qinf.f90 -llapack -lblas -o cbp_ej9_qinf.exe
!       gfortran -O3 -march=native cbp_ej9_qinf.f90 -o cbp_ej9_qinf.exe
!       gfortran -g -O0 -fcheck=all -fbacktrace -Wall -Wextra cbp_ej9_qinf.f90 -o cbp_ej9_qinf.exe
!       ./cbp_ej9_qinf.exe
!------------------------------------------------------------------------------------------------

program cbp_ej9_qinf
    implicit none

    real(8) :: init_vector(3,3)
    real(8) :: gamma_vector(3)
    real(8) :: T_final
    real(8) :: dt = 0.01d0

    integer :: i,j
    character(len=2) :: identifier

    ! Variables para medir el tiempo
    integer(8) :: start_count, count_rate

    ! Hora de inicio de la ejecución
    call print_actual_time(0)

    ! Obtenemos la frecuencia del reloj (ticks por segundo) y el conteo inicial
    call system_clock(start_count, count_rate)


    init_vector(1,:) = [0.0d0, 0.0d0, 1.0d0]  ! Polo Norte |0>
    init_vector(2,:) = [-1.0d0, 0.0d0, 0.0d0] ! Estado |-> (Destino del decaimiento)
    init_vector(3,:) = [1.0d0, 0.0d0, 0.0d0]  ! Estado |+> (Origen del decaimiento)

    gamma_vector = [0.1d0, 1.0d0, 10.0d0]

    do i = 1, 3
        do j = 1, 3
            write(identifier, '(I2.2)') i*10 + j
            T_final = 10.0d0 / gamma_vector(i)

            ! Primera parte del ejercicio: resolver la ecuación de Lindblad
            call main_lindblad(init_vector(j,:), gamma_vector(i), T_final, dt, identifier)

            ! Segunda parte del ejercicio: resolver la ecuación en Monte Carlo
            call main_monte_carlo(init_vector(j,:), gamma_vector(i), T_final, dt, 500, identifier)
        end do
    end do


    call print_actual_time(1)
    call elapsed_time(start_count, count_rate)


    write(*,*) 
    write(*,*)'########################################################'
    write(*,*)'            -> Programa terminado <-'
    write(*,*)'########################################################'
    write(*,*)    
    call elapsed_time(start_count, count_rate)


    contains


    !-------------------------------------------------------------------------------------------------
    ! Subrutina para resolver la ecuación de Lindblad para un sistema cuántico de 1 qubit sometido
    ! a un Hamiltoniano H = Gamma * sigma_z, donde sigma_z es la matriz de Pauli Z.
    !       Entrada:
    !           init_vector: Vector de estado inicial.
    !           gamma: Constante de decaimiento.
    !           T_final: Tiempo final de evolución temporal.
    !           dt: Paso de evolución temporal.
    !           identifier: Identificador para el archivo de salida.
    !-------------------------------------------------------------------------------------------------
    subroutine main_lindblad(init_vector, gamma, T_final, dt, identifier)
        implicit none

        real(8), intent(in) :: init_vector(3)
        real(8), intent(in) :: gamma
        real(8), intent(in) :: T_final
        real(8), intent(in) :: dt
        character(len=*), intent(in) :: identifier

        real(8) :: vector(3)
        integer :: Ntime_steps
        integer :: t
        real(8) :: time

        Ntime_steps = ceiling(T_final / dt)

        open(unit=10, file='out_lindblad/cbp_ej9_lindblad_'//trim(adjustl(identifier))//'.dat', status='replace')
        write(10,*) '# gamma = ', gamma, ', || T_final = ', T_final, ', || dt = ', dt
        write(10,*) '# t,             x,             y,             z'

        vector = init_vector

        ! Evolución temporal
        do t = 0, Ntime_steps
            time = t * dt
            write(10,*) time, vector(:)

            call lindblad_rk4_step(vector, gamma, dt)
        end do

        close(10)
    
    end subroutine main_lindblad


    !-------------------------------------------------------------------------------------------------
    ! Subrutina obtener el siguiente paso con el algoritmo de Runge-Kutta de orden 4 del vector a 
    ! evolucionar.
    !       Entrada:
    !           vector: Vector de estado en el punto t.
    !           gamma: Constante de decaimiento.
    !           dt: Paso de evolución temporal.
    !       Salida:
    !           vector: Vector de estado en el punto t + dt.
    !-------------------------------------------------------------------------------------------------
    subroutine lindblad_rk4_step(vector, gamma, dt)
        implicit none

        real(8), intent(inout) :: vector(3)
        real(8), intent(in) :: gamma
        real(8), intent(in) :: dt

        real(8) :: k1(3), k2(3), k3(3), k4(3), vector_temp(3)

        ! Algoritmo de Runge-Kutta de orden 4
        call drdt_rk4(vector, gamma, k1)
        vector_temp = vector + k1 * dt / 2.0d0

        call drdt_rk4(vector_temp, gamma, k2)
        vector_temp = vector + k2 * dt / 2.0d0

        call drdt_rk4(vector_temp, gamma, k3)
        vector_temp = vector + k3 * dt

        call drdt_rk4(vector_temp, gamma, k4)

        vector = vector + (k1 + k2 * 2.0d0 + k3 * 2.0d0 + k4) * dt / 6.0d0
    
    end subroutine lindblad_rk4_step


    !-------------------------------------------------------------------------------------------------
    ! Subrutina para calcular el vector derivado con el algoritmo de Runge-Kutta de orden 4.
    !       Entrada:
    !           vector: Vector de estado en el punto t.
    !           gamma: Constante de decaimiento.
    !       Salida:
    !           drdt: Vector derivado en el punto t.
    !-------------------------------------------------------------------------------------------------
    subroutine drdt_rk4(vector, gamma, drdt)
        implicit none

        real(8), intent(in) :: vector(3)
        real(8), intent(in) :: gamma
        real(8), intent(out) :: drdt(3)

        drdt(1) = - gamma * vector(1) - 2.0d0 * vector(2) - gamma
        drdt(2) = 2.0d0 * vector(1) - gamma/2.0d0 * vector(2)
        drdt(3) = - gamma/2.0d0 * vector(3)

    end subroutine drdt_rk4


!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------------------------
    ! Subrutina principal para la simulación de Monte Carlo.
    !       Entrada:
    !           init_vector: Vector de estado inicial.
    !           gamma: Constante de decaimiento.
    !           T_final: Tiempo final de evolución temporal.
    !           dt: Paso de evolución temporal.
    !           Nsim: Número de simulaciones.
    !           identifier: Identificador para el archivo de salida.
    !-------------------------------------------------------------------------------------------------
    subroutine main_monte_carlo(init_vector, gamma, T_final, dt, Nsim, identifier)
        implicit none
        
        real(8), intent(in) :: init_vector(3)
        real(8), intent(in) :: gamma
        real(8), intent(in) :: T_final
        real(8), intent(in) :: dt
        integer, intent(in) :: Nsim
        character(len=*), intent(in) :: identifier

        integer :: Ntime_steps
        integer :: t
        real(8) :: time

        real(8), allocatable :: x_mean(:), y_mean(:), z_mean(:)
        

        complex(8) :: init_psi(2), psi(2)
        real(8) :: x, y, z, dp, epsilon, norm
        integer :: m

        Ntime_steps = ceiling(T_final / dt)

        ! Pasar el estado inicial de Bloch a vector de estado cuántico
        call bloch_to_state_vector(init_vector, init_psi)
        
        allocate(x_mean(0:Ntime_steps), y_mean(0:Ntime_steps), z_mean(0:Ntime_steps))

        x_mean = 0.0d0
        y_mean = 0.0d0
        z_mean = 0.0d0

        call random_seed()
        

        ! Trayectorias independientes
        do m = 1, Nsim
            ! Reiniciamos el estado cuántico al valor inicial
            psi = init_psi

            ! Evolución temporal de la trayectoria actual
            do t = 0, Ntime_steps
                
                ! Extraemos y acumulamos las coordenadas de Bloch de la trayectoria m
                x = 2.0d0 * dble(psi(1) * conjg(psi(2)))
                y = 2.0d0 * aimag(psi(2) * conjg(psi(1)))
                z = dble(psi(1) * conjg(psi(1))) - dble(psi(2) * conjg(psi(2)))

                x_mean(t) = x_mean(t) + x
                y_mean(t) = y_mean(t) + y
                z_mean(t) = z_mean(t) + z

                ! Probabilidad de salto cuántico en este dt
                dp = gamma * dt * 0.5d0 * cdabs(psi(1) + psi(2))**2

                call random_number(epsilon)

                ! Si epsilon es menor a dp, entonces hay un salto cuántico
                if (epsilon < dp) then
                    ! El espín colapsa instantáneamente al estado |->
                    psi(1) = cmplx(1.0d0, 0.0d0, 8) / dsqrt(2.0d0)
                    psi(2) = cmplx(-1.0d0, 0.0d0, 8) / dsqrt(2.0d0)
                else
                    ! En caso contrario se evoluciona con el paso de RK4
                    call monte_carlo_rk4_step(psi, gamma, dt)
                    
                    ! Normalización
                    norm = cdabs(psi(1))**2 + cdabs(psi(2))**2
                    psi = psi / dsqrt(norm)
                end if

            end do
        end do

        ! Promedio de las trayectorias
        x_mean = x_mean / real(Nsim, 8)
        y_mean = y_mean / real(Nsim, 8)
        z_mean = z_mean / real(Nsim, 8)


        open(unit=10, file='out_monte_carlo/cbp_ej9_monte_carlo_'//trim(adjustl(identifier))//'.dat', status='replace')
        write(10,*) '# gamma = ', gamma, ', || T_final = ', T_final, ', || dt = ', dt
        write(10,*) '# t,             x,             y,             z'

        do t = 0, Ntime_steps
            time = t * dt
            write(10,*) time, x_mean(t), y_mean(t), z_mean(t)
        end do

        close(10)

    end subroutine main_monte_carlo


    !-------------------------------------------------------------------------------------------------
    ! Subrutina para avanzar un paso temporal dt en el método Monte Carlo utilizando el algoritmo
    ! de Runge-Kutta de orden 4 para vectores de estado complejos (evolución continua no-Hermítica).
    !       Entrada:
    !           psi: Vector de estado complejo (alpha, beta) en el punto t.
    !           gamma: Constante de decaimiento.
    !           dt: Paso de evolución temporal.
    !       Salida:
    !           psi: Vector de estado complejo evolucionado en el punto t + dt (sin normalizar).
    !-------------------------------------------------------------------------------------------------
    subroutine monte_carlo_rk4_step(psi, gamma, dt)
        implicit none

        complex(8), intent(inout) :: psi(2)
        real(8), intent(in) :: gamma
        real(8), intent(in) :: dt

        complex(8) :: k1(2), k2(2), k3(2), k4(2)
        complex(8) :: psi_temp(2)

        call dpsidt_rk4(psi, gamma, k1)
        psi_temp = psi + k1 * dt / 2.0d0

        call dpsidt_rk4(psi_temp, gamma, k2)
        psi_temp = psi + k2 * dt / 2.0d0

        call dpsidt_rk4(psi_temp, gamma, k3)
        psi_temp = psi + k3 * dt

        call dpsidt_rk4(psi_temp, gamma, k4)

        psi = psi + (k1 + k2 * 2.0d0 + k3 * 2.0d0 + k4) * dt / 6.0d0
    
    end subroutine monte_carlo_rk4_step


    !-------------------------------------------------------------------------------------------------
    ! Subrutina para calcular las derivadas de Heff para el algoritmo de Runge-Kutta de orden 4
    !       Entrada:
    !           psi: Vector de estado complejo (alpha, beta) en el punto t.
    !           gamma: Constante de decaimiento.
    !       Salida:
    !           dpsidt: Derivadas de Heff en el punto t.
    !-------------------------------------------------------------------------------------------------
    subroutine dpsidt_rk4(psi, gamma, dpsidt)
        implicit none
        complex(8), intent(in) :: psi(2)
        real(8), intent(in) :: gamma
        complex(8), intent(out) :: dpsidt(2)

        dpsidt(1) = - cmplx(gamma/4.0d0,  1.0d0, 8) * psi(1) - cmplx(gamma/4.0d0, 0.0d0, 8) * psi(2)
        dpsidt(2) = - cmplx(gamma/4.0d0,  0.0d0, 8) * psi(1) - cmplx(gamma/4.0d0, -1.0d0, 8) * psi(2)

    end subroutine dpsidt_rk4


    !-------------------------------------------------------------------------------------------------
    ! Subrutina para transformar un vector de Bloch real (x,y,z) a un vector de estado complejo (alpha, beta)
    !       Entrada:
    !           r: Vector de Bloch: r(1)=x, r(2)=y, r(3)=z
    !       Salida:
    !           psi: Vector de estado complejo (alpha, beta)
    !-------------------------------------------------------------------------------------------------
    subroutine bloch_to_state_vector(r, psi)
        implicit none
        real(8), intent(in) :: r(3)        ! Vector de Bloch: r(1)=x, r(2)=y, r(3)=z
        complex(8), intent(out) :: psi(2)  ! Estado cuántico: psi(1)=alpha, psi(2)=beta
        
        real(8) :: x, y, z, cos_half

        x = r(1)
        y = r(2)
        z = r(3)

        ! Control del Polo Sur (z = -1) para evitar división por cero
        if (z <= -1.0d0 + 1.0d-12) then
            psi(1) = cmplx(0.0d0, 0.0d0, 8)
            psi(2) = cmplx(1.0d0, 0.0d0, 8)
        else
            cos_half = dsqrt((1.0d0 + z) / 2.0d0)
            psi(1) = cmplx(cos_half, 0.0d0, 8)
            
            psi(2) = cmplx(x, y, 8) / (2.0d0 * cos_half)
        end if
    end subroutine bloch_to_state_vector
    




    !-------------------------------------------------------------------------------------------------
    ! Subrutina para medir el tiempo real de ejecución del programa en el momento de su llamada
    !-------------------------------------------------------------------------------------------------
    subroutine elapsed_time(start_count, count_rate)
        implicit none

        ! Variables para medir el tiempo
        integer(8), intent(in) :: start_count, count_rate
        integer(8) :: end_count
        real(8) :: elapsed_total
        integer :: hours, minutes, seconds, milliseconds

        ! Obtenemos el conteo final de tiempo real
        call system_clock(end_count)
        elapsed_total = real(end_count - start_count, 8) / real(count_rate, 8)

        ! Desglose del tiempo
        hours = int(elapsed_total / 3600.0d0)
        minutes = int(mod(elapsed_total, 3600.0d0) / 60.0d0)
        seconds = int(mod(elapsed_total, 60.0d0))
        milliseconds = int((elapsed_total - int(elapsed_total)) * 1000.0d0)

        write(*,'(A, I2.2, A, I2.2, A, I2.2, A, I3.3, A)') '  -> Tiempo real de ejecucion: ', hours, ' h ', minutes, ' min ', seconds, ' s ', milliseconds, ' ms'
        write(*,*)'!--------------------------------------------------------'
        write(*,*)

    end subroutine elapsed_time


    !-------------------------------------------------------------------------------------------------
    ! Subrutina para imprimir la hora exacta de inicio
    !-------------------------------------------------------------------------------------------------
    subroutine print_actual_time(message_display)
        implicit none
        integer, intent(in) :: message_display
        integer :: values(8)
        character(len=8)  :: date
        character(len=10) :: time
        character(len=5)  :: zone

        ! date_and_time devuelve: año, mes, día, diferencia horaria, hora, min, seg, miliseg
        call date_and_time(date, time, zone, values)

        
        if (message_display == 0) then
            write(*,*) '--------------------------------------------------------'
            write(*, '(A, I2.2, A, I2.2, A, I2.2)') &
                '  -> Inicio de la ejecucion: ', values(5), ':', values(6), ':', values(7)
            write(*,*) '--------------------------------------------------------'
            write(*,*)
        else if(message_display == 1) then
            write(*, '(A, I2.2, A, I2.2, A, I2.2)') &
                '[-]> Hora de salida del mensaje: ', values(5), ':', values(6), ':', values(7)
        end if


    end subroutine print_actual_time

end program cbp_ej9_qinf