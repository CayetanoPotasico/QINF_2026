!
!   cbp_ej10_qinf.f90
!------------------------------------------------------------------------------------------------
!   3 - Examen QINF, Ejercicio 10
!   Introducción a la Información y Computación Cuánticas. Máster en Física Avanzada UNED.
!
!   Autor: Cayetano Bayona Pacheco, Junio 2026.
!------------------------------------------------------------------------------------------------
!
!   Programa principal para obtener estado fundamental del Hamiltoniano del problema propuesto
!   (ver archivo cbp_ej10_subs_qinf.f90) de un sistema de 8 qubits en función del parámetro Gamma 
!   con el fin de calcular la entropía media y su desviación estandard de los primeros 4 qubits 
!   mediante simulaciones Monte Carlo.
!       - Requiere las librerías LAPACK y BLAS para álgebra lineal.
!
!------------------------------------------------------------------------------------------------
!   Compilación y ejecución (opciones):
!       gfortran -g -fbacktrace -fcheck=all -fbounds-check cbp_ej10_subs_qinf.f90 cbp_ej10_qinf.f90 -llapack -lblas -fopenmp -o cbp_ej10_qinf.exe
!       gfortran cbp_ej10_subs_qinf.f90 cbp_ej10_qinf.f90 -llapack -lblas -fopenmp -o cbp_ej10_qinf.exe
!       gfortran -O3 -march=native -fopenmp -fno-stack-arrays cbp_ej10_subs_qinf.f90 cbp_ej10_qinf.f90 -llapack -lblas -o cbp_ej10_qinf.exe
!       ./cbp_ej10_qinf.exe
!------------------------------------------------------------------------------------------------

program cbp_ej10_qinf
    use mod_cbp_ej10_subs_qinf

    implicit none

    integer :: N = 8                               ! Número de qubits

    ! Variables para medir el tiempo
    integer(8) :: start_count, count_rate


    ! Hora de inicio de la ejecución
    call print_actual_time(0)

    ! Obtenemos la frecuencia del reloj (ticks por segundo) y el conteo inicial
    call system_clock(start_count, count_rate)

    ! Inicializamos el generador de números aleatorios
    call random_seed()



    call main(N, 2**N, 0.0d0, 12.0d0, 0.1d0, name='cbp_ej10_entropias_S15', Nsimulations=1000)



    write(*,*) 
    write(*,*)'########################################################'
    write(*,*)'            -> Programa terminado <-'
    write(*,*)'########################################################'
    write(*,*)    
    call elapsed_time(start_count, count_rate)

    contains

    ! --------------------------------------------------------------------------------------------------
    !   SUBRUTINAS
    ! --------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------
    ! Subrutina main para obtener el conjunto de valores obtenidos de entropía de los primeros 4 qubits
    ! del estado fundamental del Hamiltoniano para diferentes valores de parámetro Gamma.
    !       Entrada:
    !           N: Número de qubits.
    !           dim: Dimensión del espacio de Hilbert y del Hamiltoniano, dim = 2^N.
    !           Gamma_init: Parámetro Gamma inicial.
    !           Gamma_end: Parámetro Gamma final.
    !           dGamma: Incremento del parámetro Gamma.
    !           name: Nombre añadido del archivo de salida.
    !           Nsimulations: Número de simulaciones de Monte Carlo.
    !-------------------------------------------------------------------------------------------------
    subroutine main(N, dim, Gamma_init, Gamma_end, dGamma, name, Nsimulations)
        use omp_lib
        implicit none

        integer, intent(in) :: N
        integer, intent(in) :: dim
        real(8), intent(in) :: Gamma_init, Gamma_end, dGamma
        character(len=*), intent(in) :: name
        integer, intent(in) :: Nsimulations

        real(8) :: J_ij(N,N), H(dim,dim)        ! Matrices de interacciones y Hamiltoniano
        real(8) :: W(dim)                       ! Array para autovalores (W(1)=E0, W(2)=E1)

        ! Variables de trabajo para DSYEV
        real(8), allocatable :: work(:)
        real(8) :: dummy_work(1)
        integer :: lwork
        integer :: info

        real(8) ::  S_vector(0:2**N-1)          ! Vector de entropías de Von Neumann
        complex(8) :: ground_state(2**N)        ! Vector fundamental del Hamiltoniano

        ! Parámetros extra
        integer :: i, g
        real(8) :: Gamma
        integer :: NGamma_steps
        

        NGamma_steps = ceiling((Gamma_end - Gamma_init)/dGamma)


        ! Consulta previa a LAPACK para obtener el tamaño óptimo de WORK
        lwork = -1
        call dsyev('N', 'U', dim, H, dim, W, dummy_work, lwork, info)
        lwork = int(dummy_work(1))

        ! allocate(work(lwork))

        open(unit=10, file=name//'.dat', status='replace')
        write(10,*) '# Gamma                Simulacion                Entropia_S15'

        do g = 0, NGamma_steps ! Bucle de parámetros Gamma

            Gamma = Gamma_init + g*dGamma

            !$omp parallel default(none) &
            !$omp private(i, J_ij, H, W, work, info, ground_state, S_vector) &
            !$omp shared(N, dim, Gamma, Nsimulations, lwork)

            ! Cada hilo debe reservar su propia memoria para 'work'
            allocate(work(lwork))

            ! Cada hilo debe inicializar su semilla aleatoria para no repetir secuencias
            call init_random_seed_omp()

            !$omp do
            do i = 1, Nsimulations ! Bucle de simulaciones de Monte Carlo

                ! Generación de la matriz de interacciones
                call generate_J_matrix_uniform(J_ij, N)

                ! Generación del Hamiltoniano
                call generate_total_Hamiltonian(J_ij, N, H, dim, Gamma)

                ! Diagonalización del Hamiltoniano
                call dsyev('V', 'U', dim, H, dim, W, work, lwork, info)
                if (info /= 0) then
                    write(*,*) "Error en DSYEV PAR:", info
                    stop
                end if

                ! Obtención del autovector del estado fundamental
                ground_state(:) = H(:, 1)

                ! Cálculo del vector de entropías de Von Neumann
                call calculate_all_entropies(N, ground_state, S_vector)

                ! Guardar todos los datos en el archivo para estadística e histogramas
                ! La entropía de los 4 primeros qubits es la correspondiente al estado |00001111>, es decir, el índice 15
                !$omp critical
                write(10,*) Gamma, i, S_vector(15)
                !$omp end critical

            end do ! Fin bucle de simulaciones
    
            deallocate(work)

            !$omp end parallel

            write(*,*) '[+] Gamma = ', Gamma, 'Porcentaje = ', real(g)/real(NGamma_steps)*100, '%'

            
        end do ! Fin bucle de parámetros Gamma
        
        close(10)

        ! deallocate(work)


    end subroutine main

    
    !-------------------------------------------------------------------------------------------------
    ! Subrutina para inicializar la semilla de aleatoriedad de cada hilo
    !-------------------------------------------------------------------------------------------------
    subroutine init_random_seed_omp()
        use omp_lib
        integer :: n, i, thread_id
        integer, allocatable :: seed(:)
        call random_seed(size=n)
        allocate(seed(n))
        thread_id = omp_get_thread_num()
        do i = 1, n
            seed(i) = 12345 + thread_id * 6789 + i*10  ! Semilla distinta por hilo
        end do
        call random_seed(put=seed)
    end subroutine init_random_seed_omp


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

end program cbp_ej10_qinf
