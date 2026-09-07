!
!   cbp_main_qinf.f90
!------------------------------------------------------------------------------------------------
!   Ejercicio 3 - Adaptación del Ejercicio 1 al Ejercicio 3 (Examen). 
!   Introducción a la Información y Computación Cuánticas. Máster en Física Avanzada UNED.
!
!   Autor: Cayetano Bayona Pacheco, Marzo-Junio 2026.
!------------------------------------------------------------------------------------------------
!   
!   Programa principal para obtener la matriz de adyacencia J_ij de un sistema cuántico de 
!   N qubits. Para ello se calculan las entropías de las diferentes particiones del sistema a
!   partir del vector de estado psi y se resuelve el sistema de ecuaciones lineales de solución
!   J_ij. El sistema proviene de un artículo científico, y se puede encontrar en arXiv:1906.05146.
!
!   ---> ES MUY IMPORTANTE CREAR LA CARPETA /out ANTES DE EJECUTAR EL PROGRAMA PARA QUE FORTRAN
!   ---> EXPORTE LOS DATOS SIN DAR ERRORES.
!
!------------------------------------------------------------------------------------------------
!   Compilación y ejecución (opciones):
!       gfortran cbp_subs_qinf.f90 cbp_main_qinf.f90 -llapack -lblas -fopenmp -o cbp_main_qinf.exe
!       gfortran -O3 -march=native -fopenmp -fno-stack-arrays cbp_subs_qinf.f90 cbp_main_qinf.f90 -llapack -lblas -o cbp_main_qinf.exe
!       ./cbp_main_qinf.exe
!------------------------------------------------------------------------------------------------

program cbp_main_qinf
    
    implicit none
    integer(kind=4) :: Nqubits                           ! Número de qubits del sistema
    integer(kind=4) :: N                                 ! Número de coeficientes: dimensión del espacio de Hilbert
    integer(kind=4) :: Np                                ! Número de enlaces / términos independientes de J_ij

    ! Variables para medir el tiempo
    integer(8) :: start_count, count_rate

    ! Hora de inicio de la ejecución
    call print_actual_time(0)

    ! Obtenemos la frecuencia del reloj (ticks por segundo) y el conteo inicial
    call system_clock(start_count, count_rate)

    ! Inicializamos el generador de números aleatorios
    call random_seed()

    ! Cálculos para los estados pedidos en el enunciado y un aleatorio de Nqubits = 8
    Nqubits = 12
    N = 2**Nqubits
    Np = Nqubits*(Nqubits-1)/2


    ! Cálculos para el estado numérico 1
    call main(Nqubits, N, Np, "wavefunction.txt", "")


    call print_actual_time(1)
    call elapsed_time(start_count, count_rate)


    write(*,*) 
    write(*,*)'########################################################'
    write(*,*)'            -> Programa terminado <-'
    write(*,*)'########################################################'
    write(*,*)    
    call elapsed_time(start_count, count_rate)

contains

    !------------------------------------------------------------------------------------------------
    ! Subrutina principal de ejecución para el análisis de un estado cuántico específico.
    ! Coordina la lectura, el cálculo de entropías, la resolución del sistema y el guardado.
    !       Entrada:
    !           Nqubits:           Número de qubits del sistema.
    !           N:                 Dimensión del espacio de Hilbert (N = 2^Nqubits).
    !           Np:                Número de pares de acoplamiento (Np = Nqubits*(Nqubits-1)/2).
    !           psi_file:          Ruta del archivo con el vector de estado.
    !           psi_savefile_indx: Índice numérico para identificar el archivo de salida.
    !-------------------------------------------------------------------------------------------------
    subroutine main(Nqubits, N, Np, psi_file, psi_savefile_type)
        use mod_cbp_subs_qinf
        implicit none

        integer(kind=4), intent(in) :: Nqubits, N, Np
        character(len=*), intent(in) :: psi_file
        character(len=*), intent(in) :: psi_savefile_type
        character(len=3) :: Nqubits_str

        ! Vectores de trabajo
        complex(kind=8), dimension(0:N-1) :: psi
        real(kind=8), dimension(0:N-1) :: S_vector
        real(kind=8), dimension(Np) :: b
        real(kind=8), dimension(Np) :: J_vector
        integer(kind=4) :: i, j

        ! Valores obtenidos del vector phi tal que |psi> = \bigo_i=1^12 |phi>
        ! |phi> = alpha |0> + beta |1>
        complex(8) :: gamma, alpha, beta

        real(8) :: norm

        write(Nqubits_str,'(i3)') Nqubits

        write(*,*)
        write(*,*) '########################################################################'
        write(*,*) '- Iniciando calculos para ' // psi_file
        
        ! Lectura de coeficientes del estado psi
        call state_vector(psi_file, N, psi)
        write(*,*) '[+] Estado psi leido correctamente.'

        norm = sqrt(sum(abs(psi)**2))
        write(*,*) '[+] Norma del estado psi: ', norm

        gamma = psi(1) / psi(0)
        alpha = 1.0d0 / sqrt(1 + gamma**2)
        beta = alpha * gamma
        open(unit=10, file='out/gamma_alpha_beta_'//psi_savefile_type//'_'//trim(adjustl(Nqubits_str))//'.dat', status='replace')
        write(10,*) '# gamma,         gamma^2,         alpha,         beta'
        write(10,*) gamma, gamma**2, alpha, beta
        close(10)

        ! Cálculo de entropías de Von Neumann para todas las particiones (2^N)
        call calculate_all_entropies(Nqubits, psi, S_vector)
        write(*,*) '[+] Entropias calculadas correctamente.'

        open(unit=10, file='out/entropias_'//psi_savefile_type//'_'//trim(adjustl(Nqubits_str))//'.dat', status='replace')
        write(10,*) '# i    S_i'
        do i = 0, N - 1
            write(10,*) i, S_vector(i)
        end do
        close(10)

        ! Verificación de la Ley de Áreas
        call analyze_area_law(Nqubits, S_vector, psi_savefile_type)

        ! Construcción del término independiente b = A^T * S
        call independent_term_calc(Nqubits, Np, S_vector, b)
        write(*,*) '[+] Vector b calculado correctamente.'

        ! Resolución del sistema de ecuaciones lineales
        write(*,*) '[-] Resolviendo el sistema lineal con LAPACK (DGESV)...'
        call solve_linear_system(Nqubits, Np, b, J_vector)

        ! Guardado de la matriz J_ij
        open(unit=20, file='out/matriz_J_'//psi_savefile_type//'_'//trim(adjustl(Nqubits_str))//'.dat', status='replace')
        
        write(20,*) "# Matriz de acoplamientos J(i,j) para "//psi_savefile_type//"con Nqubits = "//trim(adjustl(Nqubits_str))
        do i = 0, Nqubits - 1
            do j = 0, Nqubits - 1
                ! Rellenamos la matriz completa usando la simetría de los acoplamientos
                write(20, '(F12.6)', advance='no') get_J_value(i, j, Nqubits, J_vector)
            end do
            write(20,*) ! Salto de línea para formato de matriz
        end do
        close(20)
        
        write(*,*) "[+] Matriz J_"//psi_savefile_type//" guardada con exito."
        write(*,*)

    end subroutine main


    !------------------------------------------------------------------------------------------------
    ! Función auxiliar para mapear el vector unidimensional J_vec a una matriz bidimensional.
    ! Recupera el valor J_ij basándose en un ordenamiento triangular (i < j).
    !       Entrada:
    !           i, j:  Índices de los qubits de la matriz (0 a Nqubits-1).
    !           n:     Número total de qubits.
    !           J_vec: Vector de Np = Nqubits*(Nqubits-1)/2 componentes con los resultados de LAPACK.
    !       Salida:
    !           val:  Valor del acoplamiento J entre el qubit i y el j.
    !-------------------------------------------------------------------------------------------------
    function get_J_value(i, j, n, J_vec) result(val)
        implicit none
        integer, intent(in) :: i, j, n
        real(kind=8), intent(in) :: J_vec(:)
        real(kind=8) :: val

        integer :: p, r, c, k
        
        ! J_ii = 0
        if (i == j) then
            val = 0.0d0
        else
            ! Garantizamos que r sea el índice menor para cumplir la lógica triangular
            r = min(i, j)
            c = max(i, j)
            
            ! Cálculo del índice lineal p correspondiente al par (r, c)
            ! p = sumatoria de longitudes de filas anteriores + desplazamiento en fila actual
            p = 0
            do k = 0, r - 1
                p = p + (n - 1 - k)
            end do
            p = p + (c - r)
            val = J_vec(p)
        end if
    end function get_J_value



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

end program cbp_main_qinf