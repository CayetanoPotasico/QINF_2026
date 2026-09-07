!
!   cbp_subs_qinf.f90
!------------------------------------------------------------------------------------------------
!   Ejercicio 1 - Discusión de Artículo Científico. 
!   Introducción a la Información y Computación Cuánticas. Máster en Física Avanzada UNED.
!
!   Autor: Cayetano Bayona Pacheco, Marzo-Junio 2026.
!------------------------------------------------------------------------------------------------
!   
!   Archivo (módulo) de subrutinas y funciones utilizadas en el programa principal cbp_1_0_qinf.f90
!   - Requiere las librerías LAPACK y BLAS para álgebra lineal.
!
!   CONTENIDO DEL MÓDULO (mod_cbp_subs_qinf):
!       - analyze_area_law:         Estudio del escalado de la entropía (Ley de Áreas).
!       - compute_reduced_rho:      Cálculo de rho_A mediante traza parcial.
!       - calculate_all_entropies:  Gestión del barrido de las 256 particiones.
!       - solve_linear_system:      Solver lineal mediante rutina DGESV (LAPACK).
!       - get_von_neumann_entropy:  Cálculo de S mediante rutina ZHEEV (LAPACK).
!       - analyze_area_law:         Estudio del escalado de la entropía (Ley de Áreas).
!       - count_set_bits:           Función auxiliar de conteo de qubits activos.
!       - scatter_bits:             Función auxiliar de mapeo de índices globales.
!       - state_vector:             Lectura de coeficientes del estado psi.
!
!------------------------------------------------------------------------------------------------

module mod_cbp_subs_qinf
    implicit none
contains

    !------------------------------------------------------------------------------------------------
    ! Subrutina para realizar el estudio del escalado de la entropía (Ley de Áreas) mediante el
    ! promedio de las entropías de Von Neumann para distintos tamaños de bloque L.
    !       Entrada:
    !           Nqubits:  Número total de qubits del sistema.
    !           S_vector: Vector con las entropias de Von Neumann para todas las particiones.
    !           psi_savefile_type: Identificador del tipo de archivo para guardar los datos.
    !-------------------------------------------------------------------------------------------------
    subroutine analyze_area_law(Nqubits, S_vector, psi_savefile_type)
        implicit none
        integer(kind=4), intent(in) :: Nqubits
        real(kind=8), intent(in) :: S_vector(0:2**Nqubits-1)
        
        integer(kind=4) :: L, start_bit, mask, i
        real(kind=8) :: avg_S
        integer(kind=8) :: count

        character(len=*), intent(in) :: psi_savefile_type
        character(len=3) :: Nqubits_str

        write(Nqubits_str,'(i3)') Nqubits

        ! write(*,*) '    ----------------------------------------------------------------'
        ! write(*,*) '    Analisis de Ley de Areas (Entropia por tamanyo de bloque L)'
        ! write(*,*) '      L (Tamanyo)    <S_L> (Entropia Media)'
        ! write(*,*) '    ----------------------------------------------------------------'

        
        open(unit=10, file='out/mean_S_per_L_'//psi_savefile_type//'_'//trim(adjustl(Nqubits_str))//'.dat', status='replace')
        write(10, *) '# L (Tamanyo)    <S_L> (Entropia Media)'
        do L = 1, Nqubits - 1
            
            avg_S = 0.0d0
            count = 0
            
            ! Buscamos bloques de tamaño L que sean contiguos (0,1,...,L-1), (1,2,...,L), etc.
            do start_bit = 0, Nqubits - L
                mask = 0
                do i = 0, L - 1
                    mask = ibset(mask, start_bit + i)
                end do
                
                avg_S = avg_S + S_vector(mask)
                count = count + 1
            end do
            
            avg_S = avg_S / real(count, 8)
            write(10, *) L, avg_S
        end do
        close(10)

    end subroutine analyze_area_law


    !------------------------------------------------------------------------------------------------
    ! Subrutina para resolver el sistema lineal M * J = b mediante LAPACK. M = A^T * A y b = A^T * S
    ! Obtiene las constantes de acoplamiento J_ij minimizando el error cuadrático.
    !       Entrada:
    !           Nqubits:  Numero total de qubits del sistema.
    !           Np:       Numero de pares únicos o enlaces Np = Nqubits*(Nqubits-1)/2.
    !           b:        Vector de términos independientes (A^T * S).
    !       Salida:
    !           J_vector: Elementos de J_ij en forma de vector.
    !-------------------------------------------------------------------------------------------------
    subroutine solve_linear_system(Nqubits, Np, b, J_vector)
        implicit none
        integer(kind=4), intent(in) :: Nqubits
        integer(kind=4), intent(in) :: Np
        real(kind=8), intent(in)  :: b(Np)
        real(kind=8), intent(out) :: J_vector(Np)

        ! Variables internas para el solver
        real(kind=8) :: M(Np, Np)               ! Matriz de diseño analitica (A^T * A)
        integer(kind=4) :: ipiv(Np)             ! Vector de pivotes para factorizacion LU
        integer(kind=4) :: info, i, j

        real(kind=8) :: coef_diag, coef_offdiag

        ! Declaracion de la subrutina externa de LAPACK
        external dgesv

        ! Cálculo de los coeficientes de la matriz M según la combinatoria de particiones
        ! coef_diag: particiones que rompen un enlace i-j
        ! coef_offdiag: particiones que rompen simultaneamente dos enlaces
        coef_diag = 2.0d0**(Nqubits - 1)
        coef_offdiag = 2.0d0**(Nqubits - 2)

        ! Construcción de la matriz de coeficientes M
        do i = 1, Np
            do j = 1, Np
                if (i == j) then
                    M(i, j) = coef_diag
                else
                    M(i, j) = coef_offdiag
                end if
            end do
        end do

        ! Preparación del vector de salida
        J_vector = b

        ! Llamada a LAPACK (DGESV)
        ! Resuelve el sistema M * J = b
        call dgesv(Np, 1, M, Np, ipiv, J_vector, Np, info)

        ! Verificación
        if (info == 0) then
            write(*,*) "    -> Sistema lineal resuelto con exito."
        else
            write(*,*) "    -> Error en DGESV: El sistema es singular o mal condicionado. INFO =", info
        end if

    end subroutine solve_linear_system


    !------------------------------------------------------------------------------------------------
    ! Subrutina para calcular la entropía de Von Neumann S = -sum_i (\lambda_i * log(\lambda_i)).
    ! Utiliza la subrutina ZHEEV de LAPACK para diagonalizar la matriz de densidad reducida.
    !       Entrada:
    !           dim: Dimensión de la matriz reducida rho_A (2^n_A).
    !           matrix: Matriz de densidad reducida rho_A (Hermítica).
    !       Salida:
    !           S: Valor escalar de la entropía de Von Neumann resultante.
    !-------------------------------------------------------------------------------------------------
    subroutine get_von_neumann_entropy(dim, matrix, S)
        implicit none
        integer(kind=4), intent(in) :: dim
        complex(kind=8), intent(in) :: matrix(0:dim-1, 0:dim-1)
        real(kind=8), intent(out) :: S

        ! Parámetros internos para la comunicación con LAPACK (ZHEEV)
        complex(kind=8) :: A(dim, dim)            ! Copia de trabajo de la matriz
        real(kind=8)    :: W(dim)                 ! Vector donde LAPACK devuelve los autovalores
        complex(kind=8), allocatable :: WORK(:)
        real(kind=8), allocatable    :: RWORK(:)
        integer(kind=4) :: INFO, LWORK, i

        ! Inicialización y copia de seguridad
        ! ZHEEV sobrescribe la matriz de entrada con vectores propios, por eso usamos 'A'
        A = matrix
        S = 0.0d0

        ! Configuración del espacio de trabajo
        LWORK = max(1, 2*dim-1)
        allocate(WORK(LWORK), RWORK(max(1, 3*dim-2)))

        ! 'N' = No calcular autovectores
        ! 'U' = Usar triángulo superior
        call zheev('N', 'U', dim, A, dim, W, WORK, LWORK, RWORK, INFO)

        ! Cálculo de entropía
        if (INFO == 0) then
            do i = 1, dim
                ! Se ignoran autovalores nulos o negativos por ruido
                if (W(i) > 1.0d-14) then
                    S = S - W(i) * log(W(i))
                end if
            end do
        else
            write(*,*) "¡Error en ZHEEV! No se pudo diagonalizar la matriz. INFO =", INFO
        end if

        deallocate(WORK, RWORK)
        
    end subroutine get_von_neumann_entropy


    !------------------------------------------------------------------------------------------------
    ! Subrutina para calcular el vector completo de entropías de Von Neumann.
    ! Gestiona el bucle sobre las particiones y aprovecha la simetría S_A = S_B.
    !       Entrada:
    !           Nqubits: Número total de qubits del sistema.
    !           psi:      Vector de estado global.
    !       Salida:
    !           S_vector: Vector con las entropías de las 2^Nqubits posibles particiones.
    !-------------------------------------------------------------------------------------------------
    subroutine calculate_all_entropies(Nqubits, psi, S_vector)
        implicit none
        integer(kind=4), intent(in) :: Nqubits
        complex(kind=8), intent(in) :: psi(0:2**Nqubits-1)
        real(kind=8), intent(out) :: S_vector(0:2**Nqubits-1)

        integer(kind=4) :: I, complement, dim_A, n_A, b
        complex(kind=8), allocatable :: rho_A(:,:)
        real(kind=8) :: entropy

        ! Recorremos solo la mitad de las configuraciones (0 a 2^(Nqubits-1) - 1) ya que la entropía de una partición es idéntica a la de su complemento.
        !$omp parallel do private(I, n_A, b, dim_A, rho_A, entropy, complement)
        do I = 0, 2**(Nqubits-1) - 1

            ! Determinamos el número de qubits en el subsistema A y su dimensión asociada
            n_A = 0
            do b = 0, Nqubits - 1
                if (btest(I, b)) n_A = n_A + 1
            end do

            dim_A = 2**n_A

            ! Caso trivial: si tenemos el conjunto vacío o el sistema completo entonces la entropía es cero
            if (n_A == 0) then
                entropy = 0.0d0
            else

                ! Inicialización de la matriz de densidad reducida
                allocate(rho_A(0:dim_A-1, 0:dim_A-1))

                ! Tomamos la traza parcial sobre el subsistema B
                call compute_reduced_rho(Nqubits, psi, I, rho_A, dim_A)

                ! Obtención de la entropía mediante diagonalización (LAPACK) ya que S = -Tr(rho*log(rho)) = - sum_i \lambda_i * log(\lambda_i)
                call get_von_neumann_entropy(dim_A, rho_A, entropy)

                deallocate(rho_A)
            end if

            ! Aplicación de la simetría S(A) = S(B)
            S_vector(I) = entropy
            complement = (2**Nqubits - 1) - I
            S_vector(complement) = entropy

        end do
        !$omp end parallel do

    end subroutine calculate_all_entropies


    !------------------------------------------------------------------------------------------------
    ! Subrutina para calcular la matriz de densidad reducida mediante traza parcial.
    ! Implementa rho_A = Tr_B(|ψ><ψ|) = sum_b (<bxa| |ψ><ψ| |bxa>) donde b \in B y a \in A
    !       Entrada:
    !           Nqubits:  Número total de qubits del sistema.
    !           psi:      Vector de estado global (2^Nqubits componentes).
    !           mask_A:   Máscara de bits que define el subsistema A. (Será 0,1,...,2^Nqubits-1)
    !           dim_A:    Dimensión de la matriz reducida (2^n_A).
    !       Salida:
    !           rho_A:    Matriz de densidad reducida para el subsistema A.
    !-------------------------------------------------------------------------------------------------
    subroutine compute_reduced_rho(Nqubits, psi, mask_A, rho_A, dim_A)
        implicit none
        integer(kind=4), intent(in) :: Nqubits
        complex(kind=8), intent(in) :: psi(0:2**Nqubits-1)
        integer(kind=4), intent(in) :: mask_A
        integer(kind=4), intent(in) :: dim_A
        complex(kind=8), intent(out) :: rho_A(0:dim_A-1, 0:dim_A-1)

        integer(kind=4) :: i, j, k, n_A, n_B, b
        integer(kind=4) :: row_idx_global, col_idx_global
        integer(kind=4) :: map_A(30), map_B(30)

        integer(kind=4) :: scatter_bits_i_A, scatter_bits_j_A, scatter_bits_k_B

        ! Clasificación de qubits: separamos qué índices (de 0 a Nqubits-1) pertenecen a A y cuáles a B
        ! allocate(map_A(Nqubits), map_B(Nqubits))

        n_A = 0
        n_B = 0

        ! Crea los conjuntos disjuntos A y B. p.ej: A = {0, 2, 4}, B = {1, 3, 5} para Nqubits = 6
        do b = 0, Nqubits - 1
            if (btest(mask_A, b)) then
                n_A = n_A + 1
                map_A(n_A) = b
            else
                n_B = n_B + 1
                map_B(n_B) = b
            end if
        end do

        ! Obtención de la matriz de densidad reducida rho_A
        rho_A = 0.0d0

        do j = 0, dim_A - 1
            scatter_bits_j_A = scatter_bits(j, n_A, map_A)
            
            do i = 0, dim_A - 1
                scatter_bits_i_A = scatter_bits(i, n_A, map_A)
                
                do k = 0, 2**n_B - 1
                    scatter_bits_k_B = scatter_bits(k, n_B, map_B)
                    
                    row_idx_global = scatter_bits_i_A + scatter_bits_k_B
                    col_idx_global = scatter_bits_j_A + scatter_bits_k_B

                    rho_A(i, j) = rho_A(i, j) + psi(row_idx_global) * conjg(psi(col_idx_global))
                end do
            end do
        end do

        ! deallocate(map_A, map_B)
    end subroutine compute_reduced_rho


    !------------------------------------------------------------------------------------------------
    ! Función para dispersar los bits de un índice local a sus posiciones globales.
    ! Mapea índices de subsistemas reducidos al espacio de Hilbert global de 2^N.
    !       Entrada:
    !           short_idx: Índice en el espacio reducido (ej. del subsistema B).
    !           n_bits:    Número de qubits en dicho subsistema.
    !           mapping:   Vector con las posiciones reales de los qubits en la cadena.
    !       Salida:
    !           long_idx:  Índice con los bits colocados en sus posiciones originales.
    !-------------------------------------------------------------------------------------------------
    function scatter_bits(short_idx, n_bits, mapping) result(long_idx)
        implicit none
        integer(kind=4), intent(in) :: short_idx, n_bits
        integer(kind=4), intent(in) :: mapping(n_bits)
        integer(kind=4) :: long_idx, b

        long_idx = 0

        ! Recorremos cada bit del índice local (0, 1, 2, ..., n_bits)
        do b = 0, n_bits - 1
            ! Si el bit 'b' está encendido en el índice corto
            if (btest(short_idx, b)) then
                ! Lo encendemos en la posición real del qubit (mapping)
                long_idx = ibset(long_idx, mapping(b+1))
            end if
        end do
    end function scatter_bits


    !------------------------------------------------------------------------------------------------
    ! Subrutina para calcular el vector b = A^T * S completo (término independiente).
    ! Recorre todos los pares únicos de qubits (enlaces) y calcula su componente b_p.
    !       Entrada:
    !           Nqubits: Número total de qubits del sistema.
    !           Np: Número total de pares únicos o enlaces.
    !           S_vector: Vector con las 2^Nqubits entropías de todas las particiones.
    !       Salida:
    !           b: Vector de dimensión Np con los términos independientes para el sistema.
    !------------------------------------------------------------------------------------------------
    subroutine independent_term_calc(Nqubits, Np, S_vector, b)
        implicit none
        integer(kind=4), intent(in) :: Nqubits 
        integer(kind=4), intent(in) :: Np      
        real(kind=8), dimension(0:2**Nqubits-1), intent(in) :: S_vector 
        real(kind=8), dimension(Np), intent(out) :: b

        integer(kind=4) :: i, j, p

        ! Inicializamos el vector b a cero
        b = 0.0d0
        ! p: índice lineal que identifica cada enlace (de 1 a Np)
        p = 1

        ! Recorremos todos los pares únicos (i, j), 0 <= i < j < Nqubits
        do i = 0, Nqubits - 2
            do j = i + 1, Nqubits - 1
                
                ! Cálculo de la componente b_p para el par p
                call A_T_times_S_single_pair(i, j, Nqubits, S_vector, b(p))

                p = p + 1
            end do
        end do

    end subroutine independent_term_calc


    !------------------------------------------------------------------------------------------------
    ! Subrutina para calcular b_p = sum_I (A_p,I * S_I) para un par específico (i, j).
    ! Determina qué particiones rompen el enlace entre el qubit i y el j.
    !       Entrada:
    !           i, j: Índices de los qubits que forman el par (0 a Nqubits-1).
    !           Nqubits: Número total de qubits del sistema.
    !           S_vector: Vector con las 2^Nqubits entropías calculadas previamente. (Una por cada partición)
    !       Salida:
    !           single_b: Valor acumulado de la componente del vector independiente b.
    !------------------------------------------------------------------------------------------------
    subroutine A_T_times_S_single_pair(i, j, Nqubits, S_vector, single_b)
        implicit none
        integer(kind=4), intent(in) :: i, j, Nqubits
        real(kind=8), dimension(0:2**Nqubits-1), intent(in) :: S_vector
        real(kind=8), intent(out) :: single_b

        integer(kind=4) :: I_indx
        logical :: bit_i, bit_j
        real(kind=8) :: local_sum       ! Añadido para paralelización

        ! Inicializamos la suma para este par (i, j)
        local_sum = 0.0d0
        ! single_b = 0.0d0

        ! Recorremos todas las 2^Nqubits particiones posibles
        !$omp parallel do reduction(+:local_sum) private(I_indx, bit_i, bit_j)
        do I_indx = 0, 2**Nqubits - 1
            
            ! btest(I_indx, k) es .true. si el qubit k está en el subsistema A
            bit_i = btest(I_indx, i)
            bit_j = btest(I_indx, j)

            ! Si los qubits están en lados distintos (.neqv.), el enlace se rompe y por tanto A(I, p) = 1. En otro caso A(I, p) = 0
            if (bit_i .neqv. bit_j) then
                ! single_b = single_b + S_vector(I_indx)
                local_sum = local_sum + S_vector(I_indx)
            end if
        end do
        !$omp end parallel do

        single_b = local_sum        ! Asignamos el resultado a single_b
    end subroutine A_T_times_S_single_pair


    !------------------------------------------------------------------------------------------------
    ! Subrutina para obtener el vector de estado |psi⟩ a partir de los archivos.
    !       Entrada:
    !           file_name: Nombre del archivo con los coeficientes del estado |psi⟩
    !           N: Número de coeficientes del estado |psi⟩
    !       Salida:
    !           psi: Coeficientes (complejos) del estado |psi⟩ en la base computacional
    !------------------------------------------------------------------------------------------------
    subroutine state_vector(file_name, N, psi)
        implicit none
        character(len=*), intent(in) :: file_name
        integer(kind=4), intent(in) :: N
        complex(kind=8), intent(out) :: psi(N)

        integer(kind=4) :: i
        character(len=9) :: dummy
        real(kind=8) :: psi_re, psi_im

        ! Leemos los valores numéricos de los coeficientes del estado |psi⟩
        open(unit=10, file=file_name, status='old', action='read')

        ! Saltamos la cabecera
        read(10,*)
        read(10,*)

        do i = 1, N
            read(10,*) dummy, psi_re, psi_im

            psi(i) = complex(psi_re, psi_im)

            ! write(*,*) i, psi(i)
        end do
        close(10)
        
    end subroutine state_vector


end module mod_cbp_subs_qinf
