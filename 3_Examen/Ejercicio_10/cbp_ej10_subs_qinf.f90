!
!   cbp_ej10_subs_qinf.f90
!--------------------------------------------------------------------------------------------------
!   3 - Examen QINF, Ejercicio 10
!   Introducción a la Información y Computación Cuánticas. Máster en Física Avanzada UNED.
!
!   Autor: Cayetano Bayona Pacheco, Junio 2026.
!--------------------------------------------------------------------------------------------------
!   
!   Archivo (módulo) de subrutinas y funciones utilizadas en el programa principal cbp_ej10_qinf.f90.
!   - Requiere las librerías LAPACK y BLAS para álgebra lineal.
!
!   CONTENIDO DEL MÓDULO (mod_cbp_subs_qinf):
!       - Generación de la matriz de enlaces J_ij
!       - Generación de la matriz total del Hamiltoniano H
!       - Subrutinas para el cálculo del vector de entropías de Von Neumann
!
!--------------------------------------------------------------------------------------------------

module mod_cbp_ej10_subs_qinf
    implicit none
contains

    !-------------------------------------------------------------------------------------------------
    ! Subrutina para generar la matriz simétrica J_ij con valores aleatorios de una distribución 
    ! uniforme entre 0 y 1. Ponemos a 0 los valores J_ii.
    !       Entrada:
    !           N: Dimensión de la matriz cuadrada.
    !       Salida:
    !           J_ij: Matriz cuadrada J_ij con valores aleatorios de una distribución uniforme.
    !-------------------------------------------------------------------------------------------------
    subroutine generate_J_matrix_uniform(J_ij, N)
        implicit none

        integer, intent(in) :: N
        real(8), intent(out) :: J_ij(N,N)

        integer :: i, j
        real(8) :: rnd1


        J_ij = 0.0d0

        do i = 1, N
            do j = i, N
                call random_number(rnd1)

                J_ij(i,j) = rnd1

                J_ij(j,i) = J_ij(i,j)
            end do
        end do


        do i = 1, N
            J_ij(i,i) = 0.0d0
        end do

    end subroutine generate_J_matrix_uniform


    !-------------------------------------------------------------------------------------------------
    ! Subrutina para generar la matriz total del Hamiltoniano H de dimensión 2^N x 2^N. Este se compone
    ! de: H = Gamma * H0 - H1 donde:
    !       H0 = Sum_i^N sigma_i^x
    !       H1 = Sum_ij^N J_ij * sigma_i^z * sigma_j^z
    ! siendo 
    !       sigma_i^x = (I x ... x |0><1| + |1><0| x ... x I)
    !       sigma_i^z = (I x ... x |0><0| - |1><1| x ... x I) 
    ! las matrices de Pauli de dimensión 2 para el qubit i.
    !
    !       Entrada:
    !           J_ij: Matriz cuadrada J_ij de conexiones, dimensión N x N.
    !           N: Número de qubits y dimensión de J_ij.
    !           dim: Dimensión del espacio de Hilbert y del Hamiltoniano, dim = 2^N.
    !           Gamma: Parámetro Gamma de control.
    !       Salida:
    !           H: Matriz cuadrada H de dimensión 2^N x 2^N.
    !-------------------------------------------------------------------------------------------------
    subroutine generate_total_Hamiltonian(J_ij, N, H, dim, Gamma)
        implicit none

        integer, intent(in) :: N
        real(8), intent(in) :: J_ij(N,N)
        integer, intent(in) :: dim
        real(8), intent(out) :: H(dim,dim)
        real(8), intent(in) :: Gamma

        integer :: i, j, k
        real(8) :: energy_k!, s_i, s_j
        integer :: neighbour

        real(8) :: spins(N)
        real(8) :: frac1, frac0


        frac1 = -1.0d0
        frac0 = Gamma

        H = 0.0d0

        ! Bucle de 2^N estados |k>
        do k = 0, dim-1

            ! Precalcular los espines
            do i = 0, N-1
                spins(i+1) = merge(-1.0d0, 1.0d0, btest(k, i))
            end do

            
            energy_k = 0.0d0
            do i = 0, N-1

                ! Aprovechamos que J_ij es simétrica y que J_ii = 0
                do j = i + 1, N-1

                    energy_k = energy_k + J_ij(j+1,i+1) * spins(i+1) * spins(j+1)

                end do
                
                ! --- Parte Off-Diagonal (H0)
				! neighbour = ieor(k, 2**i)
                neighbour = ieor(k, ishft(1, i))   ! Para mayor eficiencia, calculamos 2^i con un desplazamiento de bits
                H(k+1, neighbour+1) = frac0
                
            end do

			! --- Parte Diagonal (H1)
            H(k+1, k+1) = energy_k * frac1


        end do


    end subroutine generate_total_Hamiltonian







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

        ! Recorremos solo la mitad de las configuraciones (0 a 2^(Nqubits-1) - 1) ya que la entropía de una partición es idéntica a la de su complemento
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
                ! Lo encendemos en la posición real del qubit
                long_idx = ibset(long_idx, mapping(b+1))
            end if
        end do
    end function scatter_bits



end module mod_cbp_ej10_subs_qinf
