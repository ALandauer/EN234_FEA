!     Subroutines for basic 3D hyperelastic elements

!==========================SUBROUTINE el_linelast_hyperelastic ==============================
subroutine el_linelast_hyperelastic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Mesh, only : extract_node_data
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
!    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D

    implicit none

    integer, intent( in )         :: lmn                                        ! Element number
    integer, intent( in )         :: element_identifier                         ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                    ! # nodes on the element
    integer, intent( in )         :: n_properties                               ! # properties for the element
    integer, intent( in )         :: length_coord_array                         ! Total # coords
    integer, intent( in )         :: length_dof_array                           ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint,ii,jj,kk

    !    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  mu1, K1                           ! Material properties
    real (prec)  ::  F(3,3),u(3,n_nodes),dNdy(n_nodes,3),G(6,9),tau(3,3),B_star(9,3*n_nodes)
    real (prec)  ::  B_def(3,3),Binv(3,3),Finv(3,3),S(3,length_dof_array/3)
    real (prec)  ::  Pmat(3*n_nodes,3*n_nodes),Smat(3*n_nodes,3*n_nodes),Sigma(3*n_nodes,3*n_nodes)
    real (prec)  ::  eye(3,3),eye6(6,6),Ivec(6),Bvec_inv(6),Bvec(6),Svec(3*n_nodes),Pvec(3*n_nodes)
    real (prec)  ::  J,Bdet,Bkk
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         shear modulus
    !     element_properties(2)         bulk modulus

    fail = .false.
    
    x = 0.d0
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    !    nn = 1
    !    call extract_node_data(nn,flag,n_coords,nodal_coords,n_dof,nodal_dof_increment,nodal_dof_total)
    call initialize_integration_points(n_points, n_nodes, xi, w)

    Ivec = 0.d0
    Ivec(1) = 1.d0
    Ivec(2) = 1.d0
    Ivec(3) = 1.d0

    element_residual = 0.d0
    element_stiffness = 0.d0
    G = 0.d0
	
    D = 0.d0
    mu1 = element_properties(1)
    K1 = element_properties(2)
    u = 0.d0
    u = reshape(dof_increment+dof_total,(/3,n_nodes/))

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        F = 0.d0
        F(1,1) = 1.d0
        F(2,2) = 1.d0
        F(3,3) = 1.d0
        do ii = 1,3
            do jj = 1,3

                !                if (ii == jj) then
                !                    F(ii,jj) = 1.d0
                !                end if
                do kk = 1,n_nodes
                    F(ii,jj) = F(ii,jj) + u(ii,kk)*dNdx(kk,jj)
                end do
            end do
        end do

        !write(6,*) F
        !J = determinate(F)
        B_def = matmul(F,transpose(F))
        call invert_small(F,Finv,J)
        call invert_small(B_def,Binv,Bdet)

        !       J = sqrt(Bdet)
        dNdy = 0.d0
        do ii = 1,n_nodes
            do jj = 1,3
                do kk = 1,3
                    dNdy(ii,jj) = dNdy(ii,jj) + dNdx(ii,kk)*Finv(kk,jj)
                end do
            end do
        end do

        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
!write(6,*) B
        !Find the G matrix from the left cauchy def. tensor
        G = 0.d0
        G(1,1) = 2*B_def(1,1)
        G(1,4) = 2*B_def(1,2)
        G(1,6) = 2*B_def(1,3)
        G(2,2) = 2*B_def(2,2)
        G(2,5) = 2*B_def(1,2)
        G(2,7) = 2*B_def(2,3)
        G(3,3) = 2*B_def(3,3)
        G(3,6) = 2*B_def(1,3)
        G(3,9) = 2*B_def(2,3)
        G(4,1) = 2*B_def(1,2)
        G(4,2) = 2*B_def(1,2)
        G(4,4) = 2*B_def(2,2)
        G(4,5) = 2*B_def(1,1)
        G(4,6) = 2*B_def(2,3)
        G(4,8) = 2*B_def(1,3)
        G(5,1) = 2*B_def(1,3)
        G(5,3) = 2*B_def(1,3)
        G(5,4) = 2*B_def(2,3)
        G(5,6) = 2*B_def(3,3)
        G(5,7) = 2*B_def(1,1)
        G(5,9) = 2*B_def(1,2)
        G(6,2) = 2*B_def(2,3)
        G(6,3) = 2*B_def(2,3)
        G(6,5) = 2*B_def(1,3)
        G(6,7) = 2*B_def(1,2)
        G(6,8) = 2*B_def(3,3)
        G(6,9) = 2*B_def(2,2)

        Bvec = 0.d0
        Bvec(1) = B_def(1,1)
        Bvec(2) = B_def(2,2)
        Bvec(3) = B_def(3,3)
        Bvec(4) = B_def(1,2)
        Bvec(5) = B_def(1,3)
        Bvec(6) = B_def(2,3)
!        write(6,*) Bvec
        Bvec_inv = 0.d0
        Bvec_inv(1) = Binv(1,1)
        Bvec_inv(2) = Binv(2,2)
        Bvec_inv(3) = Binv(3,3)
        Bvec_inv(4) = Binv(1,2)
        Bvec_inv(5) = Binv(1,3)
        Bvec_inv(6) = Binv(2,3)

        B_star = 0.d0
        B_star(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B_star(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B_star(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B_star(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B_star(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B_star(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B_star(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B_star(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B_star(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        Bkk = B_def(1,1)+B_def(2,2)+B_def(3,3)

        eye = 0.d0
        eye(1,1) = 1.d0
        eye(2,2) = 1.d0
        eye(3,3) = 1.d0

        eye6 = 0.d0
        eye6(1,1) = 1.d0
        eye6(2,2) = 1.d0
        eye6(3,3) = 1.d0
        eye6(4,4) = 0.5
        eye6(5,5) = 0.5
        eye6(6,6) = 0.5

        !Find Kirchoff stress
        do ii = 1,3
            do jj = 1,3
                tau(ii,jj)  = mu1*(B_def(ii,jj) - Bkk*eye(ii,jj)/3.d0)/(J**(2.d0/3.d0)) + K1*J*(J-1.d0)*eye(ii,jj)
            end do
        end do

        !Find D
        D = mu1/(J**(2.d0/3.d0))*eye6 + mu1/(3.d0*J**(2.d0/3.d0))*(Bkk/3.d0*spread(Ivec,dim=2,ncopies=6)* &
            spread(Bvec_inv,dim=1,ncopies=6) - spread(Ivec,dim=2,ncopies=6)*spread(Ivec,dim=1,ncopies=6) &
            - spread(Bvec,dim=2,ncopies=6)*spread(Bvec_inv,dim=1,ncopies=6)) +&
            K1*J*(J-1.d0/2.d0)*spread(Ivec,dim=2,ncopies=6)*spread(Bvec_inv,dim=1,ncopies=6)

        stress(1) = tau(1,1)
        stress(2) = tau(2,2)
        stress(3) = tau(3,3)
        stress(4) = tau(1,2)
        stress(5) = tau(1,3)
        stress(6) = tau(2,3)

        S = reshape(matmul(transpose(B),stress),(/3,length_dof_array/3/))
        !write(6,*) S
        do ii = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdy(ii:ii,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*ii-2:3*ii,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
            Svec = reshape(spread(S(1:3,ii:ii),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Smat(3*ii-2:3*ii,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
        end do
        Sigma = Pmat*transpose(Smat)

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant
!write(6,*) B


        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + (matmul(transpose(B),matmul(D,matmul(G,B_star)))-&
            Sigma)*w(kint)*determinant

    end do

    return
end subroutine el_linelast_hyperelastic

!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_linelast_hyperelastic(lmn, element_identifier, n_nodes, node_property_list, &         ! Input variables
    n_properties, element_properties,element_coords, length_coord_array,  &                     ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                              ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                       ! Input variables
    n_field_variables,field_variable_names, &                                                   ! Field variable definition
    nodal_fieldvariables)     ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Mesh, only : extract_node_data
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only : dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
!    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D

    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( inout )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables

  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,ii,jj,kk,k

    !    real (prec)  ::  strain(6), dstrain(6)          ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  symvec(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                        ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  mu1, K1                           ! Material properties
    real (prec)  ::  F(3,3),u(3,n_nodes),dNdy(n_nodes,3),G(6,9),tau(3,3),B_star(9,3*n_nodes)
    real (prec)  ::  B_def(3,3),Binv(3,3),Finv(3,3),S(3,length_dof_array/3)
    real (prec)  ::  Pmat(3*n_nodes,3*n_nodes),Smat(3*n_nodes,3*n_nodes),Sigma(3*n_nodes,3*n_nodes)
    real (prec)  ::  eye(3,3),eye6(6,6),Ivec(6),Bvec_inv(6),Bvec(6),Svec(3*n_nodes),Pvec(3*n_nodes)
    real (prec)  ::  J,Bdet,Bkk
    real (prec)  ::  p, smises                     ! Pressure and Mises stress

!        real (prec), intent(in)  ::  symvec(6)   ! symvec contains (A(1,1),A(2,2),A(3,3),A(1,2),A(1,3),A(2,3))
        real (prec) :: principalvals33(3)
        real (prec) :: B_princ(6)
        real (prec) :: p1,p2
        real (prec) :: pp,qq,rr
        real (prec) :: phi

    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Shear mod
    !     element_properties(2)         bulk mod

    nodal_fieldvariables = 0.d0

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    !    nn = 1
    !    call extract_node_data(nn,flag,n_coords,nodal_coords,n_dof,nodal_dof_increment,nodal_dof_total)
    call initialize_integration_points(n_points, n_nodes, xi, w)

     Ivec = 0.d0
    Ivec(1) = 1.d0
    Ivec(2) = 1.d0
    Ivec(3) = 1.d0

    G = 0.d0

    D = 0.d0
    mu1 = element_properties(1)
    K1 = element_properties(2)
    u = 0.d0
    u = reshape(dof_increment+dof_total,(/3,n_nodes/))

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        F = 0.d0
        F(1,1) = 1.d0
        F(2,2) = 1.d0
        F(3,3) = 1.d0
        do ii = 1,3
            do jj = 1,3

                !                if (ii == jj) then
                !                    F(ii,jj) = 1.d0
                !                end if
                do kk = 1,n_nodes
                    F(ii,jj) = F(ii,jj) + u(ii,kk)*dNdx(kk,jj)
                end do
            end do
        end do

        !write(6,*) F
        !J = determinate(F)
        B_def = matmul(F,transpose(F))
        call invert_small(F,Finv,J)
        call invert_small(B_def,Binv,Bdet)

        !       J = sqrt(Bdet)
        dNdy = 0.d0
        do ii = 1,n_nodes
            do jj = 1,3
                do kk = 1,3
                    dNdy(ii,jj) = dNdy(ii,jj) + dNdx(ii,kk)*Finv(kk,jj)
                end do
            end do
        end do

        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
!write(6,*) B
        !Find the G matrix from the left cauchy def. tensor
        G = 0.d0
        G(1,1) = 2*B_def(1,1)
        G(1,4) = 2*B_def(1,2)
        G(1,6) = 2*B_def(1,3)
        G(2,2) = 2*B_def(2,2)
        G(2,5) = 2*B_def(1,2)
        G(2,7) = 2*B_def(2,3)
        G(3,3) = 2*B_def(3,3)
        G(3,6) = 2*B_def(1,3)
        G(3,9) = 2*B_def(2,3)
        G(4,1) = 2*B_def(1,2)
        G(4,2) = 2*B_def(1,2)
        G(4,4) = 2*B_def(2,2)
        G(4,5) = 2*B_def(1,1)
        G(4,6) = 2*B_def(2,3)
        G(4,8) = 2*B_def(1,3)
        G(5,1) = 2*B_def(1,3)
        G(5,3) = 2*B_def(1,3)
        G(5,4) = 2*B_def(2,3)
        G(5,6) = 2*B_def(3,3)
        G(5,7) = 2*B_def(1,1)
        G(5,9) = 2*B_def(1,2)
        G(6,2) = 2*B_def(2,3)
        G(6,3) = 2*B_def(2,3)
        G(6,5) = 2*B_def(1,3)
        G(6,7) = 2*B_def(1,2)
        G(6,8) = 2*B_def(3,3)
        G(6,9) = 2*B_def(2,2)

        Bvec = 0.d0
        Bvec(1) = B_def(1,1)
        Bvec(2) = B_def(2,2)
        Bvec(3) = B_def(3,3)
        Bvec(4) = B_def(1,2)
        Bvec(5) = B_def(1,3)
        Bvec(6) = B_def(2,3)
!        write(6,*) Bvec
        Bvec_inv = 0.d0
        Bvec_inv(1) = Binv(1,1)
        Bvec_inv(2) = Binv(2,2)
        Bvec_inv(3) = Binv(3,3)
        Bvec_inv(4) = Binv(1,2)
        Bvec_inv(5) = Binv(1,3)
        Bvec_inv(6) = Binv(2,3)

        B_star = 0.d0
        B_star(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B_star(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B_star(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B_star(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B_star(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B_star(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B_star(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B_star(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B_star(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        Bkk = B_def(1,1)+B_def(2,2)+B_def(3,3)

        eye = 0.d0
        eye(1,1) = 1.d0
        eye(2,2) = 1.d0
        eye(3,3) = 1.d0

        eye6 = 0.d0
        eye6(1,1) = 1.d0
        eye6(2,2) = 1.d0
        eye6(3,3) = 1.d0
        eye6(4,4) = 0.5
        eye6(5,5) = 0.5
        eye6(6,6) = 0.5

        !Find Kirchoff stress
        do ii = 1,3
            do jj = 1,3
                tau(ii,jj)  = mu1*(B_def(ii,jj) - Bkk*eye(ii,jj)/3.d0)/(J**(2.d0/3.d0)) + K1*J*(J-1.d0)*eye(ii,jj)
            end do
        end do

        !Find D
        D = mu1/(J**(2.d0/3.d0))*eye6 + mu1/(3.d0*J**(2.d0/3.d0))*(Bkk/3.d0*spread(Ivec,dim=2,ncopies=6)* &
            spread(Bvec_inv,dim=1,ncopies=6) - spread(Ivec,dim=2,ncopies=6)*spread(Ivec,dim=1,ncopies=6) &
            - spread(Bvec,dim=2,ncopies=6)*spread(Bvec_inv,dim=1,ncopies=6)) +&
            K1*J*(J-1.d0/2.d0)*spread(Ivec,dim=2,ncopies=6)*spread(Bvec_inv,dim=1,ncopies=6)

        stress(1) = tau(1,1)
        stress(2) = tau(2,2)
        stress(3) = tau(3,3)
        stress(4) = tau(1,2)
        stress(5) = tau(1,3)
        stress(6) = tau(2,3)

        S = reshape(matmul(transpose(B),stress),(/3,length_dof_array/3/))
        !write(6,*) S
        do ii = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdy(ii:ii,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*ii-2:3*ii,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
            Svec = reshape(spread(S(1:3,ii:ii),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Smat(3*ii-2:3*ii,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
        end do
        Sigma = Pmat*transpose(Smat)

        !        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant
        !
        !        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
        !            + (matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,matmul(G,B(1:6,1:3*n_nodes))))-Sigma)*w(kint)*determinant

!        write(6,*) symvec

        symvec = 0.d0
        symvec(1:6) = stress(1:6)

        principalvals33 = 0.d0
        p1 = symvec(4)*symvec(4)+symvec(5)*symvec(5)+symvec(6)*symvec(6)
        if (p1 == 0.d0) then    ! Vector is diagonal - sort the eigenvalues
           principalvals33(1) = maxval(symvec(1:3))
           principalvals33(3) = minval(symvec(1:3))
           principalvals33(2) = sum(symvec(1:3)) - principalvals33(1) - principalvals33(2)
        else
          qq = sum(symvec(1:3))/3.d0
          p2 = (symvec(1) - qq)*(symvec(1)-qq) + (symvec(2) - qq)*(symvec(2) - qq) + (symvec(3) - qq)*(symvec(3) - qq) + 2.d0 * p1
          pp = dsqrt(p2 / 6.d0)
          B_princ = symvec/pp
          B_princ(1:3) = B_princ(1:3) - qq/pp
          rr = 0.5d0*(  B_princ(1)*B_princ(2)*B_princ(3)  &
               - B_princ(1)*B_princ(6)*B_princ(6)  &
               - B_princ(4)*B_princ(4)*B_princ(3)  &
               + B_princ(4)*B_princ(6)*B_princ(5)  &
               + B_princ(5)*B_princ(4)*B_princ(6)  &
               - B_princ(5)*B_princ(2)*B_princ(5) )

           if (rr < -1.d0) then
              phi = PI_D / 3.d0
           else if (rr > 1.d0) then
              phi = 0.d0
           else
              phi = dacos(rr) / 3.d0
           end if

           ! Principal values ordered from largest to smallest.
           principalvals33(1) = qq + 2.d0 * pp * cos(phi)
           principalvals33(3) = qq + 2.d0 * pp * cos(phi + (2.d0*PI_D/3.d0))
           principalvals33(2) = 3.d0 * qq - principalvals33(1) - principalvals33(3)
        endif
!write(6,*) principalvals33
        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)

        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                    stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                    stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                    stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                    stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                    stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                    stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                    smises*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'p1',2) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                    principalvals33(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'p2',2) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                    principalvals33(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'p3',2) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                    principalvals33(3)*N(1:n_nodes)*determinant*w(kint)

            endif
        end do
    end do
    return
end subroutine fieldvars_linelast_hyperelastic
