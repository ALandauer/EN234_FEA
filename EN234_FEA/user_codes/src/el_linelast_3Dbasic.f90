!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D

    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
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
    integer      :: n_points,kint,intvol

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6), dstress(6)             ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, D44, D11, D12, el_vol              ! Material properties
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    !    nn = 1
    !    call extract_node_data(nn,flag,n_coords,nodal_coords,n_dof,nodal_dof_increment,nodal_dof_total)
    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    !    E = element_properties(1)
    !    xnu = element_properties(2)
    !    d44 = 0.5D0*E/(1+xnu)
    !    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    !    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    !    D(1:3,1:3) = d12
    !    D(1,1) = d11
    !    D(2,2) = d11
    !    D(3,3) = d11
    !    D(4,4) = d44
    !    D(5,5) = d44
    !    D(6,6) = d44


    if (element_identifier == 1001) then

            !     --  Loop over integration points
        do kint = 1, n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
            B = 0.d0
            B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
            B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
            B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
            B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
            B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
            B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
            B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
            B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
            B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

            strain = matmul(B,dof_total)
            dstrain = matmul(B,dof_increment)


!            write(6,*) strain,dstrain

            if (size(element_properties) .eq. 4) then

                call calc_S_and_D(lmn, element_identifier, n_nodes, &  ! Input variables
                    n_properties,element_properties,length_coord_array, &         ! Input variables
                    dof_increment, dof_total, length_dof_array,D,stress,strain,dstrain)
            else
                !                strain = matmul(B,dof_total)
                !                dstrain = matmul(B,dof_increment)
                stress = matmul(D,strain+dstrain)
            end if

            element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

            element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
                + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

        end do
      
    elseif (element_identifier == 1002) then
        ! For B_bar elements
        el_vol = 0.d0
        dNbardx = 0.d0
        do intvol = 1,n_points
            call calculate_shapefunctions(xi(1:3,intvol),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            !            write(6,*) w(intvol)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
            el_vol = el_vol + w(intvol)*determinant
            dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3) + dNdx(1:n_nodes,1:3)*w(intvol)*determinant
        end do

        dNbardx = dNbardx/el_vol

        !         write(6,*) dNbardx

        !     --  Loop over integration points
        do kint = 1, n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

            B = 0.d0
            B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
            B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
            B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)

            B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
            B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)

            B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
            B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)

            B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
            B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

            B(1,1:3*n_nodes-2:3) = B(1,1:3*n_nodes-2:3) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3
            B(1,2:3*n_nodes-1:3) = B(1,2:3*n_nodes-1:3) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3
            B(1,3:3*n_nodes:3) = B(1,3:3*n_nodes:3) + (dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3

            B(2,1:3*n_nodes-2:3) = B(2,1:3*n_nodes-2:3) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3
            B(2,2:3*n_nodes-1:3) = B(2,2:3*n_nodes-1:3) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3
            B(2,3:3*n_nodes:3) = B(2,3:3*n_nodes:3) + (dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3

            B(3,1:3*n_nodes-2:3) = B(3,1:3*n_nodes-2:3) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3
            B(3,2:3*n_nodes-1:3) = B(3,2:3*n_nodes-1:3) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3
            B(3,3:3*n_nodes:3) = B(3,3:3*n_nodes:3) + (dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3

            strain = matmul(B,dof_total)
            dstrain = matmul(B,dof_increment)

            !            write(6,*) strain

            if (size(element_properties) .eq. 4) then

                call calc_S_and_D(lmn, element_identifier, n_nodes, &  ! Input variables
                    n_properties,element_properties,length_coord_array, &         ! Input variables
                    dof_increment, dof_total, length_dof_array,D,stress,dstress,strain,dstrain)
            else
                !                strain = matmul(B,dof_total)
                !                dstrain = matmul(B,dof_increment)
                stress = matmul(D,strain+dstrain)
            end if

            element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

            element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
                + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

        end do
    else
        write(6,*) 'no vaild element identifier'

    end if
  
    return
end subroutine el_linelast_3dbasic

!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_linelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none
!!!!!!!!!!!!!!!!!!!!
!
!
!
!  Not currently updated for b-bar and hypoelasticity
!
!
!
!
!!!!!!!!!!!!


    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
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

    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)

    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0

    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu)
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)

        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do

    return
end subroutine el_linelast_3dbasic_dynamic

!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Mesh, only : extract_node_data
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D

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
  
    integer      :: n_points,k,kint,intvol

    real (prec)  ::  strain(6), dstrain(6)          ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6), dstress(6)          ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                        ! Deviatoric stress
    real (prec)  ::  D(6,6)                         ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)          ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant        ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)      ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    !    real (prec)  :: E, xnu, D44, D11, D12           ! Material properties
    real (prec)  :: p, smises,el_vol                ! Pressure and Mises stress
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    !    nn = 1
    !    call extract_node_data(nn,flag,n_coords,nodal_coords,n_dof,nodal_dof_increment,nodal_dof_total)
    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0
	
    el_vol = 0.d0
    D = 0.d0
    !    E = element_properties(1)
    !    xnu = element_properties(2)
    !    d44 = 0.5D0*E/(1+xnu)
    !    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    !    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    !    D(1:3,1:3) = d12
    !    D(1,1) = d11
    !    D(2,2) = d11
    !    D(3,3) = d11
    !    D(4,4) = d44
    !    D(5,5) = d44
    !    D(6,6) = d44
  
    if (element_identifier == 1001) then

        !     --  Loop over integration points
        do kint = 1, n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
            B = 0.d0
            B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
            B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
            B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
            B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
            B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
            B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
            B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
            B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
            B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

            strain = matmul(B,dof_total)
            dstrain = matmul(B,dof_increment)

            !            write(6,*) strain,dstrain

            if (size(element_properties) .eq. 4) then

                call calc_S_and_D(lmn, element_identifier, n_nodes, &  ! Input variables
                    n_properties,element_properties,length_coord_array, &         ! Input variables
                    dof_increment, dof_total, length_dof_array,D,stress,dstress,strain,dstrain)

            else

                stress = matmul(D,strain+dstrain)

            end if

            !            write(6,*) strain,dstrain,stress

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
                !            else if (strcmp(field_variable_names(k),'e11',3) ) then
                !                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + (strain(1)+dstrain(1))
                !            else if (strcmp(field_variable_names(k),'e22',3) ) then
                !                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + (strain(2)+dstrain(2))
                !            else if (strcmp(field_variable_names(k),'e33',3) ) then
                !                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + (strain(3)+dstrain(3))
                !            else if (strcmp(field_variable_names(k),'e12',3) ) then
                !                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + (strain(4)+dstrain(4))/2
                !            else if (strcmp(field_variable_names(k),'e13',3) ) then
                !                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + (strain(5)+dstrain(5))/2
                !            else if (strcmp(field_variable_names(k),'e23',3) ) then
                !                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + (strain(6)+dstrain(6))/2
                endif
            end do
        end do

    elseif (element_identifier == 1002) then
        ! For B_bar elements
        el_vol = 0.d0
        dNbardx = 0.d0
        do intvol = 1,n_points

            call calculate_shapefunctions(xi(1:3,intvol),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
            el_vol = el_vol + w(intvol)*determinant
            dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3)&
                + dNdx(1:n_nodes,1:3)*w(intvol)*determinant

        end do

        dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3)/el_vol

        !        write(6,*) dNbardx

        !     --  Loop over integration points
        do kint = 1, n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

            !        write(6,*) dNbardx(1:n_nodes,1:3)

            B = 0.d0
            B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
            B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
            B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)

            B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
            B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)

            B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
            B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)

            B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
            B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

            B(1,1:3*n_nodes-2:3) = B(1,1:3*n_nodes-2:3) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3
            B(1,2:3*n_nodes-1:3) = B(1,2:3*n_nodes-1:3) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3
            B(1,3:3*n_nodes:3) = B(1,3:3*n_nodes:3) + (dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3

            B(2,1:3*n_nodes-2:3) = B(2,1:3*n_nodes-2:3) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3
            B(2,2:3*n_nodes-1:3) = B(2,2:3*n_nodes-1:3) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3
            B(2,3:3*n_nodes:3) = B(2,3:3*n_nodes:3) + (dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3

            B(3,1:3*n_nodes-2:3) = B(3,1:3*n_nodes-2:3) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3
            B(3,2:3*n_nodes-1:3) = B(3,2:3*n_nodes-1:3) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3
            B(3,3:3*n_nodes:3) = B(3,3:3*n_nodes:3) + (dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3


            stress = matmul(D,strain+dstrain)
            dstress = matmul(D,dstrain)

            if (size(element_properties) .eq. 4) then

                call calc_S_and_D(lmn, element_identifier, n_nodes, &  ! Input variables
                    n_properties,element_properties,length_coord_array, &         ! Input variables
                    dof_increment, dof_total, length_dof_array,D,stress,dstress,strain,dstrain)
            else
                !                strain = matmul(B,dof_total)
                !                dstrain = matmul(B,dof_increment)
                stress = matmul(D,strain+dstrain)
            end if

            p = sum(stress(1:3))/3.d0
            sdev = stress
            sdev(1:3) = sdev(1:3)-p
            smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)))*dsqrt(1.5d0)
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

                endif

            end do
        end do
    else
        write(6,*) 'no vaild element identifier'

    end if

    return
end subroutine fieldvars_linelast_3dbasic

subroutine calc_S_and_D(lmn, element_identifier, n_nodes, &  ! Input variables
    n_properties,element_properties,length_coord_array, &         ! Input variables
    dof_increment, dof_total, length_dof_array,D,stress,strain,dstrain)!,  & ! Input variables
    !    n_state_variables, initial_state_variables,updated_state_variables, &        ! Input variables
    !    n_field_variables,field_variable_names, &                                    ! Field variable definition
    !    nodal_fieldvariables)      ! Output variables

    use Types
    use ParamIO
    use Mesh, only : node

    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    !    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    !    integer, intent( in )         :: n_field_variables                                      ! # field variables

    !    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    !    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    !    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    !    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( in )    :: strain(6)             ! total strain
    real( prec ), intent( inout ) :: stress(6)             ! total strain
    real( prec ), intent( in )    :: dstrain(6)            ! step of strain
!    real( prec ), intent( in )    :: dstress(6)           ! step of stress
    !    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
    !   real( prec ), intent( inout ) :: D

    ! Local Variables
    real (prec), dimension(6,6)  ::  D1
    real (prec), dimension(6,6)  ::  D2
    !    real (prec)  ::  total_strain(6)
    real (prec)  ::  e_dyadic_e(6,6)
    real (prec)  ::  e(6)
    real (prec)  ::  eps(6)
    real (prec), intent(out) ::  D(6,6)
    real (prec)  ::  K,E_s,E_t,eps_e,sum_eps,stress_e,stress_0,n,dstress_e,&
        dstrain_e,strain_0
    !    integer :: ij


    !mat props
    D1 = 0.d0
    D1(1,1) = 2.d0
    D1(2,2) = 2.d0
    D1(3,3) = 2.d0
    D1(4,4) = 1.d0
    D1(5,5) = 1.d0
    D1(6,6) = 1.d0

    D2 = 0.d0
    D2(1:3,1:3) = 1.d0

    dstrain_e = 0.d0
    stress = 0.d0
    dstress_e = 0.d0
    stress_e = 0.d0

    stress_0 = element_properties(1)
    strain_0 = element_properties(2)
    n = element_properties(3)
    K = element_properties(4)

    eps = 0.d0
    eps_e = 0.d0

    !Calculate D matrix coeff variable values
    eps = strain+dstrain
    eps(4:6) = eps(4:6)/2.d0
    eps_e = sqrt(2.d0/3.d0*(eps(1)**2.d0+eps(2)**2.d0+eps(3)**2.d0+eps(4)**2.d0+&
        eps(5)**2.d0+eps(6)**2.d0))

    dstrain_e = sqrt(2.d0/3.d0*(dstrain(1)**2.d0+dstrain(2)**2.d0+dstrain(3)**2.d0+&
        dstrain(4)**2.d0+dstrain(5)**2.d0+dstrain(6)**2.d0))

!write(6,*) eps_e

    e = 0.d0
    e_dyadic_e = 0.d0
    sum_eps = 0.d0
    e = eps - 1.d0/3.d0*(eps(1)+eps(2)+eps(3)+eps(4)+eps(5)+eps(6))
    e_dyadic_e = spread(e,dim=2,ncopies=6)*spread(e,dim=1,ncopies=6)
    sum_eps = sum(eps**2.d0)

    !calculate elastic stress
    if (eps_e >= strain_0) then

        stress_e = (stress_0)*&
            (eps_e/strain_0)**1.d0/n

        dstress_e = (stress_0)*&
            (dstrain_e/strain_0)**1.d0/n

    else

        stress_e = (stress_0)*&
            (sqrt((1+n**2)/(n-1)**2-(n/(n-1)-eps_e/n)**2.d0)-1.d0/(n-1))

        dstress_e = (stress_0)*&
            (sqrt((1.d0+n**2.d0)/(n-1.d0)**2.d0-&
            (n/(n-1.d0)-dstrain_e/strain_0)**2.d0)-1.d0/(n-1.d0))

    end if

    !    write(6,*) eps_e,stress_e,dstress_e

    !stress vector
    stress(1) = 2.d0/3.d0*stress_e*e(1)/eps_e + K*(eps(1)+eps(2)+eps(3)+&
        eps(4)+eps(5)+eps(6))
    stress(2) = 2.d0/3.d0*stress_e*e(2)/eps_e + K*(eps(1)+eps(2)+eps(3)+&
        eps(4)+eps(5)+eps(6))
    stress(3) = 2.d0/3.d0*stress_e*e(3)/eps_e + K*(eps(1)+eps(2)+eps(3)+&
        eps(4)+eps(5)+eps(6))
    stress(4) = 2.d0/3.d0*stress_e*e(4)/eps_e + K*(eps(1)+eps(2)+eps(3)+&
        eps(4)+eps(5)+eps(6))
    stress(5) = 2.d0/3.d0*stress_e*e(5)/eps_e + K*(eps(1)+eps(2)+eps(3)+&
        eps(4)+eps(5)+eps(6))
    stress(6) = 2.d0/3.d0*stress_e*e(6)/eps_e + K*(eps(1)+eps(2)+eps(3)+&
        eps(4)+eps(5)+eps(6))


    !check for the zero strain case
    if (sum_eps == 0) then

        stress = 0.d0
        E_s = dstress_e/dstrain_e
        E_t = E_s

        !calculate D matrix
        D(1:6,1:6) = E_s/3.d0*D1 + (K-2.d0*E_s/9.d0)*D2

    else

        E_s = stress_e/eps_e
        E_t = dstress_e/dstrain_e

        !calculate D matrix
        D(1:6,1:6) = 4.d0/9.d0*(eps_e**2.d0)*(E_t - E_s)*e_dyadic_e + E_s/2*D1+&
            (K-2.d0*E_s/9.d0)*D2

    end if


write(6,*) D


end subroutine calc_S_and_D
