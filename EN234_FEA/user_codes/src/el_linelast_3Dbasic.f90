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
    real (prec)  ::  E, nu, D44, D11, D12, el_vol              ! Material properties
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
	
    stress = 0.d0
    D = 0.d0

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

            !write(6,*) dstrain

            if (size(element_properties) .eq. 4) then

                call calc_S_and_D(n_properties,element_properties,length_coord_array, &         ! Input variables
                    dof_increment, dof_total, length_dof_array,D,stress,strain,dstrain)

            else
                !                strain = matmul(B,dof_total)
                !                dstrain = matmul(B,dof_increment)
                E = element_properties(1)
                nu = element_properties(2)
                d44 = 0.5D0*E/(1+nu)
                d11 = (1.D0-nu)*E/( (1+nu)*(1-2.D0*nu) )
                d12 = nu*E/( (1+nu)*(1-2.D0*nu) )
                D(1:3,1:3) = d12
                D(1,1) = d11
                D(2,2) = d11
                D(3,3) = d11
                D(4,4) = d44
                D(5,5) = d44
                D(6,6) = d44
                stress = matmul(D,strain+dstrain)
            end if

            element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes)&
                - matmul(transpose(B),stress)*w(kint)*determinant

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

                call calc_S_and_D(n_properties,element_properties,length_coord_array, &         ! Input variables
                    dof_increment, dof_total, length_dof_array,D,stress,strain,dstrain)
            else
                !                strain = matmul(B,dof_total)
                !                dstrain = matmul(B,dof_increment)
                stress = matmul(D,strain+dstrain)
            end if

            element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) + matmul(transpose(B),stress)*w(kint)*determinant

            element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
                + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

        end do
    else
        write(6,*) 'no vaild element identifier'

    end if
  
    return
end subroutine el_linelast_3dbasic







!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_linelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &
    n_properties, element_properties,element_coords, length_coord_array, &
    dof_increment, dof_total, length_dof_array, &
    n_state_variables, initial_state_variables, &
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

    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)

    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint,i,j
    integer      :: intvol1,intvol2,m_count
    integer      :: state_cnt, elDelCnt

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  dof_i_reshape(3,n_nodes),dof_t_reshape(3,n_nodes)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  B_bar(6,length_dof_array)          ! strain = Bbar*(dof_total+dof_increment), B-bar method
    real (prec)  ::  temp_bar(6,length_dof_array)      ! temp matrix when constructing B-bar matrix
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node

    real (prec)  :: d_ij(3,3),midF(3,3),midF_inv(3,3),deltaF(3,3),delR_inv_temp(3,3),detF
    real (prec)  :: dNbardy(n_nodes,3),dNdy(n_nodes,3),dNdy_bar(n_nodes,3),del_L(3,3)
    real (prec)  :: deltaL(3,3),delR(3,3),delR_inv(3,3),del_eps(3,3),del_W(3,3),del_L_bar(3,3)
    real (prec)  :: eta,el_vol,del_eta,detR,el_vol_bar,fF,Vf

    !write(6,*)dof_increment
    x = reshape(element_coords,(/3,length_coord_array/3/))
    dof_i_reshape = reshape(dof_increment,(/3,length_dof_array/3/))
    dof_t_reshape = reshape(dof_total,(/3,length_dof_array/3/))


    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    midF = 0.d0
    midF_inv = 0.d0
    detF = 0.d0
    deltaF = 0.d0
    deltaL = 0.d0
    detR = 0.d0
    fF = 0.d0
    Vf = 0.d0

    d_ij = 0.d0
    d_ij(1,1) = 1.d0
    d_ij(2,2) = 1.d0
    d_ij(3,3) = 1.d0

    el_vol = 0.d0
    eta = 0.d0
    del_eta = 0.d0
    dNdy = 0.d0
    dNdy_bar = 0.d0
    dNbardy = 0.d0

    el_vol_bar = 0.d0
    !    dNbardx = 0.d0

    !    do intvol = 1,n_points
    !        call calculate_shapefunctions(xi(1:3,intvol),n_nodes,N,dNdxi)
    !        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
    !        call invert_small(dxdxi,dxidx,determinant)
    !        !            write(6,*) w(intvol)
    !        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
    !        el_vol_bar = el_vol_bar + w(intvol)*determinant
    !        dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3) + dNdx(1:n_nodes,1:3)*w(intvol)*determinant
    !    end do
    !write(6,*) 'dNbardx'
    !write(6,*) dNbardx
    !write(6,*) 'el_vol_bar'
    !write(6,*) el_vol_bar
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        midF = 0.d0;
        deltaF = 0.d0
        !        do i = 1,3
        !            do j = 1,3
        !                do a = 1,n_nodes
        deltaF = matmul(dof_i_reshape(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
        midF = matmul(dof_t_reshape(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
        midF = midF + d_ij + deltaF/2.d0

        !write(6,*) midF

        call invert_small(midF,midF_inv,detF)

        deltaL = matmul(deltaF(1:3,1:3),midF_inv(1:3,1:3))
!        write(6,*) dof_i_reshape
        !write(6,*) deltaF

        el_vol = el_vol + w(kint)*determinant
        eta = eta + detF*w(kint)*determinant
        del_eta = del_eta + detF*(deltaL(1,1)+deltaL(2,2)+deltaL(3,3))*w(kint)*determinant

        dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),midF_inv(1:3,1:3))

        do intvol1 = 1,n_nodes
            do intvol2 = 1,3
                dNdy_bar(intvol1,intvol2) = dNdy_bar(intvol1,intvol2) + &
                    detF*dNdy(intvol1,intvol2)*w(kint)*determinant
            end do
        end do

        el_vol = el_vol + w(kint)*determinant

    end do

    eta = eta/el_vol
    del_eta = del_eta/(el_vol*eta)
    dNdy_bar = dNdy_bar/(el_vol*eta)


!    stress = initial_state_variables(1:6)

    !0.d0 loop-local vars...

    state_cnt = 0
    elDelCnt = 0

    !     --  Loop over integration points a second time
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        !        do i = 1,3
        !            do j = 1,3
        !                do a = 1,n_nodes
        deltaF = matmul(dof_i_reshape(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
        midF = matmul(dof_t_reshape(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
        midF = midF + d_ij + deltaF/2.d0
        !write(6,*) midF

        call invert_small(midF,midF_inv,detF)
        !write(6,*) detF
        del_L = matmul(deltaF(1:3,1:3),midF_inv(1:3,1:3))
        del_L_bar = del_L+(del_eta - (del_L(1,1)+del_L(2,2)+del_L(3,3)))*d_ij/3.d0

        do i = 1,3
            do j = 1,3
                !                del_L_bar(i,j) = deltaF(i,j)*midF_inv(i,j) + (del_eta - deltaF(i,j)*midF_inv(i,j)*d_ij(i,j))/3.d0
                del_eps(i,j) = del_L_bar(i,j)/2.d0 + del_L_bar(j,i)/2.d0
                del_W(i,j) = del_L_bar(i,j)/2.d0 - del_L_bar(j,i)/2.d0

            end do
        end do



        delR_inv_temp = d_ij - del_W/2.d0
        !write(6,*) delR_inv

        call invert_small(delR_inv_temp,delR_inv,detR)

        delR(1:3,1:3) = matmul(delR_inv(1:3,1:3),(d_ij(1:3,1:3)+del_W(1:3,1:3)/2.d0))

        call gurson(element_properties,n_properties,n_state_variables,state_cnt,elDelCnt,initial_state_variables, &
            updated_state_variables,del_eps,delR,stress)

!write(6,*) stress
        !        dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),midF_inv(1:3,1:3))
        !        dNbardy(1:n_nodes,1:3) = matmul(dNbardx(1:n_nodes,1:3),midF_inv(1:3,1:3))
        !        B = 0.d0
        !        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        !        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        !        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        !
        !        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        !        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        !
        !        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        !        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        !
        !        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        !        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
        !
        !        B(1,1:3*n_nodes-2:3) = B(1,1:3*n_nodes-2:3) + (dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))/3
        !        B(1,2:3*n_nodes-1:3) = B(1,2:3*n_nodes-1:3) + (dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))/3
        !        B(1,3:3*n_nodes:3) = B(1,3:3*n_nodes:3) + (dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))/3
        !
        !        B(2,1:3*n_nodes-2:3) = B(2,1:3*n_nodes-2:3) + (dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))/3
        !        B(2,2:3*n_nodes-1:3) = B(2,2:3*n_nodes-1:3) + (dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))/3
        !        B(2,3:3*n_nodes:3) = B(2,3:3*n_nodes:3) + (dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))/3
        !
        !        B(3,1:3*n_nodes-2:3) = B(3,1:3*n_nodes-2:3) + (dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))/3
        !        B(3,2:3*n_nodes-1:3) = B(3,2:3*n_nodes-1:3) + (dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))/3
        !        B(3,3:3*n_nodes:3) = B(3,3:3*n_nodes:3) + (dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))/3


        dNdy(1:n_nodes,1:3)= matmul(dNdx(1:n_nodes,1:3),midF_inv)
        B=0.d0
        B_bar = 0.d0
        temp_bar = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
        do m_count = 1,n_nodes
            temp_bar(1:3:1,3*m_count-2) = dNbardy(m_count,1) - dNdy(m_count,1)
            temp_bar(1:3:1,3*m_count-1) = dNbardy(m_count,2) - dNdy(m_count,2)
            temp_bar(1:3:1,3*m_count) = dNbardy(m_count,3) - dNdy(m_count,3)
        end do
        B_bar = B + temp_bar/3.d0

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) -&
            matmul(transpose(B_bar),stress)*w(kint)*determinant

    end do

    !write(6,*) stress
    fF = element_properties(13)
    Vf = updated_state_variables(8)

    if (elDelCnt .eq. n_points) then
        element_deleted = .true.
    else
        element_deleted = .false.
    end if


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
  
    integer      :: n_points,k,kint,intvol,intvol1,intvol2,m_count,i,j,state_cnt,elDelCnt

    real (prec)  ::  strain(6), dstrain(6)          ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6), dstress(6)          ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                        ! Deviatoric stress
    real (prec)  ::  D(6,6)                         ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)          ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  B_bar(6,length_dof_array)          ! strain = Bbar*(dof_total+dof_increment), B-bar method
    real (prec)  ::  temp_bar(6,length_dof_array)      ! temp matrix when constructing B-bar matrix

    real (prec)  ::  dxidx(3,3), determinant        ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)      ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  dof_i_reshape(3,n_nodes),dof_t_reshape(3,n_nodes)
    !    real (prec)  :: E, nu, D44, D11, D12           ! Material properties
    real (prec)  :: p, smises,el_vol                ! Pressure and Mises stress

    real (prec)  :: d_ij(3,3),midF(3,3),midF_inv(3,3),deltaF(3,3),delR_inv_temp(3,3),detF
    real (prec)  :: dNbardy(n_nodes,3),dNdy(n_nodes,3),dNdy_bar(n_nodes,3),del_L(3,3)
    real (prec)  :: deltaL(3,3),delR(3,3),delR_inv(3,3),del_eps(3,3),del_W(3,3),del_L_bar(3,3)
    real (prec)  :: eta,del_eta,detR,el_vol_bar,fF,Vf
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))
    dof_i_reshape = reshape(dof_increment,(/3,length_dof_array/3/))
    dof_t_reshape = reshape(dof_total,(/3,length_dof_array/3/))

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
    !    nu = element_properties(2)
    !    d44 = 0.5D0*E/(1+nu)
    !    d11 = (1.D0-nu)*E/( (1+nu)*(1-2.D0*nu) )
    !    d12 = nu*E/( (1+nu)*(1-2.D0*nu) )
    !    D(1:3,1:3) = d12
    !    D(1,1) = d11
    !    D(2,2) = d11
    !    D(3,3) = d11
    !    D(4,4) = d44
    !    D(5,5) = d44
    !    D(6,6) = d44
  
    if (element_identifier == 1001) then

        !       write(6,*) ' Entered fieldvars'

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
            !write(6,*) dstrain

            !            write(6,*) strain,dstrain

            if (size(element_properties) .eq. 4) then
                !                         write(6,*) ' Calculating stress and D'
                call calc_S_and_D(n_properties,element_properties,length_coord_array, &         ! Input variables
                    dof_increment, dof_total, length_dof_array,D,stress,strain,dstrain)
            ! write(6,*) stress
            !               stress = 0.d0
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
        !            write(6,*) ' Field vars completed'
        end do

        !        write(6,*) 'Returning '
        return
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


            strain = matmul(B,dof_total)
            dstrain = matmul(B,dof_increment)

            if (size(element_properties) .eq. 4) then

                call calc_S_and_D(n_properties,element_properties,length_coord_array, &         ! Input variables
                    dof_increment, dof_total, length_dof_array,D,stress,strain,dstrain)
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

    elseif (element_identifier == 3001) then



        state_cnt = 0

        do kint = 1,n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

            stress = updated_state_variables(1+state_cnt:6+state_cnt)
            Vf = updated_state_variables(8+state_cnt)

            state_cnt = state_cnt+8
            !        element_residual = element_residual - matmul(transpose(B),stress)*w(kint)*determinant
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
                else if (strcmp(field_variable_names(k),'E11',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                        strain(1)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'E22',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                        strain(2)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'E33',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                        strain(3)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'E12',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                        strain(4)/0.5d0*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'E13',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                        strain(5)/0.5d0*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'E23',3) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                        strain(6)/0.5d0*N(1:n_nodes)*determinant*w(kint)
                !                else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                !                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                !                        smises*N(1:n_nodes)*determinant*w(kint)
                !                else if (strcmp(field_variable_names(k),'ematrix',6) ) then
                !                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                !                        updated_state_variables(7)*N(1:n_nodes)*determinant*w(kint)
                else if (strcmp(field_variable_names(k),'Vf',2) ) then
                    nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) +&
                        Vf*N(1:n_nodes)*determinant*w(kint)

                endif

            end do

        end do

    else
        write(6,*) 'no vaild element identifier'

    end if

    return
end subroutine fieldvars_linelast_3dbasic



subroutine calc_S_and_D(n_properties,element_properties,length_coord_array, &         ! Input variables
    dof_increment, dof_total, length_dof_array,D,stress,strain,dstrain)

    use Types
    use ParamIO
    use Mesh, only : node

    implicit none

    !    integer, intent( in )         :: lmn                                                    ! Element number
    !    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    !    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
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
    !        real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    !        real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( in )    :: strain(6)             ! total strain
    real( prec ), intent( inout ) :: stress(6)
    real( prec ), intent( in )    :: dstrain(6)            ! step of strain
    !        real( prec ), intent( in )    :: dstress(6)           ! step of stress
    !    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables

    ! Local Variables
    real (prec), dimension(6,6)  ::  D1(6,6)
    real (prec), dimension(6,6)  ::  D2(6,6)
    !    real (prec)  ::  total_strain(6)
    real (prec)  ::  e_dyadic_e(6,6),e_mat(6,6)
    real (prec)  ::  e(6),e_05(6)!,deps(6)!,deps_e,de_05(6),de(6)
    real (prec)  ::  eps(6),eps_05(6)!,deps_05(6)
    real (prec), intent(inout) ::  D(6,6)
    real (prec)  ::  K,E_s,E_t,eps_e,sum_strain,stress_e,stress_0,n,dstress_e,&
        dstrain_e,strain_0
    integer :: i,j

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

    !mat props


    stress_0 = element_properties(1)
    strain_0 = element_properties(2)
    n = element_properties(3)
    K = element_properties(4)

    eps = 0.d0
    eps_e = 0.d0
    !    deps = 0.d0
    !    deps_e = 0.d0
    !    de = 0.d0
    !    de_05 = 0.d0
    !    deps = 0.d0
    !    deps_e = 0.d0

    !Calculate D matrix coeff variable values
    eps(1:6) = strain(1:6)+dstrain(1:6)
    !    deps(1:6) = dstrain(1:6)
    eps_05(1:6) = eps(1:6)
    eps_05(4:6) = eps_05(4:6)/2.d0

    !write(6,*) eps_e

    e = 0.d0
    e_05 = 0.d0
    e_dyadic_e = 0.d0
    sum_strain = 0.d0

    e(1:3) = eps(1:3) - 1.d0/3.d0*(eps(1)+eps(2)+eps(3))
    e(4:6) = eps(4:6)
    e_05(1:3) = eps_05(1:3) - 1.d0/3.d0*(eps(1)+eps(2)+eps(3))
    e_05(4:6) = eps_05(4:6)
    e_dyadic_e(1:6,1:6) = &
        spread(e_05(1:6),dim=2,ncopies=6)*spread(e_05(1:6),dim=1,ncopies=6)
    e = e_05
    sum_strain = sqrt(sum(eps**2.d0))

    e_mat(1,1) = e(1)
    e_mat(2,2) = e(2)
    e_mat(3,3) = e(3)
    e_mat(1,2) = e(4)
    e_mat(2,1) = e_mat(1,2)
    e_mat(1,3) = e(5)
    e_mat(3,1) = e_mat(1,3)
    e_mat(2,3) = e(6)
    e_mat(3,2) = e_mat(2,3)


    do i = 1,3
        do j = 1,3
            eps_e = eps_e + e_mat(i,j)*e_mat(i,j)
        end do
    end do
    eps_e = sqrt(2.d0/3.d0*eps_e)

    !sqrt(2.d0/3.d0*(e(1)**2.d0+e(2)**2.d0+e(3)**2.d0))!+2.d0*e(4)**2.d0+&
    !        2.d0*e(5)**2.d0+2.d0*e(6)**2.d0))

    !write(6,*) deps_e

    !calculate elastic stress
    if (eps_e >= strain_0) then
        stress_e = stress_0*(eps_e/strain_0)**(1.d0/n)
        dstress_e = stress_0*((eps_e/strain_0)**(1.d0/n))/(n*eps_e)
    else
        stress_e = stress_0*&
            (sqrt((1.d0+n**2.d0)/((n-1.d0)**2.d0) -&
            (n/(n-1.d0)-eps_e/strain_0)**2.d0)-1.d0/(n-1.d0))

        dstress_e = stress_0*(n/(n-1)-eps_e/strain_0)/(strain_0*&
            sqrt((1.d0+n**2.d0)/((n-1.d0)**2.d0) -&
            (n/(n-1.d0)-eps_e/strain_0)**2.d0))

        !=(stress_0*(n/(n-1.d0)-eps_e/strain_0)/strain_0)/&
         !   sqrt((1.d0+n**2.d0)/(n-1.d0)**2.d0-(n/(n-1.d0)-eps_e/strain_0)**2.d0)

    end if
    !write(6,*) stress_e!,dstress_e


    !write(6,*) stress!,dstress_e

    !check for the 0.d0 strain case
    if (sum_strain == 0.d0) then

        stress(1:6) = 0.d0
        E_s = dstress_e

        !calculate D matrix
        D(1:6,1:6) = E_s/3.d0*D1 + (K-2.d0*E_s/9.d0)*D2
    !write(6,*) D

    else
        E_s = stress_e/eps_e
        E_t = dstress_e

        !stress vector
        stress(1) = 2.d0/3.d0*stress_e*e(1)/eps_e + K*(eps(1)+eps(2)+eps(3))
        stress(2) = 2.d0/3.d0*stress_e*e(2)/eps_e + K*(eps(1)+eps(2)+eps(3))
        stress(3) = 2.d0/3.d0*stress_e*e(3)/eps_e + K*(eps(1)+eps(2)+eps(3))
        stress(4) = 2.d0/3.d0*stress_e*e(4)/eps_e
        stress(5) = 2.d0/3.d0*stress_e*e(5)/eps_e
        stress(6) = 2.d0/3.d0*stress_e*e(6)/eps_e

        !calculate D matrix
        D(1:6,1:6) = 4.d0/(9.d0*eps_e*eps_e)*(E_t - E_s)*e_dyadic_e + E_s/3.d0*D1(1:6,1:6)+&
            (K-2.d0*E_s/9.d0)*D2(1:6,1:6)

    !write(6,*) D
    end if
!write(6,*) stress_e

end subroutine calc_S_and_D

subroutine gurson(element_properties,n_properties,n_state_variables,state_cnt,elDelCnt,initial_state_variables, &
    updated_state_variables,dstrain,dRot,stress1)
    use Types
    use ParamIO
    use Globals, only : TIME, DTIME

    use Element_Utilities, only : rotatesymvec

    implicit none

    integer, intent( in ) :: n_properties
    integer, intent( in ) :: n_state_variables
    real (prec), intent( in ) :: element_properties(n_properties)
    real (prec), intent( in ) :: initial_state_variables(n_state_variables)
    real (prec), intent( in ) :: dstrain(3,3)
    real (prec), intent( in ) :: dRot(3,3)

    real (prec), intent( out ) :: stress1(6)
    real (prec), intent( out ) :: updated_state_variables(n_state_variables)

    real (prec) :: Ss(3,3),tauDn(3,3),taun(3,3),de(3,3),stressR33(3,3),taun1(3,3)
    real (prec) :: stress0(6),stressR(6),d_ij(3,3)

    real (prec) :: E,nu,Y,e_dot,m,q1,q2,q3,fN,eN,sN,fc,fF,Vf,ematrix,fF_bar,phi_2
    real (prec) :: p_star,f_star,p,Se,Se_star,deps_barmat,Vf1,eN1matrix
    real (prec) :: dphidSe,dphidp,dev,dee,phi,temp1_F
    integer i,j,state_cnt,elDelCnt
!write(6,*) n_state_variables
    E = element_properties(1)
    nu = element_properties(2)
    Y = element_properties(3)
    e_dot = element_properties(4)
    m = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fN = element_properties(9)
    eN = element_properties(11)
    sN = element_properties(10)
    fc = element_properties(12)
    fF = element_properties(13)

    stress0 = initial_state_variables(1+state_cnt:6+state_cnt) ! Stress at start of increment
    ematrix = initial_state_variables(7+state_cnt)
    Vf = initial_state_variables(8+state_cnt)
!   write(6,*) stress0
    taun = 0.d0
    stressR = 0.d0

    d_ij = 0.d0
    d_ij(1,1) = 1.d0
    d_ij(2,2) = 1.d0
    d_ij(3,3) = 1.d0

    taun(1,1) = stress0(1)
    taun(2,2) = stress0(2)
    taun(3,3) = stress0(3)
    taun(1,2) = stress0(4)
    taun(1,3) = stress0(5)
    taun(2,3) = stress0(6)
    taun(3,2) = stress0(6)
    taun(2,1) = stress0(4)
    taun(3,1) = stress0(5)

    stressR = rotatesymvec(stress0,dRot)
    stressR33(1,1) = stressR(1)
    stressR33(2,2) = stressR(2)
    stressR33(3,3) = stressR(3)
    stressR33(1,2) = stressR(4)
    stressR33(1,3) = stressR(5)
    stressR33(2,3) = stressR(6)
    stressR33(3,2) = stressR(6)
    stressR33(2,1) = stressR(4)
    stressR33(3,1) = stressR(5)

    f_star = 0.d0
    fF_bar = q1 + sqrt(q1**2.d0 - q3)/q3

    p_star = 1.d0/3.d0*(taun(1,1)+taun(2,2)+taun(3,3))+&
        E/(3.d0*(1.d0-2.d0*nu))*(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))

    if (Vf .lt. fc) then
        f_star = Vf
    else if (Vf .gt. fF) then
        f_star = fF_bar
    else
        f_star = fc+(fF_bar-fc)/(fF-fc)*(Vf-fc)
    end if

    dee = 1.d-15
    dev = 1.d-15
    Se_star = 0.d0
    de = 0.d0
    Ss = 0.d0
    tauDn = 0.d0
    Se_star = 0.d0
    phi_2 = 0.d0

    tauDn = taun - (taun(1,1)+taun(2,2)+taun(3,3))*d_ij/3.d0
    de = dstrain - (dstrain(1,1)+dstrain(2,2)+dstrain(3,3))*d_ij/3.d0

    Ss = E*de/(1.d0+nu)+ matmul(dRot,matmul(tauDn,transpose(dRot)))!+tauDn(i,j)

    do i = 1,3
        do j = 1,3
            Se_star = Se_star + 1.5d0*Ss(i,j)*Ss(i,j)
        end do
    end do

    Se_star = sqrt(Se_star)
!write(6,*) Se_star

    phi_2 = (Se_star**2.d0)/(Y**2.d0) + 2.d0*q1*f_star*cosh(3.d0/2.d0*q2*p_star/Y) -&
        (1.d0+q3*f_star**2.d0)
!write(6,*) phi_2

!    phi = sqrt(phi_2)

    if (phi_2 .le. 0.0000000001) then
!        write(6,*) 'Elastic Step'
        do i = 1,3
            do j = 1,3
                if (i==j) then
                    taun1(i,j) = Ss(i,j) + p_star
                else
                    taun1(i,j) = Ss(i,j)
                end if

            end do
        end do

        deps_barmat = 0
        Vf1 = Vf
        eN1matrix = ematrix

    else
!        write(6,*) 'Plastic Step'
        call NewtonLoopF(E,nu,Y,e_dot,m,q1,q2,q3,fN,eN,sN,fc,fF,f_star,fF_bar,Se_star,&
            p_star,dee,dev,Se,p,phi,dphidSe,dphidp,temp1_F)

        do i = 1,3
            do j = 1,3
                if (i==j) then
                    taun1(i,j) = Ss(i,j) - dee*E*3.d0*Ss(i,j)/((1.d0+nu)*2.d0*Se_star)+&
                        p_star-E*dev/(3.d0*(1.d0-2.d0*nu))
                else
                    taun1(i,j) = Ss(i,j) - dee*E*3.d0*Ss(i,j)/((1.d0+nu)*2.d0*Se_star)
                end if

            end do
        end do

        deps_barmat = e_dot*DTIME/(1.d0-Vf)*phi**m*temp1_F**(-0.5d0)*(dphidSe*Se+1.d0/3.d0*dphidp*p)
        Vf1 = 1.d0 + (Vf-1.d0)*exp(-dev) + fN*deps_barmat/(sN*sqrt(2.d0*pi))*exp(-1.d0/2.d0*((ematrix-eN)/sN)**2.d0)
        eN1matrix = ematrix + deps_barmat

!        Se = Se_star+3.d0/2.d0*E/(1+nu)*dee
!        p = p_star+1.d0/2.d0*E/(1-2.d0*nu)*dev

    end if

!    dphidSe = (2.d0*Se_star/(Y**2.d0)-3.d0*E*dee/((Y**2.d0)*(1.d0+nu)))/(sqrt((Se**2.d0)/(Y**2.d0) +&
!        2.d0*q1*f_star*cosh(3.d0/2.d0*q2/Y*(p-E*dev/(3.d0*(1.d0-2.d0*nu))))))
!    dphidp = 3.d0*q1*q2*f_star/Y*sinh(3*q2/(2.d0*Y)*(p-E*dev/(3.d0*(1.d0-2.d0*nu))))/(sqrt((Se**2.d0)/(Y**2.d0) +&
!        2.d0*q1*f_star*cosh(3.d0/2.d0*q2/Y*(p-E*dev/(3.d0*(1.d0-2.d0*nu))))))

    if (Vf1 .gt. fF) then
        elDelCnt = elDelCnt + 1;
    end if

    stress1 = 0.d0
    stress1(1) = taun1(1,1)
    stress1(2) = taun1(2,2)
    stress1(3) = taun1(3,3)
    stress1(4) = taun1(1,2)
    stress1(5) = taun1(1,3)
    stress1(6) = taun1(2,3)

    updated_state_variables(1+state_cnt:6+state_cnt) = stress1(1:6)
    updated_state_variables(7+state_cnt) = eN1matrix
    updated_state_variables(8+state_cnt) = Vf1

    state_cnt = state_cnt + 8


    return
end subroutine gurson

subroutine NewtonLoopF(E,nu,Y,e_dot,m,q1,q2,q3,fN,eN,sN,fc,fF,f_star,fF_bar,&
    Se_star,p_star,dee,dev,Se,p,phi,dphidSe,dphidp,temp1_F)
    use Types
    use ParamIO
    use Globals, only : DTIME

    use Element_Utilities, only : invert_small

    implicit none
    real (prec), intent( in ) :: E,nu,Y,e_dot,m,q1,q2,q3,fN,eN,sN,fc,fF,f_star,Se_star,p_star

    logical convrg
    real (prec) :: max_iter,tol,relax,p
    real (prec) :: phi,dphidSe,dphidp,ddphidSedee,ddphidSedev,ddphidpdee
    real (prec) :: ddphidpdev,dphimdee,dphimdev,dF1ddee,dF1ddev,dF2ddee,dF2ddev,F1,F2
    real (prec) :: sol(2),de_col(2),F_col(2),J_invF(2),Jac(2,2),Jac_inv(2,2)
    real (prec) :: ematrix0,ematrix,Vf0,Vf,d_ij(3,3)
    real (prec) :: fF_bar,dEmatrix,dee,dev,error,Se
    real (prec) :: dsigdEe,dpdEv,ddphidpdSe,ddphiddSe
    real (prec) :: ddphiddp,temp1_F,temp2_F,temp3_F
    real (prec) :: dF_mat(2,2),dE_vec(2),sol_new(2),dF_mat_inv(2,2),det_dF_mat,phi_2
    integer iter

    dee = 1.d-15
    dev = 1.d-15
    error = Y
    tol = 1.d-8 * Y
    sol(1) =  dee
    sol(2) =  dev
    relax = 0.3d0
    iter = 0
    max_iter = 100

    do while(error .GT. tol)
        iter = iter + 1

        if (iter .EQ. max_iter) then
            stop
        end if

        Se = Se_star -1.5d0*E/(1.d0+nu)*sol(1)
        p = P_star - E*1.d0/(3.d0*(1.d0-2.d0*nu))*sol(2)

        phi = sqrt(Se**2.d0/(Y**2.d0) + 2.d0*q1*f_star*cosh(1.5d0*q2*p/Y)-(1.d0 + q3*f_star**2.d0))

        dphidSe = Se/(Y**2.d0*phi)
        dphidp = 1.5d0*q1*q2*f_star*sinh(1.5d0*q2*p/Y)/(Y*phi)
        dsigdEe = -1.5d0*E/(1.d0+nu)
        dpdEv = -1.d0/3.d0*E/(1.d0-2.d0*nu)
        ddphidpdSe = -Se/(Y**2.d0*phi**2.d0)*dphidp;
        ddphiddSe = 1.d0/(Y**2.d0*phi) - Se/(Y*phi)**2.d0*dphidSe
        ddphiddp = 1.5d0*q1*q2*f_star/Y*(cosh(1.5d0*q2*p/Y)*1.5d0*q2/(Y*phi)&
            - sinh(1.5d0*q2*p/Y)/phi**2.d0*dphidp)

        temp1_F = dphidSe**2.d0 +2.d0/9.d0*dphidp**2.d0
        temp2_F = 2.d0*dphidSe*ddphiddSe +4.d0/9.d0*dphidp*ddphidpdSe
        temp3_F = 2.d0*dphidSe*ddphidpdSe + 4.d0/9.d0*dphidp*ddphiddp
        dF1ddee = 0.5d0*temp1_F**(-0.5d0)*temp2_F*dsigdEe*sol(1)/(dtime*e_dot) + temp1_F**0.5d0/(dtime*e_dot)&
            -ddphiddSe*dsigdEe*phi**m - dphidSe*m*phi**(m-1.d0)*dphidSe*dsigdEe
        dF1ddev = 0.5*temp1_F**(-0.5d0)*temp3_F*dpdEv*sol(1)/(dtime*e_dot) -ddphidpdSe*dpdEv*phi**m -dphidSe&
            *m*phi**(m-1.d0)*dphidp*dpdEv
        dF2ddee = 0.5d0*temp1_F**(-0.5d0)*temp2_F*dsigdEe*sol(2)/(dtime*e_dot) - ddphidpdSe*dsigdEe*phi**m&
            -dphidp*m*phi**(m-1)*dphidSe*dsigdEe
        dF2ddev = 0.5*temp1_F**(-0.5d0)*temp3_F*dpdEv*sol(2)/(dtime*e_dot) + temp1_F**0.5d0/(dtime*e_dot)&
            -ddphiddp*dpdEv*phi**m-dphidp*m*phi**(m-1)*dphidp*dpdEv


        F1 = temp1_F**0.5d0*sol(1)/(dtime*e_dot) - dphidSe*phi**m
        F2 = temp1_F**0.5d0*sol(2)/(dtime*e_dot) - dphidp*phi**m

        F_col(1) = F1
        F_col(2) = F2

        dF_mat(1,1) = dF1ddee
        dF_mat(1,2) = dF1ddev
        dF_mat(2,1) = dF2ddee
        dF_mat(2,2) = dF2ddev

        call invert_small(dF_mat,dF_mat_inv,det_dF_mat)
        dE_vec = - matmul(dF_mat_inv,F_col)
        sol = sol + relax*dE_vec
        error = dsqrt(F1**2.d0 +F2**2.d0)

    end do

    dEe = sol(1)
    dEv = sol(2)

    Se = Se_star - 1.5d0*E/(1.d0+nu)*dee
    p = P_star - E*1.d0/(3.d0*(1.d0-2.d0*nu))*dev
    phi = dsqrt(Se**2.d0/(Y**2.d0) + 2.d0*q1*f_star*cosh(1.5d0*q2*p/Y)-(1.d0 + q3*f_star**2.d0))
    dphidSe = Se/(Y**2.d0*phi)
    dphidp  = q1*f_star*sinh(1.5d0*q2*p/Y)*1.5d0*q2/(Y*phi)
    temp1_F = dphidSe**2.d0 +2.d0/9.d0*dphidp**2.d0


return
end subroutine NewtonLoopF



