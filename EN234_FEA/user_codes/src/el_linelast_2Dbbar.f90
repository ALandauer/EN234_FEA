!     Subroutines for basic 2D linear elastic elements



!==========================SUBROUTINE el_linelast_2dbasic ==============================
subroutine el_linelast_2dbbar(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
dof_increment, dof_total, length_dof_array, &                                                ! Input variables
n_state_variables, initial_state_variables, &                                                ! Input variables
updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
use Types
use ParamIO
!  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
use Mesh, only : node
use Mesh, only : extract_node_data
use Element_Utilities, only : N => shape_functions_2D
use Element_Utilities, only : dNdxi => shape_function_derivatives_2D
use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
use Element_Utilities, only : dxdxi => jacobian_2D
use Element_Utilities, only : initialize_integration_points
use Element_Utilities, only : calculate_shapefunctions
use Element_Utilities, only : invert_small
use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_2D

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

real( prec )   :: nodal_coords(length_coord_array)                                      !use to call to extract_node_data
real( prec )   :: nodal_dof_increment(length_coord_array)
real( prec )   :: nodal_dof_total(length_coord_array)

logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)


! Local Variables
integer      :: n_points,kint,nn,flag,n_dof,n_coords,intvol

real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, 2e12]
real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
real (prec)  ::  B(3,length_dof_array)             ! strain = B*(dof_total+dof_increment)
real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
real (prec)  ::  E, xnu, D33, D11, D12, el_vol              ! Material properties

!    nn = 1
!    call extract_node_data(nn,flag,n_coords,nodal_coords,n_dof,nodal_dof_increment,nodal_dof_total)
!
!     Subroutine to compute element stiffness matrix and residual force vector for 2D linear elastic elements
!     El props are:

!     element_properties(1)         Young's modulus
!     element_properties(2)         Poisson's ratio

fail = .false.

x = reshape(element_coords,(/2,length_coord_array/2/))

if (n_nodes == 3) n_points = 1 !L6, P.11
if (n_nodes == 6) n_points = 4
if (n_nodes == 4) n_points = 4
if (n_nodes == 8) n_points = 9

call initialize_integration_points(n_points, n_nodes, xi, w)

element_residual = 0.d0
element_stiffness = 0.d0

D = 0.d0
E = element_properties(1)
xnu = element_properties(2)
d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
d33 = (1.D0-2*xnu)*E/( 2*(1+xnu)*(1-2.D0*xnu) )
D(1:2,1:2) = d12
D(1,1) = d11
D(2,2) = d11
D(3,3) = d33

!    write(6,*) flag


if (flag == 1) then
!     --  Loop over integration points
do kint = 1, n_points
call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
call invert_small(dxdxi,dxidx,determinant)
dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
B = 0.d0

B(1,1:2*n_nodes-2:2) = dNdx(1:n_nodes,1)
B(2,2:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
B(3,1:2*n_nodes-2:2) = dNdx(1:n_nodes,2)
B(3,2:2*n_nodes-1:2) = dNdx(1:n_nodes,1)

strain = matmul(B,dof_total)
dstrain = matmul(B,dof_increment)

stress = matmul(D,strain+dstrain)
element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

element_stiffness(1:2*n_nodes,1:2*n_nodes) = element_stiffness(1:2*n_nodes,1:2*n_nodes) &
+ matmul(transpose(B(1:3,1:2*n_nodes)),matmul(D,B(1:3,1:2*n_nodes)))*w(kint)*determinant

end do
else if (element_identifier == 1011) then

! For B_bar elements
el_vol = 0
dNbardx=0
do intvol = 1,n_points
el_vol = el_vol + w(intvol)*determinant
dNbardx(1:n_nodes,1:2) = dNbardx(1:n_nodes,1:2) + dNdx(1:n_nodes,1:2)*w(intvol)*determinant
end do

dNbardx = dNbardx/el_vol
!     --  Loop over integration points
do kint = 1, n_points
call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
call invert_small(dxdxi,dxidx,determinant)
dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
B = 0.d0

B(1,1:2*n_nodes-2:2) = dNdx(1:n_nodes,1)
B(2,2:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
B(3,1:2*n_nodes-2:2) = dNdx(1:n_nodes,2)
B(3,2:2*n_nodes-1:2) = dNdx(1:n_nodes,1)

!set up B_bar elements
B(1,1:2*n_nodes-2:2) = B(1,1:2*n_nodes-2:2) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/2
B(1,2:2*n_nodes-1:2) = B(1,2:2*n_nodes-1:2) + (dNbardx(2:n_nodes,2)-dNdx(2:n_nodes,2))/2

B(2,1:2*n_nodes-2:2) = B(2,1:2*n_nodes-2:2) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/2
B(2,2:2*n_nodes-1:2) = B(2,2:2*n_nodes-1:2) + (dNbardx(2:n_nodes,2)-dNdx(2:n_nodes,2))/2


strain = matmul(B,dof_total)
dstrain = matmul(B,dof_increment)

stress = matmul(D,strain+dstrain)
element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - &
 matmul(transpose(B),stress)*w(kint)*determinant

element_stiffness(1:2*n_nodes,1:2*n_nodes) = element_stiffness(1:2*n_nodes,1:2*n_nodes) &
+ matmul(transpose(B(1:3,1:2*n_nodes)),matmul(D,B(1:3,1:2*n_nodes)))*w(kint)*determinant
end do

else
write(6,*) 'no vaild element identifier'

end if

return
end subroutine el_linelast_2dbbar


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_linelast_2dbbar_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
n_state_variables, initial_state_variables, &                                                        ! Input variables
updated_state_variables,element_residual,element_deleted)                                            ! Output variables
use Types
use ParamIO
use Mesh, only : node
use Mesh, only : extract_node_data
use Element_Utilities, only : N => shape_functions_2D
use Element_Utilities, only:  dNdxi => shape_function_derivatives_2D
use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
use Element_Utilities, only : dxdxi => jacobian_2D
use Element_Utilities, only : initialize_integration_points
use Element_Utilities, only : calculate_shapefunctions
use Element_Utilities, only : invert_small
use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_2D

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
real( prec )   :: nodal_coords(length_coord_array)                                      !use to call to extract_node_data
real( prec )   :: nodal_dof_increment(length_coord_array)
real( prec )   :: nodal_dof_total(length_coord_array)

real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine

real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)

logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element



! Local Variables
integer      :: n_points,kint,nn,flag,n_dof,n_coords,intvol

real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, 2e12]
real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
real (prec)  ::  B(3,length_dof_array)             ! strain = B*(dof_total+dof_increment)
real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
real (prec)  ::  E, xnu, D33, D11, D12, el_vol              ! Material properties

!nn = 1
!call extract_node_data(nn,flag,n_coords,nodal_coords,n_dof,nodal_dof_increment,nodal_dof_total)
!
!     Subroutine to compute element force vector for a linear elastodynamic problem
!     El props are:

!     element_properties(1)         Young's modulus
!     element_properties(2)         Poisson's ratio

x = reshape(element_coords,(/2,length_coord_array/2/))

if (n_nodes == 3) n_points = 1 !L6, P.11
if (n_nodes == 6) n_points = 4
if (n_nodes == 4) n_points = 4
if (n_nodes == 8) n_points = 9

call initialize_integration_points(n_points, n_nodes, xi, w)

element_residual = 0.d0

D = 0.d0
E = element_properties(1)
xnu = element_properties(2)
d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
d33 = (1.D0-2*xnu)*E/( 2*(1+xnu)*(1-2.D0*xnu) )
D(1:2,1:2) = d12
D(1,1) = d11
D(2,2) = d11
D(3,3) = d33
!    D(4,4) = d44
!    D(5,5) = d44
!    D(6,6) = d44

!     --  Loop over integration points
if (flag == 1) then

do kint = 1, n_points
call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
call invert_small(dxdxi,dxidx,determinant)
dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
B = 0.d0
B(1,1:2*n_nodes-2:2) = dNdx(1:n_nodes,1)
B(2,2:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
B(3,1:2*n_nodes-2:2) = dNdx(1:n_nodes,2)
B(3,2:2*n_nodes-1:2) = dNdx(1:n_nodes,1)


strain = matmul(B,dof_total)
dstrain = matmul(B,dof_increment)

stress = matmul(D,strain+dstrain)
element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

end do

elseif (flag == 2) then
! For B_bar elements

!     --  Loop over integration points
do kint = 1, n_points
call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
call invert_small(dxdxi,dxidx,determinant)
dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
B = 0.d0

do intvol = 1,n_points
el_vol = el_vol + w(intvol)*determinant
dNbardx(1:n_nodes,1:2) = dNbardx(1:n_nodes,1:2) + dNdx(1:n_nodes,1:2)*w(intvol)*determinant
end do

dNbardx = dNbardx/el_vol

B = 0.d0
B(1,1:2*n_nodes-2:2) = dNdx(1:n_nodes,1)
B(2,2:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
B(3,1:2*n_nodes-2:2) = dNdx(1:n_nodes,2)
B(3,2:2*n_nodes-1:2) = dNdx(1:n_nodes,1)

!set up B_bar elements
B(1,1:2*n_nodes-1:2) = B(1,1:2*n_nodes-1:2) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/2
B(1,2:2*n_nodes-1:2) = B(1,2:2*n_nodes-1:2) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/2

B(2,1:2*n_nodes-1:2) = B(2,1:2*n_nodes-1:2) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/2
B(2,2:2*n_nodes-1:2) = B(2,2:2*n_nodes-1:2) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/2

B(3,1:2*n_nodes-1:2) = B(3,1:2*n_nodes-1:2) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/2
B(3,2:2*n_nodes-1:2) = B(3,2:2*n_nodes-1:2) + (dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/2


strain = matmul(B,dof_total)
dstrain = matmul(B,dof_increment)

stress = matmul(D,strain+dstrain)
element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant
end do
else
write(6,*) 'no vaild element identifier'

end if
return
end subroutine el_linelast_2dbbar_dynamic


!==========================SUBROUTINE fieldvars_linelast_2dbasic ==============================
subroutine fieldvars_linelast_2dbbar(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
n_field_variables,field_variable_names, &                                                           ! Field variable definition
nodal_fieldvariables)      ! Output variables
use Types
use ParamIO
use Mesh, only : node
use Mesh, only : extract_node_data
use Element_Utilities, only : N => shape_functions_2D
use Element_Utilities, only: dNdxi => shape_function_derivatives_2D
use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_2D
use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
use Element_Utilities, only : dxdxi => jacobian_2D
use Element_Utilities, only : initialize_integration_points
use Element_Utilities, only : calculate_shapefunctions
use Element_Utilities, only : invert_small
use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_2D

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

real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
real( prec )   :: nodal_coords(length_coord_array)                                      !use to call to extract_node_data
real( prec )   :: nodal_dof_increment(length_coord_array)
real( prec )   :: nodal_dof_total(length_coord_array)

! Local Variables
logical      :: strcmp
integer      :: n_points,kint,nn,flag,n_dof,n_coords,k,intvol

real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, 2e12]
real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
real (prec)  ::  sdev(3)                           ! Deviatoric stress
real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
real (prec)  ::  B(3,length_dof_array)             ! strain = B*(dof_total+dof_increment)
real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
real (prec)  :: E, xnu, D33, D11, D12              ! Material properties
real (prec)  :: p, smises, el_vol                          ! Pressure and Mises stress
!
!     Subroutine to compute element contribution to project element integration point data to nodes

!     element_properties(1)         Young's modulus
!     element_properties(2)         Poisson's ratio


x = reshape(element_coords,(/2,length_coord_array/2/))

if (n_nodes == 3) n_points = 1 !L6, P.11
if (n_nodes == 6) n_points = 4
if (n_nodes == 4) n_points = 4
if (n_nodes == 8) n_points = 9

!    nn = 1
!    call extract_node_data(nn,flag,n_coords,nodal_coords,n_dof,nodal_dof_increment,nodal_dof_total)
call initialize_integration_points(n_points, n_nodes, xi, w)

nodal_fieldvariables = 0.d0

D = 0.d0
E = element_properties(1)
xnu = element_properties(2)
d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
d33 = (1.D0-2*xnu)*E/( 2*(1+xnu)*(1-2.D0*xnu) )
D(1:2,1:2) = d12
D(1,1) = d11
D(2,2) = d11
D(3,3) = d33


if (flag == 1) then

!     --  Loop over integration points
do kint = 1, n_points
call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))

call invert_small(dxdxi,dxidx,determinant)
dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
B = 0.d0
B(1,1:2*n_nodes-2:2) = dNdx(1:n_nodes,1)
B(2,2:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
B(3,1:2*n_nodes-2:2) = dNdx(1:n_nodes,2)
B(3,2:2*n_nodes-1:2) = dNdx(1:n_nodes,1)


strain = matmul(B,dof_total)
dstrain = matmul(B,dof_increment)
stress = matmul(D,strain+dstrain)
p = sum(stress(1:3))/3.d0
sdev = stress
sdev(1:3) = sdev(1:3)-p
smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)))! + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
do k = 1,n_field_variables
if (strcmp(field_variable_names(k),'S11',3) ) then
nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
else if (strcmp(field_variable_names(k),'S22',3) ) then
nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
!            else if (strcmp(field_variable_names(k),'S33',3) ) then
!                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
else if (strcmp(field_variable_names(k),'S12',3) ) then
nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
!            else if (strcmp(field_variable_names(k),'S13',3) ) then
!                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
!            else if (strcmp(field_variable_names(k),'S23',3) ) then
!                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
!            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
!                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
endif
end do

end do

elseif (element_identifier == 1011) then
! For B_bar elements
el_vol = 0
dNbardx=0
do intvol = 1,n_points
el_vol = el_vol + w(intvol)*determinant
dNbardx(1:n_nodes,1:2) = dNbardx(1:n_nodes,1:2) + dNdx(1:n_nodes,1:2)*w(intvol)*determinant
end do

dNbardx = dNbardx/el_vol
!     --  Loop over integration points
do kint = 1, n_points

call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
call invert_small(dxdxi,dxidx,determinant)
dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)

B = 0.d0
B(1,1:2*n_nodes-2:2) = dNdx(1:n_nodes,1)
B(2,2:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
B(3,1:2*n_nodes-2:2) = dNdx(1:n_nodes,2)
B(3,2:2*n_nodes-1:2) = dNdx(1:n_nodes,1)

!set up B_bar elements
B(1,1:2*n_nodes-2:2) = B(1,1:2*n_nodes-2:2) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/2
B(1,2:2*n_nodes-1:2) = B(1,2:2*n_nodes-1:2) + (dNbardx(2:n_nodes,2)-dNdx(2:n_nodes,2))/2

B(2,1:2*n_nodes-2:2) = B(2,1:2*n_nodes-2:2) + (dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/2
B(2,2:2*n_nodes-1:2) = B(2,2:2*n_nodes-1:2) + (dNbardx(2:n_nodes,2)-dNdx(2:n_nodes,2))/2

strain = matmul(B,dof_total)
dstrain = matmul(B,dof_increment)
stress = matmul(D,strain+dstrain)
p = sum(stress(1:3))/3.d0
sdev = stress
sdev(1:3) = sdev(1:3)-p
smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)))! + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
do k = 1,n_field_variables
if (strcmp(field_variable_names(k),'S11',3) ) then
nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
else if (strcmp(field_variable_names(k),'S22',3) ) then
nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
!            else if (strcmp(field_variable_names(k),'S33',3) ) then
!                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
else if (strcmp(field_variable_names(k),'S12',3) ) then
nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
!            else if (strcmp(field_variable_names(k),'S13',3) ) then
!                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
!            else if (strcmp(field_variable_names(k),'S23',3) ) then
!                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
else if (strcmp(field_variable_names(k),'SMISES',6) ) then
nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
endif
end do
end do

else
write(6,*) 'no vaild element identifier'

end if

return
end subroutine fieldvars_linelast_2dbbar

