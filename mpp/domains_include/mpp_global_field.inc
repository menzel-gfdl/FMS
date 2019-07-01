!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
subroutine mpp_global_field_2d( domain, local, global, flags, position,tile_count, default_data)
  type(domain2D),      intent(in) :: domain
  class(*), intent(in)            :: local(:,:)
  class(*), intent(out)           :: global(:,:)
  integer,  intent(in),  optional :: flags
  integer,  intent(in),  optional :: position
  integer,  intent(in),  optional :: tile_count
  class(*), intent(in),  optional :: default_data

  call check_domain2d_index(domain, "mpp_global_field_2d")
  call global_field_2d_private(domainList(domain%index), local, global, flags, position,tile_count, default_data)

end subroutine mpp_global_field_2d

subroutine global_field_2d_private( domain, local, global, flags, position,tile_count, default_data)
  type(domain2D_private), intent(in) :: domain
  class(*), intent(in)            :: local(:,:)
  class(*), intent(out)           :: global(:,:)
  integer,  intent(in),  optional :: flags
  integer,  intent(in),  optional :: position
  integer,  intent(in),  optional :: tile_count
  class(*), intent(in),  optional :: default_data
  integer           :: l_size(3), g_size(3)
  class(*), pointer :: local3D(:,:,:) => NULL()
  class(*), pointer :: global3D(:,:,:) => NULL()

  if( .NOT.module_is_initialized )call mpp_error( FATAL, 'mpp_global_field_2d: must first call mpp_domains_init.' )
  l_size(1) = size(local,1);  l_size(2) = size(local,2);  l_size(3)=1
  g_size(1) = size(global,1); g_size(2) = size(global,2); g_size(3)=1
  call point_3D_to_2D(l_size(1), l_size(2), local, local3D)
  call point_3D_to_2D(g_size(1), g_size(2), global, global3D)
  call global_field_3d_private(domain, local3D, global3d, flags, position,tile_count, default_data)

end subroutine global_field_2d_private

subroutine mpp_global_field_3d( domain, local, global, flags, position,tile_count, default_data)
  type(domain2D),     intent(in) :: domain
  class(*),          intent(in)  :: local(:,:,:)
  class(*),          intent(out) :: global(:,:,:)
  integer,  intent(in), optional :: flags
  integer,  intent(in), optional :: position
  integer,  intent(in), optional :: tile_count
  class(*), intent(in), optional :: default_data

  call check_domain2d_index(domain, "mpp_global_field_3d")
  call global_field_3d_private(domainList(domain%index), local, global, flags, position,tile_count, default_data)

end subroutine mpp_global_field_3d

subroutine global_field_3d_private( domain, local, global, flags, position,tile_count, default_data)
  type(domain2D_private), intent(in) :: domain
  class(*),          intent(in)  :: local(:,:,:)
  class(*),          intent(out) :: global(:,:,:)
  integer,  intent(in), optional :: flags
  integer,  intent(in), optional :: position
  integer,  intent(in), optional :: tile_count
  class(*), intent(in), optional :: default_data
  integer     :: l_size(3), g_size(3)

  logical :: xonly, yonly, root_only, global_on_this_pe
  integer :: tile, ishift, jshift, stackuse, ipos, jpos, ioff, joff
  integer :: local_type, global_type, default_type
  real(DOUBLE_KIND) :: default_r8
  real(FLOAT_KIND)  :: default_r4
  character(len=8) :: text

  if( .NOT.module_is_initialized )call mpp_error( FATAL, 'mpp_global_field_3d: must first call mpp_domains_init.' )
  tile = 1; if(PRESENT(tile_count)) tile = tile_count
  l_size(1) = size(local,1);  l_size(2) = size(local,2);  l_size(3)=size(local,3)
  g_size(1) = size(global,1); g_size(2) = size(global,2); g_size(3)=size(global,3)
  xonly = .FALSE.
  yonly = .FALSE.
  root_only = .FALSE.
  if( PRESENT(flags) ) then
     xonly = BTEST(flags,EAST)
     yonly = BTEST(flags,SOUTH)
     if( .NOT.xonly .AND. .NOT.yonly )call mpp_error( WARNING,  &
          'mpp_global_field_3d: you must have flags=XUPDATE, YUPDATE or XUPDATE+YUPDATE' )
     if(xonly .AND. yonly) then
        xonly = .false.; yonly = .false.
     endif
     root_only = BTEST(flags, ROOT_GLOBAL)
     if( (xonly .or. yonly) .AND. root_only ) then
        call mpp_error( WARNING, 'mpp_global_field_3d: flags = XUPDATE+GLOBAL_ROOT_ONLY or ' // &
             'flags = YUPDATE+GLOBAL_ROOT_ONLY is not supported, will ignore GLOBAL_ROOT_ONLY' )
        root_only = .FALSE.
     endif
  endif

  call mpp_get_domain_shift(domain, ishift, jshift, position)
  stackuse = (domain%x(1)%compute%size+ishift)*(domain%y(1)%compute%size+jshift)*l_size(3)
  if( stackuse.GT.mpp_domains_stack_size )then
     write( text, '(i8)' )stackuse
     call mpp_error( FATAL, &
          'mpp_global_field_3d: user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
  end if
  mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, stackuse )

  global_on_this_pe =  .NOT. root_only .OR. domain%pe == domain%tile_root_pe
  ipos = 1-domain%x(tile)%global%begin; jpos = 1-domain%y(tile)%global%begin
  if(global_on_this_pe ) then
     if(l_size(3) .NE. g_size(3) ) call mpp_error( FATAL, &
          'mpp_global_field_3d: mismatch of third dimension size of global and local')
     if( g_size(1) .NE.(domain%x(tile)%global%size+ishift) .OR. g_size(2).NE.(domain%y(tile)%global%size+jshift))then
        if(xonly) then
           if(g_size(1).NE.(domain%x(tile)%global%size+ishift) .OR. g_size(2).NE.(domain%y(tile)%compute%size+jshift)) &
                call mpp_error( FATAL, 'mpp_global_field_3d: incoming arrays do not match domain for xonly global field.' )
           jpos = domain%y(tile)%global%begin-domain%y(tile)%compute%begin 
        else if(yonly) then
           if(g_size(1).NE.(domain%x(tile)%compute%size+ishift) .OR. g_size(2).NE.(domain%y(tile)%global%size+jshift)) &
                call mpp_error( FATAL, 'mpp_global_field_3d: incoming arrays do not match domain for yonly global field.' )
           ipos = domain%x(tile)%global%begin-domain%x(tile)%compute%begin
        else
           call mpp_error( FATAL, 'mpp_global_field_3d: incoming arrays do not match domain.' )
        endif
     endif
  endif

  if(l_size(1).EQ.(domain%x(tile)%compute%size+ishift) .AND. l_size(2).EQ.(domain%y(tile)%compute%size+jshift) )then
     !local is on compute domain
     ioff =  -domain%x(tile)%compute%begin + 1
     joff =  -domain%y(tile)%compute%begin + 1
  else if(l_size(1).EQ.(domain%x(tile)%memory%size+ishift) .AND. l_size(2).EQ.(domain%y(tile)%memory%size+jshift))then
     !local is on data domain
     ioff = - domain%x(tile)%data%begin + 1
     joff = - domain%y(tile)%data%begin + 1
  else
     call mpp_error( FATAL, 'mpp_global_field_3d: incoming field array must match either compute domain or memory domain.' )
  end if

  local_type = get_field_type(local)
  global_type = get_field_type(global)

  if(local_type .NE. global_type ) call mpp_error( FATAL, 'mpp_global_field_3d: local and glboal are not the same data type')
  default_r8 = 0
  default_r4 = 0
  if(present(default_data)) then
     default_type = get_field_type(default_data)
     if(local_type .NE. default_type) call mpp_error( FATAL, &
           'mpp_global_field_3d: local and default_data are not the same data type')
     select type(default_data)
     type is(real(DOUBLE_KIND))
        default_r8 = default_data
     type is(real(FLOAT_KIND))
        default_r4 = default_data
     class default
        call mpp_error( FATAL, 'mpp_global_field_3d: optional argument default_data must be 4-byte or 8-byte real')
     end select
  endif

  select type(local)
  type is(real(DOUBLE_KIND))
     select type(global)
     type is(real(DOUBLE_KIND))
        call do_global_field_r8(domain, local, global, xonly, yonly, &
                                root_only, tile, ishift, jshift, ioff, joff, ipos, jpos, default_r8)
     end select
  type is(real(FLOAT_KIND))
     select type(global)
     type is(real(FLOAT_KIND))
        call do_global_field_r4(domain, local, global, xonly, yonly, &
                                root_only, tile, ishift, jshift, ioff, joff, ipos, jpos, default_r4)
     end select
  class default
     call mpp_error(FATAL, "mpp_global_field_3d: only support 4-byte and 8-byte real, contact developer")   
  end select

end subroutine global_field_3d_private

subroutine mpp_global_field_4d( domain, local, global, flags, position,tile_count, default_data)
  type(domain2D),     intent(in) :: domain
  class(*), intent(in)           :: local(:,:,:,:)
  class(*), intent(out)          :: global(:,:,:,:)
  integer,  intent(in), optional :: flags
  integer,  intent(in), optional :: position
  integer,  intent(in), optional :: tile_count
  class(*), intent(in), optional :: default_data

  call check_domain2d_index(domain, "mpp_global_field_4d")
  call global_field_4d_private(domainList(domain%index), local, global, flags, position,tile_count, default_data)

end subroutine mpp_global_field_4d

subroutine global_field_4d_private( domain, local, global, flags, position,tile_count, default_data)
  type(domain2D_private), intent(in) :: domain
  class(*), intent(in)           :: local(:,:,:,:)
  class(*), intent(out)          :: global(:,:,:,:)
  integer,  intent(in), optional :: flags
  integer,  intent(in), optional :: position
  integer,  intent(in), optional :: tile_count
  class(*), intent(in), optional :: default_data
  integer     :: l_size(3), g_size(3)
  class(*), pointer :: local3D(:,:,:) => NULL()
  class(*), pointer :: global3D(:,:,:) => NULL()

  if( .NOT.module_is_initialized )call mpp_error( FATAL, 'mpp_global_field_4d: must first call mpp_domains_init.' )
  l_size(1) = size(local,1); l_size(2) = size(local,2); l_size(3)=size(local,3)*size(local,4)
  g_size(1) =size(global,1); g_size(2) =size(global,2); g_size(3)=size(global,3)*size(global,4)
  call point_3D_to_4D(l_size(1), l_size(2), size(local,3), size(local,4), local, local3D)
  call point_3D_to_4D(g_size(1), g_size(2), size(global,3), size(global,4), global, global3D)

  call global_field_3d_private(domain, local3D, global3D, flags, position,tile_count, default_data)

end subroutine global_field_4d_private

subroutine mpp_global_field_5d( domain, local, global, flags, position,tile_count, default_data)
  type(domain2D),     intent(in)  :: domain
  class(*), intent(in)            :: local(:,:,:,:,:)
  class(*), intent(out)           :: global(:,:,:,:,:)
  integer,  intent(in),  optional :: flags
  integer,  intent(in),  optional :: position
  integer,  intent(in),  optional :: tile_count
  class(*), intent(in),  optional :: default_data

  call check_domain2d_index(domain, "mpp_global_field_5d")
  call global_field_5d_private(domain, local, global, flags, position,tile_count, default_data)

end subroutine mpp_global_field_5d

subroutine global_field_5d_private( domain, local, global, flags, position,tile_count, default_data)
  type(domain2D),     intent(in)  :: domain
  class(*), intent(in)            :: local(:,:,:,:,:)
  class(*), intent(out)           :: global(:,:,:,:,:)
  integer,  intent(in),  optional :: flags
  integer,  intent(in),  optional :: position
  integer,  intent(in),  optional :: tile_count
  class(*), intent(in),  optional :: default_data
  integer           :: l_size(3), g_size(3)
  class(*), pointer :: local3D(:,:,:) => NULL()
  class(*), pointer :: global3D(:,:,:) => NULL()

  if( .NOT.module_is_initialized )call mpp_error( FATAL, 'mpp_global_field_5d: must first call mpp_domains_init.' )
  l_size(1) = size(local,1);  l_size(2) = size(local,2);  l_size(3)=size(local,3)*size(local,4)*size(local,5)
  g_size(1) = size(global,1); g_size(2) = size(global,2); g_size(3)=size(global,3)*size(global,4)*size(global,5)
  call point_3D_to_5D(l_size(1), l_size(2), size(local,3), size(local,4), size(local,5), local, local3D)
  call point_3D_to_5D(g_size(1), g_size(2), size(global,3), size(global,4), size(global,5), global, global3D)
  call global_field_3d_private(domainList(domain%index), local3D, global3D, flags, position,tile_count, default_data)

end subroutine global_field_5d_private

subroutine do_global_field_r8(domain, local, global, xonly, yonly, &
     root_only, tile, ishift, jshift, ioff, joff, ipos, jpos, default_data)
  type(domain2D_private), intent(in) :: domain
  real(DOUBLE_KIND),  intent(in) :: local (:,:,:)
  real(DOUBLE_KIND), intent(out) :: global(:,:,:)
  logical,            intent(in) :: xonly, yonly, root_only
  integer,            intent(in) :: tile, ishift, jshift, ioff, joff, ipos, jpos
  real(DOUBLE_KIND),  intent(in) :: default_data
  real(DOUBLE_KIND) :: buffer(mpp_domains_stack_size)

#include "do_global_field.inc"

end subroutine do_global_field_r8

subroutine do_global_field_r4(domain, local, global, xonly, yonly, &
     root_only, tile, ishift, jshift, ioff, joff, ipos, jpos, default_data)
  type(domain2D_private), intent(in) :: domain
  real(FLOAT_KIND),   intent(in) :: local(:,:,:)
  real(FLOAT_KIND),  intent(out) :: global(:,:,:)
  logical,            intent(in) :: xonly, yonly, root_only
  integer,            intent(in) :: tile, ishift, jshift, ioff, joff, ipos, jpos
  real(FLOAT_KIND),   intent(in) :: default_data
  real(FLOAT_KIND) :: buffer(mpp_domains_stack_size)

#include "do_global_field.inc"

end subroutine do_global_field_r4




