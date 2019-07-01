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
subroutine mpp_update_nest_fine_core(nest_domain, f_cptr_in, isize_in, jsize_in, ksize_in, data_type_in, &
                                     b_cptr_in, b_size_in, flags, complete, position, extra_halo, tile_count) 
    type(nest_domain_type), intent(inout)  :: nest_domain
    type(C_PTR),       intent(in) :: f_cptr_in, b_cptr_in(:)
    integer,           intent(in) :: isize_in, jsize_in, ksize_in, data_type_in
    integer,           intent(in) :: b_size_in(:,:)
    integer,          intent(in), optional :: flags
    logical,          intent(in), optional :: complete
    integer,          intent(in), optional :: position
    integer,          intent(in), optional :: extra_halo
    integer,          intent(in), optional :: tile_count

   type(nestSpec), pointer :: update=>NULL()
   type(holder_r8_type), dimension(MAX_DOMAIN_FIELDS),save :: f_r8
   type(holder_r4_type), dimension(MAX_DOMAIN_FIELDS),save :: f_r4
   type(holder_r8_type), dimension(4,MAX_DOMAIN_FIELDS),save :: b_r8
   type(holder_r4_type), dimension(4,MAX_DOMAIN_FIELDS),save :: b_r4
   character(len=3) :: text
   logical          :: is_complete, set_mismatch
   integer          :: tile, i
   integer          :: add_halo, update_flags, update_position
   integer, save    :: isize=0, jsize=0, ksize=0, l_size=0, data_type=0
   integer, save    :: b_size(4,2) = 0
   integer, save    :: add_halo_save=0, update_flags_save=0, update_position_save=0
   integer, save    :: list=0 

   add_halo = 0
   if(present(extra_halo)) add_halo = add_halo
   update_position = CENTER
   if(present(position)) update_position = position
   update_flags = XUPDATE+YUPDATE   !default
   if( PRESENT(flags) )update_flags = flags

   is_complete = .true.
   if(PRESENT(complete)) then
      is_complete = complete
   end if
   tile = 1
   if(present(tile_count)) tile = tile_count
   if( tile > 1 ) then
      call mpp_error(FATAL,' mpp_update_nest_fine_core: currently do not support multiple tile per pe')
   endif

   list = list+1
   if(list > MAX_DOMAIN_FIELDS)then
      write( text,'(i2)' ) MAX_DOMAIN_FIELDS
      call mpp_error(FATAL,' mpp_update_nest_fine_core: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
   endif

   if(data_type_in == MPP_DOUBLE) then
      call C_F_POINTER(f_cptr_in, f_r8(list)%data, (/isize_in,jsize_in,ksize_in/))
      do i = 1, 4
         call C_F_POINTER(b_cptr_in(i), b_r8(i,list)%data, (/b_size_in(i,1),b_size_in(i,2),ksize_in/))
      enddo
   else if(data_type_in == MPP_FLOAT) then
      call C_F_POINTER(f_cptr_in, f_r4(list)%data, (/isize_in,jsize_in,ksize_in/))
      do i = 1, 4
         call C_F_POINTER(b_cptr_in(i), b_r4(i,list)%data, (/b_size_in(i,1),b_size_in(i,2),ksize_in/))
      enddo
   endif

   if(list == 1)then
      isize=isize_in; jsize=jsize_in; ksize = ksize_in; data_type = data_type_in
      b_size = b_size_in
      update_position_save = update_position
      update_flags_save    = update_flags
      add_halo_save = add_halo
   else
      set_mismatch = .false.
      set_mismatch = set_mismatch .OR. (isize /= isize_in)
      set_mismatch = set_mismatch .OR. (jsize /= jsize_in)
      set_mismatch = set_mismatch .OR. (ksize /= ksize_in)
      set_mismatch = set_mismatch .OR. (update_position_save /= update_position)
      set_mismatch = set_mismatch .OR. ANY( b_size .NE. b_size_in )
      set_mismatch = set_mismatch .OR. (data_type .NE. data_type_in)
      set_mismatch = set_mismatch .OR. (update_flags_save /= update_flags)
      set_mismatch = set_mismatch .OR. (add_halo_save /= add_halo)

      if(set_mismatch)then
         write( text,'(i2)' ) list
         call mpp_error(FATAL,' mpp_update_nest_fine_core: Incompatible field at count '//text//' for group update.' )
      endif
   endif

   if(is_complete) then
      l_size = list
      update => search_C2F_nest_overlap(nest_domain, add_halo, update_position)
      if(data_type == MPP_DOUBLE) then
         call do_update_nest_fine_r8(f_r8(1:l_size), update, b_r8(:, 1:l_size), ksize, update_flags )
      else if(data_type == MPP_FLOAT) then
         call do_update_nest_fine_r4(f_r4(1:l_size), update, b_r4(:, 1:l_size), ksize, update_flags )
      else
         call mpp_error(FATAL,"mpp_update_nest_fine_core: only support 4-byte and 8-byte real, contact developer")
      endif
      list = 0
   endif


end subroutine mpp_update_nest_fine_core

!#######################################################################
subroutine mpp_update_nest_fine_2d(field, nest_domain, wbuffer, ebuffer, sbuffer, nbuffer, &
     flags, complete, position, extra_halo, tile_count) 
  class(*),              intent(in)      :: field(:,:)
  type(nest_domain_type), intent(inout)  :: nest_domain
  class(*),              intent(inout)   :: wbuffer(:,:)
  class(*),              intent(inout)   :: ebuffer(:,:)
  class(*),              intent(inout)   :: sbuffer(:,:)
  class(*),              intent(inout)   :: nbuffer(:,:)
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: extra_halo
  integer,          intent(in), optional :: tile_count
  integer     :: isize, jsize, ksize, data_type, b_size(4,2)
  type(C_PTR) :: f_cptr, b_cptr(4)

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_update_nest_fine_2d: field is not contiguous")
  f_cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = 1
  data_type = get_field_type(field)
  b_size(1,1)  = size(ebuffer,1)
  b_size(2,1)  = size(sbuffer,1)
  b_size(3,1)  = size(wbuffer,1)
  b_size(4,1)  = size(nbuffer,1)
  b_size(1,2)  = size(ebuffer,2)
  b_size(2,2)  = size(sbuffer,2)
  b_size(3,2)  = size(wbuffer,2)
  b_size(4,2)  = size(nbuffer,2)
  b_cptr(1) = C_LOC(ebuffer)
  b_cptr(2) = C_LOC(sbuffer)
  b_cptr(3) = C_LOC(wbuffer)
  b_cptr(4) = C_LOC(nbuffer)

  call mpp_update_nest_fine_core(nest_domain, f_cptr, isize, jsize, ksize, data_type, &
       b_cptr, b_size, flags, complete, position, extra_halo, tile_count) 

  return

end subroutine mpp_update_nest_fine_2d

!#######################################################################
subroutine mpp_update_nest_fine_3d(field, nest_domain, wbuffer, ebuffer, sbuffer, nbuffer, &
     flags, complete, position, extra_halo, tile_count) 
  class(*),              intent(in)      :: field(:,:,:)
  type(nest_domain_type), intent(inout)  :: nest_domain
  class(*),              intent(inout)   :: wbuffer(:,:,:)
  class(*),              intent(inout)   :: ebuffer(:,:,:)
  class(*),              intent(inout)   :: sbuffer(:,:,:)
  class(*),              intent(inout)   :: nbuffer(:,:,:)
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: extra_halo
  integer,          intent(in), optional :: tile_count
  integer     :: isize, jsize, ksize, data_type, b_size(4,2)
  type(C_PTR) :: f_cptr, b_cptr(4)

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_update_nest_fine_3d: field is not contiguous")
  f_cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = size(field,3)
  data_type = get_field_type(field)
  b_size(1,1)  = size(ebuffer,1)
  b_size(2,1)  = size(sbuffer,1)
  b_size(3,1)  = size(wbuffer,1)
  b_size(4,1)  = size(nbuffer,1)
  b_size(1,2)  = size(ebuffer,2)
  b_size(2,2)  = size(sbuffer,2)
  b_size(3,2)  = size(wbuffer,2)
  b_size(4,2)  = size(nbuffer,2)
  b_cptr(1) = C_LOC(ebuffer)
  b_cptr(2) = C_LOC(sbuffer)
  b_cptr(3) = C_LOC(wbuffer)
  b_cptr(4) = C_LOC(nbuffer)

  call mpp_update_nest_fine_core(nest_domain, f_cptr, isize, jsize, ksize, data_type, &
       b_cptr, b_size, flags, complete, position, extra_halo, tile_count) 

  return

end subroutine mpp_update_nest_fine_3d

!#######################################################################
subroutine mpp_update_nest_fine_4d(field, nest_domain, wbuffer, ebuffer, sbuffer, nbuffer, &
     flags, complete, position, extra_halo, tile_count) 
  class(*),              intent(in)      :: field(:,:,:,:)
  type(nest_domain_type), intent(inout)  :: nest_domain
  class(*),              intent(inout)   :: wbuffer(:,:,:,:)
  class(*),              intent(inout)   :: ebuffer(:,:,:,:)
  class(*),              intent(inout)   :: sbuffer(:,:,:,:)
  class(*),              intent(inout)   :: nbuffer(:,:,:,:)
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: extra_halo
  integer,          intent(in), optional :: tile_count
  integer     :: isize, jsize, ksize, data_type, b_size(4,2)
  type(C_PTR) :: f_cptr, b_cptr(4)

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_update_nest_fine_4d: field is not contiguous")
  f_cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = size(field,3)*size(field,4)
  data_type = get_field_type(field)
  b_size(1,1)  = size(ebuffer,1)
  b_size(2,1)  = size(sbuffer,1)
  b_size(3,1)  = size(wbuffer,1)
  b_size(4,1)  = size(nbuffer,1)
  b_size(1,2)  = size(ebuffer,2)
  b_size(2,2)  = size(sbuffer,2)
  b_size(3,2)  = size(wbuffer,2)
  b_size(4,2)  = size(nbuffer,2)
  b_cptr(1) = C_LOC(ebuffer)
  b_cptr(2) = C_LOC(sbuffer)
  b_cptr(3) = C_LOC(wbuffer)
  b_cptr(4) = C_LOC(nbuffer)

  call mpp_update_nest_fine_core(nest_domain, f_cptr, isize, jsize, ksize, data_type, &
       b_cptr, b_size, flags, complete, position, extra_halo, tile_count) 

  return

end subroutine mpp_update_nest_fine_4d

!#######################################################################
subroutine do_update_nest_fine_r8(f_hold, update, b_hold, ke, flags)
  type(holder_r8_type),    intent(in) :: f_hold(:)
  type(nestSpec),          intent(in) :: update
  type(holder_r8_type), intent(inout) :: b_hold(:,:)
  integer,                 intent(in) :: ke
  integer,                 intent(in) :: flags
  real(DOUBLE_KIND) :: buffer(mpp_domains_stack_size)

#include "do_update_nest_fine.inc"

end subroutine do_update_nest_fine_r8

!#######################################################################
subroutine do_update_nest_fine_r4(f_hold, update, b_hold, ke, flags)
  type(holder_r4_type),    intent(in) :: f_hold(:)
  type(nestSpec),          intent(in) :: update
  type(holder_r4_type), intent(inout) :: b_hold(:,:)
  integer,                 intent(in) :: ke
  integer,                 intent(in) :: flags
  real(FLOAT_KIND) :: buffer(mpp_domains_stack_size)

#include "do_update_nest_fine.inc"

end subroutine do_update_nest_fine_r4

!#######################################################################
subroutine mpp_update_nest_coarse_core(nest_domain, f_cptr_in, isize_in, jsize_in, ksize_in, data_type_in, &
                                       b_cptr_in, b_size_in,  complete, position, tile_count)
   type(nest_domain_type), intent(inout)  :: nest_domain
   type(C_PTR),                intent(in) :: f_cptr_in, b_cptr_in
   integer,                    intent(in) :: isize_in, jsize_in, ksize_in, data_type_in
   integer,                    intent(in) :: b_size_in(:)
   logical,          intent(in), optional :: complete
   integer,          intent(in), optional :: position
   integer,          intent(in), optional :: tile_count

   type(nestSpec), pointer :: update=>NULL()
   type(holder_r8_type), dimension(MAX_DOMAIN_FIELDS),save :: f_r8, b_r8
   type(holder_r4_type), dimension(MAX_DOMAIN_FIELDS),save :: f_r4, b_r4
   character(len=3) :: text
   logical          :: is_complete, set_mismatch
   integer          :: tile, i
   integer          :: update_position
   integer, save    :: update_position_save=0
   integer, save    :: isize=0, jsize=0, ksize=0, l_size=0, data_type=0
   integer, save    :: b_size(2) = 0
   integer, save    :: list=0

   update_position = CENTER
   if(present(position)) update_position = position

   is_complete = .true.
   if(PRESENT(complete)) then
      is_complete = complete
   end if
   tile = 1
   if(present(tile_count)) tile = tile_count
   if( tile > 1 ) then
      call mpp_error(FATAL,'mpp_update_nest_coarse_core: currently do not support multiple tile per pe')
   endif

   list = list+1
   if(list > MAX_DOMAIN_FIELDS)then
      write( text,'(i2)' ) MAX_DOMAIN_FIELDS
      call mpp_error(FATAL,'mpp_update_nest_coarse_core: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
   endif
   if(data_type_in == MPP_DOUBLE) then
      call C_F_POINTER(f_cptr_in, f_r8(list)%data, (/isize_in,jsize_in,ksize_in/))
      call C_F_POINTER(b_cptr_in, b_r8(list)%data, (/b_size_in(1),b_size_in(2),ksize_in/))
   else if(data_type_in == MPP_FLOAT) then
      call C_F_POINTER(f_cptr_in, f_r4(list)%data, (/isize_in,jsize_in,ksize_in/))
      call C_F_POINTER(b_cptr_in, b_r4(list)%data, (/b_size_in(1),b_size_in(2),ksize_in/))
   endif

   if(list == 1)then
      isize=isize_in; jsize=jsize_in; ksize = ksize_in; data_type = data_type_in
      b_size = b_size_in
      update_position_save = update_position
   else
      set_mismatch = .false.
      set_mismatch = set_mismatch .OR. (isize /= isize_in)
      set_mismatch = set_mismatch .OR. (jsize /= jsize_in)
      set_mismatch = set_mismatch .OR. (ksize /= ksize_in)
      set_mismatch = set_mismatch .OR. ANY(b_size .NE. b_size_in) 
      set_mismatch = set_mismatch .OR. (data_type .NE. data_type_in)
      set_mismatch = set_mismatch .OR. (update_position_save /= update_position)

      if(set_mismatch)then
         write( text,'(i2)' ) list
         call mpp_error(FATAL,'mpp_update_nest_coarse_core: Incompatible field at count '//text//' for group update.' )
      endif
   endif

   if(is_complete) then
      l_size = list
      list = 0
   end if

   if(is_complete)then
      update => search_F2C_nest_overlap(nest_domain, update_position)
      if(data_type == MPP_DOUBLE) then
         call do_update_nest_coarse_r8(f_r8(1:l_size), update, b_r8(1:l_size), ksize)
      else if(data_type == MPP_FLOAT) then
         call do_update_nest_coarse_r4(f_r4(1:l_size), update, b_r4(1:l_size), ksize)
      else
         call mpp_error(FATAL,"mpp_update_nest_coarse_core: only support 4-byte and 8-byte real, contact developer")
      endif
   endif

end subroutine mpp_update_nest_coarse_core

!#######################################################################
subroutine mpp_update_nest_coarse_2d(field, nest_domain, buffer, complete, position,  tile_count) 
  class(*),              intent(in)      :: field(:,:)
  type(nest_domain_type), intent(inout)  :: nest_domain
  class(*),              intent(inout)   :: buffer(:,:)
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: tile_count
  integer     :: isize, jsize, ksize, data_type, b_size(2)
  type(C_PTR) :: f_cptr, b_cptr

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_update_nest_coarse_2d: field is not contiguous")
  f_cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = 1
  data_type = get_field_type(field)
  b_size(1)  = size(buffer,1)
  b_size(2)  = size(buffer,2)
  b_cptr = C_LOC(buffer)

  call mpp_update_nest_coarse_core(nest_domain, f_cptr, isize, jsize, ksize, data_type, &
       b_cptr, b_size, complete, position, tile_count) 

  return

end subroutine mpp_update_nest_coarse_2d

!#######################################################################
subroutine mpp_update_nest_coarse_3d(field, nest_domain, buffer, complete, position, tile_count) 
  class(*),              intent(in)      :: field(:,:,:)
  type(nest_domain_type), intent(inout)  :: nest_domain
  class(*),              intent(inout)   :: buffer(:,:,:)
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: tile_count
  integer     :: isize, jsize, ksize, data_type, b_size(2)
  type(C_PTR) :: f_cptr, b_cptr

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_update_nest_coarse_3d: field is not contiguous")
  f_cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = size(field,3)
  data_type = get_field_type(field)
  b_size(1)  = size(buffer,1)
  b_size(2)  = size(buffer,2)
  b_cptr = C_LOC(buffer)

  call mpp_update_nest_coarse_core(nest_domain, f_cptr, isize, jsize, ksize, data_type, &
       b_cptr, b_size, complete, position, tile_count) 

  return

end subroutine mpp_update_nest_coarse_3d

!#######################################################################
subroutine mpp_update_nest_coarse_4d(field, nest_domain, buffer, complete, position, tile_count) 
  class(*),              intent(in)      :: field(:,:,:,:)
  type(nest_domain_type), intent(inout)  :: nest_domain
  class(*),              intent(inout)   :: buffer(:,:,:,:)
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: tile_count
  integer     :: isize, jsize, ksize, data_type, b_size(2)
  type(C_PTR) :: f_cptr, b_cptr

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_update_nest_coarse_4d: field is not contiguous")
  f_cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = size(field,3)*size(field,4)
  data_type = get_field_type(field)
  b_size(1)  = size(buffer,1)
  b_size(2)  = size(buffer,2)
  b_cptr = C_LOC(buffer)

  call mpp_update_nest_coarse_core(nest_domain, f_cptr, isize, jsize, ksize, data_type, &
       b_cptr, b_size, complete, position, tile_count) 

  return

end subroutine mpp_update_nest_coarse_4d

!#######################################################################
subroutine do_update_nest_coarse_r8(f_hold, update, b_hold, ke)
  type(holder_r8_type),    intent(in) :: f_hold(:)
  type(nestSpec),          intent(in) :: update
  type(holder_r8_type), intent(inout) :: b_hold(:)
  integer,                 intent(in) :: ke
  real(DOUBLE_KIND) :: buffer(mpp_domains_stack_size)

#include "do_update_nest_coarse.inc"

end subroutine do_update_nest_coarse_r8

!#######################################################################
subroutine do_update_nest_coarse_r4(f_hold, update, b_hold, ke)
  type(holder_r4_type),    intent(in) :: f_hold(:)
  type(nestSpec),          intent(in) :: update
  type(holder_r4_type), intent(inout) :: b_hold(:)
  integer,                 intent(in) :: ke
  real(FLOAT_KIND) :: buffer(mpp_domains_stack_size)

#include "do_update_nest_coarse.inc"

end subroutine do_update_nest_coarse_r4
