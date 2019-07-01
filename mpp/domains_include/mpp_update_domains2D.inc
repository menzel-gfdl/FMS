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
subroutine mpp_update_domain_core(domainObj, cptr, isize_in, jsize_in, ksize_in, data_type_in, flags, complete, position, &
     whalo, ehalo, shalo, nhalo, name, tile_count)
  !updates data domain of 3D field whose computational domain have been computed
  type(domain2D),   intent(inout)        :: domainObj
  type(C_PTR),         intent(in)        :: cptr
  integer,          intent(in)           :: isize_in, jsize_in, ksize_in, data_type_in
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count

  integer                 :: update_position, update_whalo, update_ehalo, update_shalo, update_nhalo, ntile
  type(holder_r8_type), dimension(MAX_DOMAIN_FIELDS, MAX_TILES),save :: s_r8
  type(holder_r4_type), dimension(MAX_DOMAIN_FIELDS, MAX_TILES),save :: s_r4
  type(holder_i4_type), dimension(MAX_DOMAIN_FIELDS, MAX_TILES),save :: i_r4
  integer          :: tile, max_ntile, update_flags
  character(len=3) :: text
  logical          :: set_mismatch, is_complete
  logical          :: do_update, send(8), recv(8)
  integer, save    :: isize=0, jsize=0, ke=0, l_size=0, list=0, data_type=0
  integer, save    :: pos, whalosz, ehalosz, shalosz, nhalosz, flags_saved
  type(overlapSpec), pointer :: update => NULL()
  type(overlapSpec), pointer :: check  => NULL()
  type(domain2d_private), pointer :: domain => NULL()

  call check_domain2d_index(domainObj, "mpp_update_domain_core")

  domain => domainList(domainObj%index)
  call get_domain_halos(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, whalo, ehalo, shalo, nhalo)

  !--- when there is NINETY or MINUS_NINETY rotation for some contact, the salar data can not be on E or N-cell,
  if(present(position)) then
     if(domain%rotated_ninety .AND. ( position == EAST .OR. position == NORTH ) )  &
          call mpp_error(FATAL, 'mpp_update_domain_core: when there is NINETY or MINUS_NINETY rotation, ' // &
          'can not use scalar version update_domain for data on E or N-cell' )
  end if

  max_ntile = domain%max_ntile_pe
  ntile = size(domain%x(:))
  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if
  tile = 1

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'mpp_update_domain_core: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "mpp_update_domain_core: "// &
          "optional argument tile_count should be present when number of tiles on this pe is more than 1")
     tile = tile_count
  end if
  do_update = (tile == ntile) .AND. is_complete
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'mpp_update_domain_core: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif
  if(data_type_in == MPP_DOUBLE) then
    call C_F_POINTER(cptr, s_r8(list,tile)%data, (/isize_in,jsize_in,ksize_in/))
  else if(data_type_in == MPP_FLOAT) then
    call C_F_POINTER(cptr, s_r4(list,tile)%data, (/isize_in,jsize_in,ksize_in/))
  else if(data_type_in == MPP_INT) then
    call C_F_POINTER(cptr, i_r4(list,tile)%data, (/isize_in,jsize_in,ksize_in/))
  endif

  update_position = CENTER
  if(present(position)) update_position = position
  update_flags = XUPDATE+YUPDATE   !default
  if( PRESENT(flags) )update_flags = flags

  if(list == 1 .AND. tile == 1 )then
     isize=isize_in; jsize=jsize_in; ke = ksize_in; pos = update_position; data_type = data_type_in
     whalosz = update_whalo; ehalosz = update_ehalo; shalosz = update_shalo; nhalosz = update_nhalo
     flags_saved = update_flags
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize /= isize_in)
     set_mismatch = set_mismatch .OR. (jsize /= jsize_in)
     set_mismatch = set_mismatch .OR. (ke /= ksize_in)
     set_mismatch = set_mismatch .OR. (data_type /= data_type_in)
     set_mismatch = set_mismatch .OR. (update_position /= pos)
     set_mismatch = set_mismatch .OR. (update_whalo /= whalosz)
     set_mismatch = set_mismatch .OR. (update_ehalo /= ehalosz)
     set_mismatch = set_mismatch .OR. (update_shalo /= shalosz)
     set_mismatch = set_mismatch .OR. (update_nhalo /= nhalosz)
     set_mismatch = set_mismatch .OR. (update_flags /= flags_saved)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'mpp_update_domain_core: Incompatible field at count '//text//' for group update.' )
     endif
  endif
  if(is_complete) then
     l_size = list
     list = 0
  end if
  if(do_update )then
     if( domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) )then
        if(debug_update_level .NE. NO_CHECK) then
           check => search_check_overlap(domain, update_position) 
           if(ASSOCIATED(check) ) then
              if( data_type == MPP_DOUBLE) then
                 call do_check_r8(s_r8(1:l_size,1:ntile), check, ke, name )
              else
                 call do_check_r4(s_r4(1:l_size,1:ntile), check, ke, name ) 
              endif
           endif
        endif
        update => search_update_overlap(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, update_position)

        if(debug_message_passing) then
           call get_send_recv_value(update_flags, send, recv) 
           call check_message_size(domain, update, send, recv, name=name)
        endif
        if( data_type == MPP_DOUBLE) then
           call do_update_r8(s_r8(1:l_size,1:ntile), update, ke, update_flags )
        else if(data_type == MPP_FLOAT) then
           call do_update_r4(s_r4(1:l_size,1:ntile), update, ke, update_flags )                 
        else if(data_type == MPP_INT) then
           call do_update_i4(i_r4(1:l_size,1:ntile), update, ke, update_flags )
        else
           call mpp_error(FATAL,"mpp_update_domain_core: only support 4-byte real, 8-byte real and 4-byte int, contact developer")
        endif
     end if
     l_size=0; isize=0;  jsize=0;  ke=0
  endif
  return

end subroutine mpp_update_domain_core

subroutine mpp_update_domain_2d( field, domain, flags, complete, position, &
     whalo, ehalo, shalo, nhalo, name, tile_count)
  !updates data domain of 2D field whose computational domains have been computed
  class(*),         intent(inout)        :: field(:,:)
  type(domain2D),   intent(inout)        :: domain  
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer     :: isize, jsize, ksize, data_type
  type(C_PTR) :: cptr

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_update_domain_2d: field is not contiguous")
  cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = 1
  data_type = get_field_type(field)
  call mpp_update_domain_core( domain, cptr, isize, jsize, ksize, data_type, flags, complete, position, &
       whalo, ehalo, shalo, nhalo, name, tile_count )
  return
end subroutine mpp_update_domain_2d

subroutine mpp_update_domain_3d( field, domain, flags, complete, position, &
     whalo, ehalo, shalo, nhalo, name, tile_count)
  !updates data domain of 2D field whose computational domains have been computed
  class(*),         intent(inout)        :: field(:,:,:)
  type(domain2D),   intent(inout)        :: domain  
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer     :: isize, jsize, ksize, data_type
  type(C_PTR) :: cptr
  
  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_update_domain_3d: field is not contiguous")
  cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = size(field,3)
  data_type = get_field_type(field)
  call mpp_update_domain_core( domain, cptr, isize, jsize, ksize, data_type, flags, complete, position, &
       whalo, ehalo, shalo, nhalo, name, tile_count )
  return
end subroutine mpp_update_domain_3d

subroutine mpp_update_domain_4d( field, domain, flags, complete, position, &
     whalo, ehalo, shalo, nhalo, name, tile_count)
  !updates data domain of 2D field whose computational domains have been computed
  class(*),         intent(inout)        :: field(:,:,:,:)
  type(domain2D),   intent(inout)        :: domain  
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer     :: isize, jsize, ksize, data_type
  type(C_PTR) :: cptr

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_update_domain_4d: field is not contiguous")
  cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = size(field,3)*size(field,4)
  data_type = get_field_type(field)
  call mpp_update_domain_core( domain, cptr, isize, jsize, ksize, data_type, flags, complete, position, &
       whalo, ehalo, shalo, nhalo, name, tile_count )
  return
end subroutine mpp_update_domain_4d

subroutine mpp_update_domain_5d( field, domain, flags, complete, position, &
     whalo, ehalo, shalo, nhalo, name, tile_count)
  !updates data domain of 2D field whose computational domains have been computed
  class(*),         intent(inout)        :: field(:,:,:,:,:)
  type(domain2D),   intent(inout)        :: domain  
  integer,          intent(in), optional :: flags
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer     :: isize, jsize, ksize, data_type
  type(C_PTR) :: cptr

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_update_domain_5d: field is not contiguous")
  cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = size(field,3)*size(field,4)*size(field,5)
  data_type = get_field_type(field)
  call mpp_update_domain_core( domain, cptr, isize, jsize, ksize, data_type, flags, complete, position, &
       whalo, ehalo, shalo, nhalo, name, tile_count )
  return
end subroutine mpp_update_domain_5d

subroutine do_update_r8(sHold, update, ke, update_flags)
  !updates data domain of 3D field whose computational domains have been computed
  type(holder_r8_type),    intent(inout) :: sHold(:,:)
  type(overlapSpec),          intent(in) :: update
  integer,                    intent(in) :: ke
  integer,                    intent(in) :: update_flags
  real(DOUBLE_KIND) :: buffer(mpp_domains_stack_size)

#include "do_update.inc"

  return
end subroutine do_update_r8

subroutine do_update_r4( sHold, update, ke, update_flags)
  !updates data domain of 3D field whose computational domains have been computed
  type(holder_r4_type),    intent(inout) :: sHold(:,:)
  type(overlapSpec),          intent(in) :: update
  integer,                    intent(in) :: ke
  integer,                    intent(in) :: update_flags
  real(FLOAT_KIND) :: buffer(mpp_domains_stack_size)

#include "do_update.inc"

  return
end subroutine do_update_r4

subroutine do_update_i4( sHold, update, ke, update_flags)
  !updates data domain of 3D field whose computational domains have been computed
  type(holder_i4_type),    intent(inout) :: sHold(:,:)
  type(overlapSpec),          intent(in) :: update
  integer,                    intent(in) :: ke
  integer,                    intent(in) :: update_flags
  integer(INT_KIND) :: buffer(mpp_domains_stack_size)

#include "do_update.inc"

  return
end subroutine do_update_i4

subroutine mpp_update_domainV_core( domainObj, cptr_x, cptr_y, isizex_in, jsizex_in, isizey_in, jsizey_in, ksizex_in, ksizey_in, &
                                    data_typex_in, data_typey_in, flags, gridtype, complete,  &
                                    whalo, ehalo, shalo, nhalo, name, tile_count)
  !updates data domain of 3D field whose computational domains have been computed
  type(domain2D),   intent(inout)        :: domainObj
   type(C_PTR),        intent(in)        :: cptr_x, cptr_y
  integer,          intent(in)           :: isizex_in, jsizex_in, ksizex_in, data_typex_in
  integer,          intent(in)           :: isizey_in, jsizey_in, ksizey_in, data_typey_in
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  
  integer                                :: update_whalo, update_ehalo, update_shalo, update_nhalo, ntile    
  integer                                :: grid_offset_type
  logical                                :: exchange_uv
  type(holder_r8_type), dimension(MAX_DOMAIN_FIELDS, MAX_TILES),save :: x_r8, y_r8
  type(holder_r4_type), dimension(MAX_DOMAIN_FIELDS, MAX_TILES),save :: x_r4, y_r4
  logical          :: do_update, is_complete
  integer, save    :: isize(2)=0,jsize(2)=0,ke=0,l_size=0, offset_type=0, list=0
  integer, save    :: whalosz, ehalosz, shalosz, nhalosz, flags_saved, data_type=0
  integer          :: tile, max_ntile, update_flags
  integer          :: position_x, position_y
  logical          :: set_mismatch, send(8), recv(8)
  character(len=3) :: text
  type(overlapSpec),  pointer :: updatex => NULL()
  type(overlapSpec),  pointer :: updatey => NULL()
  type(overlapSpec),  pointer :: checkx  => NULL()
  type(overlapSpec),  pointer :: checky  => NULL()
  type(domain2d_private), pointer :: domain=>NULL()

  call check_domain2d_index(domainObj, "mpp_update_domainV_core")

  domain=>domainList(domainObj%index)

  call get_domain_halos(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, whalo, ehalo, shalo, nhalo)

  grid_offset_type = AGRID
  if( PRESENT(gridtype) ) grid_offset_type = gridtype

  exchange_uv = .false.
  if(grid_offset_type == DGRID_NE) then
     exchange_uv = .true.
     grid_offset_type = CGRID_NE
  else if( grid_offset_type == DGRID_SW ) then
     exchange_uv = .true.
     grid_offset_type = CGRID_SW
  end if

  max_ntile = domain%max_ntile_pe
  ntile = size(domain%x(:))

  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if
  tile = 1

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'mpp_update_domainV_core: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "mpp_update_domainV_core: "// &
          "optional argument tile_count should be present when number of tiles on some pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'mpp_update_domainV_core: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif

  update_flags = XUPDATE+YUPDATE   !default
  if( PRESENT(flags) ) then
     update_flags = flags
     if(BTEST(update_flags,SCALAR_BIT)) then
        if (.NOT.(BTEST(update_flags,WEST) .OR. BTEST(update_flags,EAST) &
              .OR. BTEST(update_flags,NORTH) .OR. BTEST(update_flags,SOUTH))) &
           update_flags = update_flags + XUPDATE+YUPDATE   !default with SCALAR_PAIR
     endif
  endif
  if(data_typex_in == MPP_DOUBLE) then
    call C_F_POINTER(cptr_x, x_r8(list, tile)%data, (/isizex_in,jsizex_in,ksizex_in/))
    call C_F_POINTER(cptr_y, y_r8(list, tile)%data, (/isizey_in,jsizey_in,ksizey_in/))
  else if(data_typey_in == MPP_FLOAT) then
    call C_F_POINTER(cptr_x, x_r4(list, tile)%data, (/isizex_in,jsizex_in,ksizex_in/))
    call C_F_POINTER(cptr_y, y_r4(list, tile)%data, (/isizey_in,jsizey_in,ksizey_in/))
  endif

  if(list == 1 .AND. tile == 1)then
     if( BTEST(update_flags,NORTH) .AND. BTEST(domain%fold,NORTH) .AND. BTEST(grid_offset_type,SOUTH) ) &
           call mpp_error( FATAL, 'mpp_update_domainV_core: Incompatible grid offset and fold.' )
     if(ksizex_in .NE. ksizey_in) call mpp_error(FATAL,  &
          'mpp_update_domainV_core: mismatch of ksize between fieldx and fieldy')
     if(data_typex_in .NE. data_typey_in) call mpp_error(FATAL,  &
          'mpp_update_domainV_core: mismatch of data type between fieldx and fieldy')
     isize(1)=isizex_in; jsize(1)=jsizex_in; ke = ksizex_in
     isize(2)=isizey_in; jsize(2)=jsizey_in; data_type = data_typex_in
     offset_type = grid_offset_type
     flags_saved = update_flags
     whalosz = update_whalo; ehalosz = update_ehalo; shalosz = update_shalo; nhalosz = update_nhalo
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize(1) /= isizex_in)
     set_mismatch = set_mismatch .OR. (jsize(1) /= jsizex_in)
     set_mismatch = set_mismatch .OR. (ke /= ksizex_in)
     set_mismatch = set_mismatch .OR. (isize(2) /= isizey_in)
     set_mismatch = set_mismatch .OR. (jsize(2) /= jsizey_in)
     set_mismatch = set_mismatch .OR. (ke /= ksizey_in)
     set_mismatch = set_mismatch .OR. (data_type /= data_typex_in)
     set_mismatch = set_mismatch .OR. (data_type /= data_typey_in)
     set_mismatch = set_mismatch .OR. (grid_offset_type /= offset_type)
     set_mismatch = set_mismatch .OR. (update_whalo /= whalosz)
     set_mismatch = set_mismatch .OR. (update_ehalo /= ehalosz)
     set_mismatch = set_mismatch .OR. (update_shalo /= shalosz)
     set_mismatch = set_mismatch .OR. (update_nhalo /= nhalosz)
     set_mismatch = set_mismatch .OR. (update_flags /= flags_saved)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'MPP_UPDATE_3D_V: Incompatible field at count '//text//' for group vector update.' )
     end if
  end if
  if(is_complete) then
     l_size = list
     list = 0
  end if
  if(do_update)then
     if( domain_update_is_needed(domain, update_whalo, update_ehalo, update_shalo, update_nhalo) )then
        select case(grid_offset_type)
        case (AGRID)
           position_x = CENTER
           position_y = CENTER
        case (BGRID_NE, BGRID_SW)
           position_x = CORNER
           position_y = CORNER
        case (CGRID_NE, CGRID_SW)
           position_x = EAST
           position_y = NORTH
        case default
           call mpp_error(FATAL, "mpp_update_domain2D.h: invalid value of grid_offset_type")
        end select

        if(debug_update_level .NE. NO_CHECK) then
           checkx => search_check_overlap(domain, position_x)
           checky => search_check_overlap(domain, position_y)
           if(ASSOCIATED(checkx)) then
              if(exchange_uv) then
                 if( data_type == MPP_DOUBLE) then
                    call do_checkV_r8(x_r8(1:l_size,1:ntile),y_r8(1:l_size,1:ntile), domain,        &
                         checky, checkx, ke, update_flags, name)
                 else
                    call do_checkV_r4(x_r4(1:l_size,1:ntile),y_r4(1:l_size,1:ntile), domain,        &
                         checky, checkx, ke, update_flags, name)
                 endif
              else
                 if( data_type == MPP_DOUBLE) then
                    call do_checkV_r8(x_r8(1:l_size,1:ntile),y_r8(1:l_size,1:ntile), domain,        &
                         checkx, checky, ke, update_flags, name)
                 else
                    call do_checkV_r4(x_r4(1:l_size,1:ntile),y_r4(1:l_size,1:ntile), domain,        &
                         checkx, checky, ke, update_flags, name)
                 endif
              end if
           endif
        endif
        updatex => search_update_overlap(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, position_x)
        updatey => search_update_overlap(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, position_y)
        if(debug_message_passing) then
           call get_send_recv_value(update_flags, send, recv)
           call check_message_size_vector(domain, updatex, updatey, send, recv)
        endif
        if(exchange_uv) then
           if( data_type == MPP_DOUBLE) then
              call do_updateV_r8(x_r8(1:l_size,1:ntile),y_r8(1:l_size,1:ntile), domain, updatey, updatex, &
                   ke, grid_offset_type, update_flags)
           else
              call do_updateV_r4(x_r4(1:l_size,1:ntile),y_r4(1:l_size,1:ntile), domain, updatey, updatex, &
                   ke, grid_offset_type, update_flags)
           endif
        else
           if( data_type == MPP_DOUBLE) then
              call do_updateV_r8(x_r8(1:l_size,1:ntile),y_r8(1:l_size,1:ntile), domain, updatex, updatey, &
                   ke,grid_offset_type, update_flags)
           else
              call do_updateV_r4(x_r4(1:l_size,1:ntile),y_r4(1:l_size,1:ntile), domain, updatex, updatey, &
                   ke,grid_offset_type, update_flags)
           endif
        end if
     end if
     l_size=0; isize=0;  jsize=0;  ke=0
  end if

  return
end subroutine mpp_update_domainV_core


subroutine mpp_update_domain_2dv( fieldx, fieldy, domain, flags, gridtype, complete, &
     whalo, ehalo, shalo, nhalo, name, tile_count)
  !updates data domain of 2D field whose computational domains have been computed
  class(*),         intent(inout)        :: fieldx(:,:), fieldy(:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer     :: isizex,jsizex,isizey,jsizey,ksizex,ksizey,data_typex,data_typey
  type(C_PTR) :: cptrx, cptry

  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_update_domain_2dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_update_domain_2dv: fieldy is not contiguous")
  cptrx = C_LOC(fieldx)
  cptry = C_LOC(fieldy)
  isizex = size(fieldx,1)
  jsizex = size(fieldx,2)
  isizey = size(fieldy,1)
  jsizey = size(fieldy,2)
  ksizex = 1
  ksizey = 1
  data_typex = get_field_type(fieldx)
  data_typey = get_field_type(fieldy)

  call mpp_update_domainV_core(domain, cptrx, cptry, isizex, jsizex, isizey, jsizey, ksizex, ksizey, &
                data_typex, data_typey, flags, gridtype, complete, whalo, ehalo, shalo, nhalo, name, tile_count)
  return
end subroutine mpp_update_domain_2dv

subroutine mpp_update_domain_3dv( fieldx, fieldy, domain, flags, gridtype, complete, &
     whalo, ehalo, shalo, nhalo, name, tile_count)
  !updates data domain of 2D field whose computational domains have been computed
  class(*),         intent(inout)        :: fieldx(:,:,:), fieldy(:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer     :: isizex,jsizex,isizey,jsizey,ksizex,ksizey,data_typex,data_typey 
  type(C_PTR) :: cptrx, cptry
  
  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_update_domain_3dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_update_domain_3dv: fieldy is not contiguous")
  cptrx = C_LOC(fieldx)
  cptry = C_LOC(fieldy)
  isizex = size(fieldx,1)
  jsizex = size(fieldx,2)
  isizey = size(fieldy,1)
  jsizey = size(fieldy,2)
  ksizex = size(fieldx,3)
  ksizey = size(fieldy,3)
  data_typex = get_field_type(fieldx)
  data_typey = get_field_type(fieldy)

  call mpp_update_domainV_core(domain, cptrx, cptry, isizex, jsizex, isizey, jsizey, ksizex, ksizey, &
                data_typex, data_typey, flags, gridtype, complete, whalo, ehalo, shalo, nhalo, name, tile_count)
  return
end subroutine mpp_update_domain_3dv

subroutine mpp_update_domain_4dv( fieldx, fieldy, domain, flags, gridtype, complete, &
     whalo, ehalo, shalo, nhalo, name, tile_count)
  !updates data domain of 2D field whose computational domains have been computed
  class(*),         intent(inout)        :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer     :: isizex,jsizex,isizey,jsizey,ksizex,ksizey,data_typex,data_typey
  type(C_PTR) :: cptrx, cptry
  
  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_update_domain_4dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_update_domain_4dv: fieldy is not contiguous")
  cptrx = C_LOC(fieldx)
  cptry = C_LOC(fieldy)
  isizex = size(fieldx,1)
  jsizex = size(fieldx,2)
  isizey = size(fieldy,1)
  jsizey = size(fieldy,2)
  ksizex = size(fieldx,3)*size(fieldx,4)
  ksizey = size(fieldy,3)*size(fieldy,4)
  data_typex = get_field_type(fieldx)
  data_typey = get_field_type(fieldy)

  call mpp_update_domainV_core(domain, cptrx, cptry, isizex, jsizex, isizey, jsizey, ksizex, ksizey, &
                data_typex, data_typey, flags, gridtype, complete, whalo, ehalo, shalo, nhalo, name, tile_count)
  return
end subroutine mpp_update_domain_4dv

subroutine mpp_update_domain_5dv( fieldx, fieldy, domain, flags, gridtype, complete, &
     whalo, ehalo, shalo, nhalo, name, tile_count)
  !updates data domain of 2D field whose computational domains have been computed
  class(*),         intent(inout)        :: fieldx(:,:,:,:,:), fieldy(:,:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  logical,          intent(in), optional :: complete
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  integer     :: isizex,jsizex,isizey,jsizey,ksizex,ksizey,data_typex,data_typey
  type(C_PTR) :: cptrx, cptry
  
  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_update_domain_5dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_update_domain_5dv: fieldy is not contiguous")
  cptrx = C_LOC(fieldx)
  cptry = C_LOC(fieldy)
  isizex = size(fieldx,1)
  jsizex = size(fieldx,2)
  isizey = size(fieldy,1)
  jsizey = size(fieldy,2)
  ksizex = size(fieldx,3)*size(fieldx,4)*size(fieldx,5)
  ksizey = size(fieldy,3)*size(fieldy,4)*size(fieldy,5)
  data_typex = get_field_type(fieldx)
  data_typey = get_field_type(fieldy)

  call mpp_update_domainV_core(domain, cptrx, cptry, isizex, jsizex, isizey, jsizey, ksizex, ksizey, &
                data_typex, data_typey, flags, gridtype, complete, whalo, ehalo, shalo, nhalo, name, tile_count)
  return
end subroutine mpp_update_domain_5dv

subroutine do_updateV_r8(xHold, yHold, domain, update_x, update_y, ke, gridtype, update_flags)
!updates data domain of 3D field whose computational domains have been computed
  type(holder_r8_type), intent(inout) :: xHold(:,:), yHold(:,:)
  type(domain2d_private),  intent(in) :: domain
  type(overlapSpec),       intent(in) :: update_x, update_y
  integer,                 intent(in) :: ke
  integer,                 intent(in) :: gridtype
  integer,                 intent(in) :: update_flags
  integer :: ioff, joff
  real(DOUBLE_KIND) :: buffer(mpp_domains_stack_size)

#include "do_updateV.inc"

end subroutine do_updateV_r8

subroutine do_updateV_r4(xHold, yHold, domain, update_x, update_y, ke, gridtype, update_flags)
!updates data domain of 3D field whose computational domains have been computed
  type(holder_r4_type),    intent(inout) :: xHold(:,:), yHold(:,:)
  type(domain2d_private),  intent(in)    :: domain
  type(overlapSpec),   intent(in)        :: update_x, update_y
  integer,             intent(in)        :: ke
  integer, intent(in)                    :: gridtype
  integer, intent(in)                    :: update_flags
  integer :: ioff, joff
  real(FLOAT_KIND) :: buffer(mpp_domains_stack_size)

#include "do_updateV.inc"

end subroutine do_updateV_r4

!######################################################################
subroutine do_check_r8(sHold, check, ke, name)
   type(holder_r8_type),       intent(in) :: sHold(:,:)
   type(overlapSpec),          intent(in) :: check
   integer,                    intent(in) :: ke
   character(len=*), optional, intent(in) :: name
   real(DOUBLE_KIND) :: buffer(mpp_domains_stack_size)

#include "do_check.inc"

end subroutine do_check_r8

!######################################################################
subroutine do_check_r4( sHold, check, ke, name)
   type(holder_r4_type),       intent(in) :: sHold(:,:)
   type(overlapSpec),          intent(in) :: check
   integer,                    intent(in) :: ke
   character(len=*), optional, intent(in) :: name
   real(FLOAT_KIND) :: buffer(mpp_domains_stack_size)

#include "do_check.inc"

end subroutine do_check_r4

!######################################################################
subroutine do_checkV_r8(xHold, yHold, domain, check_x, check_y, ke, flags, name)
   type(holder_r8_type), intent(inout)    :: xHold(:,:), yHold(:,:)
   type(domain2d_private),  intent(in)    :: domain
   type(overlapSpec),   intent(in)        :: check_x, check_y
   integer,             intent(in)        :: ke
   integer,             intent(in)        :: flags
   character(len=*), intent(in), optional :: name
   real(DOUBLE_KIND) :: buffer(mpp_domains_stack_size)

#include "do_checkV.inc"

end subroutine do_checkV_r8

subroutine do_checkV_r4(xHold, yHold, domain, check_x, check_y, ke, flags, name)
   type(holder_r4_type), intent(inout)    :: xHold(:,:), yHold(:,:)
   type(domain2d_private),  intent(in)    :: domain
   type(overlapSpec),   intent(in)        :: check_x, check_y
   integer,             intent(in)        :: ke
   integer,             intent(in)        :: flags
   character(len=*), intent(in), optional :: name
   real(FLOAT_KIND) :: buffer(mpp_domains_stack_size)

#include "do_checkV.inc"

end subroutine do_checkV_r4

