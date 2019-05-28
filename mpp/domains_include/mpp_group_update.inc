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
subroutine create_group_update_core(group, domain, cptr, fsize, data_type, &
                                    flags, position, whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D_private),      intent(inout) :: domain  
  type(C_PTR),          intent(in)           :: cptr
  integer,              intent(in)           :: fsize(3), data_type
  integer,              intent(in), optional :: flags
  integer,              intent(in), optional :: position
  integer,              intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.

  integer          :: update_position, update_whalo, update_ehalo, update_shalo, update_nhalo
  integer          :: update_flags
  integer          :: nscalar
  character(len=3) :: text
  logical          :: set_mismatch
  logical          :: recv(8), send(8)

  if(group%initialized) then
      call mpp_error(FATAL, "create_group_update_core: group is already initialized")
  endif
  call get_domain_halos(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, whalo, ehalo, shalo, nhalo)
  update_position = CENTER
  !--- when there is NINETY or MINUS_NINETY rotation for some contact, the salar data can not be on E or N-cell,
  if(present(position)) then
     update_position = position 
     if(domain%rotated_ninety .AND. ( position == EAST .OR. position == NORTH ) )  &
          call mpp_error(FATAL, 'create_group_update_core: hen there is NINETY or MINUS_NINETY rotation, ' // &
          'can not use scalar version update_domain for data on E or N-cell' )
  end if

  if( domain%max_ntile_pe > 1 ) then
     call mpp_error(FATAL,'create_group_update_core: do not support multiple tile per processor')
  endif

  if(group%data_type >0) then
     if(data_type .NE. group%data_type) call mpp_error(FATAL,  &
           'create_group_update_core: mismatch of data_type between scalar and or vector fields')
  else
     group%data_type = data_type
     if(data_type .NE. MPP_DOUBLE .and. data_type .NE. MPP_FLOAT) &
        call mpp_error(FATAL,'create_group_update_core: only support 4-byte and 8-byte real, contact developer')
  endif

  update_flags = XUPDATE+YUPDATE 
  if(present(flags)) update_flags = flags
  group%nscalar = group%nscalar + 1
  nscalar = group%nscalar
  if( nscalar > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'mpp_create_group_update_core: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif
  if(group%data_type == MPP_DOUBLE) then
    call C_F_POINTER(cptr, group%s_r8(nscalar,1)%data, (/fsize(1),fsize(2),fsize(3)/))
  else if(group%data_type == MPP_FLOAT) then
    call C_F_POINTER(cptr, group%s_r4(nscalar,1)%data, (/fsize(1),fsize(2),fsize(3)/))
  endif

  if( group%ksize_s == 0 ) then
     group%flags_s  = update_flags
     group%whalo_s  = update_whalo
     group%ehalo_s  = update_ehalo
     group%shalo_s  = update_shalo
     group%nhalo_s  = update_nhalo
     group%position = update_position
     group%isize_s  = fsize(1)
     group%jsize_s  = fsize(2)
     group%ksize_s  = fsize(3)
     call get_send_recv_value(update_flags, send, recv)
     group%recv_s = recv
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (group%flags_s  .NE. update_flags)
     set_mismatch = set_mismatch .OR. (group%whalo_s  .NE. update_whalo)
     set_mismatch = set_mismatch .OR. (group%ehalo_s  .NE. update_ehalo)
     set_mismatch = set_mismatch .OR. (group%shalo_s  .NE. update_shalo)
     set_mismatch = set_mismatch .OR. (group%nhalo_s  .NE. update_nhalo)
     set_mismatch = set_mismatch .OR. (group%position .NE. update_position)
     set_mismatch = set_mismatch .OR. (group%isize_s  .NE. fsize(1))
     set_mismatch = set_mismatch .OR. (group%jsize_s  .NE. fsize(2))
     set_mismatch = set_mismatch .OR. (group%ksize_s  .NE. fsize(3))

     if(set_mismatch)then
        write( text,'(i2)' ) nscalar
        call mpp_error(FATAL,'create_group_update_core: Incompatible field at count '//text//' for group update.' )
     endif
  endif

  return

end subroutine create_group_update_core

subroutine mpp_create_group_update_2d(group, field, domainObj, flags, position, whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                    intent(inout) :: field(:,:)
  type(domain2D),       intent(inout)        :: domainObj
  integer,              intent(in), optional :: flags
  integer,              intent(in), optional :: position
  integer,              intent(in), optional :: whalo, ehalo, shalo, nhalo
  integer        :: fsize(3), data_type, nscalar
  type(C_PTR)    :: cptr
  type(domain2d_private), pointer :: domain=>NULL()

  call check_domain2d_index(domainObj, "mpp_create_group_update_2d")
  domain => domainList(domainObj%index)

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_create_group_update_2d: field is not contiguous")        
  fsize(1) = size(field,1)
  fsize(2) = size(field,2)
  fsize(3) = 1
  nscalar = group%nscalar
  data_type = get_field_type(field)  
  cptr = C_LOC(field)
  call create_group_update_core(group, domain, cptr, fsize, data_type, flags, position, whalo, ehalo, shalo, nhalo)
  return

end subroutine mpp_create_group_update_2d

subroutine mpp_create_group_update_3d(group, field, domainObj, flags, position, whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                    intent(inout) :: field(:,:,:)
  type(domain2D),              intent(inout) :: domainObj
  integer,              intent(in), optional :: flags
  integer,              intent(in), optional :: position
  integer,              intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  integer     :: fsize(3), data_type
  type(C_PTR)       :: cptr
  type(domain2d_private), pointer :: domain=>NULL()

  call check_domain2d_index(domainObj, "mpp_create_group_update_3d")
  domain => domainList(domainObj%index)

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_create_group_update_3d: field is not contiguous")
  fsize(1) = size(field,1)
  fsize(2) = size(field,2)
  fsize(3) = size(field,3)
  data_type = get_field_type(field)
  cptr = C_LOC(field)
  call create_group_update_core(group, domain, cptr, fsize, data_type, flags, position, whalo, ehalo, shalo, nhalo)
  return

end subroutine mpp_create_group_update_3d

subroutine mpp_create_group_update_4d(group, field, domainObj, flags, position, whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                    intent(inout) :: field(:,:,:,:)
  type(domain2D),              intent(inout) :: domainObj
  integer,              intent(in), optional :: flags
  integer,              intent(in), optional :: position
  integer,              intent(in), optional :: whalo, ehalo, shalo, nhalo
  integer     :: fsize(3), data_type
  type(C_PTR)       :: cptr  
  type(domain2d_private), pointer :: domain=>NULL()

  call check_domain2d_index(domainObj, "mpp_create_group_update_4d")
  domain => domainList(domainObj%index)

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_create_group_update_4d: field is not contiguous")
  fsize(1) = size(field,1)
  fsize(2) = size(field,2)
  fsize(3) = size(field,3)*size(field,4)
  data_type = get_field_type(field)
  cptr = C_LOC(field)
  call create_group_update_core(group, domain, cptr, fsize, data_type, flags, position, whalo, ehalo, shalo, nhalo)

  return

end subroutine mpp_create_group_update_4d

subroutine create_group_updateV_core( group, domain, cptrx, cptry, fsizex, fsizey, &
                                      data_typex, data_typey, flags, gridtype, whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D_private),      intent(inout) :: domain
  type(C_PTR),      intent(in)               :: cptrx, cptry
  integer,                        intent(in) :: fsizex(3), data_typex
  integer,                        intent(in) :: fsizey(3), data_typey
  integer,              intent(in), optional :: flags, gridtype
  integer,              intent(in), optional :: whalo, ehalo, shalo, nhalo

  integer          :: update_whalo, update_ehalo, update_shalo, update_nhalo
  integer          :: update_flags
  integer          :: nvector, update_gridtype, position_x, position_y
  character(len=3) :: text
  logical          :: set_mismatch
  logical          :: recv_x(8), send_x(8)
  logical          :: recv_y(8), send_y(8)

  if(group%initialized) then
     call mpp_error(FATAL, "create_group_updateV_core: group is already initialized")
  endif
  call get_domain_halos(domain, update_whalo, update_ehalo, update_shalo, update_nhalo, whalo, ehalo, shalo, nhalo)

  update_gridtype = AGRID
  if(PRESENT(gridtype)) update_gridtype = gridtype

  if( domain%max_ntile_pe > 1 ) then
     call mpp_error(FATAL,'create_group_updateV_core: do not support multiple tile per processor')
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

  group%nvector = group%nvector + 1
  nvector = group%nvector
  if( nvector > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'create_group_updateV_core: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif

  if(fsizex(3) .NE. fsizey(3)) call mpp_error(FATAL,  &
           'create_group_updateV_core: mismatch of ksize between fieldx and fieldy')
  if(data_typex .NE. data_typey) call mpp_error(FATAL,  &
           'create_group_updateV_core: mismatch of data_type between fieldx and fieldy')
  if(group%data_type >0) then
     if(data_typex .NE. group%data_type) call mpp_error(FATAL,  &
           'create_group_updateV_core: mismatch of data_type between scalar and or vector fields')
  else
     group%data_type = data_typex
     if(data_typex .NE. MPP_DOUBLE .and. data_typex .NE. MPP_FLOAT) &
        call mpp_error(FATAL,'create_group_updateV_core: only support 4-byte and 8-byte real, contact developer')
  endif

  if(group%data_type == MPP_DOUBLE) then
    call C_F_POINTER(cptrx, group%x_r8(nvector,1)%data, (/fsizex(1),fsizex(2),fsizex(3)/))
    call C_F_POINTER(cptry, group%y_r8(nvector,1)%data, (/fsizey(1),fsizey(2),fsizey(3)/))
  else if(data_typey == MPP_FLOAT) then
    call C_F_POINTER(cptrx, group%x_r4(nvector,1)%data, (/fsizex(1),fsizex(2),fsizex(3)/))
    call C_F_POINTER(cptry, group%y_r4(nvector,1)%data, (/fsizey(1),fsizey(2),fsizey(3)/))
  endif

  if( group%ksize_v == 0 ) then
     group%flags_v  = update_flags
     group%whalo_v  = update_whalo
     group%ehalo_v  = update_ehalo
     group%shalo_v  = update_shalo
     group%nhalo_v  = update_nhalo
     group%gridtype = update_gridtype
     group%isize_x  = fsizex(1)
     group%jsize_x  = fsizex(2)
     group%isize_y  = fsizey(1)
     group%jsize_y  = fsizey(2)
     group%ksize_v  = fsizex(3)

     call get_vector_send_recv_value(update_flags, update_gridtype, domain%symmetry, group%nonsym_edge, &
                                     send_x, send_y, recv_x, recv_y)
     group%recv_x = recv_x
     group%recv_y = recv_y
     select case(group%gridtype)
     case (AGRID)
        position_x = CENTER
        position_y = CENTER
     case (BGRID_NE, BGRID_SW)
        position_x = CORNER
        position_y = CORNER
     case (CGRID_NE, CGRID_SW)
        position_x = EAST
        position_y = NORTH
     case (DGRID_NE, DGRID_SW)
        position_x = NORTH
        position_y = EAST
     case default
        call mpp_error(FATAL, "create_group_updateV_core: invalid value of gridtype")
     end select

  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (group%flags_v  .NE. update_flags)
     set_mismatch = set_mismatch .OR. (group%whalo_v  .NE. update_whalo)
     set_mismatch = set_mismatch .OR. (group%ehalo_v  .NE. update_ehalo)
     set_mismatch = set_mismatch .OR. (group%shalo_v  .NE. update_shalo)
     set_mismatch = set_mismatch .OR. (group%nhalo_v  .NE. update_nhalo)
     set_mismatch = set_mismatch .OR. (group%gridtype .NE. update_gridtype)
     set_mismatch = set_mismatch .OR. (group%isize_x  .NE. fsizex(1))
     set_mismatch = set_mismatch .OR. (group%jsize_x  .NE. fsizex(2))
     set_mismatch = set_mismatch .OR. (group%isize_y  .NE. fsizey(1))
     set_mismatch = set_mismatch .OR. (group%jsize_y  .NE. fsizey(2))
     set_mismatch = set_mismatch .OR. (group%ksize_v  .NE. fsizex(3))

     if(set_mismatch)then
        write( text,'(i2)' ) nvector
        call mpp_error(FATAL,'create_group_updateV_core: Incompatible field at count '//text//' for group update.' )
     endif
  endif
  
  return

end subroutine create_group_updateV_core

subroutine mpp_create_group_update_2dv( group, fieldx, fieldy, domainObj, flags, gridtype, &
                                        whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                    intent(inout) :: fieldx(:,:), fieldy(:,:)
  type(domain2D),              intent(inout) :: domainObj
  integer,          intent(in),     optional :: flags, gridtype
  integer,          intent(in),     optional :: whalo, ehalo, shalo, nhalo
  integer :: fsizex(3), fsizey(3), data_typex, data_typey
  type(C_PTR)       :: cptrx, cptry
  type(domain2d_private), pointer :: domain=>NULL()

  call check_domain2d_index(domainObj, "mpp_create_group_update_2dv")
  domain => domainList(domainObj%index)

  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_create_group_update_2dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_create_group_update_2dv: fieldy is not contiguous")
  fsizex(1) = size(fieldx,1)
  fsizey(1) = size(fieldy,1)
  fsizex(2) = size(fieldx,2)
  fsizey(2) = size(fieldy,2)
  fsizex(3) = 1
  fsizey(3) = 1
  data_typex = get_field_type(fieldx)
  data_typey = get_field_type(fieldy)
  cptrx = C_LOC(fieldx)
  cptry = C_LOC(fieldy)

  call  create_group_updateV_core(group, domain, cptrx, cptry, fsizex, fsizey, &
                data_typex, data_typey,flags, gridtype, whalo, ehalo, shalo, nhalo)
  return

end subroutine mpp_create_group_update_2dv

subroutine mpp_create_group_update_3dv( group, fieldx, fieldy, domainObj, flags, gridtype, &
                                        whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                    intent(inout) :: fieldx(:,:,:), fieldy(:,:,:)
  type(domain2D),              intent(inout) :: domainObj
  integer,          intent(in),     optional :: flags, gridtype
  integer,          intent(in),     optional :: whalo, ehalo, shalo, nhalo
  integer     :: fsizex(3), fsizey(3), data_typex, data_typey
  type(C_PTR)       :: cptrx, cptry
  type(domain2d_private), pointer :: domain=>NULL()

  call check_domain2d_index(domainObj, "mpp_create_group_update_3dv")
  domain => domainList(domainObj%index)  

  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_create_group_update_3dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_create_group_update_3dv: fieldy is not contiguous")
  fsizex(1) = size(fieldx,1)
  fsizey(1) = size(fieldy,1)
  fsizex(2) = size(fieldx,2)
  fsizey(2) = size(fieldy,2)
  fsizex(3) = size(fieldx,3)
  fsizey(3) = size(fieldy,3)
  data_typex = get_field_type(fieldx)
  data_typey = get_field_type(fieldy)
  cptrx = C_LOC(fieldx)
  cptry = C_LOC(fieldy)

  call  create_group_updateV_core(group, domain, cptrx, cptry, fsizex, fsizey, &
                data_typex, data_typey,flags, gridtype, whalo, ehalo, shalo, nhalo)

  return

end subroutine mpp_create_group_update_3dv

subroutine mpp_create_group_update_4dv( group, fieldx, fieldy, domainObj, flags, gridtype, &
                                         whalo, ehalo, shalo, nhalo)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                    intent(inout) :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  type(domain2D),              intent(inout) :: domainObj
  integer,          intent(in),     optional :: flags, gridtype
  integer,          intent(in),     optional :: whalo, ehalo, shalo, nhalo
  integer     :: fsizex(3), fsizey(3), data_typex, data_typey
  type(C_PTR)       :: cptrx, cptry
  type(domain2d_private), pointer :: domain=>NULL()

  call check_domain2d_index(domainObj, "mpp_create_group_update_4dv")
  domain => domainList(domainObj%index)  

  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_create_group_update_4dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_create_group_update_4dv: fieldy is not contiguous")
  fsizex(1) = size(fieldx,1)
  fsizey(1) = size(fieldy,1)
  fsizex(2) = size(fieldx,2)
  fsizey(2) = size(fieldy,2)
  fsizex(3) = size(fieldx,3)*size(fieldx,4)
  fsizey(3) = size(fieldy,3)*size(fieldy,4)
  data_typex = get_field_type(fieldx)
  data_typey = get_field_type(fieldy)
  cptrx = C_LOC(fieldx)
  cptry = C_LOC(fieldy)

  call  create_group_updateV_core(group, domain, cptrx, cptry, fsizex, fsizey, &
                data_typex, data_typey,flags, gridtype, whalo, ehalo, shalo, nhalo)
  return

end subroutine mpp_create_group_update_4dv

subroutine mpp_do_group_update(group, domainObj)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D),              intent(inout) :: domainObj 
  real(FLOAT_KIND), pointer :: buffer_r4(:) => NULL()
  type(domain2d_private), pointer :: domain=>NULL()

  call check_domain2d_index(domainObj, "mpp_do_group_update")
  domain => domainList(domainObj%index)  

  if(num_nonblock_group_update .NE. 0) then
     call mpp_error(FATAL,"mpp_do_group_update: can not be called between mpp_start_group_update/mpp_complete_group_update")
  endif

  if(.not. group%initialized) call set_group_update(group,domain)
  group%buffer_start_pos = 0
  if( group%data_type == MPP_DOUBLE) then
     call start_group_update_r8(group, domain, group%s_r8, group%x_r8, group%y_r8, mpp_domains_stack)
     call complete_group_update_r8(group, domain, group%s_r8, group%x_r8, group%y_r8, mpp_domains_stack)
  else if(  group%data_type == MPP_FLOAT) then
     call C_F_POINTER(C_LOC(mpp_domains_stack), buffer_r4, (/mpp_domains_stack_size/) )
     call start_group_update_r4(group, domain, group%s_r4, group%x_r4, group%y_r4, buffer_r4)
     call complete_group_update_r4(group, domain, group%s_r4, group%x_r4, group%y_r4, buffer_r4)
  else
     call mpp_error(FATAL,"mpp_do_group_update: only support 4-byte and 8-byte real, contact developer")
  endif

end subroutine mpp_do_group_update

subroutine mpp_start_group_update(group, domainObj, reuse_buffer)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D),              intent(inout) :: domainObj  
  logical,  optional,          intent(in)    :: reuse_buffer
  logical :: reuse_buf_pos
  character(len=8)          :: text
  real(FLOAT_KIND), pointer :: buffer_r4(:) => NULL()
  type(domain2d_private), pointer :: domain=>NULL()

  call check_domain2d_index(domainObj, "mpp_start_group_update")
  domain => domainList(domainObj%index)

  if(.not. group%initialized) call set_group_update(group,domain)
  reuse_buf_pos = .FALSE.
  if (PRESENT(reuse_buffer)) reuse_buf_pos = reuse_buffer
  if (.not. reuse_buf_pos) then
     group%buffer_start_pos = nonblock_group_buffer_pos
     nonblock_group_buffer_pos = nonblock_group_buffer_pos + group%tot_msgsize
     mpp_domains_stack_hwm = nonblock_group_buffer_pos + 1
     if( mpp_domains_stack_hwm .GT. mpp_domains_stack_size )then
        write( text,'(i8)' )mpp_domains_stack_hwm
        call mpp_error( FATAL, 'mpp_start_group_update: mpp_domains_stack overflow, '// &
                     'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
     end if

  else if( group%buffer_start_pos < 0 ) then
     call mpp_error(FATAL, "mpp_start_group_update: group%buffer_start_pos is not set")
  endif

  num_nonblock_group_update = num_nonblock_group_update + 1

  if( group%data_type == MPP_DOUBLE) then
     call start_group_update_r8(group, domain, group%s_r8, group%x_r8, group%y_r8, mpp_domains_stack)
  else if(  group%data_type == MPP_FLOAT) then
     call C_F_POINTER(C_LOC(mpp_domains_stack), buffer_r4, (/mpp_domains_stack_size/) )
     call start_group_update_r4(group, domain, group%s_r4, group%x_r4, group%y_r4, buffer_r4)
  else
     call mpp_error(FATAL,"mpp_start_group_update: only support 4-byte and 8-byte real, contact developer")
  endif

end subroutine mpp_start_group_update

subroutine mpp_complete_group_update(group, domainObj)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D),              intent(inout) :: domainObj  
  real(FLOAT_KIND), pointer :: buffer_r4(:) => NULL()
  type(domain2d_private), pointer :: domain=>NULL()

  call check_domain2d_index(domainObj, "mpp_complete_group_update")
  domain => domainList(domainObj%index)

  if(num_nonblock_group_update < 1) call mpp_error(FATAL, &
    'mpp_start_group_update must be called before calling mpp_end_group_update')
  num_nonblock_group_update = num_nonblock_group_update - 1

  if( group%data_type == MPP_DOUBLE) then
     call complete_group_update_r8(group, domain, group%s_r8, group%x_r8, group%y_r8, mpp_domains_stack)
  else if(  group%data_type == MPP_FLOAT) then
     call C_F_POINTER(C_LOC(mpp_domains_stack), buffer_r4, (/mpp_domains_stack_size/) )
     call complete_group_update_r4(group, domain, group%s_r4, group%x_r4, group%y_r4, buffer_r4)
  else
     call mpp_error(FATAL,"mpp_complete_group_update: only support 4-byte and 8-byte real, contact developer")
  endif

  if( num_nonblock_group_update == 0) then
     nonblock_group_buffer_pos   = 0
  endif

end subroutine mpp_complete_group_update

subroutine start_group_update_r8(group, domain, sHold, xHold, yHold, buffer)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D_private),      intent(inout) :: domain
  type(holder_r8_type),        intent(inout) :: sHold(:,:), xHold(:,:), yHold(:,:)
  real(DOUBLE_KIND),           intent(inout) :: buffer(mpp_domains_stack_size)  

  integer   :: nvector, nlist
  integer   :: nsend, nrecv, update_flags
  integer   :: msgsize
  integer   :: from_pe, to_pe, buffer_pos, pos
  integer   :: is, ie, js, je, ke
  integer   :: n, l, m, buffer_start_pos
  integer   :: gridtype
  integer   :: npack, rotation
  integer   :: mpitype

#ifdef use_libMPI
  mpitype = MPI_REAL8
#else
  mpitype = 0
#endif

#include "start_group_update.inc"

end subroutine start_group_update_r8

subroutine complete_group_update_r8(group, domain, sHold, xHold, yHold, buffer)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D_private),      intent(inout) :: domain  
  type(holder_r8_type),        intent(inout) :: sHold(:,:), xHold(:,:), yHold(:,:)
  real(DOUBLE_KIND),           intent(inout) :: buffer(mpp_domains_stack_size)

  integer   :: nscalar, nvector, l_size, tMe
  logical   :: recv_x(8), recv_y(8)
  integer   :: nsend, nrecv, update_flags
  integer   :: buffer_pos, pos
  integer   :: is, ie, js, je, ke
  integer   :: isd, ied, jsd, jed
  integer   :: n, l, m, i, j, k, buffer_start_pos
  integer   :: shift, gridtype, midpoint, ioff, joff
  integer   :: nunpack, rotation, whalo_v, ehalo_v, shalo_v, nhalo_v

#include "complete_group_update.inc"

end subroutine complete_group_update_r8

subroutine start_group_update_r4(group, domain, sHold, xHold, yHold, buffer)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D_private),      intent(inout) :: domain  
  type(holder_r4_type),        intent(inout) :: sHold(:,:), xHold(:,:), yHold(:,:)
  real(FLOAT_KIND),           intent(inout) :: buffer(mpp_domains_stack_size)

  integer   :: nvector, nlist
  integer   :: nsend, nrecv, update_flags
  integer   :: msgsize
  integer   :: from_pe, to_pe, buffer_pos, pos
  integer   :: is, ie, js, je, ke
  integer   :: n, l, m, buffer_start_pos
  integer   :: gridtype
  integer   :: npack, rotation
  integer   :: mpitype

#ifdef use_libMPI
  mpitype = MPI_REAL4
#else
  mpitype = 0
#endif

#include "start_group_update.inc"

end subroutine start_group_update_r4

subroutine complete_group_update_r4(group, domain, sHold, xHold, yHold, buffer)
  type(mpp_group_update_type), intent(inout) :: group
  type(domain2D_private),      intent(inout) :: domain  
  type(holder_r4_type),        intent(inout) :: sHold(:,:), xHold(:,:), yHold(:,:)
  real(FLOAT_KIND),           intent(inout) :: buffer(mpp_domains_stack_size)

  integer   :: nscalar, nvector, l_size, tMe
  logical   :: recv_x(8), recv_y(8)
  integer   :: nsend, nrecv, update_flags
  integer   :: buffer_pos, pos
  integer   :: is, ie, js, je, ke
  integer   :: isd, ied, jsd, jed
  integer   :: n, l, m, i, j, k, buffer_start_pos
  integer   :: shift, gridtype, midpoint, ioff, joff
  integer   :: nunpack, rotation, whalo_v, ehalo_v, shalo_v, nhalo_v

#include "complete_group_update.inc"

end subroutine complete_group_update_r4

subroutine mpp_reset_group_update_field_2d(group, field)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                       intent(in) :: field(:,:)
  integer :: isize, jsize, ksize, data_type
  type(C_PTR)       :: cptr

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_reset_group_update_field_2d: fieldx is not contiguous")

  group%reset_index_s = group%reset_index_s + 1
  if(group%reset_index_s > group%nscalar) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_2d: group%reset_index_s > group%nscalar")
  isize = size(field,1)
  jsize = size(field,2)
  ksize = 1
  data_type = get_field_type(field)
  if(isize .NE. group%isize_s .OR. jsize .NE. group%jsize_s .OR. group%ksize_s .NE. ksize) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_2d: size of field does not match the size stored in group")
  if(data_type .NE. group%data_type)  &
     call mpp_error(FATAL, "mpp_reset_group_update_field_2d: type of field does not match the type stored in group")
  cptr = C_LOC(field)
  select type(field)
  type is (real(DOUBLE_KIND))
    call C_F_POINTER(cptr, group%s_r8(group%reset_index_s,1)%data, (/isize,jsize,ksize/))
  type is (real(FLOAT_KIND))
    call C_F_POINTER(cptr, group%s_r4(group%reset_index_s,1)%data, (/isize,jsize,ksize/))
  class default
    call mpp_error(FATAL,"mpp_reset_group_update_2d: only support 4-byte and 8-byte real, contact developer")
  end select

end subroutine mpp_reset_group_update_field_2d

subroutine mpp_reset_group_update_field_3d(group, field)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                       intent(in) :: field(:,:,:)
  integer :: isize, jsize, ksize, data_type
  type(C_PTR)       :: cptr

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_reset_group_update_field_3d: fieldx is not contiguous")
  group%reset_index_s = group%reset_index_s + 1
  if(group%reset_index_s > group%nscalar) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_3d: group%reset_index_s > group%nscalar")
  isize = size(field,1)
  jsize = size(field,2)
  ksize = size(field,3)
  data_type = get_field_type(field)
  if(isize .NE. group%isize_s .OR. jsize .NE. group%jsize_s .OR. group%ksize_s .NE. ksize) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_3d: size of field does not match the size stored in group")
  if(data_type .NE. group%data_type)  &
     call mpp_error(FATAL, "mpp_reset_group_update_field_3d: type of field does not match the type stored in group")
  cptr = C_LOC(field)
  select type(field)
  type is (real(DOUBLE_KIND))
    call C_F_POINTER(cptr, group%s_r8(group%reset_index_s,1)%data, (/isize,jsize,ksize/))
  type is (real(FLOAT_KIND))
    call C_F_POINTER(cptr, group%s_r4(group%reset_index_s,1)%data, (/isize,jsize,ksize/))
  class default
    call mpp_error(FATAL,"mpp_reset_group_update_3d: only support 4-byte and 8-byte real, contact developer")
  end select

end subroutine mpp_reset_group_update_field_3d

subroutine mpp_reset_group_update_field_4d(group, field)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                       intent(in) :: field(:,:,:,:)
  integer :: isize, jsize, ksize, data_type
  type(C_PTR)       :: cptr
  group%reset_index_s = group%reset_index_s + 1

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_reset_group_update_field_4d: fieldx is not contiguous")
  if(group%reset_index_s > group%nscalar) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_4d: group%reset_index_s > group%nscalar")
  isize = size(field,1)
  jsize = size(field,2)
  ksize = size(field,3)*size(field,4)
  data_type = get_field_type(field)
  if(isize .NE. group%isize_s .OR. jsize .NE. group%jsize_s .OR. group%ksize_s .NE. ksize) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_4d: size of field does not match the size stored in group")
  if(data_type .NE. group%data_type)  &
     call mpp_error(FATAL, "mpp_reset_group_update_field_4d: type of field does not match the type stored in group")
  cptr = C_LOC(field)
  select type(field)
  type is (real(DOUBLE_KIND))
    call C_F_POINTER(cptr, group%s_r8(group%reset_index_s,1)%data, (/isize,jsize,ksize/))
  type is (real(FLOAT_KIND))
    call C_F_POINTER(cptr, group%s_r4(group%reset_index_s,1)%data, (/isize,jsize,ksize/))
  end select

end subroutine mpp_reset_group_update_field_4d

subroutine mpp_reset_group_update_field_2dv(group, fieldx, fieldy)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                       intent(in) :: fieldx(:,:), fieldy(:,:)
  integer :: isizex, jsizex, ksizex, data_typex
  integer :: isizey, jsizey, ksizey, data_typey
  type(C_PTR)       :: cptr

  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_reset_group_update_field_2dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_reset_group_update_field_2dv: fieldy is not contiguous")
  group%reset_index_v = group%reset_index_v + 1
  if(group%reset_index_v > group%nvector) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_2dv: group%reset_index_v > group%nvector")
  isizex = size(fieldx,1)
  jsizex = size(fieldx,2)
  ksizex = 1
  data_typex = get_field_type(fieldx)
  isizey = size(fieldy,1)
  jsizey = size(fieldy,2)
  ksizey = 1
  data_typey = get_field_type(fieldy)
  if(isizex .NE. group%isize_x .OR. jsizex .NE. group%jsize_x .OR. ksizex .NE. group%ksize_v) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_2dv: size of fieldx does not match the size stored in group")
  if(isizey .NE. group%isize_y .OR. jsizey .NE. group%jsize_y .OR. ksizey .NE. group%ksize_v) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_2dv: size of fieldy does not match the size stored in group")
  if(data_typex .NE. group%data_type .or. data_typey .NE. group%data_type)  &
     call mpp_error(FATAL, "mpp_reset_group_update_field_2dv: type of field does not match the type stored in group")
  cptr = C_LOC(fieldx)
  select type(fieldx)
  type is (real(DOUBLE_KIND))
    call C_F_POINTER(cptr, group%s_r8(group%reset_index_v,1)%data, (/isizex,jsizex,ksizex/))
  type is (real(FLOAT_KIND))
    call C_F_POINTER(cptr, group%s_r4(group%reset_index_v,1)%data, (/isizex,jsizex,ksizex/))
  end select
  cptr = C_LOC(fieldy)
  select type(fieldy)
  type is (real(DOUBLE_KIND))
    call C_F_POINTER(cptr, group%s_r8(group%reset_index_v,1)%data, (/isizey,jsizey,ksizey/))
  type is (real(FLOAT_KIND))
    call C_F_POINTER(cptr, group%s_r4(group%reset_index_v,1)%data, (/isizey,jsizey,ksizey/))
  end select

end subroutine mpp_reset_group_update_field_2dv

subroutine mpp_reset_group_update_field_3dv(group, fieldx, fieldy)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                       intent(in) :: fieldx(:,:,:), fieldy(:,:,:)
  integer :: isizex, jsizex, ksizex, data_typex
  integer :: isizey, jsizey, ksizey, data_typey
  type(C_PTR)       :: cptr

  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_reset_group_update_field_3dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_reset_group_update_field_3dv: fieldy is not contiguous")
  group%reset_index_v = group%reset_index_v + 1
  if(group%reset_index_v > group%nvector) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_3dv: group%reset_index_v > group%nvector")
  isizex = size(fieldx,1)
  jsizex = size(fieldx,2)
  ksizex = size(fieldx,3)
  data_typex = get_field_type(fieldx)
  isizey = size(fieldy,1)
  jsizey = size(fieldy,2)
  ksizey = size(fieldy,3)
  data_typey = get_field_type(fieldy)
  if(isizex .NE. group%isize_x .OR. jsizex .NE. group%jsize_x .OR. ksizex .NE. group%ksize_v) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_3dv: size of fieldx does not match the size stored in group")
  if(isizey .NE. group%isize_y .OR. jsizey .NE. group%jsize_y .OR. ksizey .NE. group%ksize_v) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_3dv: size of fieldy does not match the size stored in group")
  if(data_typex .NE. group%data_type .or. data_typey .NE. group%data_type)  &
     call mpp_error(FATAL, "mpp_reset_group_update_field_3dv: type of field does not match the type stored in group")
  cptr = C_LOC(fieldx)
  select type(fieldx)
  type is (real(DOUBLE_KIND))
    call C_F_POINTER(cptr, group%s_r8(group%reset_index_v,1)%data, (/isizex,jsizex,ksizex/))
  type is (real(FLOAT_KIND))
    call C_F_POINTER(cptr, group%s_r4(group%reset_index_v,1)%data, (/isizex,jsizex,ksizex/))
  end select
  cptr = C_LOC(fieldy)
  select type(fieldy)
  type is (real(DOUBLE_KIND))
    call C_F_POINTER(cptr, group%s_r8(group%reset_index_v,1)%data, (/isizey,jsizey,ksizey/))
  type is (real(FLOAT_KIND))
    call C_F_POINTER(cptr, group%s_r4(group%reset_index_v,1)%data, (/isizey,jsizey,ksizey/))
  end select

end subroutine mpp_reset_group_update_field_3dv

subroutine mpp_reset_group_update_field_4dv(group, fieldx, fieldy)
  type(mpp_group_update_type), intent(inout) :: group
  class(*),                       intent(in) :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  integer :: isizex, jsizex, ksizex, data_typex
  integer :: isizey, jsizey, ksizey, data_typey 
  type(C_PTR)       :: cptr

  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_reset_group_update_field_4dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_reset_group_update_field_4dv: fieldy is not contiguous")
  group%reset_index_v = group%reset_index_v + 1
  if(group%reset_index_v > group%nvector) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_4dv: group%reset_index_v > group%nvector")
  isizex = size(fieldx,1)
  jsizex = size(fieldx,2)
  ksizex = size(fieldx,3)*size(fieldx,4)
  data_typex = get_field_type(fieldx)
  isizey = size(fieldy,1)
  jsizey = size(fieldy,2)
  ksizey = size(fieldy,3)*size(fieldy,4)
  data_typey = get_field_type(fieldy)
  if(isizex .NE. group%isize_x .OR. jsizex .NE. group%jsize_x .OR. ksizex .NE. group%ksize_v) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_3dv: size of fieldx does not match the size stored in group")
  if(isizey .NE. group%isize_y .OR. jsizey .NE. group%jsize_y .OR. ksizey .NE. group%ksize_v) &
     call mpp_error(FATAL, "mpp_reset_group_update_field_3dv: size of fieldy does not match the size stored in group")
  if(data_typex .NE. group%data_type .or. data_typey .NE. group%data_type)  &
     call mpp_error(FATAL, "mpp_reset_group_update_field_3dv: type of field does not match the type stored in group")
  cptr = C_LOC(fieldx)
  select type(fieldx)
  type is (real(DOUBLE_KIND))
    call C_F_POINTER(cptr, group%s_r8(group%reset_index_v,1)%data, (/isizex,jsizex,ksizex/))
  type is (real(FLOAT_KIND))
    call C_F_POINTER(cptr, group%s_r4(group%reset_index_v,1)%data, (/isizex,jsizex,ksizex/))
  end select
  cptr = C_LOC(fieldy)
  select type(fieldy)
  type is (real(DOUBLE_KIND))
    call C_F_POINTER(cptr, group%s_r8(group%reset_index_v,1)%data, (/isizey,jsizey,ksizey/))
  type is (real(FLOAT_KIND))
    call C_F_POINTER(cptr, group%s_r4(group%reset_index_v,1)%data, (/isizey,jsizey,ksizey/))
  end select

end subroutine mpp_reset_group_update_field_4dv

!######################################################################
function mpp_start_update_domain_2d( field, domain, flags, position, &
                                     whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  type(domain2D),   intent(inout)        :: domain
  class(*),         intent(inout)        :: field(:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  integer                                :: mpp_start_update_domain_2d

  !--- local variables
  character(len=128)         :: text
  logical                    :: is_complete
  logical, save              :: new_update = .true.

  call check_domain2d_index(domain, "mpp_start_update_domain_2d")

  if(domainList(domain%index)%max_ntile_pe>1) &
     call mpp_error(FATAL,'mpp_start_update_domain_2d: do not support multiple tile per pe, contact developer' )
  if(new_update) then
     num_update = num_update + 1
     current_id_update = current_id_update + 1
     if( current_id_update > MAX_NONBLOCK_UPDATE ) then
        write( text,'(a,i8,a,i8)' ) 'num_fields =', current_id_update, &
              ' greater than MAX_NONBLOCK_UPDATE =', MAX_NONBLOCK_UPDATE
        call mpp_error(FATAL,'mpp_start_update_domain_2d: '//trim(text))
     endif
     new_update = .false.
  endif

  call mpp_create_group_update_2d(nonblock_data(current_id_update), field, domain, flags, position, whalo, ehalo, shalo, nhalo)

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete

  mpp_start_update_domain_2d = 0
  if(is_complete) then
     new_update = .true.
     mpp_start_update_domain_2d = current_id_update
     call mpp_start_group_update(nonblock_data(current_id_update), domain)
  endif

end function mpp_start_update_domain_2d

function mpp_start_update_domain_3d( field, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  type(domain2D),   intent(inout)        :: domain
  class(*),         intent(inout)        :: field(:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  integer                                :: mpp_start_update_domain_3D

  !--- local variables
  character(len=128)         :: text
  logical                    :: is_complete
  logical, save              :: new_update = .true.

  call check_domain2d_index(domain, "mpp_start_update_domain_3d")

  if(domainList(domain%index)%max_ntile_pe>1) &
     call mpp_error(FATAL,'mpp_start_update_domain_3d: do not support multiple tile per pe, contact developer' )
  if(new_update) then
     num_update = num_update + 1
     current_id_update = current_id_update + 1
     if( current_id_update > MAX_NONBLOCK_UPDATE ) then
        write( text,'(a,i8,a,i8)' ) 'num_fields =', current_id_update, &
              ' greater than MAX_NONBLOCK_UPDATE =', MAX_NONBLOCK_UPDATE
        call mpp_error(FATAL,'mpp_start_update_domain_3D: '//trim(text))
     endif
     new_update = .false.
  endif
  call mpp_create_group_update_3d(nonblock_data(current_id_update), field, domain, flags, position, whalo, ehalo, shalo, nhalo)

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete

  mpp_start_update_domain_3d = 0
  if(is_complete) then
     new_update = .true.
     mpp_start_update_domain_3d = current_id_update
     call mpp_start_group_update(nonblock_data(current_id_update), domain)
  endif

end function mpp_start_update_domain_3d

function mpp_start_update_domain_4D( field, domain, flags, position, &
                                       whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  type(domain2D),   intent(inout)        :: domain
  class(*),         intent(inout)        :: field(:,:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  integer                                :: mpp_start_update_domain_4D

  !--- local variables
  character(len=128)         :: text
  logical                    :: is_complete
  logical, save              :: new_update = .true.

  call check_domain2d_index(domain, "mpp_start_update_domain_4D")

  if(domainList(domain%index)%max_ntile_pe>1) &
     call mpp_error(FATAL,'mpp_start_update_domain_3d: do not support multiple tile per pe, contact developer' )
  if(new_update) then
     current_id_update = current_id_update + 1
     num_update = num_update+1
     if( current_id_update > MAX_NONBLOCK_UPDATE ) then
        write( text,'(a,i8,a,i8)' ) 'num_fields =', current_id_update, &
              ' greater than MAX_NONBLOCK_UPDATE =', MAX_NONBLOCK_UPDATE
        call mpp_error(FATAL,'mpp_start_update_domain_4D: '//trim(text))
     endif
     new_update = .false.
  endif

  call mpp_create_group_update_4d(nonblock_data(current_id_update), field, domain, flags, position, whalo, ehalo, shalo, nhalo)

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete

  mpp_start_update_domain_4d = 0
  if(is_complete) then
     new_update = .true.
     mpp_start_update_domain_4d = current_id_update
     call mpp_start_group_update(nonblock_data(current_id_update), domain)
  endif

end function mpp_start_update_domain_4d

!##################################################################################
subroutine mpp_complete_update_domain_2d( id_update, field, domain, flags, position, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  integer,          intent(in)           :: id_update
  type(domain2D),   intent(inout)        :: domain  
  class(*),         intent(inout)        :: field(:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  logical :: is_complete
  integer :: n

  call check_domain2d_index(domain, "mpp_complete_update_domain_2d")

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete
  
  if(is_complete) then
     call mpp_complete_group_update(nonblock_data(id_update), domain)
     call mpp_clear_group_update(nonblock_data(id_update))
     num_update = num_update - 1
     if(num_update == 0) then
        do n = 1, num_update
           call mpp_clear_group_update(nonblock_data(n))
        enddo
        current_id_update = 0
     endif
  endif

end subroutine mpp_complete_update_domain_2d

subroutine mpp_complete_update_domain_3d( id_update, field, domain, flags, position, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  integer,          intent(in)           :: id_update
  type(domain2D),   intent(inout)        :: domain
  class(*),         intent(inout)        :: field(:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  logical :: is_complete
  integer :: n

  call check_domain2d_index(domain, "mpp_complete_update_domain_3d")

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete

  if(is_complete) then
     call mpp_complete_group_update(nonblock_data(id_update), domain)
     call mpp_clear_group_update(nonblock_data(id_update))
     num_update = num_update - 1
     if(num_update == 0) then
        do n = 1, num_update
           call mpp_clear_group_update(nonblock_data(n))
        enddo
        current_id_update = 0
     endif
  endif
end subroutine mpp_complete_update_domain_3d

subroutine mpp_complete_update_domain_4d( id_update, field, domain, flags, position, &
                                            whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  integer,          intent(in)           :: id_update
  type(domain2D),   intent(inout)        :: domain
  class(*),         intent(inout)        :: field(:,:,:,:)
  integer,          intent(in), optional :: flags
  integer,          intent(in), optional :: position
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo ! specify halo region to be updated.
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  logical :: is_complete
  integer :: n

  call check_domain2d_index(domain, "mpp_complete_update_domain_4d")

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete

  if(is_complete) then
     call mpp_complete_group_update(nonblock_data(id_update), domain)
     call mpp_clear_group_update(nonblock_data(id_update))
     num_update = num_update - 1
     if(num_update == 0) then
        do n = 1, num_update
           call mpp_clear_group_update(nonblock_data(n))
        enddo
        current_id_update = 0
     endif
  endif

end subroutine mpp_complete_update_domain_4d

!######################################################################
function mpp_start_update_domain_2dv( fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  !updates data domain of 3D field whose computational domains have been computed
  class(*),         intent(inout)        :: fieldx(:,:), fieldy(:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  integer                                :: mpp_start_update_domain_2dv
  character(len=128)         :: text
  logical                    :: is_complete
  logical, save              :: new_update = .true.

  call check_domain2d_index(domain, "mpp_start_update_domain_2dv")

  if(domainList(domain%index)%max_ntile_pe>1) &
     call mpp_error(FATAL,'mpp_start_update_domain_2dv: do not support multiple tile per pe, contact developer' )
  if(new_update) then
     num_update = num_update + 1
     current_id_update = current_id_update + 1
     if( current_id_update > MAX_NONBLOCK_UPDATE ) then
        write( text,'(a,i8,a,i8)' ) 'num_fields =', current_id_update, &
              ' greater than MAX_NONBLOCK_UPDATE =', MAX_NONBLOCK_UPDATE
        call mpp_error(FATAL,'mpp_start_update_domain_2D_v: '//trim(text))
     endif
     new_update = .false.
  endif
  call mpp_create_group_update_2dv(nonblock_data(current_id_update), fieldx, fieldy, &
                                    domain, flags, gridtype, whalo, ehalo, shalo, nhalo)

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete

  mpp_start_update_domain_2dv = 0
  if(is_complete) then
     new_update = .true.
     mpp_start_update_domain_2dv = current_id_update
     call mpp_start_group_update(nonblock_data(current_id_update), domain)
  endif

  return

end function mpp_start_update_domain_2dv

function mpp_start_update_domain_3dv( fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  !updates data domain of 3D field whose computational domains have been computed
  class(*),         intent(inout)        :: fieldx(:,:,:), fieldy(:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  integer                                :: mpp_start_update_domain_3dv
  character(len=128)         :: text
  logical                    :: is_complete
  logical, save              :: new_update = .true.

  call check_domain2d_index(domain, "mpp_start_update_domain_3dv")

  if(domainList(domain%index)%max_ntile_pe>1) &
     call mpp_error(FATAL,'mpp_start_update_domain_3dv: do not support multiple tile per pe, contact developer' )
  if(new_update) then
     num_update = num_update + 1
     current_id_update = current_id_update + 1
     if( current_id_update > MAX_NONBLOCK_UPDATE ) then
        write( text,'(a,i8,a,i8)' ) 'num_fields =', current_id_update, &
              ' greater than MAX_NONBLOCK_UPDATE =', MAX_NONBLOCK_UPDATE
        call mpp_error(FATAL,'mpp_start_update_domain_3D_v: '//trim(text))
     endif
     new_update = .false.
  endif
  call mpp_create_group_update_3dv(nonblock_data(current_id_update), fieldx, fieldy, &
                                    domain, flags, gridtype, whalo, ehalo, shalo, nhalo)

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete

  mpp_start_update_domain_3dv = 0
  if(is_complete) then
     new_update = .true.
     mpp_start_update_domain_3dv = current_id_update
     call mpp_start_group_update(nonblock_data(current_id_update), domain)
  endif

  return

end function mpp_start_update_domain_3dv

function mpp_start_update_domain_4dv( fieldx, fieldy, domain, flags, gridtype, &
     whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  !updates data domain of 3D field whose computational domains have been computed
  class(*),         intent(inout)        :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  integer                                :: mpp_start_update_domain_4dv
  character(len=128)         :: text
  logical                    :: is_complete
  logical, save              :: new_update = .true.

  call check_domain2d_index(domain, "mpp_start_update_domain_4dv")

  if(domainList(domain%index)%max_ntile_pe>1) &
     call mpp_error(FATAL,'mpp_start_update_domain_4dv: do not support multiple tile per pe, contact developer' )
  if(new_update) then
     num_update = num_update + 1
     current_id_update = current_id_update + 1
     if( current_id_update > MAX_NONBLOCK_UPDATE ) then
        write( text,'(a,i8,a,i8)' ) 'num_fields =', current_id_update, &
              ' greater than MAX_NONBLOCK_UPDATE =', MAX_NONBLOCK_UPDATE
        call mpp_error(FATAL,'mpp_start_update_domain_4D_v: '//trim(text))
     endif
     new_update = .false.
  endif
  call mpp_create_group_update_4dv(nonblock_data(current_id_update), fieldx, fieldy, &
                                    domain, flags, gridtype, whalo, ehalo, shalo, nhalo)

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete

  mpp_start_update_domain_4dv = 0
  if(is_complete) then
     new_update = .true.
     mpp_start_update_domain_4dv = current_id_update
     call mpp_start_group_update(nonblock_data(current_id_update), domain)
  endif

  return

end function mpp_start_update_domain_4dv

!####################################################################################
subroutine mpp_complete_update_domain_2dv( id_update, fieldx, fieldy, domain, flags, gridtype, &
                whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  !updates data domain of 3D field whose computational domains have been computed
  integer,          intent(in)           :: id_update
  class(*),         intent(inout)        :: fieldx(:,:), fieldy(:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  integer :: n
  logical :: is_complete

  call check_domain2d_index(domain, "mpp_complete_update_domain_2dv")

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete

  if(is_complete) then
     current_id_update = current_id_update - 1
     call mpp_complete_group_update(nonblock_data(id_update), domain)
     call mpp_clear_group_update(nonblock_data(id_update))
     num_update = num_update - 1
     if(num_update == 0) then
        do n = 1, num_update
           call mpp_clear_group_update(nonblock_data(n))
        enddo
        current_id_update = 0
     endif
  endif

end subroutine mpp_complete_update_domain_2dv

subroutine mpp_complete_update_domain_3dv( id_update, fieldx, fieldy, domain, flags, gridtype, &
                whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  !updates data domain of 3D field whose computational domains have been computed
  integer,          intent(in)           :: id_update
  class(*),         intent(inout)        :: fieldx(:,:,:), fieldy(:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  integer :: n
  logical :: is_complete

  call check_domain2d_index(domain, "mpp_complete_update_domain_3dv")

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete

  if(is_complete) then
     call mpp_complete_group_update(nonblock_data(id_update), domain)
     call mpp_clear_group_update(nonblock_data(id_update))
     num_update = num_update - 1
     if(num_update == 0) then
        do n = 1, num_update
           call mpp_clear_group_update(nonblock_data(n))
        enddo
        current_id_update = 0
     endif
  endif

end subroutine mpp_complete_update_domain_3dv

subroutine mpp_complete_update_domain_4dv( id_update, fieldx, fieldy, domain, flags, gridtype, &
                whalo, ehalo, shalo, nhalo, name, tile_count, complete )
  !updates data domain of 3D field whose computational domains have been computed
  integer,          intent(in)           :: id_update
  class(*),         intent(inout)        :: fieldx(:,:,:,:), fieldy(:,:,:,:)
  type(domain2D),   intent(inout)        :: domain
  integer,          intent(in), optional :: flags, gridtype
  integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
  character(len=*), intent(in), optional :: name
  integer,          intent(in), optional :: tile_count
  logical,          intent(in), optional :: complete
  integer :: n
  logical :: is_complete

  call check_domain2d_index(domain, "mpp_complete_update_domain_4dv")

  is_complete = .true.
  if(PRESENT(complete)) is_complete = complete

  if(is_complete) then
     call mpp_complete_group_update(nonblock_data(id_update), domain)
     call mpp_clear_group_update(nonblock_data(id_update))
     num_update = num_update - 1
     if(num_update == 0) then
        do n = 1, num_update
           call mpp_clear_group_update(nonblock_data(n))
        enddo
        current_id_update = 0
     endif
  endif

end subroutine mpp_complete_update_domain_4dv

