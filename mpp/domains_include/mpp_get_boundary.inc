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
!###############################################################################################
subroutine get_boundary_core(domain, f_cptr_in, isize_in, jsize_in, ksize_in, data_type_in, &
                                 b_cptr_in, b_size_in, position, complete, tile_count)
  type(domain2D_private),    intent(in) :: domain
  type(C_PTR),       intent(in) :: f_cptr_in, b_cptr_in(:)
  integer,           intent(in) :: isize_in, jsize_in, ksize_in, data_type_in
  integer,           intent(in) :: b_size_in(:) 
  integer, intent(in), optional :: position, tile_count
  logical, intent(in), optional :: complete

  integer                  :: ntile
  logical                  :: need_ebuffer, need_sbuffer, need_wbuffer, need_nbuffer
  type(holder_r8_type), dimension(MAX_DOMAIN_FIELDS, MAX_TILES),save :: f_r8
  type(holder_r4_type), dimension(MAX_DOMAIN_FIELDS, MAX_TILES),save :: f_r4
  type(bound_r8_type), dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: b_r8
  type(bound_r4_type), dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: b_r4
  integer, save    :: b_size(4)=0, isize=0, jsize=0, ksize=0, pos, list=0, l_size=0
  integer, save    :: data_type=0
  integer          :: max_ntile, tile, update_position, ishift, jshift, i
  logical          :: do_update, is_complete, set_mismatch
  logical          :: send(4), recv(4)
  character(len=3) :: text
  type(overlapSpec), pointer :: bound => NULL()

  ntile = size(domain%x(:))

  update_position = CENTER
  if(present(position)) update_position = position

  !--- check if the suitable buffer are present
  need_ebuffer=.false.; need_sbuffer=.false.; need_wbuffer=.false.; need_nbuffer=.false.
  if( domain%symmetry .AND. PRESENT(position) ) then
     select case(position)
     case(CORNER)
       need_ebuffer=.true.; need_sbuffer=.true.; need_wbuffer=.true.; need_nbuffer=.true.
     case(NORTH)
       need_sbuffer=.true.; need_nbuffer=.true.
     case(EAST)
       need_ebuffer=.true.; need_wbuffer=.true.
     end select
  end if
   
  tile = 1
  max_ntile = domain%max_ntile_pe
  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'get_boundary_core: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "get_boundary_core: "// &
          "optional argument tile_count should be present when number of tiles on this pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete        
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'get_boundary_core: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif
  if(b_size_in(1) > 0) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
           'mpp_get_boundary_core: ebuffer should not be present when north is folded')
     if(.not. need_ebuffer) call mpp_error(FATAL,'get_boundary_core: ebuffer should not be present')
  endif
  if(b_size_in(2) > 0) then
     if(.not. need_sbuffer) call mpp_error(FATAL,'get_boundary_core: sbuffer should not be present')
  end if
  if(b_size_in(3) > 0) then
     if(.not. need_wbuffer) call mpp_error(FATAL,'get_boundary_core: wbuffer should not be present')
  end if
  if(b_size_in(4) > 0) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'mpp_get_boundary_core: nbuffer should not be present when north is folded')
     if(.not. need_nbuffer) call mpp_error(FATAL,'get_boundary_core: nbuffer should not be present')
  end if

  if(data_type_in == MPP_DOUBLE) then
     call C_F_POINTER(f_cptr_in, f_r8(list,tile)%data, (/isize_in,jsize_in,ksize_in/))
     do i = 1, 4
        if(b_size_in(i) > 0) then
           call C_F_POINTER(b_cptr_in(i), b_r8(i,list,tile)%data, (/b_size_in(i),ksize_in/))
        endif
     enddo
  else if(data_type_in == MPP_FLOAT) then
     call C_F_POINTER(f_cptr_in, f_r4(list,tile)%data, (/isize_in,jsize_in,ksize_in/))
     do i = 1, 4
        if(b_size_in(i) > 0) &
           call C_F_POINTER(b_cptr_in(i), b_r4(i,list,tile)%data, (/b_size_in(i),ksize_in/))
     enddo
  endif

  if(list == 1 .AND. tile == 1 )then
     isize=isize_in; jsize=jsize_in; ksize = ksize_in; pos = update_position
     b_size = b_size_in; data_type = data_type_in
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize .NE. isize_in)
     set_mismatch = set_mismatch .OR. (jsize .NE. jsize_in)
     set_mismatch = set_mismatch .OR. (ksize .NE. ksize_in)
     set_mismatch = set_mismatch .OR. ANY( b_size .NE. b_size_in )
     set_mismatch = set_mismatch .OR. (update_position .NE. pos)
     set_mismatch = set_mismatch .OR. (data_type .NE. data_type_in)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'get_boundary_core: Incompatible field at count '//text//' for group update.' )
     endif
  endif
  if(is_complete) then
     l_size = list
     list = 0
  end if

  if(do_update )then 
     !--- only non-center data in symmetry domain will be retrieved.
     if(position == CENTER .OR. (.NOT. domain%symmetry) ) return 
     bound => search_bound_overlap(domain, update_position)
     call mpp_get_domain_shift(domain, ishift, jshift, update_position)
     if(isize .NE. domain%x(1)%memory%size+ishift .OR. jsize .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "get_boundary_core: field is not on memory domain")
     if(ASSOCIATED(bound)) then
        if(debug_message_passing) then
           recv(:) = b_size(:) > 0
           send = recv
           call check_message_size(domain, bound, send, recv)
        endif
        if(data_type == MPP_DOUBLE) then
           call do_get_boundary_r8(f_r8(1:l_size,1:ntile), domain, bound, b_r8(:,1:l_size,1:ntile), &
               b_size, ksize)
        else if(data_type == MPP_FLOAT) then
           call do_get_boundary_r4(f_r4(1:l_size,1:ntile), domain, bound, b_r4(:,1:l_size,1:ntile), &
               b_size, ksize)
        else
           call mpp_error(FATAL, "get_oundary_core: only support 4-byte and 8-byte real, contact developer")
        endif
     endif
     l_size=0; b_size=0; isize=0;  jsize=0;  ksize=0
  end if

end subroutine get_boundary_core

!#######################################################################
subroutine mpp_get_boundary_2d(field, domain, ebuffer, sbuffer, wbuffer, nbuffer, &
                                position, complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  class(*),             intent(in)   :: field(:,:)
  class(*),  intent(inout), optional :: ebuffer(:), sbuffer(:), wbuffer(:), nbuffer(:)
  integer,      intent(in), optional :: position, tile_count
  logical,      intent(in), optional :: complete  

  integer     :: isize, jsize, ksize, data_type
  integer     :: b_size(4)
  type(C_PTR) :: f_cptr, b_cptr(4)

  call check_domain2d_index(domain, "mpp_get_boundary_2d")

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_get_boundary_2d: field is not contiguous")
  f_cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = 1
  data_type = get_field_type(field)
  b_size = 0
  if(present(ebuffer)) then
     b_size(1)  = 1
     b_cptr(1) = C_LOC(ebuffer)
  endif
  if(present(sbuffer)) then
     b_size(2)  = 1
     b_cptr(2) = C_LOC(sbuffer)
  endif
  if(present(wbuffer)) then
     b_size(3)  = 1
     b_cptr(3) = C_LOC(wbuffer)
  endif
  if(present(ebuffer)) then
     b_size(4)  = 1
     b_cptr(4) = C_LOC(nbuffer)
  endif
  call get_boundary_core(domainList(domain%index), f_cptr, isize, jsize, ksize, data_type, &
                             b_cptr, b_size, position, complete, tile_count)

end subroutine mpp_get_boundary_2d

!#######################################################################
subroutine mpp_get_boundary_3d(field, domain, ebuffer, sbuffer, wbuffer, nbuffer, &
                               position, complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  class(*),             intent(in)   :: field(:,:,:)
  class(*),  intent(inout), optional :: ebuffer(:,:), sbuffer(:,:), wbuffer(:,:), nbuffer(:,:)
  integer,      intent(in), optional :: position, tile_count
  logical,      intent(in), optional :: complete  

  integer     :: isize, jsize, ksize, data_type
  integer     :: b_size(4)
  type(C_PTR) :: f_cptr, b_cptr(4)

  call check_domain2d_index(domain, "mpp_get_boundary_3d")

  if(.not. is_contiguous(field)) call mpp_error(FATAL, "mpp_get_boundary_3d: field is not contiguous")
  f_cptr = C_LOC(field)
  isize = size(field,1)
  jsize = size(field,2)
  ksize = size(field,3)
  data_type = get_field_type(field)
  b_size = 0
  if(present(ebuffer)) then
     b_size(1)  = size(ebuffer,1)
     b_cptr(1) = C_LOC(ebuffer)
  endif
  if(present(sbuffer)) then
     b_size(2)  = size(sbuffer,1)
     b_cptr(2) = C_LOC(sbuffer)
  endif
  if(present(wbuffer)) then
     b_size(3)  = size(wbuffer,1)
     b_cptr(3) = C_LOC(wbuffer)
  endif
  if(present(ebuffer)) then
     b_size(4)  = size(nbuffer,1)
     b_cptr(4) = C_LOC(nbuffer)
  endif
  call get_boundary_core(domainList(domain%index), f_cptr, isize, jsize, ksize, data_type, &
                             b_cptr, b_size, position, complete, tile_count)

end subroutine mpp_get_boundary_3d

!#######################################################################
subroutine do_get_boundary_r8(f_hold, domain, bound, b_hold, b_size, ke)
  type(holder_r8_type),   intent(in) :: f_hold(:,:)
  type(domain2D_private), intent(in) :: domain
  type(overlapSpec),      intent(in) :: bound
  type(bound_r8_type), intent(inout) :: b_hold(:,:,:)
  integer,                intent(in) :: b_size(:), ke
  real(DOUBLE_KIND) :: buffer(mpp_domains_stack_size)

#include "do_get_boundary.inc"

end subroutine do_get_boundary_r8

!#######################################################################
subroutine do_get_boundary_r4(f_hold, domain, bound, b_hold, b_size, ke)
  type(holder_r4_type),   intent(in) :: f_hold(:,:)
  type(domain2D_private), intent(in) :: domain
  type(overlapSpec),      intent(in) :: bound
  type(bound_r4_type), intent(inout) :: b_hold(:,:,:)
  integer,                intent(in) :: b_size(:), ke
  real(FLOAT_KIND) :: buffer(mpp_domains_stack_size)

#include "do_get_boundary.inc"

end subroutine do_get_boundary_r4

!###############################################################################################
subroutine get_boundaryV_core(domain, f_cptrx_in, f_cptry_in, isizex_in, jsizex_in, ksizex_in, &
                                  isizey_in, jsizey_in,  ksizey_in, data_typex_in, data_typey_in,     &
                                  b_cptrx_in, b_cptry_in, b_sizex_in, b_sizey_in, flags, gridtype, complete, tile_count)
  type(domain2D_private),    intent(in) :: domain
  type(C_PTR),       intent(in) :: f_cptrx_in, f_cptry_in 
  type(C_PTR),       intent(in) :: b_cptrx_in(:), b_cptry_in(:)
  integer,           intent(in) :: isizex_in, jsizex_in, ksizex_in, data_typex_in
  integer,           intent(in) :: isizey_in, jsizey_in, ksizey_in, data_typey_in
  integer,           intent(in) :: b_sizex_in(:), b_sizey_in(:) 
  integer, intent(in), optional :: flags, gridtype, tile_count
  logical, intent(in), optional :: complete

  integer                 :: ntile, update_flags
  logical                 :: need_ebufferx, need_sbufferx, need_wbufferx, need_nbufferx
  logical                 :: need_ebuffery, need_sbuffery, need_wbuffery, need_nbuffery

  type(holder_r8_type), dimension(MAX_DOMAIN_FIELDS, MAX_TILES),save :: fx_r8, fy_r8
  type(holder_r4_type), dimension(MAX_DOMAIN_FIELDS, MAX_TILES),save :: fx_r4, fy_r4
  type(bound_r8_type), dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: bx_r8, by_r8
  type(bound_r4_type), dimension(4,MAX_DOMAIN_FIELDS, MAX_TILES),save :: bx_r4, by_r4
  integer, save    :: b_sizex(4)=0, b_sizey(4)=0, isize(2)=0, jsize(2)=0, ksize=0, l_size=0, list=0
  integer, save    :: offset_type, upflags
  integer          :: max_ntile, tile, grid_offset_type, i
  integer, save    :: data_typex=0, data_typey=0
  logical          :: do_update, is_complete, set_mismatch
  logical          :: sendx(4), recvx(4), sendy(4), recvy(4)
  character(len=3) :: text
  type(overlapSpec), pointer :: boundx=>NULL()
  type(overlapSpec), pointer :: boundy=>NULL()
  integer                     :: position_x, position_y, ishift, jshift

  ntile = size(domain%x(:))
  update_flags = 0
  if( PRESENT(flags) ) then 
     update_flags = flags
  end if

  !--- check if the suitable buffer are present
  need_ebufferx=.FALSE.; need_sbufferx=.FALSE.
  need_wbufferx=.FALSE.; need_nbufferx=.FALSE.
  need_ebuffery=.FALSE.; need_sbuffery=.FALSE.
  need_wbuffery=.FALSE.; need_nbuffery=.FALSE.
  if( domain%symmetry .AND. PRESENT(gridtype) ) then
     select case(gridtype)
     case(BGRID_NE, BGRID_SW)
       need_ebufferx=.true.; need_sbufferx=.true.; need_wbufferx=.true.; need_nbufferx=.true.
       need_ebuffery=.true.; need_sbuffery=.true.; need_wbuffery=.true.; need_nbuffery=.true.
     case(CGRID_NE, CGRID_SW)
       need_ebufferx=.true.; need_wbufferx=.true.; need_sbuffery=.true.; need_nbuffery=.true.
     case(DGRID_NE, DGRID_SW)
       need_ebuffery=.true.; need_wbuffery=.true.; need_sbufferx=.true.; need_nbufferx=.true.
     end select
  end if

  tile = 1
  max_ntile = domain%max_ntile_pe
  is_complete = .true.
  if(PRESENT(complete)) then
     is_complete = complete
  end if

  if(max_ntile>1) then
     if(ntile>MAX_TILES) then
        write( text,'(i2)' ) MAX_TILES
        call mpp_error(FATAL,'mpp_get_boundaryV_core: MAX_TILES='//text//' is less than number of tiles on this pe.' )
     endif
     if(.NOT. present(tile_count) ) call mpp_error(FATAL, "get_boundaryV_core: "// &
          "optional argument tile_count should be present when number of tiles on this pe is more than 1")
     tile = tile_count
  end if

  do_update = (tile == ntile) .AND. is_complete        
  list = list+1
  if(list > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'get_boundaryV_core: MAX_DOMAIN_FIELDS='//text//' exceeded for group update.' )
  endif

  if(b_sizex_in(1) >0) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'get_boundaryV_core: ebufferx should not be present when north is folded')
     if(.not. need_ebufferx) call mpp_error(FATAL,'get_boundaryV_core: ebufferx should not be present')
  end if
  if(b_sizex_in(2) >0) then
     if(.not. need_sbufferx) call mpp_error(FATAL,'get_boundaryV_core: sbufferx should not be present')
  end if
  if(b_sizex_in(3) >0) then
     if(.not. need_wbufferx) call mpp_error(FATAL,'get_boundaryV_core: wbufferx should not be present')
  end if
  if(b_sizex_in(4) >0) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'get_boundaryV_core: nbufferx should not be present when north is folded')
     if(.not. need_nbufferx) call mpp_error(FATAL,'get_boundaryV_core: nbufferx should not be present')
  end if

  if(b_sizey_in(1) >0) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'get_boundaryV_core: ebuffery should not be present when north is folded')
     if(.not. need_ebuffery) call mpp_error(FATAL,'get_boundaryV_core: ebuffery should not be present')
  end if
  if(b_sizey_in(2) >0) then
     if(.not. need_sbuffery) call mpp_error(FATAL,'get_boundaryV_core: sbuffery should not be present')
  end if
  if(b_sizey_in(3) >0) then
     if(.not. need_wbuffery) call mpp_error(FATAL,'get_boundaryV_core: wbuffery should not be present')
  end if
  if(b_sizey_in(4) >0) then
     if(BTEST(domain%fold,NORTH)) call mpp_error(FATAL, &
             'get_boundaryV_core: nbuffery should not be present when north is folded')
     if(.not. need_nbuffery) call mpp_error(FATAL,'get_boundaryV_core: nbuffery should not be present')
  end if

  if(data_typex_in == MPP_DOUBLE) then
     call C_F_POINTER(f_cptrx_in, fx_r8(list,tile)%data, (/isizex_in,jsizex_in,ksizex_in/))
     call C_F_POINTER(f_cptry_in, fy_r8(list,tile)%data, (/isizey_in,jsizey_in,ksizey_in/))
     do i = 1, 4
        if(b_sizex_in(i) > 0) &
           call C_F_POINTER(b_cptrx_in(i), bx_r8(i,list,tile)%data, (/b_sizex_in(i),ksizex_in/))
        if(b_sizey_in(i) > 0) &
           call C_F_POINTER(b_cptry_in(i), by_r8(i,list,tile)%data, (/b_sizey_in(i),ksizey_in/))
     enddo
  else if(data_typex_in == MPP_FLOAT) then
     call C_F_POINTER(f_cptrx_in, fx_r4(list,tile)%data, (/isizex_in,jsizex_in,ksizex_in/))
     call C_F_POINTER(f_cptry_in, fy_r4(list,tile)%data, (/isizey_in,jsizey_in,ksizey_in/))
     do i = 1, 4
        if(b_sizex_in(i) > 0) &
           call C_F_POINTER(b_cptrx_in(i), bx_r4(i,list,tile)%data, (/b_sizex_in(i),ksizex_in/))
        if(b_sizey_in(i) > 0) &
           call C_F_POINTER(b_cptry_in(i), by_r4(i,list,tile)%data, (/b_sizey_in(i),ksizey_in/))
     enddo
  endif

  grid_offset_type = AGRID
  if(present(gridtype)) grid_offset_type = gridtype
  if(list == 1 .AND. tile == 1 )then
     if(ksizex_in .NE. ksizey_in) call mpp_error(FATAL,  &
          'get_boundaryV_core: mismatch of ksize between fieldx and fieldy')
     if(data_typex_in .NE. data_typey_in) call mpp_error(FATAL,  &
          'get_boundaryV_core: mismatch of data type between fieldx and fieldy')
     isize(1)=isizex_in; jsize(1)=jsizex_in; isize(2)=isizey_in; jsize(2)=jsizey_in
     ksize = ksizex_in; offset_type = grid_offset_type
     data_typex = data_typex_in; data_typey = data_typey_in
     b_sizex = b_sizex_in; b_sizey = b_sizey_in; upflags = update_flags
  else
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (isize(1) .NE. isizex_in)
     set_mismatch = set_mismatch .OR. (jsize(1) .NE. jsizex_in)
     set_mismatch = set_mismatch .OR. (ksize    .NE. ksizex_in)
     set_mismatch = set_mismatch .OR. (isize(2) .NE. isizey_in)
     set_mismatch = set_mismatch .OR. (jsize(2) .NE. jsizey_in)
     set_mismatch = set_mismatch .OR. (ksize    .NE. ksizey_in)
     set_mismatch = set_mismatch .OR. ANY( b_sizex .NE. b_sizex_in )
     set_mismatch = set_mismatch .OR. ANY( b_sizey .NE. b_sizey_in )
     set_mismatch = set_mismatch .OR. (data_typex .NE. data_typex_in )
     set_mismatch = set_mismatch .OR. (data_typey .NE. data_typey_in )
     set_mismatch = set_mismatch .OR. (offset_type .NE. grid_offset_type)
     set_mismatch = set_mismatch .OR. (upflags .NE. update_flags)
     if(set_mismatch)then
        write( text,'(i2)' ) list
        call mpp_error(FATAL,'get_boundaryV_core: Incompatible field at count '//text//' for group update.' )
     endif
  endif
  if(is_complete) then
     l_size = list
     list = 0
  end if

  if(do_update )then
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
     case (DGRID_NE, DGRID_SW)
        position_x = NORTH
        position_y = EAST
     case default
        call mpp_error(FATAL, "get_boundaryV_core: invalid value of grid_offset_type")
     end select

     boundx => search_bound_overlap(domain, position_x)
     boundy => search_bound_overlap(domain, position_y)  

     call mpp_get_domain_shift(domain, ishift, jshift, position_x)
     if(isize(1) .NE. domain%x(1)%memory%size+ishift .OR. jsize(1) .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "get_boundaryV_core: fieldx is not on memory domain")
     call mpp_get_domain_shift(domain, ishift, jshift, position_y)
     if(isize(2) .NE. domain%x(1)%memory%size+ishift .OR. jsize(2) .NE. domain%y(1)%memory%size+jshift ) &
          call mpp_error(FATAL, "get_boundaryV_core: fieldy is not on memory domain")
     if(ASSOCIATED(boundx) ) then
        if(debug_message_passing) then
           recvx(:) = b_sizex(:) > 0
           recvy(:) = b_sizey(:) > 0
           sendx = recvx
           sendy = recvy
           call check_message_size_vector(domain, boundx, boundy, sendx, recvx, sendy, recvy)
        endif
        if(data_typex == MPP_DOUBLE) then
           call do_get_boundary_r8v(fx_r8(1:l_size,1:ntile), fy_r8(1:l_size,1:ntile), domain, boundx, boundy, &
                bx_r8(:,1:l_size,1:ntile), by_r8(:,1:l_size,1:ntile), b_sizex, &
                b_sizey, ksize, update_flags, grid_offset_type)
        else if(data_typex == MPP_FLOAT) then
           call do_get_boundary_r4v(fx_r4(1:l_size,1:ntile), fy_r4(1:l_size,1:ntile), domain, boundx, boundy, &
                bx_r4(:,1:l_size,1:ntile), by_r4(:,1:l_size,1:ntile), b_sizex, &
                b_sizey, ksize, update_flags, grid_offset_type)
        else
           call mpp_error(FATAL, "get_boundaryV_core: only support 4-byte and 8-byte real, contact developer")
        endif
     endif
     l_size=0; isize=0;  jsize=0;  ksize=0; b_sizex=0; b_sizey=0
  end if

end subroutine  get_boundaryV_core

!#######################################################################
subroutine mpp_get_boundary_2dv(fieldx, fieldy, domain, ebufferx, sbufferx, wbufferx, nbufferx, &
                                ebuffery, sbuffery, wbuffery, nbuffery, flags, gridtype, &
                                complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  class(*),             intent(in)   :: fieldx(:,:), fieldy(:,:)
  class(*),  intent(inout), optional :: ebufferx(:), sbufferx(:), wbufferx(:), nbufferx(:)
  class(*),  intent(inout), optional :: ebuffery(:), sbuffery(:), wbuffery(:), nbuffery(:)
  integer,      intent(in), optional :: flags, gridtype, tile_count
  logical,      intent(in), optional :: complete  
  integer     :: isizex,jsizex,isizey,jsizey,ksizex,ksizey,data_typex,data_typey
  integer     :: b_sizex(4), b_sizey(4)
  type(C_PTR) :: f_cptrx, f_cptry, b_cptrx(4), b_cptry(4)

  call check_domain2d_index(domain, "mpp_get_boundary_2dv")

  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_get_boundary_2dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_get_boundary_2dv: fieldy is not contiguous")
  f_cptrx = C_LOC(fieldx)
  f_cptry = C_LOC(fieldy)
  isizex = size(fieldx,1)
  jsizex = size(fieldx,2)
  isizey = size(fieldy,1)
  jsizey = size(fieldy,2)
  ksizex = 1
  ksizey = 1
  data_typex = get_field_type(fieldx)
  data_typey = get_field_type(fieldy)
  b_sizex = 0
  b_sizey = 0
  if(present(ebufferx)) then
     b_sizex(1)  = size(ebufferx(:))
     b_cptrx(1) = C_LOC(ebufferx)
  endif
  if(present(sbufferx)) then
     b_sizex(2)  = size(sbufferx(:))
     b_cptrx(2) = C_LOC(sbufferx)
  endif
  if(present(wbufferx)) then
     b_sizex(3)  = size(wbufferx(:))
     b_cptrx(3) = C_LOC(wbufferx)
  endif
  if(present(nbufferx)) then
     b_sizex(4)  = size(nbufferx(:))
     b_cptrx(4) = C_LOC(nbufferx)
  endif
  if(present(ebuffery)) then
     b_sizey(1)  = size(ebuffery(:))
     b_cptry(1) = C_LOC(ebuffery)
  endif
  if(present(sbuffery)) then
     b_sizey(2)  = size(sbuffery(:))
     b_cptry(2) = C_LOC(sbuffery)
  endif
  if(present(wbuffery)) then
     b_sizey(3)  = size(wbuffery(:))
     b_cptry(3) = C_LOC(wbuffery)
  endif
  if(present(nbuffery)) then
     b_sizey(4)  = size(nbuffery(:))
     b_cptry(4) = C_LOC(nbuffery)
  endif  

  call get_boundaryV_core(domainList(domain%index), f_cptrx, f_cptry, isizex, jsizex, ksizex, isizey, jsizey,  ksizey, &
                              data_typex, data_typey, b_cptrx, b_cptry, b_sizex, b_sizey, flags, gridtype, complete, tile_count)

  return

end subroutine mpp_get_boundary_2dv

!#######################################################################
subroutine mpp_get_boundary_3dv(fieldx, fieldy, domain, ebufferx, sbufferx, wbufferx, nbufferx, &
                                ebuffery, sbuffery, wbuffery, nbuffery, flags, gridtype, &
                                complete, tile_count)
  type(domain2D),       intent(in)   :: domain
  class(*),             intent(in)   :: fieldx(:,:,:), fieldy(:,:,:)
  class(*),  intent(inout), optional :: ebufferx(:,:), sbufferx(:,:), wbufferx(:,:), nbufferx(:,:)
  class(*),  intent(inout), optional :: ebuffery(:,:), sbuffery(:,:), wbuffery(:,:), nbuffery(:,:)
  integer,      intent(in), optional :: flags, gridtype, tile_count
  logical,      intent(in), optional :: complete  
  integer     :: isizex,jsizex,isizey,jsizey,ksizex,ksizey,data_typex,data_typey
  integer     :: b_sizex(4), b_sizey(4)
  type(C_PTR) :: f_cptrx, f_cptry, b_cptrx(4), b_cptry(4)

  call check_domain2d_index(domain, "mpp_get_boundary_3dv")

  if(.not. is_contiguous(fieldx)) call mpp_error(FATAL, "mpp_get_boundary_3dv: fieldx is not contiguous")
  if(.not. is_contiguous(fieldy)) call mpp_error(FATAL, "mpp_get_boundary_3dv: fieldy is not contiguous")
  f_cptrx = C_LOC(fieldx)
  f_cptry = C_LOC(fieldy)
  isizex = size(fieldx,1)
  jsizex = size(fieldx,2)
  isizey = size(fieldy,1)
  jsizey = size(fieldy,2)
  ksizex = size(fieldx,3)
  ksizey = size(fieldy,3)
  data_typex = get_field_type(fieldx)
  data_typey = get_field_type(fieldy)
  b_sizex = 0
  b_sizey = 0
  if(present(ebufferx)) then
     b_sizex(1)  = size(ebufferx,1)
     b_cptrx(1) = C_LOC(ebufferx)
  endif
  if(present(sbufferx)) then
     b_sizex(2)  = size(sbufferx,1)
     b_cptrx(2) = C_LOC(sbufferx)
  endif
  if(present(wbufferx)) then
     b_sizex(3)  = size(wbufferx,1)
     b_cptrx(3) = C_LOC(wbufferx)
  endif
  if(present(nbufferx)) then
     b_sizex(4)  = size(nbufferx,1)
     b_cptrx(4) = C_LOC(nbufferx)
  endif
  if(present(ebuffery)) then
     b_sizey(1)  = size(ebuffery,1)
     b_cptry(1) = C_LOC(ebuffery)
  endif
  if(present(sbuffery)) then
     b_sizey(2)  = size(sbuffery,1)
     b_cptry(2) = C_LOC(sbuffery)
  endif
  if(present(wbuffery)) then
     b_sizey(3)  = size(wbuffery,1)
     b_cptry(3) = C_LOC(wbuffery)
  endif
  if(present(nbuffery)) then
     b_sizey(4)  = size(nbuffery,1)
     b_cptry(4) = C_LOC(nbuffery)
  endif  

  call get_boundaryV_core(domainList(domain%index), f_cptrx, f_cptry, isizex, jsizex, ksizex, isizey, jsizey, &
           ksizey, data_typex, data_typey, b_cptrx, b_cptry, b_sizex, b_sizey, flags, gridtype, complete, tile_count)

  return

end subroutine mpp_get_boundary_3dv

!#######################################################################
subroutine do_get_boundary_r8v(fx_hold, fy_hold, domain, boundx, boundy, bx_hold, by_hold,  &
                                   b_sizex, b_sizey, ke, flags, gridtype)
  type(holder_r8_type),   intent(in) :: fx_hold(:,:), fy_hold(:,:)
  type(domain2D_private), intent(in) :: domain
  type(overlapSpec),      intent(in) :: boundx, boundy
  type(bound_r8_type), intent(inout) :: bx_hold(:,:,:), by_hold(:,:,:)
  integer,                intent(in) :: b_sizex(:), b_sizey(:), ke
  integer, intent(in)                :: flags
  integer, intent(in)                :: gridtype
  real(DOUBLE_KIND) :: buffer(mpp_domains_stack_size)

#include "do_get_boundaryV.inc"

end subroutine do_get_boundary_r8v

!#######################################################################
subroutine do_get_boundary_r4v(fx_hold, fy_hold, domain, boundx, boundy, bx_hold, by_hold,  &
                                   b_sizex, b_sizey, ke, flags, gridtype)
  type(holder_r4_type),   intent(in) :: fx_hold(:,:), fy_hold(:,:)
  type(domain2D_private), intent(in) :: domain
  type(overlapSpec),      intent(in) :: boundx, boundy
  type(bound_r4_type), intent(inout) :: bx_hold(:,:,:), by_hold(:,:,:)
  integer,                intent(in) :: b_sizex(:), b_sizey(:), ke
  integer,                intent(in) :: flags
  integer,                intent(in) :: gridtype
  real(FLOAT_KIND) :: buffer(mpp_domains_stack_size)

#include "do_get_boundaryV.inc"

end subroutine do_get_boundary_r4v





