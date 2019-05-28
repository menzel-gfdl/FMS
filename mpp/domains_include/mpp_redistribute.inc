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
subroutine redistribute_core( domain_in, cptr_in, f_loc_in, isize_in, jsize_in, ke_in, data_type_in, &
      domain_out, cptr_out, f_loc_out, isize_out, jsize_out, ke_out, data_type_out, complete, dc_handle, position )
  type(domain2D_private), intent(in) :: domain_in, domain_out
  type(C_PTR),    intent(in) :: cptr_in, cptr_out
  integer(LONG_KIND), intent(in) :: f_loc_in, f_loc_out
  integer,        intent(in) :: isize_in, jsize_in, ke_in, data_type_in
  integer,        intent(in) :: isize_out, jsize_out, ke_out, data_type_out
  logical, intent(in), optional :: complete
  integer, intent(in), optional :: position
  type(DomainCommunicator2D),pointer,optional :: dc_handle
  type(DomainCommunicator2D),pointer,save :: d_comm =>NULL()
  logical                       :: do_redist
  type(holder_r8_type), dimension(MAX_DOMAIN_FIELDS),save :: s_r8_in, s_r8_out
  type(holder_r4_type), dimension(MAX_DOMAIN_FIELDS),save :: s_r4_in, s_r4_out
  integer(LONG_KIND), dimension(MAX_DOMAIN_FIELDS),save :: l_loc_in, l_loc_out
  integer, save :: l_isize_in=0, l_jsize_in=0, l_ke_in=0,l_size=0
  integer, save :: l_isize_out=0,l_jsize_out=0,l_ke_out=0
  integer, save :: l_data_type_in=0, l_data_type_out=0;
  logical       :: set_mismatch
  character(len=2) :: text
  if(present(position)) then
     if(position .NE. CENTER) call mpp_error( FATAL,  &
          'MPP_REDISTRIBUTE_3Dold_: only position = CENTER is implemented, contact author')
  endif

  do_redist=.true.; if(PRESENT(complete))do_redist=complete
  l_size = l_size+1
  if(l_size > MAX_DOMAIN_FIELDS)then
     write( text,'(i2)' ) MAX_DOMAIN_FIELDS
     call mpp_error(FATAL,'mpp_redistribute_core: MAX_DOMAIN_FIELDS='//text//' exceeded for group redistribute.' )
  end if
  if(data_type_in == MPP_DOUBLE) then
    if(domain_in%initialized) &
       call C_F_POINTER(cptr_in, s_r8_in(l_size)%data, (/isize_in,jsize_in,ke_in/))
    if(domain_out%initialized) &
       call C_F_POINTER(cptr_out, s_r8_out(l_size)%data, (/isize_out,jsize_out,ke_out/))
  else if(data_type_in == MPP_FLOAT) then
    if(domain_in%initialized) &
       call C_F_POINTER(cptr_in, s_r4_in(l_size)%data, (/isize_in,jsize_in,ke_in/))
    if(domain_out%initialized) &
       call C_F_POINTER(cptr_out, s_r4_out(l_size)%data, (/isize_out,jsize_out,ke_out/))
  endif

  l_loc_in(l_size) = f_loc_in; l_loc_out(l_size) = f_loc_out
  if(l_size == 1)then
     if(data_type_in .NE. data_type_out) then
        call mpp_error(FATAL,'mpp_redistribute_core: mismatch of data_type between fieldx and fieldy')
     endif
     l_data_type_in = data_type_in; l_data_type_out = data_type_out
     if(domain_in%initialized)then
        l_isize_in=isize_in; l_jsize_in=jsize_in; l_ke_in = ke_in
     end if
     if(domain_out%initialized)then
        l_isize_out=isize_out; l_jsize_out=jsize_out; l_ke_out = ke_out
     endif
  else   
     set_mismatch = .false.
     set_mismatch = set_mismatch .OR. (l_data_type_in .NE. data_type_in)
     set_mismatch = set_mismatch .OR. (l_data_type_out .NE. data_type_out)
     if(domain_in%initialized)then
        set_mismatch = set_mismatch .OR. (l_isize_in /= isize_in)
        set_mismatch = set_mismatch .OR. (l_jsize_in /= jsize_in)
        set_mismatch = set_mismatch .OR. (l_ke_in /= ke_in)
     endif
     if(domain_out%initialized)then
        set_mismatch = set_mismatch .OR. (l_isize_out /= isize_out)
        set_mismatch = set_mismatch .OR. (l_jsize_out /= jsize_out)
        set_mismatch = set_mismatch .OR. (l_ke_out /= ke_out)
     endif
     if(set_mismatch)then
        write( text,'(i2)' ) l_size
        call mpp_error(FATAL,'mpp_redistribute_core: Incompatible field at count '//text//' for group redistribute.' )
     endif
  endif
  if(do_redist)then
     if(PRESENT(dc_handle))d_comm =>dc_handle  ! User has kept pointer to d_comm
     if(.not.ASSOCIATED(d_comm))then  ! d_comm needs initialization or lookup
        d_comm =>redistribute_init_comm(domain_in,l_loc_in(1:l_size),domain_out,l_loc_out(1:l_size), &
             isize_in,jsize_in,ke_in,isize_out,jsize_out,ke_out)
        if(PRESENT(dc_handle))dc_handle =>d_comm  ! User wants to keep pointer to d_comm
     endif
     if(data_type_in == MPP_DOUBLE) then 
       call mpp_do_redistribute_r8( s_r8_in(1:l_size), s_r8_out(1:l_size), d_comm)
     else if(data_type_in == MPP_FLOAT) then
       call mpp_do_redistribute_r4( s_r4_in(1:l_size), s_r4_out(1:l_size), d_comm)
     else
       call mpp_error(FATAL,"mpp_redistribute_core: only support 4-byte and 8-byte real, contact developer")
     endif
     l_size=0
     l_isize_in=0;  l_jsize_in=0;  l_ke_in=0
     l_isize_out=0; l_jsize_out=0; l_ke_out=0
     d_comm =>NULL()
  endif

end subroutine redistribute_core

subroutine mpp_redistribute_2d( domainObj_in, field_in, domainObj_out, field_out, complete, dc_handle, position )
   type(domain2D), intent(in) :: domainObj_in, domainObj_out
   class(*), intent(in)  :: field_in (:,:)
   class(*), intent(out) :: field_out(:,:)
   logical, intent(in), optional :: complete
   integer, intent(in), optional :: position
   type(DomainCommunicator2D),pointer,optional :: dc_handle
   type(domain2d_private), pointer :: domain_in=>NULL()
   type(domain2d_private), pointer :: domain_out=>NULL()

   call check_domain2d_index(domainObj_in, "mpp_redistribute_2d domainObj_in")
   call check_domain2d_index(domainObj_out, "mpp_redistribute_2d domainObj_out")
   domain_in => domainList(domainObj_in%index)
   domain_out => domainList(domainObj_out%index)

   call redistribute_2d( domain_in, field_in, domain_out, field_out, complete, dc_handle, position )

end subroutine mpp_redistribute_2d

subroutine redistribute_2d( domain_in, field_in, domain_out, field_out, complete, dc_handle, position )
   type(domain2D_private), intent(in) :: domain_in, domain_out
   class(*), intent(in)  :: field_in (:,:)
   class(*), intent(out) :: field_out(:,:)
   logical, intent(in), optional :: complete
   integer, intent(in), optional :: position
   type(DomainCommunicator2D),pointer,optional :: dc_handle
   integer     :: isize_in, jsize_in, ksize_in, data_type_in
   integer     :: isize_out, jsize_out, ksize_out, data_type_out
   type(C_PTR) :: cptr_in, cptr_out
   integer(LONG_KIND) :: f_loc_in, f_loc_out

   data_type_in = get_field_type(field_in)
   data_type_out = get_field_type(field_out)
   isize_in = 0; jsize_in = 0; ksize_in = 1
   isize_out = 0; jsize_out = 0; ksize_out = 1
   cptr_in = c_null_ptr; cptr_out = c_null_ptr
   if(domain_in%initialized) then
      if(.not. is_contiguous(field_in)) call mpp_error(FATAL, "mpp_redistribute_2d: field_in is not contiguous")
      cptr_in = C_LOC(field_in)
      isize_in = size(field_in,1)
      jsize_in = size(field_in,2)
      ksize_in = 1
   endif
   if(domain_out%initialized) then
      if(.not. is_contiguous(field_out)) call mpp_error(FATAL, "mpp_redistribute_2d: field_out is not contiguous")
      cptr_out = C_LOC(field_out)
      isize_out = size(field_out,1)
      jsize_out = size(field_out,2)
      ksize_out = 1
   endif
   
   f_loc_in = LOC(field_in)
   f_loc_out = LOC(field_out)
   call redistribute_core( domain_in, cptr_in, f_loc_in, isize_in, jsize_in, ksize_in, data_type_in, &
        domain_out, cptr_out, f_loc_out, isize_out, jsize_out, ksize_out, data_type_out, complete, dc_handle, position )
   return
end subroutine redistribute_2d

!######################################################################
subroutine mpp_redistribute_3d( domainObj_in, field_in, domainObj_out, field_out, complete, dc_handle, position )
   type(domain2D), intent(in) :: domainObj_in, domainObj_out
   class(*), intent(in)  :: field_in (:,:,:)
   class(*), intent(out) :: field_out(:,:,:)
   logical, intent(in), optional :: complete
   integer, intent(in), optional :: position
   type(DomainCommunicator2D),pointer,optional :: dc_handle
   type(domain2d_private), pointer :: domain_in=>NULL()
   type(domain2d_private), pointer :: domain_out=>NULL()

   call check_domain2d_index(domainObj_in, "mpp_redistribute_3d domainObj_in")
   call check_domain2d_index(domainObj_out, "mpp_redistribute_3d domainObj_out")

   domain_in => domainList(domainObj_in%index)
   domain_out => domainList(domainObj_out%index)

   call redistribute_3d( domain_in, field_in, domain_out, field_out, complete, dc_handle, position )

end subroutine mpp_redistribute_3d

subroutine redistribute_3d( domain_in, field_in, domain_out, field_out, complete, dc_handle, position )
   type(domain2D_private), intent(in) :: domain_in, domain_out
   class(*), intent(in)  :: field_in (:,:,:)
   class(*), intent(out) :: field_out(:,:,:)
   logical, intent(in), optional :: complete
   integer, intent(in), optional :: position
   type(DomainCommunicator2D),pointer,optional :: dc_handle
   integer     :: isize_in, jsize_in, ksize_in, data_type_in
   integer     :: isize_out, jsize_out, ksize_out, data_type_out
   type(C_PTR) :: cptr_in, cptr_out
   integer(LONG_KIND) :: f_loc_in, f_loc_out

   data_type_in = get_field_type(field_in)
   data_type_out = get_field_type(field_out)
   isize_in = 0; jsize_in = 0; ksize_in = 1
   isize_out = 0; jsize_out = 0; ksize_out = 1
   cptr_in = c_null_ptr; cptr_out = c_null_ptr
   if(domain_in%initialized) then
      if(.not. is_contiguous(field_in)) call mpp_error(FATAL, "mpp_redistribute_3d: field_in is not contiguous")
      cptr_in = C_LOC(field_in)
      isize_in = size(field_in,1)
      jsize_in = size(field_in,2)
      ksize_in = size(field_in,3)
   endif
   if(domain_out%initialized) then
      if(.not. is_contiguous(field_out)) call mpp_error(FATAL, "mpp_redistribute_3d: field_out is not contiguous")
      cptr_out = C_LOC(field_out)
      isize_out = size(field_out,1)
      jsize_out = size(field_out,2)
      ksize_out = size(field_out,3)
   endif

   f_loc_in = LOC(field_in)
   f_loc_out = LOC(field_out)
   call redistribute_core( domain_in, cptr_in, f_loc_in, isize_in, jsize_in, ksize_in, data_type_in, &
        domain_out, cptr_out, f_loc_out, isize_out, jsize_out, ksize_out, data_type_out, complete, dc_handle, position )

   return
end subroutine redistribute_3d

!######################################################################
subroutine mpp_redistribute_4d( domainObj_in, field_in, domainObj_out, field_out, complete, dc_handle, position )
   type(domain2D), intent(in) :: domainObj_in, domainObj_out
   class(*), intent(in)  :: field_in (:,:,:,:)
   class(*), intent(out) :: field_out(:,:,:,:)
   logical, intent(in), optional :: complete
   integer, intent(in), optional :: position
   type(DomainCommunicator2D),pointer,optional :: dc_handle
   type(domain2d_private), pointer :: domain_in=>NULL()
   type(domain2d_private), pointer :: domain_out=>NULL()

   call check_domain2d_index(domainObj_in, "mpp_redistribute_4d domainObj_in")
   call check_domain2d_index(domainObj_out, "mpp_redistribute_4d domainObj_out")

   domain_in => domainList(domainObj_in%index)
   domain_out => domainList(domainObj_out%index)

   call redistribute_4d( domain_in, field_in, domain_out, field_out, complete, dc_handle, position )

end subroutine mpp_redistribute_4d

subroutine redistribute_4d( domain_in, field_in, domain_out, field_out, complete, dc_handle, position )
   type(domain2D_private), intent(in) :: domain_in, domain_out
   class(*), intent(in)  :: field_in (:,:,:,:)
   class(*), intent(out) :: field_out(:,:,:,:)
   logical, intent(in), optional :: complete
   integer, intent(in), optional :: position 
   type(DomainCommunicator2D),pointer,optional :: dc_handle
   integer     :: isize_in, jsize_in, ksize_in, data_type_in
   integer     :: isize_out, jsize_out, ksize_out, data_type_out
   type(C_PTR) :: cptr_in, cptr_out
   integer(LONG_KIND) :: f_loc_in, f_loc_out

   data_type_in = get_field_type(field_in)
   data_type_out = get_field_type(field_out)
   isize_in = 0; jsize_in = 0; ksize_in = 1
   isize_out = 0; jsize_out = 0; ksize_out = 1
   cptr_in = c_null_ptr; cptr_out = c_null_ptr
   if(domain_in%initialized) then
      if(.not. is_contiguous(field_in)) call mpp_error(FATAL, "mpp_redistribute_4d: field_in is not contiguous")
      cptr_in = C_LOC(field_in)
      isize_in = size(field_in,1)
      jsize_in = size(field_in,2)
      ksize_in = size(field_in,3)*size(field_in,4)
   endif
   if(domain_out%initialized) then
      if(.not. is_contiguous(field_out)) call mpp_error(FATAL, "mpp_redistribute_4d: field_out is not contiguous")
      cptr_out = C_LOC(field_out)
      isize_out = size(field_out,1)
      jsize_out = size(field_out,2)
      ksize_out = size(field_out,3)*size(field_out,4)
   endif

   f_loc_in = LOC(field_in)
   f_loc_out = LOC(field_out)
   call redistribute_core( domain_in, cptr_in, f_loc_in, isize_in, jsize_in, ksize_in, data_type_in, &
        domain_out, cptr_out, f_loc_out, isize_out, jsize_out, ksize_out, data_type_out, complete, dc_handle, position )

   return
end subroutine redistribute_4d

!######################################################################
subroutine mpp_redistribute_5d( domainObj_in, field_in, domainObj_out, field_out, complete, dc_handle, position )
   type(domain2D),    intent(in) :: domainObj_in, domainObj_out
   class(*),         intent(in)  :: field_in (:,:,:,:,:)
   class(*),         intent(out) :: field_out(:,:,:,:,:)
   logical, intent(in), optional :: complete
   integer, intent(in), optional :: position
   type(DomainCommunicator2D),pointer,optional :: dc_handle
   type(domain2d_private), pointer :: domain_in=>NULL()
   type(domain2d_private), pointer :: domain_out=>NULL()

   call check_domain2d_index(domainObj_in, "mpp_redistribute_5d domainObj_in")
   call check_domain2d_index(domainObj_out, "mpp_redistribute_5d domainObj_out")

   domain_in => domainList(domainObj_in%index)
   domain_out => domainList(domainObj_out%index)

   call redistribute_5d( domain_in, field_in, domain_out, field_out, complete, dc_handle, position )

end subroutine mpp_redistribute_5d

subroutine redistribute_5d( domain_in, field_in, domain_out, field_out, complete, dc_handle, position )
   type(domain2D_private),    intent(in) :: domain_in, domain_out
   class(*),         intent(in)  :: field_in (:,:,:,:,:)
   class(*),         intent(out) :: field_out(:,:,:,:,:)
   logical, intent(in), optional :: complete 
   integer, intent(in), optional :: position 
   type(DomainCommunicator2D),pointer,optional :: dc_handle
   integer     :: isize_in, jsize_in, ksize_in, data_type_in
   integer     :: isize_out, jsize_out, ksize_out, data_type_out
   type(C_PTR) :: cptr_in, cptr_out
   integer(LONG_KIND) :: f_loc_in, f_loc_out

   data_type_in = get_field_type(field_in)
   data_type_out = get_field_type(field_out)
   isize_in = 0; jsize_in = 0; ksize_in = 1
   isize_out = 0; jsize_out = 0; ksize_out = 1
   cptr_in = c_null_ptr; cptr_out = c_null_ptr
   if(domain_in%initialized) then
      if(.not. is_contiguous(field_in)) call mpp_error(FATAL, "mpp_redistribute_5d: field_in is not contiguous")
      cptr_in = C_LOC(field_in)
      isize_in = size(field_in,1)
      jsize_in = size(field_in,2)
      ksize_in = size(field_in,3)*size(field_in,4)*size(field_in,5)
   endif
   if(domain_out%initialized) then
      if(.not. is_contiguous(field_out)) call mpp_error(FATAL, "mpp_redistribute_5d: field_out is not contiguous")
      cptr_out = C_LOC(field_out)
      isize_out = size(field_out,1)
      jsize_out = size(field_out,2)
      ksize_out = size(field_out,3)*size(field_out,4)*size(field_out,5)
   endif

   f_loc_in = LOC(field_in)
   f_loc_out = LOC(field_out)
   call redistribute_core( domain_in, cptr_in, f_loc_in, isize_in, jsize_in, ksize_in, data_type_in, &
        domain_out, cptr_out, f_loc_out, isize_out, jsize_out, ksize_out, data_type_out, complete, dc_handle, position )

   return
end subroutine redistribute_5d

subroutine mpp_do_redistribute_r8(dHold_in, dHold_out, d_comm)
   type(holder_r8_type),    intent(in)    :: dHold_in(:)
   type(holder_r8_type),  intent(inout)   :: dHold_out(:)
   type(DomainCommunicator2D), intent(in) :: d_comm

   real(DOUBLE_KIND) :: buffer(mpp_domains_stack_size)

#include "do_redistribute.inc"
   return

end subroutine mpp_do_redistribute_r8

subroutine mpp_do_redistribute_r4(dHold_in, dHold_out, d_comm)
   type(holder_r4_type),    intent(in)    :: dHold_in(:)
   type(holder_r4_type),  intent(inout)   :: dHold_out(:)
   type(DomainCommunicator2D), intent(in) :: d_comm

   real(FLOAT_KIND) :: buffer(mpp_domains_stack_size)

#include "do_redistribute.inc"
   return

end subroutine mpp_do_redistribute_r4

