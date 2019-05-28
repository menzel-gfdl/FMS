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

   type(domain2D_private), pointer :: domain_in, domain_out
   integer :: i, j, k, l, n, l_size
   integer :: is, ie, js, je
   integer :: ke
   integer :: list, pos, msgsize
   integer :: to_pe, from_pe
   integer :: buffer_pos, errunit
   integer :: isize_in, jsize_in, isize_out, jsize_out

   errunit = stderr()
   l_size = size(dHold_out(:))  ! equal to size(dHold_in(:))
   ke = d_comm%ke
   domain_in =>d_comm%domain_in; domain_out =>d_comm%domain_out
   isize_in = domain_in%x(1)%data%size; jsize_in = domain_in%y(1)%data%size
   isize_out = domain_out%x(1)%data%size; jsize_out = domain_out%y(1)%data%size

   buffer_pos = 0

!pre-post recv
   n = d_comm%Rlist_size
   do list = 0,n-1
      if( .NOT. d_comm%R_do_buf(list) )cycle
      from_pe = d_comm%cfrom_pe(list)
      msgsize = d_comm%R_msize(list)*l_size
      call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_1 )
      buffer_pos = buffer_pos + msgsize
   enddo

!send
   n = d_comm%Slist_size
   do list = 0,n-1
      if( .NOT. d_comm%S_do_buf(list) )cycle
      to_pe = d_comm%cto_pe(list)
      is=d_comm%sendis(1,list); ie=d_comm%sendie(1,list)
      js=d_comm%sendjs(1,list); je=d_comm%sendje(1,list)
      pos = buffer_pos
      do l=1,l_size  ! loop over number of fields
         do k = 1,ke
            do j = js,je
               do i = is,ie
                  pos = pos+1
                  buffer(pos) = dHold_in(l)%data(i,j,k)
               end do
            end do
         end do
       end do
       if( debug )write( errunit,* )'PE', pe, ' to PE ', to_pe, 'is,ie,js,je=', is, ie, js, je
       msgsize = pos - buffer_pos
       call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1  )
       buffer_pos = pos
   end do

   call mpp_sync_self(check=EVENT_RECV)

!unpack buffer
   buffer_pos = 0
   n = d_comm%Rlist_size
   do list = 0,n-1
      if( .NOT. d_comm%R_do_buf(list) )cycle
      from_pe = d_comm%cfrom_pe(list)
      is=d_comm%recvis(1,list); ie=d_comm%recvie(1,list)
      js=d_comm%recvjs(1,list); je=d_comm%recvje(1,list)
      if( debug )write( errunit,* )'PE', pe, ' from PE ', from_pe, 'is,ie,js,je=', is, ie, js, je
      pos = buffer_pos
      do l=1,l_size  ! loop over number of in/out fields
         do k = 1,ke
            do j = js,je
               do i = is,ie
                  pos = pos+1
                  dHold_out(l)%data(i,j,k) = buffer(pos)
               end do
            end do
         end do
       end do
       buffer_pos = pos
   end do

   call mpp_sync_self()


