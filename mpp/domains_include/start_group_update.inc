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

nlist   = size(domain%list(:))
gridtype = group%gridtype

!--- set reset_index_s and reset_index_v to 0
group%reset_index_s = 0
group%reset_index_v = 0

nrecv = group%nrecv
nsend = group%nsend
buffer_start_pos = group%buffer_start_pos

!---pre-post receive.
call mpp_clock_begin(group_recv_clock)
do m = 1, nrecv
   msgsize = group%recv_size(m)
   from_pe = group%from_pe(m)
   if( msgsize .GT. 0 )then
      buffer_pos = group%buffer_pos_recv(m) + group%buffer_start_pos
      call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.false., &
             tag=COMM_TAG_1, request=group%request_recv(m))
      group%type_recv(m) = mpitype
   end if
end do

!pack the data
call mpp_clock_end(group_recv_clock)

update_flags = group%flags_v
npack = group%npack

call mpp_clock_begin(group_pack_clock)
!pack the data

if( group%k_loop_inside ) then
!$OMP parallel do default(none) shared(npack,group,nvector,buffer_start_pos,update_flags,buffer,sHold,xHold,yHold) &
!$OMP                           private(buffer_pos,pos,m,is,ie,js,je,rotation,ke)
   do n = 1, npack
      buffer_pos = group%pack_buffer_pos(n) + buffer_start_pos
      pos  = buffer_pos
      is = group%pack_is(n); ie = group%pack_ie(n)
      js = group%pack_js(n); je = group%pack_je(n)
      rotation = group%pack_rotation(n)
      if( group%pack_type(n) == FIELD_S ) then
         ke = group%ksize_s
         do l=1, group%nscalar
            call pack_buffer_s(group%isize_s,group%jsize_s,sHold(l,1)%data,buffer, &
                                    pos, rotation, is, ie, js, je, ke)
         enddo
      else if( group%pack_type(n) == FIELD_X  ) then
         ke = group%ksize_v
         do l=1, group%nvector
            call pack_buffer_x(group%isize_x,group%jsize_x,xHold(l,1)%data,group%isize_y,group%jsize_y,yHold(l,1)%data, &
                    buffer,pos, rotation, update_flags, is, ie, js, je, ke)
         enddo
      else if( group%pack_type(n) == FIELD_Y  ) then
         ke = group%ksize_v
         do l=1, group%nvector
            call pack_buffer_y(group%isize_x,group%jsize_x,xHold(l,1)%data,group%isize_y,group%jsize_y,yHold(l,1)%data, &
                    buffer,pos, rotation, update_flags, is, ie, js, je, ke)
         enddo
      endif
   enddo
else
   call mpp_error(FATAL, "mpp_do_group_update: k_loop_inside = false is not implemented, contact developer")
endif
call mpp_clock_end(group_pack_clock)

call mpp_clock_begin(group_send_clock)
do n = 1, nsend
   msgsize = group%send_size(n)
   if( msgsize .GT. 0 )then
      buffer_pos = group%buffer_pos_send(n) + group%buffer_start_pos
      to_pe = group%to_pe(n)
      call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1, &
                       request=group%request_send(n))
   endif
enddo
call mpp_clock_end(group_send_clock)
