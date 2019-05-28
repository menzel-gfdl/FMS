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

      character(len=8)            :: text
      character(len=64)           :: field_name

!equate to mpp_domains_stack
      integer :: buffer_pos
      integer :: to_pe, from_pe, pos, msgsize
      integer :: n, l_size, l, m, i, j, k
      integer :: is, ie, js, je, tMe
      integer :: buffer_recv_size
      integer :: outunit

      outunit = stdout()
      l_size = size(sHold,1)

      !--- if debug_update_level is not NO_DEBUG, check the consistency on the bounds
      !--- (domain is symmetry or folded north edge). North bound will be checked when north edge is folded.
      !--- when domain is symmetry, For data on T-cell, no check is needed; for data on E-cell,
      !--- data on East and West boundary will be checked ; For data on N-cell, data on North and South
      !--- boundary will be checked; For data on C-cell, data on West, East, South, North will be checked.
      !--- The check will be done in the following way: Western boundary data sent to Eastern boundary to check
      !--- and Southern boundary to check

      if(present(name)) then
         field_name = name
      else
         field_name = "un-named"
      end if

      buffer_pos = 0
      !--- pre-post recv the data
      do m = 1, check%nrecv
         msgsize = 0
         do n = 1, check%recv(m)%count
            is = check%recv(m)%is(n); ie = check%recv(m)%ie(n)
            js = check%recv(m)%js(n); je = check%recv(m)%je(n)
            msgsize = msgsize + (ie-is+1)*(je-js+1)
         end do
         msgsize = msgsize*ke*l_size

         if( msgsize.GT.0 )then
            from_pe = check%recv(m)%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_CHECK: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
            end if
            call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_2 )
            buffer_pos = buffer_pos + msgsize
         end if
      end do
      buffer_recv_size = buffer_pos

      !--- send the data
      do m = 1, check%nsend
         pos = buffer_pos
         do n = 1, check%recv(m)%count
            is = check%recv(m)%is(n); ie = check%recv(m)%ie(n)
            js = check%recv(m)%js(n); je = check%recv(m)%je(n)
            tMe = check%recv(m)%tileMe(n)
            do l=1,l_size  ! loop over number of fields
               call pack_buffer_s(check%xsize, check%ysize,sHold(l,tMe)%data,buffer, &
                                 pos, check%recv(m)%rotation(n), is, ie, js, je, ke)
            end do
         end do
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            to_pe = check%recv(m)%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos)
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_CHECK: mpp_domains_stack overflow, ' // &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
            end if
            call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_2 )
            buffer_pos = pos
         end if
      end do ! end do list = 0,nlist-1

      call mpp_sync_self(check=EVENT_RECV) ! To ensure recv is completed.
      buffer_pos = buffer_recv_size
      !--- compare the data in reverse order
      CHECK_LOOP: do m = check%nrecv, 1, -1
         do n = check%recv(m)%count, 1, -1
            is = check%recv(m)%is(n); ie = check%recv(m)%ie(n)
            js = check%recv(m)%js(n); je = check%recv(m)%je(n)
            msgsize = (ie-is+1)*(je-js+1)*ke*l_size
            pos = buffer_pos - msgsize
            buffer_pos = pos
            tMe = check%recv(m)%tileMe(n)
            do l=1, l_size  ! loop over number of fields
               do k = 1,ke
                  do j = js, je
                     do i = is, ie
                        pos = pos + 1
                        if( sHold(l,tMe)%data(i,j,k) .NE. buffer(pos) ) then
                           print*,"Error from MPP_DO_CHECK on pe = ", mpp_pe(), ": field ", &
                                trim(field_name), " at point (", i, ",", j, ",", k, ") = ", sHold(l,tMe)%data(i,j,k), &
                                " does not equal to the value = ", buffer(pos), " on pe ", check%recv(m)%pe
                           call mpp_error(debug_update_level, "MPP_DO_CHECK: mismatch on the boundary for symmetry point")
                           exit CHECK_LOOP
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do CHECK_LOOP ! end do list = nlist-1,0,-1
      call mpp_sync_self()
