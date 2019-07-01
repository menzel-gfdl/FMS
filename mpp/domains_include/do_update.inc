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

      integer :: buffer_pos
      type(overlap_type), pointer :: overPtr => NULL()      
      character(len=8)            :: text     
      logical :: send(8), recv(8)
      integer :: to_pe, from_pe, pos, msgsize
      integer :: n, l_size, l, m, i, j, k
      integer :: is, ie, js, je, tMe, dir
      integer :: buffer_recv_size
      integer :: send_start_pos
      integer :: send_msgsize(MAXLIST)

      if(num_nonblock_group_update .NE. 0) then
         call mpp_error(FATAL,"mpp_do_update: can not be called between mpp_start_group_update/mpp_complete_group_update")
      endif

      l_size = size(sHold,1)
  
      call get_send_recv_value(update_flags, send, recv)
      !recv
      buffer_pos = 0
      call mpp_clock_begin(recv_clock)
      do m = 1, update%nrecv
         overPtr => update%recv(m)
         if( overPtr%count == 0 )cycle
         msgsize = 0
         do n = 1, overPtr%count
            dir = overPtr%dir(n)
            if(recv(dir)) then
               is = overPtr%is(n); ie = overPtr%ie(n)
               js = overPtr%js(n); je = overPtr%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end if
         end do

         msgsize = msgsize*ke*l_size
         if( msgsize.GT.0 )then
            from_pe = overPtr%pe
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_UPDATE: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
            end if
            call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_2 )
            buffer_pos = buffer_pos + msgsize
         end if
      end do ! end do m = 1, update%nrecv
      call mpp_clock_end(recv_clock)

      buffer_recv_size = buffer_pos
      send_start_pos = buffer_pos
      ! pack
      call  mpp_clock_begin(pack_clock)
      do m = 1, update%nsend
         send_msgsize(m) = 0
         overPtr => update%send(m)
         if( overPtr%count == 0 )cycle
         pos = buffer_pos
         msgsize = 0
         do n = 1, overPtr%count
            dir = overPtr%dir(n)
            if( send(dir) )  msgsize = msgsize + overPtr%msgsize(n)
         enddo
         if( msgsize.GT.0 )then
            msgsize = msgsize*ke*l_size
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_UPDATE: mpp_domains_stack overflow, ' // &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
            end if
         end if

         do n = 1, overPtr%count
            dir = overPtr%dir(n)
            if( send(dir) ) then
               tMe = overPtr%tileMe(n)
               is = overPtr%is(n); ie = overPtr%ie(n)
               js = overPtr%js(n); je = overPtr%je(n)
               do l=1,l_size  ! loop over number of fields
                  call pack_buffer_s(update%xsize, update%ysize,sHold(l,tMe)%data,buffer, &
                                    pos, overPtr%rotation(n), is, ie, js, je, ke)
               end do
            endif
         end do ! do n = 1, overPtr%count
         send_msgsize(m) = pos-buffer_pos
         buffer_pos = pos
      end do ! end do m = 1, nsend
      call mpp_clock_end(pack_clock)

      buffer_pos = send_start_pos
      call mpp_clock_begin(send_clock)
      do m = 1, update%nsend
         msgsize = send_msgsize(m)
         if(msgsize == 0) cycle
         to_pe = update%send(m)%pe
         call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_2 )
         buffer_pos = buffer_pos + msgsize
      end do ! end do ist = 0,nlist-1
      call mpp_clock_end(send_clock)

      !unpack halos in reverse order
      call mpp_clock_begin(wait_clock)
      call mpp_sync_self(check=EVENT_RECV)
      call mpp_clock_end(wait_clock)
      buffer_pos = buffer_recv_size      

      call mpp_clock_begin(unpk_clock)
      do m = update%nrecv, 1, -1
         overPtr => update%recv(m)
         if( overPtr%count == 0 )cycle
         pos = buffer_pos
         do n = overPtr%count, 1, -1
            dir = overPtr%dir(n)
            if( recv(dir) ) then
               tMe = overPtr%tileMe(n)
               is = overPtr%is(n); ie = overPtr%ie(n)
               js = overPtr%js(n); je = overPtr%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke*l_size
               pos = buffer_pos - msgsize
               buffer_pos = pos
               do l=1,l_size  ! loop over number of fields
                  call unpack_buffer(update%xsize, update%ysize, sHold(l,tMe)%data, buffer, pos, is, ie, js, je, ke)
               end do
            endif
         end do ! do n = 1, overPtr%count
      end do
      call mpp_clock_end(unpk_clock)

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)

