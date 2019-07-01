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
  type(overlap_type), pointer :: overPtr => NULL()
  logical   :: send(8), recv(8)
  integer   :: from_pe, to_pe, dir
  integer   :: m, n, l, i, j, k
  integer   :: is, ie, js, je, l_size
  integer   :: buffer_pos, msgsize
  integer   :: buffer_recv_size, pos

  call get_send_recv_value(flags, send, recv)
  l_size = size(f_hold(:))

  !--- pre-post receiving
  buffer_pos = 0  
  do m = 1, update%nrecv
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(nest_recv_clock)
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
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_FINE_3D_: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_1 )
        buffer_pos = buffer_pos + msgsize
     end if
     call mpp_clock_end(nest_recv_clock)
  end do ! end do m = 1, update%nrecv
  buffer_recv_size = buffer_pos

  !--- pack and send the data
  do m = 1, update%nsend
     overPtr => update%send(m)
     if( overPtr%count == 0 )cycle
     call mpp_clock_begin(nest_pack_clock)
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
           call mpp_error( FATAL, 'MPP_DO_UPDATE_NEST_FINE_3D_: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
     end if

     do n = 1, overPtr%count
        dir = overPtr%dir(n)
        if( send(dir) ) then
           is = overPtr%is(n); ie = overPtr%ie(n)
           js = overPtr%js(n); je = overPtr%je(n)
           do l=1,l_size  ! loop over number of fields
              do k = 1,ke  
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = f_hold(l)%data(i,j,k)
                    end do
                 end do
              end do
           end do
        endif
     end do ! do n = 1, overPtr%count

     call mpp_clock_end(nest_pack_clock)
     call mpp_clock_begin(nest_send_clock)
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then
        to_pe = overPtr%pe
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1 )
        buffer_pos = pos
     end if
     call mpp_clock_end(nest_send_clock)
  end do ! end do list = 0,nlist-1

  !unpack buffer
  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self(check=EVENT_RECV)
  call mpp_clock_end(nest_wait_clock)

  buffer_pos = buffer_recv_size      

  call mpp_clock_begin(nest_unpk_clock)
  do m = update%nrecv, 1, -1
     overPtr => update%recv(m)
     if( overPtr%count == 0 )cycle

     pos = buffer_pos
     do n = overPtr%count, 1, -1
        dir = overPtr%dir(n)
        if( recv(dir) ) then
           is = overPtr%is(n); ie = overPtr%ie(n)
           js = overPtr%js(n); je = overPtr%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke*l_size
           pos = buffer_pos - msgsize
           buffer_pos = pos
           dir = (dir+1)/2
           do l=1,l_size  ! loop over number of fields
              do k = 1,ke
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       b_hold(dir,l)%data(i,j,k) = buffer(pos)
                    end do
                 end do
              end do
           end do
        endif
     end do ! do n = 1, overPtr%count
  end do
  call mpp_clock_end(nest_unpk_clock)

  call mpp_clock_begin(nest_wait_clock)
  call mpp_sync_self( )
  call mpp_clock_end(nest_wait_clock)
