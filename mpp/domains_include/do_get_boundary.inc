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

  logical                 :: recv(4), send(4)
  integer                 :: nlist, buffer_pos, pos, tMe
  integer                 :: i, j, k, l, m, n, index, buffer_recv_size
  integer                 :: is, ie, js, je, msgsize, l_size, dir
  character(len=8)        :: text
  integer                 :: outunit

  outunit = stdout()
  l_size = size(f_hold,1)

  recv = .false.
  do i = 1, 4
     if( b_size(i) > 0 ) recv(i) = .true.
  enddo
  send = recv
  nlist = size(domain%list(:))  

  !recv
  buffer_pos = 0     
  do m = 1, bound%nrecv
     msgsize = 0
     do n = 1, bound%recv(m)%count
        if(recv(bound%recv(m)%dir(n))) then
           is = bound%recv(m)%is(n); ie = bound%recv(m)%ie(n)
           js = bound%recv(m)%js(n); je = bound%recv(m)%je(n)
           msgsize = msgsize + (ie-is+1)*(je-js+1)
        end if
     end do
     msgsize = msgsize*ke*l_size
     if( msgsize.GT.0 )then
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_OLD: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=bound%recv(m)%pe, block=.false., tag=COMM_TAG_2 )
        buffer_pos = buffer_pos + msgsize
     end if
  end do
  buffer_recv_size = buffer_pos

  ! send
  do m = 1, bound%nsend
     pos = buffer_pos
     do n = 1, bound%send(m)%count
        if(send(bound%send(m)%dir(n))) then
           is = bound%send(m)%is(n); ie = bound%send(m)%ie(n)
           js = bound%send(m)%js(n); je = bound%send(m)%je(n)
           tMe = bound%send(m)%tileMe(n)
           do l=1,l_size
              call pack_buffer_s(bound%xsize, bound%ysize,f_hold(l,tMe)%data,buffer, &
                                    pos, bound%send(m)%rotation(n), is, ie, js, je, ke)
           end do
        end if ! if(send(bound%dir(n)))
     end do ! do n = 1, bound%count
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then  
        !--- maybe we do not need the following stack size check.
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_OLD: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=bound%send(m)%pe, tag=COMM_TAG_2 )
        buffer_pos = pos
     end if
  end do

  call mpp_clock_begin(wait_clock)
  call mpp_sync_self(check=EVENT_RECV)
  call mpp_clock_end(wait_clock)
  buffer_pos = buffer_recv_size  

  !unpack recv
  !unpack buffer in reverse order.
  do m = bound%nrecv, 1, -1
     do n = bound%recv(m)%count, 1, -1
        if(recv(bound%recv(m)%dir(n))) then
           is = bound%recv(m)%is(n); ie = bound%recv(m)%ie(n)
           js = bound%recv(m)%js(n); je = bound%recv(m)%je(n)
           msgsize = (ie-is+1)*(je-js+1)*ke*l_size
           pos = buffer_pos - msgsize
           buffer_pos = pos
           tMe = bound%recv(m)%tileMe(n)
           dir = bound%recv(m)%dir(n)

           do l=1,l_size
              do k = 1, ke
                 index = bound%recv(m)%index(n)
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       b_hold(dir,l,tMe)%data(index,k) = buffer(pos)
                       index = index + 1
                    end do
                 end do
              end do
           end do
        end if
     end do
  end do

  call mpp_sync_self( )

