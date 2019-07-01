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

  ke = size(field_SG,3)
  buffer_pos = 0
  !---pre-post receive
  do m = 1, UG_domain%SG2UG%nrecv
     msgsize = UG_domain%SG2UG%recv(m)%count*ke
     if( msgsize.GT.0 )then
        from_pe = UG_domain%SG2UG%recv(m)%pe
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'mpp_pass_SG_to_UG_3D: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_1 )
        buffer_pos = buffer_pos + msgsize
     end if
  end do

  !---pack and send data
  do m = 1, UG_domain%SG2UG%nsend
     pos = buffer_pos
     msgsize = UG_domain%SG2UG%send(m)%count * ke
     if( msgsize.GT.0 )then
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'mpp_pass_SG_to_UG_3D: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if

        do k = 1, ke
           do l = 1, UG_domain%SG2UG%send(m)%count
              i = UG_domain%SG2UG%send(m)%i(l)+ioff
              j = UG_domain%SG2UG%send(m)%j(l)+joff
              pos = pos+1
              buffer(pos) = field_SG(i,j,k)
           end do
        end do
        to_pe = UG_domain%SG2UG%send(m)%pe
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1 )
        buffer_pos = buffer_pos + msgsize
     end if
  end do

  call mpp_sync_self(check=EVENT_RECV)
  !--- unpack the buffer
  buffer_pos = 0
  do m = 1, UG_domain%SG2UG%nrecv
     pos = buffer_pos
     do k = 1, ke
        do l = 1, UG_domain%SG2UG%recv(m)%count
           pos = pos+1
           i = UG_domain%SG2UG%recv(m)%i(l)
           field_UG(i,k) = buffer(pos)
        enddo
     enddo
     buffer_pos = pos
  enddo

  call mpp_sync_self( )


