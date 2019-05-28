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

      integer :: l_size, l, i, j, k, is, ie, js, je, n, m
      integer :: pos, nlist, msgsize
      integer :: to_pe, from_pe, midpoint
      integer :: tMe, dir

      integer :: send_start_pos, nsend
      integer :: send_msgsize(2*MAXLIST)
      integer :: send_pe(2*MAXLIST)
      logical :: send_x(8), recv_x(8)
      logical :: send_y(8), recv_y(8)
      logical :: nonsym_edge
      integer :: buffer_pos, isd
      character(len=8) :: text
      integer :: buffer_recv_size, shift
      integer :: nsend_x, nsend_y, nrecv_x, nrecv_y
      integer :: cur_rank, rank_x, rank_y, ind_x, ind_y
      integer :: whalo_v, ehalo_v, shalo_v, nhalo_v
      integer :: ied, jsd, jed

      if(num_nonblock_group_update .NE. 0) then
         call mpp_error(FATAL,"mpp_do_updateV: can not be called between mpp_start_group_update/mpp_complete_group_update")
      endif

      if( BTEST(update_flags,NORTH) .AND. BTEST(domain%fold,NORTH) .AND. BTEST(gridtype,SOUTH) ) &
           call mpp_error( FATAL, 'MPP_DO_UPDATE_V: Incompatible grid offset and fold.' )

      call get_vector_send_recv_value(update_flags, gridtype, domain%symmetry, nonsym_edge, send_x, send_y, recv_x, recv_y)

      l_size = size(xHold,1)
      nlist = size(domain%list(:))
      nsend_x = update_x%nsend
      nsend_y = update_y%nsend
      nrecv_x = update_x%nrecv
      nrecv_y = update_y%nrecv
      whalo_v = update_x%whalo
      ehalo_v = update_x%ehalo
      shalo_v = update_y%shalo
      nhalo_v = update_y%nhalo

      !--- recv
      buffer_pos = 0
      cur_rank = get_rank_recv(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y) 
      call mpp_clock_begin(recv_clock)
      do while (ind_x .LE. nrecv_x .OR. ind_y .LE. nrecv_y)
         msgsize = 0
         if(cur_rank == rank_x) then
            from_pe = update_x%recv(ind_x)%pe
            do n = 1, update_x%recv(ind_x)%count
               dir = update_x%recv(ind_x)%dir(n)
               if(recv_x(dir)) then
                  is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                  js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end if
            end do
            ind_x = ind_x+1
            if(ind_x .LE. nrecv_x) then
               rank_x = update_x%recv(ind_x)%pe - domain%pe 
               if(rank_x .LE.0) rank_x = rank_x + nlist
            else
               rank_x = -1
            endif
         endif
         if(cur_rank == rank_y) then
            from_pe = update_y%recv(ind_y)%pe
            do n = 1, update_y%recv(ind_y)%count
               dir = update_y%recv(ind_y)%dir(n)
               if(recv_y(dir)) then
                  is = update_y%recv(ind_y)%is(n); ie = update_y%recv(ind_y)%ie(n)
                  js = update_y%recv(ind_y)%js(n); je = update_y%recv(ind_y)%je(n)
                  msgsize = msgsize + (ie-is+1)*(je-js+1)
               end if
            end do
            ind_y = ind_y+1
            if(ind_y .LE. nrecv_y) then
               rank_y = update_y%recv(ind_y)%pe - domain%pe 
               if(rank_y .LE.0) rank_y = rank_y + nlist
            else
               rank_y = -1
            endif
         endif
         cur_rank = max(rank_x, rank_y)
         msgsize = msgsize*ke*l_size
   
         if( msgsize.GT.0 )then
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, buffer_pos+msgsize )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                write( text,'(i8)' )mpp_domains_stack_hwm
                call mpp_error( FATAL, 'MPP_DO_UPDATE_V: mpp_domains_stack overflow, '// &
                     'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.false., tag=COMM_TAG_2 )
             buffer_pos = buffer_pos + msgsize
         end if
      end do
      call mpp_clock_end(recv_clock)
      buffer_recv_size = buffer_pos
      send_start_pos = buffer_pos

      !--- send
      cur_rank = get_rank_send(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y) 
      nsend = 0
      call mpp_clock_begin(pack_clock)
      do while (ind_x .LE. nsend_x .OR. ind_y .LE. nsend_y)
         pos = buffer_pos
         !--- make sure the domain stack size is big enough
         msgsize = 0
         if(cur_rank == rank_x) then
            do n = 1, update_x%send(ind_x)%count
               dir = update_x%send(ind_x)%dir(n)
               if( send_x(dir) ) msgsize = msgsize +  update_x%send(ind_x)%msgsize(n)
            enddo
         endif
         if(cur_rank == rank_y) then
            do n = 1, update_y%send(ind_y)%count
               dir = update_y%send(ind_y)%dir(n)
               if( send_y(dir) ) msgsize = msgsize +  update_y%send(ind_y)%msgsize(n)
            enddo
         endif

         if( msgsize.GT.0 )then
            msgsize = msgsize*ke*l_size
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_UPDATE_V: mpp_domains_stack overflow, ' // &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
            end if
         end if
         if(cur_rank == rank_x) then
            to_pe = update_x%send(ind_x)%pe
            do n = 1, update_x%send(ind_x)%count
               dir = update_x%send(ind_x)%dir(n)
               if( send_x(dir) ) then
                  tMe = update_x%send(ind_x)%tileMe(n)
                  is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
                  js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)
                  do l=1,l_size  ! loop over number of fields
                     call pack_buffer_x(update_x%xsize,update_x%ysize,xHold(l,tMe)%data,update_y%xsize,update_y%ysize, &
                          yHold(l,tMe)%data, buffer,pos, update_x%send(ind_x)%rotation(n), update_flags, is, ie, js, je, ke)
                  end do
               endif
            enddo
            ind_x = ind_x+1
            if(ind_x .LE. nsend_x) then
               rank_x = update_x%send(ind_x)%pe - domain%pe
               if(rank_x .LT.0) rank_x = rank_x + nlist
            else
               rank_x = nlist+1
            endif
         endif
         if(cur_rank == rank_y) then
            to_pe = update_y%send(ind_y)%pe
            do n = 1, update_y%send(ind_y)%count
               dir = update_y%send(ind_y)%dir(n)
               if( send_y(dir) ) then
                  tMe = update_y%send(ind_y)%tileMe(n)
                  is = update_y%send(ind_y)%is(n); ie = update_y%send(ind_y)%ie(n)
                  js = update_y%send(ind_y)%js(n); je = update_y%send(ind_y)%je(n)
                  do l=1,l_size  ! loop over number of fields
                     call pack_buffer_y(update_x%xsize,update_x%ysize,xHold(l,tMe)%data,update_y%xsize,update_y%ysize, &
                          yHold(l,tMe)%data, buffer,pos, update_y%send(ind_y)%rotation(n), update_flags, is, ie, js, je, ke)
                  end do
               endif
            enddo
            ind_y = ind_y+1
            if(ind_y .LE. nsend_y) then
               rank_y = update_y%send(ind_y)%pe - domain%pe 
               if(rank_y .LT.0) rank_y = rank_y + nlist
            else
               rank_y = nlist+1
            endif
         endif
         cur_rank = min(rank_x, rank_y)
         nsend = nsend + 1
         send_pe(nsend) = to_pe
         send_msgsize(nsend) = pos - buffer_pos
         buffer_pos = pos
      end do

      buffer_pos = send_start_pos
      call mpp_clock_end(pack_clock)
      call mpp_clock_begin(send_clock)
      do m = 1, nsend
         msgsize = send_msgsize(m)
         if( msgsize.GT.0 )then
            call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=send_pe(m), tag=COMM_TAG_2 )
            buffer_pos = buffer_pos + msgsize
         end if
      end do 
      call mpp_clock_end(send_clock)

!unpack recv
!unpack halos in reverse order
      call mpp_clock_begin(wait_clock)
      call mpp_sync_self(check=EVENT_RECV)
      call mpp_clock_end(wait_clock)
      buffer_pos = buffer_recv_size      
      cur_rank = get_rank_unpack(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y) 

      call mpp_clock_begin(unpk_clock)
      do while (ind_x > 0 .OR. ind_y > 0)
         pos = buffer_pos
         if(cur_rank == rank_y) then
            do n = update_y%recv(ind_y)%count, 1, -1
               dir = update_y%recv(ind_y)%dir(n)
               if( recv_y(dir) ) then
                  tMe = update_y%recv(ind_y)%tileMe(n)
                  is = update_y%recv(ind_y)%is(n); ie = update_y%recv(ind_y)%ie(n)
                  js = update_y%recv(ind_y)%js(n); je = update_y%recv(ind_y)%je(n)
                  msgsize = (ie-is+1)*(je-js+1)*ke*l_size
                  pos = buffer_pos - msgsize
                  buffer_pos = pos
                  do l=1,l_size  ! loop over number of fields
                     call unpack_buffer(update_y%xsize, update_y%ysize, yHold(l,tMe)%data, buffer, pos, is, ie, js, je, ke)
                  enddo
               end if
            end do
            ind_y = ind_y-1
            if(ind_y .GT. 0) then
               rank_y = update_y%recv(ind_y)%pe - domain%pe
               if(rank_y .LE.0) rank_y = rank_y + nlist
            else
               rank_y = nlist+1
             endif
         endif
         if(cur_rank == rank_x) then
            do n = update_x%recv(ind_x)%count, 1, -1
               dir = update_x%recv(ind_x)%dir(n)
               if( recv_x(dir) ) then
                  tMe = update_x%recv(ind_x)%tileMe(n)
                  is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                  js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n) 
                  msgsize = (ie-is+1)*(je-js+1)*ke*l_size
                  pos = buffer_pos - msgsize
                  buffer_pos = pos
                  do l=1,l_size  ! loop over number of fields
                  call unpack_buffer(update_x%xsize, update_x%ysize, xHold(l,tMe)%data, buffer, pos, is, ie, js, je, ke)
                  end do
               end if
            end do
            ind_x = ind_x-1
            if(ind_x .GT. 0) then
               rank_x = update_x%recv(ind_x)%pe - domain%pe 
               if(rank_x .LE.0) rank_x = rank_x + nlist
            else
               rank_x = nlist+1
            endif
         endif
         cur_rank = min(rank_x, rank_y)
      end do
      call mpp_clock_end(unpk_clock)

#include "update_fold.inc"

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)


