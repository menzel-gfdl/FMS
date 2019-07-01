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

      integer :: l_size, l, i, j, k, is, ie, js, je, n
      integer :: pos, nlist, msgsize
      integer :: to_pe, from_pe
      integer :: tMe
      integer :: buffer_pos
      character(len=8) :: text
      character(len=64) :: field_name      
      integer :: buffer_recv_size
      integer :: rank_x, rank_y, ind_x, ind_y, cur_rank
      integer :: nsend_x, nsend_y, nrecv_x, nrecv_y
      integer :: outunit

      outunit = stdout()
      buffer_pos = 0        !this initialization goes away if update_domains becomes non-blocking
      l_size = size(xHold,1)
      nlist = size(domain%list(:))
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

      nsend_x = check_x%nsend
      nsend_y = check_y%nsend
      nrecv_x = check_x%nrecv
      nrecv_y = check_y%nrecv

      !--- recv the data       
      cur_rank = get_rank_recv(domain, check_x, check_y, rank_x, rank_y, ind_x, ind_y) 

      do while ( ind_x .LE. nrecv_x .OR. ind_y .LE. nrecv_y )
         msgsize = 0
         if(cur_rank == rank_x) then
            from_pe = check_x%recv(ind_x)%pe
            do n = 1, check_x%recv(ind_x)%count
               is = check_x%recv(ind_x)%is(n); ie = check_x%recv(ind_x)%ie(n)
               js = check_x%recv(ind_x)%js(n); je = check_x%recv(ind_x)%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end do
            ind_x = ind_x+1
            if(ind_x .LE. nrecv_x) then
               rank_x = check_x%recv(ind_x)%pe - domain%pe 
               if(rank_x .LE.0) rank_x = rank_x + nlist
            else
               rank_x = -1
            endif
         endif
         if(cur_rank == rank_y) then
            from_pe = check_y%recv(ind_y)%pe
            do n = 1, check_y%recv(ind_y)%count
               is = check_y%recv(ind_y)%is(n); ie = check_y%recv(ind_y)%ie(n)
               js = check_y%recv(ind_y)%js(n); je = check_y%recv(ind_y)%je(n)
               msgsize = msgsize + (ie-is+1)*(je-js+1)
            end do
            ind_y = ind_y+1
            if(ind_y .LE. nrecv_y) then
               rank_y = check_y%recv(ind_y)%pe - domain%pe 
               if(rank_y .LE.0) rank_y = rank_y + nlist
            else
               rank_y = -1
            endif
         endif
         cur_rank = max(rank_x, rank_y)
         msgsize = msgsize*ke*l_size
         if( msgsize.GT.0 )then
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_CHECK_V: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
            end if
            call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.false., tag=COMM_TAG_2 )
            buffer_pos = buffer_pos + msgsize
         end if
      enddo
      buffer_recv_size = buffer_pos

      !--- send the data
      cur_rank = get_rank_send(domain, check_x, check_y, rank_x, rank_y, ind_x, ind_y) 

      do while (ind_x .LE. nsend_x .OR. ind_y .LE. nsend_y)
         pos = buffer_pos
         if(cur_rank == rank_x) then
            to_pe = check_x%send(ind_x)%pe
            do n = 1, check_x%send(ind_x)%count
               is = check_x%send(ind_x)%is(n); ie = check_x%send(ind_x)%ie(n)
               js = check_x%send(ind_x)%js(n); je = check_x%send(ind_x)%je(n)
               tMe = check_x%send(ind_x)%tileMe(n)
               do l=1,l_size  ! loop over number of fields
                  call pack_buffer_x(check_x%xsize,check_x%ysize,xHold(l,tMe)%data,check_y%xsize,check_y%ysize, &
                       yHold(l,tMe)%data, buffer,pos, check_x%send(ind_x)%rotation(n), flags, is, ie, js, je, ke)
               end do
            end do
            ind_x = ind_x+1
            if(ind_x .LE. nsend_x) then
               rank_x = check_x%send(ind_x)%pe - domain%pe 
               if(rank_x .LT.0) rank_x = rank_x + nlist
            else
               rank_x = nlist+1
            endif
         endif

         if(cur_rank == rank_y) then
            to_pe = check_y%send(ind_y)%pe            
            do n = 1, check_y%send(ind_y)%count
               is = check_y%send(ind_y)%is(n); ie = check_y%send(ind_y)%ie(n)
               js = check_y%send(ind_y)%js(n); je = check_y%send(ind_y)%je(n)
               tMe = check_y%send(ind_y)%tileMe(n)
               do l=1,l_size  ! loop over number of fields
                  call pack_buffer_y(check_x%xsize,check_x%ysize,xHold(l,tMe)%data,check_y%xsize,check_y%ysize, &
                       yHold(l,tMe)%data, buffer,pos, check_y%send(ind_x)%rotation(n), flags, is, ie, js, je, ke)
               end do
            end do
            ind_y = ind_y+1
            if(ind_y .LE. nsend_y) then
               rank_y = check_y%send(ind_y)%pe - domain%pe 
               if(rank_y .LT.0) rank_y = rank_y + nlist
            else
               rank_y = nlist+1
            endif
         endif
         cur_rank = min(rank_x, rank_y)
         msgsize = pos - buffer_pos
         if( msgsize.GT.0 )then
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos)
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_CHECK_V: mpp_domains_stack overflow, ' // &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
            end if
            call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_2 )
            buffer_pos = pos
         end if
      end do ! end do list = 0,nlist-1

      call mpp_sync_self(check=EVENT_RECV) ! To ensure recv is completed.
      buffer_pos = buffer_recv_size

      !--- compare the data in reverse order
      cur_rank = get_rank_unpack(domain, check_x, check_y, rank_x, rank_y, ind_x, ind_y) 

      CHECK_LOOP: do while(ind_x >0 .OR. ind_y >0)
         if(cur_rank == rank_y) then
            do n = check_y%recv(ind_y)%count, 1, -1
               is = check_y%recv(ind_y)%is(n); ie = check_y%recv(ind_y)%ie(n)
               js = check_y%recv(ind_y)%js(n); je = check_y%recv(ind_y)%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke*l_size
               pos = buffer_pos - msgsize
               buffer_pos = pos
               tMe = check_y%recv(ind_y)%tileMe(n)
               do l=1,l_size  ! loop over number of fields
                  do k = 1,ke
                     do j = js, je
                        do i = is, ie
                           pos = pos + 1
                           if( yHold(l,tMe)%data(i,j,k) .NE. buffer(pos) ) then
                              print*,"Error from MPP_DO_CHECK_V on pe = ", mpp_pe(), ": y component of vector ", &
                                   trim(field_name), " at point (", i, ",", j, ",", k, ") = ", yHold(l,tMe)%data(i,j,k), &
                                   " does not equal to the value = ", buffer(pos), " on pe ", check_y%recv(ind_y)%pe
                              call mpp_error(debug_update_level, "MPP_DO_CHECK_V: mismatch on the boundary for symmetry point")
                              exit CHECK_LOOP
                           end if
                        end do
                     end do
                  end do
               end do
            end do
            ind_y = ind_y-1
            if(ind_y .GT. 0) then
               rank_y = check_y%recv(ind_y)%pe - domain%pe 
               if(rank_y .LE.0) rank_y = rank_y + nlist
            else
               rank_y = nlist+1
            endif
         endif

         if(cur_rank == rank_x) then
            do n = check_x%recv(ind_x)%count, 1, -1
               is = check_x%recv(ind_x)%is(n); ie = check_x%recv(ind_x)%ie(n)
               js = check_x%recv(ind_x)%js(n); je = check_x%recv(ind_x)%je(n)
               msgsize = (ie-is+1)*(je-js+1)*ke*l_size
               pos = buffer_pos - msgsize
               buffer_pos = pos
               tMe = check_x%recv(ind_x)%tileMe(n)
               do l=1,l_size  ! loop over number of fields
                  do k = 1,ke
                     do j = js, je
                        do i = is, ie
                           pos = pos + 1
                           if( xHold(l,tMe)%data(i,j,k) .NE. buffer(pos) ) then
                              print*,"Error from MPP_DO_CHECK_V on pe = ", mpp_pe(), ": x-component of vector ", &
                                   trim(field_name), " at point (", i, ",", j, ",", k, ") = ", xHold(l,tMe)%data(i,j,k), &
                                   " does not equal to the value = ", buffer(pos), " on pe ", check_x%recv(ind_x)%pe
                              call mpp_error(debug_update_level, "MPP_DO_CHECK_V: mismatch on the boundary for symmetry point")
                              exit CHECK_LOOP
                           end if
                        end do
                     end do
                  end do
               end do
            end do
            ind_x = ind_x-1
            if(ind_x .GT. 0) then
               rank_x = check_x%recv(ind_x)%pe - domain%pe 
               if(rank_x .LE.0) rank_x = rank_x + nlist
            else
               rank_x = nlist+1
            endif
         endif
         cur_rank = min(rank_x, rank_y)
      end do CHECK_LOOP ! end do list = nlist-1,0,-1
      call mpp_sync_self()
