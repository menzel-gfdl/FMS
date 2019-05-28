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

  logical                 :: recvx(4), sendx(4)
  logical                 :: recvy(4), sendy(4)
  integer                 :: nlist, buffer_pos, pos, tMe
  integer                 :: is, ie, js, je, msgsize, l_size, buffer_recv_size
  integer                 :: i, j, k, l, n, index, to_pe, from_pe, dir
  integer                 :: rank_x, rank_y, cur_rank, ind_x, ind_y
  integer                 :: nsend_x, nsend_y, nrecv_x, nrecv_y
  character(len=8)        :: text
  integer                 :: outunit, shift, midpoint

  outunit = stdout()
  l_size = size(fx_hold,1)
  !---- determine recv(1) based on b_addrs ( east boundary )
  recvx = .false.; recvy = .false.
  do i = 1, 4
     if( b_sizex(i) > 0 ) recvx(i) = .true.
     if( b_sizey(i) > 0 ) recvy(i) = .true.
  enddo

  sendx = recvx
  sendy = recvy

  nlist = size(domain%list(:))  

  nsend_x = boundx%nsend
  nsend_y = boundy%nsend
  nrecv_x = boundx%nrecv
  nrecv_y = boundy%nrecv

  !recv
  buffer_pos = 0     
  cur_rank = get_rank_recv(domain, boundx, boundy, rank_x, rank_y, ind_x, ind_y)   

  do while ( ind_x .LE. nrecv_x .OR. ind_y .LE. nrecv_y )
     msgsize = 0
     if(cur_rank == rank_x) then
        from_pe = boundx%recv(ind_x)%pe
        do n = 1, boundx%recv(ind_x)%count
           if(recvx(boundx%recv(ind_x)%dir(n))) then
              is = boundx%recv(ind_x)%is(n); ie = boundx%recv(ind_x)%ie(n)
              js = boundx%recv(ind_x)%js(n); je = boundx%recv(ind_x)%je(n)
              msgsize = msgsize + (ie-is+1)*(je-js+1)
           end if
        end do
        ind_x = ind_x+1
        if(ind_x .LE. nrecv_x) then
           rank_x = boundx%recv(ind_x)%pe - domain%pe 
           if(rank_x .LE.0) rank_x = rank_x + nlist
        else
           rank_x = -1
        endif
     endif

     if(cur_rank == rank_y) then
        from_pe = boundy%recv(ind_y)%pe
        do n = 1, boundy%recv(ind_y)%count
           if(recvy(boundy%recv(ind_y)%dir(n))) then
              is = boundy%recv(ind_y)%is(n); ie = boundy%recv(ind_y)%ie(n)
              js = boundy%recv(ind_y)%js(n); je = boundy%recv(ind_y)%je(n)
              msgsize = msgsize + (ie-is+1)*(je-js+1)
           end if
        end do
        ind_y = ind_y+1
        if(ind_y .LE. nrecv_y) then
           rank_y = boundy%recv(ind_y)%pe - domain%pe 
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
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_V_: mpp_domains_stack overflow, '// &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_4 )
        buffer_pos = buffer_pos + msgsize
     end if
  end do
  buffer_recv_size = buffer_pos

  ! send
  cur_rank = get_rank_send(domain, boundx, boundy, rank_x, rank_y, ind_x, ind_y) 

  do while (ind_x .LE. nsend_x .OR. ind_y .LE. nsend_y)
     pos = buffer_pos
     if(cur_rank == rank_x) then
        to_pe = boundx%send(ind_x)%pe
        do n = 1, boundx%send(ind_x)%count
           if(sendx(boundx%send(ind_x)%dir(n))) then
              is = boundx%send(ind_x)%is(n); ie = boundx%send(ind_x)%ie(n)
              js = boundx%send(ind_x)%js(n); je = boundx%send(ind_x)%je(n)
              tMe = boundx%send(ind_x)%tileMe(n)
              do l=1,l_size
                 call pack_buffer_x(boundx%xsize,boundx%ysize,fx_hold(l,tMe)%data,boundy%xsize,boundy%ysize, &
                          fy_hold(l,tMe)%data, buffer,pos, boundx%send(ind_x)%rotation(n), flags, is, ie, js, je, ke)
              end do
           end if ! if(send(boundx%dir(n)))
        end do  !do n = 1, boundx%count
        ind_x = ind_x+1
        if(ind_x .LE. nsend_x) then
           rank_x = boundx%send(ind_x)%pe - domain%pe 
           if(rank_x .LT.0) rank_x = rank_x + nlist
        else
           rank_x = nlist+1
        endif
     endif

     if(cur_rank == rank_y) then
        to_pe = boundy%send(ind_y)%pe
        do n = 1, boundy%send(ind_y)%count
           if(sendy(boundy%send(ind_y)%dir(n))) then
              is = boundy%send(ind_y)%is(n); ie = boundy%send(ind_y)%ie(n)
              js = boundy%send(ind_y)%js(n); je = boundy%send(ind_y)%je(n)
              tMe = boundy%send(ind_y)%tileMe(n)
              do l=1,l_size
                 call pack_buffer_y(boundx%xsize,boundx%ysize,fx_hold(l,tMe)%data,boundy%xsize,boundy%ysize, &
                          fy_hold(l,tMe)%data, buffer,pos, boundy%send(ind_y)%rotation(n), flags, is, ie, js, je, ke)
              end do
           end if ! if(send(boundy%dir(n)))
        end do    ! do n = 1, boundy%count
        ind_y = ind_y+1
        if(ind_y .LE. nsend_y) then
           rank_y = boundy%send(ind_y)%pe - domain%pe 
           if(rank_y .LT.0) rank_y = rank_y + nlist
        else
           rank_y = nlist+1
        endif
     endif
     cur_rank = min(rank_x, rank_y)
     msgsize = pos - buffer_pos
     if( msgsize.GT.0 )then  
        !--- maybe we do not need the following stack size check.
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'MPP_DO_GET_BOUNDARY_V_: mpp_domains_stack overflow, ' // &
                'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
        end if
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_4 )
        buffer_pos = pos
     end if

  end do       

  call mpp_sync_self(check=EVENT_RECV)

  !unpack recv
  !unpack buffer in reverse order.
  buffer_pos = buffer_recv_size  
  cur_rank = get_rank_unpack(domain, boundx, boundy, rank_x, rank_y, ind_x, ind_y) 

  do while(ind_x >0 .OR. ind_y >0)
     if(cur_rank == rank_y) then
        do n = boundy%recv(ind_y)%count, 1, -1
           if(recvy(boundy%recv(ind_y)%dir(n))) then
              is = boundy%recv(ind_y)%is(n); ie = boundy%recv(ind_y)%ie(n)
              js = boundy%recv(ind_y)%js(n); je = boundy%recv(ind_y)%je(n)
              msgsize = (ie-is+1)*(je-js+1)*ke*l_size
              pos = buffer_pos - msgsize
              buffer_pos = pos
              tMe = boundy%recv(ind_y)%tileMe(n)
              dir = boundy%recv(ind_y)%dir(n)
              do l=1,l_size
                 do k = 1, ke
                    index = boundy%recv(ind_y)%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          by_hold(dir,l,tMe)%data(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           end if
        end do
        ind_y = ind_y-1
        if(ind_y .GT. 0) then
           rank_y = boundy%recv(ind_y)%pe - domain%pe 
           if(rank_y .LE.0) rank_y = rank_y + nlist
        else
           rank_y = nlist+1
        endif
     endif

     if(cur_rank == rank_x) then
        do n = boundx%recv(ind_x)%count, 1, -1
           if(recvx(boundx%recv(ind_x)%dir(n))) then
              is = boundx%recv(ind_x)%is(n); ie = boundx%recv(ind_x)%ie(n)
              js = boundx%recv(ind_x)%js(n); je = boundx%recv(ind_x)%je(n)
              msgsize = (ie-is+1)*(je-js+1)*ke*l_size
              pos = buffer_pos - msgsize
              buffer_pos = pos
              tMe = boundx%recv(ind_x)%tileMe(n)
              dir = boundx%recv(ind_x)%dir(n)
              do l=1,l_size
                 do k = 1, ke
                    index = boundx%recv(ind_x)%index(n)
                    do j = js, je
                       do i = is, ie
                          pos = pos + 1
                          bx_hold(dir,l,tMe)%data(index,k) = buffer(pos)
                          index = index + 1
                       end do
                    end do
                 end do
              end do
           end if
        end do
        ind_x = ind_x-1
        if(ind_x .GT. 0) then
           rank_x = boundx%recv(ind_x)%pe - domain%pe 
           if(rank_x .LE.0) rank_x = rank_x + nlist
        else
           rank_x = nlist+1
        endif
     endif
     cur_rank = min(rank_x, rank_y)
  end do

  !--- domain always is symmetry
  shift = 1
  tMe = 1
  if( BTEST(domain%fold,NORTH) .AND. (.NOT.BTEST(flags,SCALAR_BIT)) )then
     j = domain%y(1)%global%end+shift
     if( domain%y(1)%data%begin.LE.j .AND. j.LE.domain%y(1)%data%end+shift )then !fold is within domain
        !poles set to 0: BGRID only
        if( gridtype.EQ.BGRID_NE )then
           midpoint = (domain%x(1)%global%begin+domain%x(1)%global%end-1+shift)/2
           j  = domain%y(1)%global%end+shift - domain%y(1)%compute%begin + 1
           is = domain%x(1)%global%begin; ie = domain%x(1)%global%end+shift
           do i = is ,ie, midpoint
              if( domain%x(1)%compute%begin == i )then
                 do l=1,l_size
                    do k = 1,ke
                       bx_hold(3,l,tMe)%data(j,k) = 0
                       by_hold(3,l,tMe)%data(j,k) = 0
                    end do
                 end do
              end if
           end do
        endif
     endif
  endif

  call mpp_sync_self( )


