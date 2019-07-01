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

  integer :: i, j, k, m, n, nd, nwords, lpos, rpos, from_pe, root_pe, tile_id
  integer :: ke, isc, iec, jsc, jec, is, ie, js, je, nword_me, send_buf_pos, recv_buf_pos
  logical :: global_on_this_pe

  ke = size(local,3)
  isc = domain%x(tile)%compute%begin; iec = domain%x(tile)%compute%end+ishift
  jsc = domain%y(tile)%compute%begin; jec = domain%y(tile)%compute%end+jshift
  nword_me = (iec-isc+1)*(jec-jsc+1)*ke

  ! make contiguous array from compute domain
  m = 0
  global_on_this_pe =  .NOT. root_only .OR. domain%pe == domain%tile_root_pe

  if(global_on_this_pe) then
     global = default_data
     do k = 1, ke
        do j = jsc, jec
           do i = isc, iec
              m = m + 1
              buffer(m) = local(i+ioff,j+joff,k)
              global(i+ipos,j+jpos,k) = buffer(m) !always fill local domain directly
           end do
        end do
     end do
  else
     do k = 1, ke
        do j = jsc, jec
           do i = isc, iec
              m = m + 1
              buffer(m) = local(i+ioff,j+joff,k)
           end do
        end do
     end do
  endif

  ! if there is more than one tile on this pe, then no decomposition for all tiles on this pe, so we can just return
  if(size(domain%x(:))>1) then
     !--- the following is needed to avoid deadlock.
     if( tile == size(domain%x(:)) ) call mpp_sync_self( )
     return
  end if

  root_pe = mpp_root_pe()
  send_buf_pos = 0
  recv_buf_pos = (domain%x(1)%compute%size+ishift)*(domain%y(1)%compute%size+jshift)*size(local,3)
  !fill off-domains (note loops begin at an offset of 1)
  if( xonly )then
     nd = size(domain%x(1)%list(:))
     do n = 1,nd-1
        lpos = mod(domain%x(1)%pos+nd-n,nd)
        rpos = mod(domain%x(1)%pos   +n,nd)
        from_pe = domain%x(1)%list(rpos)%pe
        rpos = from_pe - root_pe ! for concurrent run, root_pe may not be 0.
        if (from_pe == NULL_PE) then
           nwords = 0
        else
           nwords = (domain%list(rpos)%x(1)%compute%size+ishift) &
                  * (domain%list(rpos)%y(1)%compute%size+jshift) * ke
        endif

        call mpp_send(buffer(send_buf_pos+1), plen=nword_me, to_pe=domain%x(1)%list(lpos)%pe, tag=COMM_TAG_1)
        call mpp_recv(buffer(recv_buf_pos+1), glen=nwords, from_pe=from_pe, tag=COMM_TAG_1)
        m = recv_buf_pos
        if (from_pe /= NULL_PE) then
          is = domain%list(rpos)%x(1)%compute%begin; ie = domain%list(rpos)%x(1)%compute%end+ishift
          do k = 1, ke
             do j = jsc, jec
                do i = is, ie
                   m = m + 1
                   global(i+ipos,j+jpos,k) = buffer(m)
                end do
             end do
          end do
          call mpp_sync_self()  !-ensure MPI_ISEND is done.
        endif
     end do
  else if( yonly )then
     nd = size(domain%y(1)%list(:))
     do n = 1,nd-1
        lpos = mod(domain%y(1)%pos+nd-n,nd)
        rpos = mod(domain%y(1)%pos   +n,nd)
        from_pe = domain%y(1)%list(rpos)%pe
        rpos = from_pe - root_pe
        if (from_pe == NULL_PE) then
           nwords = 0
        else
           nwords = (domain%list(rpos)%x(1)%compute%size+ishift) &
                  * (domain%list(rpos)%y(1)%compute%size+jshift) * ke
        endif
        ! Force use of scalar, integer pointer interface
        call mpp_send(buffer(send_buf_pos+1), plen=nword_me, to_pe=domain%y(1)%list(lpos)%pe, tag=COMM_TAG_1)
        call mpp_recv(buffer(recv_buf_pos+1), glen=nwords, from_pe=from_pe, tag=COMM_TAG_1)
        m = recv_buf_pos
        if (from_pe /= NULL_PE) then
          js = domain%list(rpos)%y(1)%compute%begin; je = domain%list(rpos)%y(1)%compute%end+jshift
          do k = 1,ke
             do j = js, je
                do i = isc, iec
                   m = m + 1
                   global(i+ipos,j+jpos,k) = buffer(m)
                end do
             end do
          end do
        endif
        call mpp_sync_self()  !-ensure MPI_ISEND is done.
     end do
  else
     tile_id = domain%tile_id(1)
     nd = size(domain%list(:))
     if(root_only) then
        if(domain%pe .NE. domain%tile_root_pe) then
           call mpp_send( buffer(send_buf_pos+1), plen=nword_me, to_pe=domain%tile_root_pe, tag=COMM_TAG_1 )
        else
           do n = 1,nd-1
              rpos = mod(domain%pos+n,nd)
              if( domain%list(rpos)%tile_id(1) .NE. tile_id ) cycle
              nwords = (domain%list(rpos)%x(1)%compute%size+ishift) * (domain%list(rpos)%y(1)%compute%size+jshift) * ke
              call mpp_recv(buffer(recv_buf_pos+1), glen=nwords, from_pe=domain%list(rpos)%pe, tag=COMM_TAG_1 )
              m = recv_buf_pos
              is = domain%list(rpos)%x(1)%compute%begin; ie = domain%list(rpos)%x(1)%compute%end+ishift
              js = domain%list(rpos)%y(1)%compute%begin; je = domain%list(rpos)%y(1)%compute%end+jshift

              do k = 1,ke
                 do j = js, je
                    do i = is, ie
                       m = m + 1
                       global(i+ipos,j+jpos,k) = buffer(m)
                    end do
                 end do
              end do
           end do
        endif
     else
        do n = 1,nd-1
           lpos = mod(domain%pos+nd-n,nd)
           if( domain%list(lpos)%tile_id(1).NE. tile_id ) cycle ! global field only within tile
           call mpp_send( buffer(send_buf_pos+1), plen=nword_me, to_pe=domain%list(lpos)%pe, tag=COMM_TAG_2 )
        end do
        do n = 1,nd-1
           rpos = mod(domain%pos+n,nd)
           if( domain%list(rpos)%tile_id(1) .NE. tile_id ) cycle ! global field only within tile
           nwords = (domain%list(rpos)%x(1)%compute%size+ishift) * (domain%list(rpos)%y(1)%compute%size+jshift) * ke
           call mpp_recv( buffer(recv_buf_pos+1), glen=nwords, from_pe=domain%list(rpos)%pe, tag=COMM_TAG_2 )
           m = recv_buf_pos
           is = domain%list(rpos)%x(1)%compute%begin; ie = domain%list(rpos)%x(1)%compute%end+ishift
           js = domain%list(rpos)%y(1)%compute%begin; je = domain%list(rpos)%y(1)%compute%end+jshift

           do k = 1,ke
              do j = js, je
                 do i = is, ie
                    m = m + 1
                    global(i+ipos,j+jpos,k) = buffer(m)
                 end do
              end do
           end do
        end do
     endif
  end if

  call mpp_sync_self()
