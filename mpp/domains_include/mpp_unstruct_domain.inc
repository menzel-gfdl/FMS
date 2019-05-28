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

  !#####################################################################
  subroutine mpp_define_unstruct_domain(UG_domainObj, SG_domainObj, npts_tile, grid_nlev, ndivs, &
                                        npes_io_group, grid_index, name)
     type(domainUG),          intent(inout) :: UG_domainObj
     type(domain2d),             intent(in) :: SG_domainObj 
     integer,                    intent(in) :: npts_tile(:) ! number of unstructured points on each tile
     integer,                    intent(in) :: grid_nlev(:) ! number of levels in each unstructured grid.
     integer,                    intent(in) :: ndivs       
     integer,                    intent(in) :: npes_io_group  ! number of processors in a io group. Only pe with same tile_id 
                                                              ! in the same group
     integer,                    intent(in) :: grid_index(:)
     character(len=*), optional, intent(in) :: name
     integer, dimension(size(npts_tile(:))) :: ndivs_tile, pe_start, pe_end
     integer, dimension(0:ndivs-1)          :: ibegin, iend, costs_list
     integer :: ntiles, ntotal_pts, ndivs_used, max_npts, cur_tile, cur_npts
     integer :: n, ts, te, p, pos, tile_id, ngroup, group_id, my_pos, i
     integer :: npes_in_group, is, ie, ntotal_costs, max_cost, cur_cost, costs_left
     integer :: npts_left, ndiv_left, cur_pos, ndiv, prev_cost, ioff
     real    :: avg_cost
     integer :: costs(size(npts_tile(:)))
     type(domain2d_private), pointer :: SG_domain
     type(domainUG_private), pointer :: UG_domain
     type(domainUG_private), pointer :: io_domain     
     type(domainUG) :: io_domainObj

     if(UG_domainObj%index .NE. DEFAULT_DOMAINUG_INDEX) &
        call mpp_error(FATAL, "mpp_define_unstruct_domain: UG_domainObj%index should be DEFAULT_DOMAINUG_INDEX")
     cur_domainUG_index = cur_domainUG_index+1
     UG_domainObj%index = cur_domainUG_index
     call check_domainUG_index(UG_domainObj, "mpp_define_unstruct_domain")
     UG_domain => domainUGList(cur_domainUG_index)
     SG_domain => domainList(SG_domainObj%index)

     UG_domain%SG_domain_index = SG_domainObj%index
     ntiles = size(npts_tile(:))
     UG_domain%ntiles = ntiles

     !--- total number of points must be no less than ndivs
     if(sum(npts_tile)<ndivs) call mpp_error(FATAL, "mpp_define_unstruct_domain: total number of points is less than ndivs")
     !--- We are assuming nlev on each grid is at least one.
     do n = 1, size(grid_nlev(:))
        if(grid_nlev(n) < 1) call mpp_error(FATAL, "mpp_define_unstruct_domain: grid_nlev at some point is less than 1")
     enddo

     !-- costs for each tile.
     pos = 0
     do n = 1, ntiles
        costs(n) = 0
        do i = 1, npts_tile(n)
           pos = pos + 1
           costs(n) = costs(n) + grid_nlev(pos)
        enddo
     enddo
     ! compute number of divisions for each tile.
     ntotal_costs = sum(costs)
     !--- get the upper limit of ndivs for each tile.
     do n = 1, ntiles
        ndivs_tile(n) = CEILING(real(costs(n)*ndivs)/ntotal_costs)
     enddo

     ndivs_used = sum(ndivs_tile)
     do while (ndivs_used > ndivs)
        max_cost = 0
        cur_tile = 0
        do n = 1, ntiles
           if( ndivs_tile(n) > 1 ) then
              cur_cost = CEILING(real(costs(n))/(ndivs_tile(n)-1))
              if( max_cost == 0 .OR. cur_cost<max_cost) then
                 max_cost = cur_cost
                 cur_tile = n
              endif
           endif
        enddo
        ndivs_used = ndivs_used-1
        ndivs_tile(cur_tile) = ndivs_tile(cur_tile) - 1    
     enddo

     te = -1
     ioff = 0
     do n = 1, ntiles
        ts = te + 1
        te = ts + ndivs_tile(n) - 1
        costs_left = costs(n)
        ndiv_left = ndivs_tile(n)
        npts_left = npts_tile(n)
        cur_pos = 1
        do ndiv = 1, ndivs_tile(n)
           cur_cost = 0
           ibegin(ts+ndiv-1) = cur_pos
           avg_cost = real(costs_left)/ndiv_left
           do i = cur_pos, npts_tile(n)
              cur_cost = cur_cost + grid_nlev(i+ioff)
              costs_left = costs_left - grid_nlev(i+ioff)
              if(npts_left < ndiv_left ) then
                 call mpp_error(FATAL, "mpp_define_unstruct_domain: npts_left < ndiv_left")
              else if(npts_left == ndiv_left ) then
                 cur_pos = i + 1
                 exit
              else if(cur_cost .GE. avg_cost) then
                 prev_cost = cur_cost - grid_nlev(i+ioff)
                 if(i==cur_pos) then
                    cur_pos = i + 1
                    exit
                 else if( cur_cost - avg_cost .LE. avg_cost - prev_cost ) then
                    cur_pos = i + 1
                    exit
                 else
                    cur_pos = i
                    cur_cost = prev_cost
                    costs_left = costs_left + grid_nlev(i+ioff)
                    npts_left = npts_left+1
                    exit
                 endif
              endif
              npts_left = npts_left-1
           enddo
           iend(ts+ndiv-1) = cur_pos - 1
           costs_list(ts+ndiv-1) = cur_cost
           ndiv_left = ndiv_left-1
           npts_left = npts_left-1
        enddo
        pe_start(n) = ts
        pe_end(n) = te
        ioff = ioff+ npts_tile(n)
     enddo

     allocate(UG_domain%list(0:ndivs-1))
     do p = 0, ndivs-1
        UG_domain%list(p)%compute%begin = ibegin(p)
        UG_domain%list(p)%compute%end = iend(p)
        UG_domain%list(p)%compute%size = UG_domain%list(p)%compute%end - UG_domain%list(p)%compute%begin + 1
        UG_domain%list(p)%compute%max_size = 0
        UG_domain%list(p)%pos = p
        UG_domain%list(p)%pe = p + mpp_root_pe()
        pos = 0
        do n = 1, ntiles
           if( p .GE. pe_start(n) .AND. p .LE. pe_end(n) ) then
              UG_domain%list(p)%tile_id = n
              exit
           endif
           pos = pos + npts_tile(n)
        enddo
        is = UG_domain%list(p)%compute%begin+pos
        ie = UG_domain%list(p)%compute%end+pos
        UG_domain%list(p)%compute%begin_index = minval(grid_index(is:ie))
        UG_domain%list(p)%compute%end_index = maxval(grid_index(is:ie))
     enddo

     !--- write out domain decomposition from root pe
     if(mpp_pe() == mpp_root_pe() .and. present(name)) then
        write(stdout(),*) "unstruct domain name = ", trim(name)
        write(stdout(),*) UG_domain%list(:)%compute%size
     endif

     pos = mpp_pe() - mpp_root_pe() 
     UG_domain%pe = mpp_pe()
     UG_domain%pos = pos
     UG_domain%tile_id = UG_domain%list(pos)%tile_id
     p = pe_start(UG_domain%tile_id)
     UG_domain%tile_root_pe = UG_domain%list(p)%pe
     UG_domain%tile_npes = pe_end(UG_domain%tile_id) - pe_start(UG_domain%tile_id) + 1
     UG_domain%compute = UG_domain%list(pos)%compute
     UG_domain%compute%max_size = MAXVAL( UG_domain%list(:)%compute%size )
     UG_domain%global%begin = 1
     UG_domain%global%end   = npts_tile(UG_domain%tile_id)
     UG_domain%global%size  = UG_domain%global%end - UG_domain%global%begin + 1
     UG_domain%global%max_size = -1   ! currently this is not supposed to be used.
     pos = 0
     do n = 1, UG_domain%tile_id-1
        pos = pos + npts_tile(n)
     enddo
     UG_domain%global%begin_index = grid_index(pos+1)
     UG_domain%global%end_index = grid_index(pos+npts_tile(n))

     allocate(UG_domain%grid_index(UG_domain%compute%size))
     do n = 1, UG_domain%compute%size
        UG_domain%grid_index(n) = grid_index(pos+UG_domain%compute%begin+n-1)
     enddo

     !--- define io_domain
     if(UG_domain%io_domain_index .NE. 0) call mpp_error(FATAL, 'mpp_define_unstruct_domain: UG_domain%io_domain_index .NE. 0')
     cur_domainUG_index = cur_domainUG_index+1
     UG_domain%io_domain_index = cur_domainUG_index
     io_domainObj%index = cur_domainUG_index
     call check_domainUG_index(io_domainObj, "mpp_define_unstruct_domain io_domain")
     io_domain => domainUGList(cur_domainUG_index)

     tile_id = UG_domain%tile_id
     io_domain%pe = UG_domain%pe
     !--- figure out number groups for current tile
     if(npes_io_group == 0) then
        ngroup = 1
     else
        ngroup = CEILING(real(ndivs_tile(tile_id))/ npes_io_group)
     endif

!----------
!ug support
     UG_domain%npes_io_group = npes_io_group
     UG_domain%io_layout = ngroup
!----------

     call mpp_compute_extent(1, ndivs_tile(tile_id), ngroup, ibegin(0:ngroup-1), iend(0:ngroup-1))
     my_pos = UG_domain%pe - UG_domain%tile_root_pe + 1
     do n = 0, ngroup-1
        if( my_pos .GE. ibegin(n) .AND. my_pos .LE. iend(n) ) then
           group_id = n
           exit
        endif
     enddo
     
     io_domain%tile_id            = group_id+1
     io_domain%compute            = UG_domain%compute
     io_domain%pe                 = UG_domain%pe
     io_domain%pos                = my_pos - ibegin(group_id) + 1
     io_domain%tile_root_pe       = ibegin(group_id) + UG_domain%tile_root_pe - 1
     pos = io_domain%tile_root_pe - mpp_root_pe() 
     io_domain%global%begin       = UG_domain%list(pos)%compute%begin
     io_domain%global%begin_index = UG_domain%list(pos)%compute%begin_index
     pos = iend(group_id) + UG_domain%tile_root_pe - mpp_root_pe() - 1
     io_domain%global%end         = UG_domain%list(pos)%compute%end
     io_domain%global%end_index   = UG_domain%list(pos)%compute%end_index
     io_domain%global%size        = io_domain%global%end -  io_domain%global%begin + 1
     
     npes_in_group = iend(group_id) - ibegin(group_id) + 1
     allocate(io_domain%list(0:npes_in_group-1))
     do n = 0, npes_in_group-1
        pos = io_domain%tile_root_pe - mpp_root_pe() + n
        io_domain%list(n)%compute = UG_domain%list(pos)%compute
        io_domain%list(n)%pos = n
        io_domain%list(n)%pe = UG_domain%list(pos)%pe
        io_domain%list(n)%tile_id = group_id+1
     enddo    
 
     call compute_overlap_SG2UG(UG_domain, SG_domain)
     call compute_overlap_UG2SG(UG_domain)

     return

  end subroutine mpp_define_unstruct_domain


  !####################################################################
  subroutine compute_overlap_SG2UG(UG_domain, SG_domain)
     type(domainUG_private),   intent(inout) :: UG_domain
     type(domain2d_private),      intent(in) :: SG_domain
     integer, dimension(0:size(SG_domain%list(:))-1) :: send_cnt, recv_cnt
     integer, dimension(0:size(SG_domain%list(:))-1) :: send_buffer_pos, recv_buffer_pos
     integer, dimension(:), allocatable              :: send_buffer, recv_buffer, index_list
     integer, dimension(:), allocatable              :: buffer_pos
     integer :: tile_id, nlist, nxg, begin_index, end_index, i, j
     integer :: m, n, list, l, isc, iec, jsc, jec, ibegin, iend, grid_index
     integer :: nrecv, nsend, send_pos, recv_pos, pos

     !--- figure out the recv index information.
     tile_id = UG_domain%tile_id
     nlist = size(SG_domain%list(:))
     nxg = SG_domain%x(1)%global%size
     begin_index = UG_domain%compute%begin_index
     end_index = UG_domain%compute%end_index
     pos = 0
     recv_cnt = 0
     allocate(index_list(UG_domain%compute%size))
     allocate(send_buffer(UG_domain%compute%size))
     index_list = -1
     do n = 0, nlist-1
        if(SG_domain%list(n)%tile_id(1) .NE. tile_id) cycle
        isc = SG_domain%list(n)%x(1)%compute%begin; iec = SG_domain%list(n)%x(1)%compute%end
        jsc = SG_domain%list(n)%y(1)%compute%begin; jec = SG_domain%list(n)%y(1)%compute%end
        ibegin = (jsc-1)*nxg + isc
        iend   = (jec-1)*nxg + iec
        if(ibegin > end_index .OR. iend < begin_index) cycle
        do l = 1, UG_domain%compute%size
           grid_index = UG_domain%grid_index(l)
           i = mod((grid_index-1), nxg) + 1
           j = (grid_index-1)/nxg + 1 
           if( i .GE. isc .AND. i .LE. iec .and. j .GE. jsc .AND. j .LE. jec ) then
              recv_cnt(n) = recv_cnt(n) + 1
              pos = pos + 1
              if(pos > UG_domain%compute%size) call mpp_error(FATAL, &
                  'compute_overlap_SG2UG: pos > UG_domain%compute%size')
              index_list(pos) = l
              send_buffer(pos) = grid_index
           endif
        enddo
     enddo

     !--- make sure sum(recv_cnt) == UG_domain%compute%size
     if( UG_domain%compute%size .NE. sum(recv_cnt) ) then
        print*,"pe=", mpp_pe(), UG_domain%compute%size, sum(recv_cnt)
        call mpp_error(FATAL, &
          "compute_overlap_SG2UG: UG_domain%compute%size .NE. sum(recv_cnt)")
     endif
     allocate(buffer_pos(0:nlist-1))
     pos = 0
     do list = 0,nlist-1
        buffer_pos(list) = pos
        pos = pos + recv_cnt(list)
     enddo

     nrecv = count( recv_cnt > 0 )
     UG_domain%SG2UG%nrecv = nrecv
     allocate(UG_domain%SG2UG%recv(nrecv))
     nrecv = 0
     pos = 0
     do list = 0,nlist-1
        m = mod( SG_domain%pos+nlist-list, nlist )
        if( recv_cnt(m) > 0 ) then
           nrecv = nrecv+1
           UG_domain%SG2UG%recv(nrecv)%count = recv_cnt(m)
           UG_domain%SG2UG%recv(nrecv)%pe = UG_domain%list(m)%pe
           allocate(UG_domain%SG2UG%recv(nrecv)%i(recv_cnt(m)))
           pos = buffer_pos(m)
           do l = 1, recv_cnt(m)
              pos = pos + 1
              UG_domain%SG2UG%recv(nrecv)%i(l) = index_list(pos)
           enddo
        endif
     enddo

     !--- figure out the send index information.
     send_cnt = recv_cnt
     recv_cnt = 0
     call mpp_alltoall(send_cnt,1,recv_cnt,1)
     !--- make sure sum(send_cnt) == UG_domain%compute%size
     if( UG_domain%compute%size .NE. sum(send_cnt) ) call mpp_error(FATAL, &
          "compute_overlap_SG2UG: UG_domain%compute%size .NE. sum(send_cnt)")
     allocate(recv_buffer(sum(recv_cnt)))
     send_buffer_pos = 0; recv_buffer_pos = 0
     send_pos = 0; recv_pos = 0
     do n = 0, nlist-1
        if(send_cnt(n) > 0) then
           send_buffer_pos(n) = send_pos
           send_pos = send_pos + send_cnt(n)
        endif
        if(recv_cnt(n) > 0) then
           recv_buffer_pos(n) = recv_pos 
           recv_pos = recv_pos + recv_cnt(n)
        endif
     enddo    

     call mpp_alltoall(send_buffer, send_cnt, send_buffer_pos, &
                       recv_buffer, recv_cnt, recv_buffer_pos)

     nsend = count( recv_cnt(:) > 0 )
     UG_domain%SG2UG%nsend = nsend
     allocate(UG_domain%SG2UG%send(nsend))
     nsend = 0
     isc = SG_domain%x(1)%compute%begin
     jsc = SG_domain%y(1)%compute%begin
     do list = 0,nlist-1
        m = mod( SG_domain%pos+list, nlist )
        if( recv_cnt(m) > 0 ) then
           nsend = nsend+1
           UG_domain%SG2UG%send(nsend)%count = recv_cnt(m)
           UG_domain%SG2UG%send(nsend)%pe = UG_domain%list(m)%pe
           allocate(UG_domain%SG2UG%send(nsend)%i(recv_cnt(m)))
           allocate(UG_domain%SG2UG%send(nsend)%j(recv_cnt(m)))
           pos = recv_buffer_pos(m)
           do l = 1, recv_cnt(m)
              grid_index = recv_buffer(pos+l)
              UG_domain%SG2UG%send(nsend)%i(l) = mod(grid_index-1,nxg) + 1
              UG_domain%SG2UG%send(nsend)%j(l) = (grid_index-1)/nxg + 1
           enddo
        endif
     enddo
     deallocate(send_buffer, recv_buffer, index_list, buffer_pos)

return

  end subroutine compute_overlap_SG2UG

  !####################################################################
  subroutine compute_overlap_UG2SG(UG_domain)
     type(domainUG_private),   intent(inout) :: UG_domain

     !--- UG2SG is the reverse of SG2UG
     UG_domain%UG2SG%nsend = UG_domain%SG2UG%nrecv
     UG_domain%UG2SG%send => UG_domain%SG2UG%recv
     UG_domain%UG2SG%nrecv = UG_domain%SG2UG%nsend
     UG_domain%UG2SG%recv => UG_domain%SG2UG%send

     return

  end subroutine compute_overlap_UG2SG

  !####################################################################
  subroutine mpp_get_UG_SG_domain(UG_domain,SG_domain)
     type(domainUG), intent(inout) :: UG_domain
     type(domain2d), intent(out)   :: SG_domain

     call check_domainUG_index(UG_domain, "mpp_get_UG_SG_domain")
     SG_domain%index = domainUGList(UG_domain%index)%sg_domain_index

     return

  end subroutine mpp_get_UG_SG_domain

  !####################################################################
  function mpp_get_UG_io_domain(domain)
     type(domainUG), intent(in) :: domain
     type(domainUG)             :: mpp_get_UG_io_domain

     call check_domainUG_index(domain, "mpp_get_UG_io_domain")
     if(domainUGList(domain%index)%io_domain_index == 0) &
         call mpp_error(FATAL, "mpp_get_UG_io_domain: io_domain is not defined, contact developer")
     mpp_get_UG_io_domain%index = domainUGList(domain%index)%io_domain_index

  end function mpp_get_UG_io_domain

  !#####################################################################
  subroutine mpp_get_UG_compute_domain( domainObj, begin, end, size)
    type(domainUG),  intent(in) :: domainObj
    integer, intent(out), optional :: begin, end, size
    type(domainUG_private), pointer :: domain=> NULL()

    call check_domainUG_index(domainObj, "mpp_get_UG_compute_domain")
    domain => domainUGList(domainObj%index)
    if( PRESENT(begin)     )begin     = domain%compute%begin
    if( PRESENT(end)       )end       = domain%compute%end
    if( PRESENT(size)      )size      = domain%compute%size
    return
  end subroutine mpp_get_UG_compute_domain

  !#####################################################################
  subroutine mpp_get_UG_global_domain( domainObj, begin, end, size)
    type(domainUG),  intent(in) :: domainObj
    integer, intent(out), optional :: begin, end, size
    type(domainUG_private), pointer :: domain=> NULL()

    call check_domainUG_index(domainObj, "mpp_get_UG_global_domain")
    domain => domainUGList(domainObj%index)
    if( PRESENT(begin)     )begin     = domain%global%begin
    if( PRESENT(end)       )end       = domain%global%end
    if( PRESENT(size)      )size      = domain%global%size
    return
  end subroutine mpp_get_UG_global_domain

  !#####################################################################
  subroutine mpp_get_UG_compute_domains( domainObj, begin, end, size )
    type(domainUG),                   intent(in) :: domainObj
    integer, intent(out), optional, dimension(:) :: begin, end, size
    type(domainUG_private), pointer :: domain=> NULL()

    call check_domainUG_index(domainObj, "mpp_get_UG_compute_domains")
    domain => domainUGList(domainObj%index)

    !we use shape instead of size for error checks because size is used as an argument
    if( PRESENT(begin) )then
       if( any(shape(begin).NE.shape(domain%list)) ) &
            call mpp_error( FATAL, 'mpp_get_UG_compute_domains: begin array size does not match domain.' )
       begin(:) = domain%list(:)%compute%begin
    end if
    if( PRESENT(end) )then
       if( any(shape(end).NE.shape(domain%list)) ) &
            call mpp_error( FATAL, 'mpp_get_UG_compute_domains: end array size does not match domain.' )
            end(:) = domain%list(:)%compute%end
    end if
    if( PRESENT(size) )then
       if( any(shape(size).NE.shape(domain%list)) ) &
           call mpp_error( FATAL, 'mpp_get_UG_compute_domains: size array size does not match domain.' )
       size(:) = domain%list(:)%compute%size
    end if
    return
  end subroutine mpp_get_UG_compute_domains

  !#####################################################################
  subroutine mpp_get_UG_domains_index( domainObj, begin, end)
    type(domainUG),         intent(in) :: domainObj
    integer, intent(out), dimension(:) :: begin, end
    type(domainUG_private), pointer :: domain=> NULL()

    call check_domainUG_index(domainObj, "mpp_get_UG_domains_index")
    domain => domainUGList(domainObj%index)

    !we use shape instead of size for error checks because size is used as an argument
    if( any(shape(begin).NE.shape(domain%list)) ) &
         call mpp_error( FATAL, 'mpp_get_UG_compute_domains: begin array size does not match domain.' )
    begin(:) = domain%list(:)%compute%begin_index
    if( any(shape(end).NE.shape(domain%list)) ) &
         call mpp_error( FATAL, 'mpp_get_UG_compute_domains: end array size does not match domain.' )
         end(:) = domain%list(:)%compute%end_index
    return
  end subroutine mpp_get_UG_domains_index

  !#####################################################################
  function mpp_get_UG_domain_ntiles(domain)
    type(domainUG),  intent(in) :: domain
    integer :: mpp_get_UG_domain_ntiles

    call check_domainUG_index(domain, "mpp_get_UG_domain_ntiles")
    mpp_get_UG_domain_ntiles = domainUGList(domain%index)%ntiles
    return
  end function mpp_get_UG_domain_ntiles

  !#######################################################################
  subroutine mpp_get_ug_domain_tile_list(domainObj, tiles)
     type(domainUG), intent(in) :: domainObj
     integer,     intent(inout) :: tiles(:)
     integer                    :: i
     type(domainUG_private), pointer :: domain=> NULL()

     call check_domainUG_index(domainObj, "mpp_get_ug_domain_tile_list")
     domain => domainUGList(domainObj%index)

     if( size(tiles(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'mpp_get_ug_domain_tile_list: tiles array size does not match domain.' )
     do i = 1, size(tiles(:))
        tiles(i) = domain%list(i-1)%tile_id
     end do

  end subroutine mpp_get_ug_domain_tile_list

  !#####################################################################
  function mpp_get_UG_domain_tile_id(domain)
    type(domainUG),  intent(in) :: domain
    integer :: mpp_get_UG_domain_tile_id

    call check_domainUG_index(domain, "mpp_get_UG_domain_tile_id")
    mpp_get_UG_domain_tile_id = domainUGList(domain%index)%tile_id
    return
  end function mpp_get_UG_domain_tile_id

  !####################################################################
  function mpp_get_UG_domain_npes(domain)
     type(domainUG), intent(in) :: domain
     integer :: mpp_get_UG_domain_npes

     call check_domainUG_index(domain, "mpp_get_UG_domain_npes")
     mpp_get_UG_domain_npes = size(domainUGList(domain%index)%list(:))
     return

  end function mpp_get_UG_domain_npes


  !####################################################################
  subroutine mpp_get_UG_domain_pelist( domainObj, pelist)
     type(domainUG), intent(in) :: domainObj
     integer,              intent(out) :: pelist(:)
     type(domainUG_private), pointer :: domain=> NULL()

     call check_domainUG_index(domainObj, "mpp_get_UG_domain_pelist")
     domain => domainUGList(domainObj%index)

     if( size(pelist(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'mpp_get_UG_domain_pelist: pelist array size does not match domain.' )

     pelist(:) = domain%list(:)%pe
     return

  end subroutine mpp_get_UG_domain_pelist

  !###################################################################
  subroutine mpp_get_UG_domain_tile_pe_inf( domainObj, root_pe, npes, pelist)
     type(domainUG),     intent(in) :: domainObj
     integer, optional, intent(out) :: root_pe, npes
     integer, optional, intent(out) :: pelist(:)
     type(domainUG_private), pointer :: domain=> NULL()

     call check_domainUG_index(domainObj, "mpp_get_UG_domain_tile_pe_inf")
     domain => domainUGList(domainObj%index)
     if(present(root_pe)) root_pe = domain%tile_root_pe
     if(present(npes)) root_pe = domain%tile_npes

     if(present(pelist)) then
        if( size(pelist(:)).NE. domain%tile_npes ) &
           call mpp_error( FATAL, 'mpp_get_UG_domain_tile_pe_inf: pelist array size does not match domain.' )
        pelist(:) = domain%list(domain%pos:domain%pos+domain%tile_npes-1)%pe
     endif
     return

  end subroutine mpp_get_UG_domain_tile_pe_inf


  !####################################################################
  subroutine mpp_get_UG_domain_grid_index( domainObj, grid_index)
     type(domainUG), intent(in) :: domainObj
     integer,              intent(out) :: grid_index(:)
     type(domainUG_private), pointer :: domain=> NULL()

     call check_domainUG_index(domainObj, "mpp_get_UG_domain_grid_index")
     domain => domainUGList(domainObj%index)

     if( size(grid_index(:)).NE.size(domain%grid_index(:)) ) &
         call mpp_error( FATAL, 'mpp_get_UG_domain_grid_index: grid_index array size does not match domain.' )

     grid_index(:) = domain%grid_index(:)
     return

  end subroutine mpp_get_UG_domain_grid_index

  !###################################################################
  subroutine mpp_define_null_UG_domain(domainObj)
     type(domainUG), intent(inout) :: domainObj
     type(domainUG_private), pointer :: domain=>NULL();

     domainObj%index = DOMAINUG_INDEX_START
     domain => domainUGList(domainObj%index)
  
     domain%global%begin  = -1; domain%global%end  = -1; domain%global%size = 0
     domain%compute%begin = -1; domain%compute%end = -1; domain%compute%size = 0
     domain%pe = NULL_PE
     domain%ntiles = -1
     domain%pos = -1
     domain%tile_id = -1
     domain%tile_root_pe = -1

  end subroutine mpp_define_null_UG_domain

!##############################################################################
    subroutine mpp_broadcast_domain_ug( domainObj )
!broadcast domain (useful only outside the context of its own pelist)
      type(domainUG), intent(inout) :: domainObj
      integer, allocatable :: pes(:)
      logical :: native         !true if I'm on the pelist of this domain
      integer :: listsize, listpos
      integer :: n
      integer, dimension(7) :: msg, info         !pe and compute domain of each item in list
      integer                :: errunit
      type(domainUG_private), pointer :: domain=>NULL()

      errunit = stderr()
      if( .NOT.module_is_initialized ) &
                 call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN_ug: You must first call mpp_domains_init.' )

!get the current pelist
      allocate( pes(0:mpp_npes()-1) )
      call mpp_get_current_pelist(pes)

!am I part of this domain?
      native = domainObj%index > 0

!set local list size
      if( native )then
          call check_domainUG_index(domainObj, "mpp_broadcast_domain_ug 1")
          domain => domainUGList(domainObj%index)
          listsize = size(domain%list(:))
      else
          cur_domainUG_index = cur_domainUG_index+1
          domainObj%index = cur_domainUG_index
          call check_domainUG_index(domainObj, "mpp_broadcast_domain_ug 2")
          domain => domainUGList(domainObj%index)
          listsize = 0
      end if
      call mpp_max(listsize)

      if( .NOT.native )then
!initialize domain%list and set null values in message
          allocate( domain%list(0:listsize-1) )
          domain%pe = NULL_PE
          domain%pos = -1
          domain%ntiles = -1
          domain%compute%begin =  1
          domain%compute%end   = -1
          domain%compute%begin_index =  1
          domain%compute%end_index   = -1
          domain%global %begin = -1
          domain%global %end   = -1
          domain%tile_id       = -1
          domain%tile_root_pe  = -1
      end if
!initialize values in info
      info(1) = domain%pe
      info(2) = domain%pos
      info(3) = domain%tile_id
      call mpp_get_UG_compute_domain( domainObj, info(4), info(5))
      info(6) = domain%compute%begin_index
      info(7) = domain%compute%end_index
!broadcast your info across current pelist and unpack if needed
      listpos = 0
      do n = 0,mpp_npes()-1
         msg = info
         if( mpp_pe().EQ.pes(n) .AND. debug )write( errunit,* )'PE ', mpp_pe(), 'broadcasting msg ', msg
         call mpp_broadcast( msg, 7, pes(n) )
!no need to unpack message if native
!no need to unpack message from non-native PE
         if( .NOT.native .AND. msg(1).NE.NULL_PE )then
             domain%list(listpos)%pe            = msg(1)
             domain%list(listpos)%pos           = msg(2)
             domain%list(listpos)%tile_id       = msg(3)
             domain%list(listpos)%compute%begin = msg(4)
             domain%list(listpos)%compute%end   = msg(5)
             domain%list(listpos)%compute%begin_index = msg(6)
             domain%list(listpos)%compute%end_index   = msg(7)
             listpos = listpos + 1
             if( debug )write( errunit,* )'PE ', mpp_pe(), 'received domain from PE ', msg(1), 'ls,le=', msg(4:5)
         end if
      end do

    end subroutine mpp_broadcast_domain_ug

!------------------------------------------------------------------------------
function mpp_domain_UG_is_tile_root_pe(domain) result(is_root)

   !<Inputs/Outputs
    type(domainUG),intent(in) :: domain
    logical(INT_KIND)         :: is_root

    call check_domainUG_index(domain, "mpp_domain_UG_is_tile_root_pe")
    if (domainUGList(domain%index)%pe .eq. domainUGList(domain%index)%tile_root_pe) then
        is_root = .true.
    else
        is_root = .false.
    endif

    return
end function mpp_domain_UG_is_tile_root_pe

!------------------------------------------------------------------------------
!HELP: There needs to be a subroutine to return the "io_layout" for
!      an unstructured domain, so I made one.  Someone should check
!      to see if this is correct.
function mpp_get_io_domain_UG_layout(domain) result(io_layout)

   !<Inputs/Outputs
    type(domainUG),intent(in) :: domain
    integer(INT_KIND)         :: io_layout

    call check_domainUG_index(domain, "mpp_get_io_domain_UG_layout")
    io_layout = domainUGList(domain%index)%io_layout

    return
end function


!------------------------------------------------------------------
subroutine deallocate_unstruct_overlap_type(overlap)
  type(unstruct_overlap_type), intent(inout) :: overlap

  if(associated(overlap%i)) deallocate(overlap%i)
  if(associated(overlap%j)) deallocate(overlap%j)

end subroutine deallocate_unstruct_overlap_type

!------------------------------------------------------------------
subroutine deallocate_unstruct_pass_type(passobj)
  type(unstruct_pass_type), intent(inout) :: passobj
  integer :: n

  do n = 1, passobj%nsend
     call deallocate_unstruct_overlap_type(passobj%send(n))
  enddo
  do n = 1, passobj%nrecv
     call deallocate_unstruct_overlap_type(passobj%recv(n))
  enddo

  if(associated(passobj%send)) deallocate(passobj%send)
  if(associated(passobj%recv)) deallocate(passobj%recv)
 
end subroutine deallocate_unstruct_pass_type

!------------------------------------------------------------------
subroutine mpp_deallocate_domainUG(domainObj)

   !<Inputs/Outputs
    type(domainUG),intent(inout) :: domainObj

   !<Local variables
    integer(INT_KIND) :: i !<Loop variable.
    type(domainUG_private), pointer :: domain=>NULL()
    type(domainUG_private), pointer :: io_domain=>NULL()

    call check_domainUG_index(domainObj, "mpp_deallocate_domainUG")
    domain=>domainUGList(domainObj%index)

    if (associated(domain%list)) then
        deallocate(domain%list)
        domain%list => null()
    endif

    if (domain%io_domain_index >0) then
       io_domain => domainUGList(domain%io_domain_index)
        if (associated(io_domain%list)) then
            deallocate(io_domain%list)
            io_domain%list => null()
        endif
        domain%io_domain_index = 0
    endif

    call deallocate_unstruct_pass_type(domain%SG2UG)
    call deallocate_unstruct_pass_type(domain%UG2SG)

    if (associated(domain%grid_index)) then
        deallocate(domain%grid_index)
        domain%grid_index => null()
    endif

    domain%SG_domain_index = 0

    return
end subroutine mpp_deallocate_domainUG

  !###################################################################
  !> Overload the .eq. for UG
  function mpp_domainUG_eq( aObj, bObj )
    logical                    :: mpp_domainUG_eq
    type(domainUG), intent(in) :: aObj, bObj

    call check_domainUG_index(aObj, "mpp_domainUG_eq aObj")
    call check_domainUG_index(bObj, "mpp_domainUG_eq bObj")

    if( aObj%index == bObj%index ) then
       mpp_domainUG_eq = .true.
    else
       mpp_domainUG_eq = mpp_domainUG_eq_private( domainUGList(aObj%index), domainUGList(bObj%index) )
    endif    

    return
  end function mpp_domainUG_eq

  function mpp_domainUG_eq_private( a, b )
    logical                    :: mpp_domainUG_eq_private
    type(domainUG_private), intent(in) :: a, b
    type(domain2D) :: a_domain2d, b_domain2d

    a_domain2d%index = a%sg_domain_index
    b_domain2d%index = b%sg_domain_index

    if (a_domain2d .ne. b_domain2d) then
        mpp_domainUG_eq_private = .false.
        return
    endif

    mpp_domainUG_eq_private = (a%npes_io_group .EQ. b%npes_io_group) .AND. &
                      (a%pos .EQ. b%pos)                     .AND. &
                      (a%ntiles .EQ. b%ntiles)               .AND. &
                      (a%tile_id .EQ. b%tile_id)             .AND. &
                      (a%tile_npes .EQ. b%tile_npes)         .AND. &
                      (a%tile_root_pe .EQ. b%tile_root_pe)

    if(.not. mpp_domainUG_eq_private) return

    mpp_domainUG_eq_private = ( a%compute%begin.EQ.b%compute%begin .AND. &
         a%compute%end  .EQ.b%compute%end   .AND. &
         a%global%begin .EQ.b%global%begin  .AND. &
         a%global%end   .EQ.b%global%end    .AND. &
         a%SG2UG%nsend  .EQ.b%SG2UG%nsend   .AND. &
         a%SG2UG%nrecv  .EQ.b%SG2UG%nrecv   .AND. &
         a%UG2SG%nsend  .EQ.b%UG2SG%nsend   .AND. &
         a%UG2SG%nrecv  .EQ.b%UG2SG%nrecv         &
    )

    return
  end function mpp_domainUG_eq_private


  !> Overload the .ne. for UG
  function mpp_domainUG_ne( a, b )
    logical                    :: mpp_domainUG_ne
    type(domainUG), intent(in) :: a, b

    mpp_domainUG_ne = .NOT. ( a.EQ.b )
    return
  end function mpp_domainUG_ne

!###############################################################################
! This subroutine pass data from unstructured domain2d to domain.
! First only implement for data at grid cell center.
subroutine mpp_pass_SG_to_UG_2D(UG_domain, field_SG, field_UG)
  type(domainUG), intent(in) :: UG_domain
  class(*),    intent(inout) :: field_UG(:)
  class(*),      intent(in) :: field_SG(:,:)
  class(*), pointer :: field2D_UG(:,:) => NULL()
  class(*), pointer :: field3D_SG(:,:,:) => NULL()

  call point_3D_to_2D(size(field_SG,1), size(field_SG,2), field_SG, field3D_SG)
  call point_2D_to_1D(size(field_UG(:)), field_UG, field2D_uG)
  call mpp_pass_SG_to_UG_3D(UG_domain, field3D_SG, field2D_UG)

end subroutine mpp_pass_SG_to_UG_2D

!##############################################################################
SUBROUTINE mpp_pass_SG_to_UG_3D(UG_domainObj, field_SG, field_UG)
  type(domainUG), intent(in) :: UG_domainObj
  class(*),    intent(inout) :: field_UG(:,:)
  class(*),       intent(in) :: field_SG(:,:,:)
  integer :: ioff, joff
  type(domainUG_private), pointer :: UG_domain=>NULL()
  type(domain2D_private), pointer :: SG_domain=>NULL()
 
  call check_domainUG_index(UG_domainObj, "mpp_pass_SG_to_UG_3D") 
  UG_domain => domainUGList(UG_domainObj%index)
  SG_domain => domainList(UG_domain%sg_domain_index)

  !--- check if data is on data or computing domain
  if(size(field_SG,1) .EQ. SG_domain%x(1)%compute%size .AND.  &
     size(field_SG,2) .EQ. SG_domain%y(1)%compute%size) then
     ioff = 1 - SG_domain%x(1)%compute%begin
     joff = 1 - SG_domain%y(1)%compute%begin
  else if(size(field_SG,1) .EQ. SG_domain%x(1)%data%size .AND.  &
          size(field_SG,2) .EQ. SG_domain%y(1)%data%size) then
     ioff = 1 - SG_domain%x(1)%data%begin
     joff = 1 - SG_domain%y(1)%data%begin
  else
     call mpp_error( FATAL, 'mpp_pass_SG_to_UG_3D_: field_SG must match either compute domain or data domain.' )
  endif

  select type(field_SG)
  type is(real(DOUBLE_KIND))
     select type(field_UG)
     type is(real(DOUBLE_KIND))
        call do_pass_SG_to_UG_r8(UG_domain, field_SG, field_UG, ioff, joff)
     class default
        call mpp_error(FATAL, "mpp_pass_SG_to_UG_3D: field_SG=real(DOUBLE_KIND) is not the same type as field_UG")
     end select
  class default
        call mpp_error(FATAL, "mpp_pass_SG_to_UG_3D: only support 4-byte and 8-byte real, contact developer")
  end select

end SUBROUTINE mpp_pass_SG_to_UG_3D

!##############################################################################
subroutine do_pass_SG_to_UG_r8(UG_domain, field_SG, field_UG, ioff, joff)
  type(domainUG_private),       intent(in) :: UG_domain
  real(DOUBLE_KIND), intent(inout) :: field_UG(:,:)
  real(DOUBLE_KIND),    intent(in) :: field_SG(:,:,:)
  integer,              intent(in) :: ioff, joff
  real(DOUBLE_KIND)        :: buffer(mpp_domains_stack_size)
  character(len=8) :: text
  integer :: i, j, k, l, m, ke, from_pe, to_pe
  integer :: buffer_pos, pos, msgsize

#include "do_pass_sg_to_ug.inc"

end subroutine do_pass_SG_to_UG_r8

!##############################################################################
! This subroutine pass data from unstructured domain2d to domain.
! First only implement for data at grid cell center.
subroutine mpp_pass_UG_to_SG_2D(UG_domain, field_UG, field_SG)
  type(domainUG), intent(in) :: UG_domain
  class(*),       intent(in) :: field_UG(:)
  class(*),    intent(inout) :: field_SG(:,:)
  class(*), pointer :: field2D_UG(:,:) => NULL()
  class(*), pointer :: field3D_SG(:,:,:) => NULL()

  call point_3D_to_2D(size(field_SG,1), size(field_SG,2), field_SG, field3D_SG)
  call point_2D_to_1D(size(field_UG(:)), field_UG, field2D_uG)

  call mpp_pass_UG_to_SG_3D(UG_domain, field2D_UG, field3D_SG)

end SUBROUTINE mpp_pass_UG_to_SG_2D

!##############################################################################
SUBROUTINE mpp_pass_UG_to_SG_3D(UG_domainObj, field_UG, field_SG)
  type(domainUG), intent(in) :: UG_domainObj
  class(*),       intent(in) :: field_UG(:,:)
  class(*),    intent(inout) :: field_SG(:,:,:)
  integer :: ioff, joff
  type(domainUG_private), pointer :: UG_domain=>NULL()
  type(domain2D_private), pointer :: SG_domain=>NULL()

  call check_domainUG_index(UG_domainObj, "mpp_pass_UG_to_SG_3D")
  UG_domain => domainUGList(UG_domainObj%index)
  SG_domain => domainList(UG_domain%sg_domain_index)

  !--- check if data is on data or computing domain
  if(size(field_SG,1) .EQ. SG_domain%x(1)%compute%size .AND.  &
     size(field_SG,2) .EQ. SG_domain%y(1)%compute%size) then
     ioff = 1 - SG_domain%x(1)%compute%begin
     joff = 1 - SG_domain%y(1)%compute%begin
  else if(size(field_SG,1) .EQ. SG_domain%x(1)%data%size .AND.  &
          size(field_SG,2) .EQ. SG_domain%y(1)%data%size) then
     ioff = 1 - SG_domain%x(1)%data%begin
     joff = 1 - SG_domain%y(1)%data%begin
  else
     call mpp_error( FATAL, 'mpp_pass_UG_to_SG_3D_: field_SG must match either compute domain or data domain.' )
  endif

  select type(field_SG)
  type is(real(DOUBLE_KIND))
     select type(field_UG)
     type is(real(DOUBLE_KIND))
        call do_pass_UG_to_SG_r8(UG_domain, field_UG, field_SG, ioff, joff)
     class default
        call mpp_error(FATAL, "mpp_pass_UG_to_SG_3D: field_SG=real(DOUBLE_KIND) is not the same type as field_UG")
     end select
  class default
        call mpp_error(FATAL, "mpp_pass_UG_to_SG_3D: only support 4-byte and 8-byte real, contact developer")
  end select



end SUBROUTINE mpp_pass_UG_to_SG_3D

!##############################################################################
subroutine do_pass_UG_to_SG_r8(UG_domain, field_UG, field_SG, ioff, joff)
  type(domainUG_private), intent(in) :: UG_domain
  real(DOUBLE_KIND),    intent(in) :: field_UG(:,:)
  real(DOUBLE_KIND), intent(inout) :: field_SG(:,:,:)
  integer,              intent(in) :: ioff, joff
  real(DOUBLE_KIND)        :: buffer(mpp_domains_stack_size)
  character(len=8) :: text
  integer :: i, j, k, l, m, ke, from_pe, to_pe
  integer :: buffer_pos, pos, msgsize

#include "do_pass_ug_to_sg.inc"

end subroutine do_pass_UG_to_SG_r8

function mpp_domainUG_is_defined(domain)
   type(domainUG), intent(in) :: domain
   logical                    :: mpp_domainUG_is_defined

   mpp_domainUG_is_defined = domain%index > 0

end function mpp_domainUG_is_defined

  !####################################################################
  subroutine check_domainUG_index(domain, msg)
     type(domainUG),   intent(in) :: domain
     character(len=*), intent(in) :: msg

     if(domain%index < DOMAINUG_INDEX_START) call mpp_error(FATAL, trim(msg)//": domain%index < DOMAINUG_INDEX_START")
     if(domain%index > DOMAINUG_INDEX_END) call mpp_error(FATAL, trim(msg)//": domain%index > DOMAINUG_INDEX_END, contact developer" )

  end subroutine check_domainUG_index


