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

  integer :: isc, iec, jsc, jec, gxsize, gysize, ishift, jshift, i, j, n, m
  integer :: ioffset, joffset, ntile, nlist, list, is, ie, js, je

  call mpp_get_domain_shift(domain, ishift, jshift, position)
  call mpp_get_compute_domain( domain, is,  ie,  js,  je,  tile_count = tile )
  call mpp_get_compute_domain( domain, isc, iec, jsc, jec, tile_count = tile, position=position )
  call mpp_get_global_domain(domain, xsize = gxsize, ysize = gysize )
  ntile = size(domain%x(:))
  nlist = size(domain%list(:))
  global_sum = 0
  if( global_flag == BITWISE_EXACT_SUM )then
     !this is bitwise exact across different PE counts.
     allocate( field2D (isc:iec,jsc:jec) )
     do j = jsc, jec
        do i = isc, iec
           field2D(i,j) = sum( field(i+ioff:i+ioff,j+joff:j+joff,: ) )
        end do
     end do
     allocate( global2D( gxsize+ishift, gysize+jshift ) )
     global2D = 0
     call global_field_2d_private( domain, field2D, global2D, position=position, tile_count=tile)
      
     ioffset = domain%x(tile)%goffset*ishift; joffset = domain%y(tile)%goffset*jshift
     mygsum(tile) = sum(global2D(1:gxsize+ioffset,1:gysize+joffset))
     deallocate(global2D, field2d)
     if( tile == ntile) then 
        if(domain%ntiles == 1 ) then
           global_sum = mygsum(tile)
        else if( nlist == 1) then
           global_sum = sum(mygsum(1:ntile))
        else ! need to sum by the order of tile_count
           ! first fill the global sum on current pe.
           do n = 1, ntile
              gsum(domain%tile_id(n)) = mygsum(n)
           end do
           !--- send the data to other pe if the current pe is the root pe of any tile
           if( domain%pe == domain%tile_root_pe ) then
              do list = 1, nlist - 1
                 m = mod( domain%pos+list, nlist )
                 call mpp_send( mygsum(1), plen=ntile, to_pe=domain%list(m)%pe, tag=COMM_TAG_1 )
              end do
           end if
           call mpp_sync_self()
           !--- receive data from root_pe of each tile
           do list = 1, nlist - 1
              m = mod( domain%pos+nlist-list, nlist )
              if( domain%list(m)%pe == domain%list(m)%tile_root_pe ) then
                  call mpp_recv( nbrgsum(1), glen=size(domain%list(m)%x(:)), from_pe=domain%list(m)%pe, tag=COMM_TAG_1)
                  do n = 1, size(domain%list(m)%x(:))
                     gsum(domain%list(m)%tile_id(n)) = nbrgsum(n)
                  end do
              end if
           end do

           global_sum = sum(gsum(1:domain%ntiles))
        end if
     end if
  else if ( global_flag == BITWISE_EFP_SUM )then
     !this is bitwise across different PE counts using EFP sum
     if( ntile > 1 ) then 
        call mpp_error( FATAL, 'MPP_GLOBAL_SUM_: multiple tile per pe is not supported for BITWISE_EFP_SUM')
     endif
     allocate( field2D (isc:iec,jsc:jec) )
     do j = jsc, jec
        do i = isc, iec
           field2D(i,j) = sum( field(i+ioff:i+ioff,j+joff:j+joff,:) )
        end do
     end do
     !--- using efp sum.
     if(efp_sum_overflow_check) then
        global_sum = mpp_reproducing_sum(field2D, overflow_check=.true.)
     else
        global_sum = mpp_reproducing_sum(field2D, overflow_check=overflow_check)
     endif
  else  !this is not bitwise-exact across different PE counts
     ioffset = domain%x(tile)%loffset*ishift; joffset = domain%y(tile)%loffset*jshift
     mygsum(tile) = sum( field(is+ioff:ie+ioff+ioffset, js+joff:je+joff+joffset,:) )
     if(tile == ntile) then
        global_sum = sum(mygsum(1:ntile))
        call mpp_sum(global_sum, domain%list(:)%pe )
     end if
  end if

