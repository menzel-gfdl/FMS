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

     ! ---northern boundary fold
      shift = 0
      if(domain%symmetry) shift = 1
      ioff = 1- domain%x(1)%data%begin
      joff = 1- domain%y(1)%data%begin
      if( BTEST(domain%fold,NORTH) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then
         isd = domain%x(1)%compute%begin - whalo_v;
         ied = domain%x(1)%compute%end + ehalo_v;
         jsd = domain%y(1)%compute%begin - shalo_v;
         jed = domain%y(1)%compute%end + nhalo_v;
         j = domain%y(1)%global%end+shift
         if( jsd .LE. j .AND. j.LE.jed+shift )then !fold is within domain
            !poles set to 0: BGRID only
            if( gridtype.EQ.BGRID_NE )then
               midpoint = (domain%x(1)%global%begin+domain%x(1)%global%end-1+shift)/2
               j  = domain%y(1)%global%end+shift + joff
               is = domain%x(1)%global%begin; ie = domain%x(1)%global%end+shift
               if( .NOT. domain%symmetry ) is = is - 1
               do i = is ,ie, midpoint
                  if( isd.LE.i .AND. i.LE. ied+shift )then
                     do l=1,l_size
                        do k = 1,ke
                           xHold(l,tMe)%data(i+ioff,j,k) = 0.
                           yHold(l,tMe)%data(i+ioff,j,k) = 0.
                        end do
                     end do
                  end if
               end do
            endif

            ! the following code code block correct an error where the data in your halo coming from 
            ! other half may have the wrong sign
            !off west edge, when update north or west direction
            j = domain%y(1)%global%end+shift + joff
            if ( recv_y(7) .OR. recv_y(5) ) then
               select case(gridtype)
               case(BGRID_NE)
                  if(domain%symmetry) then
                     is = domain%x(1)%global%begin
                  else
                     is = domain%x(1)%global%begin - 1
                  end if
                  if( is.GT.isd )then
                     if( 2*is-domain%x(1)%data%begin.GT.domain%x(1)%data%end+shift ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-north BGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        do k = 1,ke
                           do i = isd,is-1
                              xHold(l,tMe)%data(i+ioff,j,k) = xHold(l,tMe)%data(2*is-i+ioff,j,k)
                              yHold(l,tMe)%data(i+ioff,j,k) = yHold(l,tMe)%data(2*is-i+ioff,j,k)
                           end do
                        end do
                     end do
                  end if
               case(CGRID_NE)
                  is = domain%x(1)%global%begin
                  if( is.GT.isd )then
                     if( 2*is-domain%x(1)%data%begin-1.GT.domain%x(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-north CGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        do k = 1,ke
                           do i = isd,is-1
                              yHold(l,tMe)%data(i+ioff,j,k) = yHold(l,tMe)%data(2*is-i-1+ioff,j,k)
                           end do
                        end do
                     end do
                  end if
               end select
            end if

            !off east edge
            is = domain%x(1)%global%end
            if(domain%x(1)%cyclic .AND. is.LT.ied )then
               ie = ied
               is = is + 1
               select case(gridtype)
               case(BGRID_NE)
                  is = is + shift
                  ie = ie + shift
                  do l=1,l_size
                     do k = 1,ke
                        do i = is,ie
                           xHold(l,tMe)%data(i+ioff,j,k) = -xHold(l,tMe)%data(i+ioff,j,k)
                           yHold(l,tMe)%data(i+ioff,j,k) = -yHold(l,tMe)%data(i+ioff,j,k)
                        end do
                     end do
                  end do
               case(CGRID_NE)
                  do l=1,l_size
                     do k = 1,ke
                        do i = is, ie
                           yHold(l,tMe)%data(i+ioff,j,k) = -yHold(l,tMe)%data(i+ioff,j,k)
                        end do
                     end do
                  end do
               end select
            end if
         end if
      else if( BTEST(domain%fold,SOUTH) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then      ! ---southern boundary fold
         ! NOTE: symmetry is assumed for fold-south boundary
         j = domain%y(1)%global%begin
         if( domain%y(1)%data%begin.LE.j .AND. j.LE.domain%y(1)%data%end+shift )then !fold is within domain
            midpoint = (domain%x(1)%global%begin+domain%x(1)%global%end-1+shift)/2
            !poles set to 0: BGRID only
            if( gridtype.EQ.BGRID_NE )then
               j  = domain%y(1)%global%begin+joff
               is = domain%x(1)%global%begin; ie = domain%x(1)%global%end+shift
               do i = is ,ie, midpoint
                  if( domain%x(1)%data%begin.LE.i .AND. i.LE. domain%x(1)%data%end+shift )then
                     do l=1,l_size
                        do k = 1,ke
                           xHold(l,tMe)%data(i+ioff,j,k) = 0.
                           yHold(l,tMe)%data(i+ioff,j,k) = 0.
                        end do
                     end do
                  end if
               end do
            endif

            ! the following code code block correct an error where the data in your halo coming from 
            ! other half may have the wrong sign
            !off west edge, when update north or west direction
            j = domain%y(1)%global%begin + joff
            if ( recv_x(3) .OR. recv_x(5) ) then
               select case(gridtype)
               case(BGRID_NE)
                  is = domain%x(1)%global%begin
                  if( is.GT.domain%x(1)%data%begin )then

                     if( 2*is-domain%x(1)%data%begin.GT.domain%x(1)%data%end+shift ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-south BGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        do k = 1,ke
                           do i = domain%x(1)%data%begin,is-1
                              xHold(l,tMe)%data(i+ioff,j,k) = xHold(l,tMe)%data(2*is-i+ioff,j,k)
                              yHold(l,tMe)%data(i+ioff,j,k) = yHold(l,tMe)%data(2*is-i+ioff,j,k)
                           end do
                        end do
                     end do
                  end if
               case(CGRID_NE)
                  is = domain%x(1)%global%begin
                  if( is.GT.domain%x(1)%data%begin )then
                     if( 2*is-domain%x(1)%data%begin-1.GT.domain%x(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-south CGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        do k = 1,ke
                           do i = domain%x(1)%data%begin,is-1
                              yHold(l,tMe)%data(i+ioff,j,k) = yHold(l,tMe)%data(2*is-i-1+ioff,j,k)
                           end do
                        end do
                     end do
                  end if
               end select
            end if

            !off east edge
            is = domain%x(1)%global%end
            if(domain%x(1)%cyclic .AND. is.LT.domain%x(1)%data%end )then
               ie = domain%x(1)%data%end
               is = is + 1
               select case(gridtype)
               case(BGRID_NE)
                  is = is + shift
                  ie = ie + shift
                  do l=1,l_size
                     do k = 1,ke
                        do i = is,ie
                           xHold(l,tMe)%data(i+ioff,j,k) = -xHold(l,tMe)%data(i+ioff,j,k)
                           yHold(l,tMe)%data(i+ioff,j,k) = -yHold(l,tMe)%data(i+ioff,j,k)
                        end do
                     end do
                  end do
               case(CGRID_NE)
                  do l=1,l_size
                     do k = 1,ke
                        do i = is, ie
                           yHold(l,tMe)%data(i+ioff,j,k) = -yHold(l,tMe)%data(i+ioff,j,k)
                        end do
                     end do
                  end do
               end select
            end if
         end if
      else if( BTEST(domain%fold,WEST) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then      ! ---eastern boundary fold
         ! NOTE: symmetry is assumed for fold-west boundary
         i = domain%x(1)%global%begin
         if( domain%x(1)%data%begin.LE.i .AND. i.LE.domain%x(1)%data%end+shift )then !fold is within domain
            midpoint = (domain%y(1)%global%begin+domain%y(1)%global%end-1+shift)/2
            !poles set to 0: BGRID only
            if( gridtype.EQ.BGRID_NE )then
               i  = domain%x(1)%global%begin+ioff
               js = domain%y(1)%global%begin; je = domain%y(1)%global%end+shift
               do j = js ,je, midpoint
                  if( domain%y(1)%data%begin.LE.j .AND. j.LE. domain%y(1)%data%end+shift )then
                     do l=1,l_size
                        do k = 1,ke
                           xHold(l,tMe)%data(i,j+joff,k) = 0.
                           yHold(l,tMe)%data(i,j+joff,k) = 0.
                        end do
                     end do
                  end if
               end do
            endif

            ! the following code code block correct an error where the data in your halo coming from 
            ! other half may have the wrong sign
            !off south edge, when update south or west direction
            i = domain%x(1)%global%begin+ioff
            if ( recv_x(3) .OR. recv_x(5) ) then
               select case(gridtype)
               case(BGRID_NE)
                  js = domain%y(1)%global%begin
                  if( js.GT.domain%y(1)%data%begin )then

                     if( 2*js-domain%y(1)%data%begin.GT.domain%y(1)%data%end+shift ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-west BGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        do k = 1,ke
                           do j = domain%y(1)%data%begin,js-1
                              xHold(l,tMe)%data(i,j+joff,k) = xHold(l,tMe)%data(i,2*js-j+joff,k)
                              yHold(l,tMe)%data(i,j+joff,k) = yHold(l,tMe)%data(i,2*js-j+joff,k)
                           end do
                        end do
                     end do
                  end if
               case(CGRID_NE)
                  js = domain%y(1)%global%begin
                  if( js.GT.domain%y(1)%data%begin )then
                     if( 2*js-domain%y(1)%data%begin-1.GT.domain%y(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-west CGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        do k = 1,ke
                           do j = domain%y(1)%data%begin,js-1
                              xHold(l,tMe)%data(i,j+joff,k) = xHold(l,tMe)%data(i, 2*js-j-1+joff,k)
                           end do
                        end do
                     end do
                  end if
               end select
            end if

            !off north edge
            js = domain%y(1)%global%end
            if(domain%y(1)%cyclic .AND. js.LT.domain%y(1)%data%end )then
               je = domain%y(1)%data%end
               js = js + 1
               select case(gridtype)
               case(BGRID_NE)
                  js = js + shift
                  je = je + shift
                  do l=1,l_size
                     do k = 1,ke
                        do j = js,je
                           xHold(l,tMe)%data(i,j+joff,k) = -xHold(l,tMe)%data(i,j+joff,k)
                           yHold(l,tMe)%data(i,j+joff,k) = -yHold(l,tMe)%data(i,j+joff,k)
                        end do
                     end do
                  end do
               case(CGRID_NE)
                  do l=1,l_size
                     do k = 1,ke
                        do j = js, je
                           xHold(l,tMe)%data(i,j+joff,k) = -xHold(l,tMe)%data(i,j+joff,k)
                        end do
                     end do
                  end do
               end select
            end if
         end if
      else if( BTEST(domain%fold,EAST) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then      ! ---eastern boundary fold
         ! NOTE: symmetry is assumed for fold-west boundary
         i = domain%x(1)%global%end+shift
         if( domain%x(1)%data%begin.LE.i .AND. i.LE.domain%x(1)%data%end+shift )then !fold is within domain
            midpoint = (domain%y(1)%global%begin+domain%y(1)%global%end-1+shift)/2
            !poles set to 0: BGRID only
            if( gridtype.EQ.BGRID_NE )then
               i  = domain%x(1)%global%end+shift+ioff
               js = domain%y(1)%global%begin; je = domain%y(1)%global%end+shift
               do j = js ,je, midpoint
                  if( domain%y(1)%data%begin.LE.j .AND. j.LE. domain%y(1)%data%end+shift )then
                     do l=1,l_size
                        do k = 1,ke
                           xHold(l,tMe)%data(i,j+joff,k) = 0.
                           yHold(l,tMe)%data(i,j+joff,k) = 0.
                        end do
                     end do
                  end if
               end do
            endif

            ! the following code code block correct an error where the data in your halo coming from 
            ! other half may have the wrong sign
            !off south edge, when update south or west direction
            i = domain%x(1)%global%end+shift+ioff
            if ( recv_x(3) .OR. recv_x(1) ) then
               select case(gridtype)
               case(BGRID_NE)
                  js = domain%y(1)%global%begin
                  if( js.GT.domain%y(1)%data%begin )then

                     if( 2*js-domain%y(1)%data%begin.GT.domain%y(1)%data%end+shift ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-east BGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        do k = 1,ke
                           do j = domain%y(1)%data%begin,js-1
                              xHold(l,tMe)%data(i,j+joff,k) = xHold(l,tMe)%data(i,2*js-j+joff,k)
                              yHold(l,tMe)%data(i,j+joff,k) = yHold(l,tMe)%data(i,2*js-j+joff,k)
                           end do
                        end do
                     end do
                  end if
               case(CGRID_NE)
                  js = domain%y(1)%global%begin
                  if( js.GT.domain%y(1)%data%begin )then
                     if( 2*js-domain%y(1)%data%begin-1.GT.domain%y(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-east CGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        do k = 1,ke
                           do j = domain%y(1)%data%begin,js-1
                              xHold(l,tMe)%data(i,j+joff,k) = xHold(l,tMe)%data(i, 2*js-j-1+joff,k)
                           end do
                        end do
                     end do
                  end if
               end select
            end if

            !off north edge
            js = domain%y(1)%global%end
            if(domain%y(1)%cyclic .AND. js.LT.domain%y(1)%data%end )then
               je = domain%y(1)%data%end
               js = js + 1
               select case(gridtype)
               case(BGRID_NE)
                  js = js + shift
                  je = je + shift
                  do l=1,l_size
                     do k = 1,ke
                        do j = js,je
                           xHold(l,tMe)%data(i,j+joff,k) = -xHold(l,tMe)%data(i,j+joff,k)
                           yHold(l,tMe)%data(i,j+joff,k) = -yHold(l,tMe)%data(i,j+joff,k)
                        end do
                     end do
                  end do
               case(CGRID_NE)
                  do l=1,l_size
                     do k = 1,ke
                        do j = js, je
                           xHold(l,tMe)%data(i,j+joff,k) = -xHold(l,tMe)%data(i,j+joff,k)
                        end do
                     end do
                  end do
               end select
            end if
         end if
      end if

