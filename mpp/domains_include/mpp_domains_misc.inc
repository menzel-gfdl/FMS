! -*-f90-*-
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                 MPP_DOMAINS: initialization and termination                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! <SUBROUTINE NAME="mpp_domains_init">
!  <OVERVIEW>
!    Initialize domain decomp package.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Called to initialize the <TT>mpp_domains_mod</TT> package.
!    
!    <TT>flags</TT> can be set to <TT>MPP_VERBOSE</TT> to have
!    <TT>mpp_domains_mod</TT> keep you informed of what it's up
!    to. <TT>MPP_DEBUG</TT> returns even more information for debugging.
!    
!    <TT>mpp_domains_init</TT> will call <TT>mpp_init</TT>, to make sure
!    <LINK SRC="mpp.html"><TT>mpp_mod</TT></LINK> is initialized. (Repeated
!    calls to <TT>mpp_init</TT> do no harm, so don't worry if you already
!    called it).
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_domains_init(flags)
!  </TEMPLATE>
!  <IN NAME="flags" TYPE="integer"></IN>
! </SUBROUTINE>
    subroutine mpp_domains_init(flags)
      integer, intent(in), optional :: flags
      integer                       :: unit_begin, unit_end, unit_nml, io_status, unit
      logical                       :: opened

      if( module_is_initialized )return
      call mpp_init(flags)           !this is a no-op if already initialized
      call mpp_pset_init             !this is a no-op if already initialized
      module_is_initialized = .TRUE.
      pe = mpp_root_pe()
      unit = stdlog()
      if( mpp_pe() .EQ.mpp_root_pe() ) write( unit,'(/a)' )'MPP_DOMAINS module '//trim(version)

      if( PRESENT(flags) )then
          debug   = flags.EQ.MPP_DEBUG
          verbose = flags.EQ.MPP_VERBOSE .OR. debug
          domain_clocks_on = flags.EQ.MPP_DOMAIN_TIME
      end if

      !--- namelist 
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, mpp_domains_nml, iostat=io_status)
#else
      unit_begin = 103
      unit_end   = 512
      do unit_nml = unit_begin, unit_end
         inquire( unit_nml,OPENED=opened )
         if( .NOT.opened )exit
      end do

      open(unit_nml,file='input.nml', iostat=io_status)
      read(unit_nml,mpp_domains_nml,iostat=io_status)
      close(unit_nml)
#endif

      if (io_status > 0) then
         call mpp_error(FATAL,'=>mpp_domains_init: Error reading input.nml')
      endif


      select case(lowercase(trim(debug_update_domain)))
      case("none")
         debug_update_level = NO_CHECK
      case("fatal")
         debug_update_level = FATAL
      case("warning")
         debug_update_level = WARNING
      case("note")
         debug_update_level = NOTe
      case default
         call mpp_error(FATAL, "mpp_domains_init: debug_update_level should be 'none', 'fatal', 'warning', or 'note'") 
      end select

      call mpp_domains_set_stack_size(32768) !default, pretty arbitrary

!NULL_DOMAIN is a domaintype that can be used to initialize to undef
      call mpp_define_null_domain(NULL_DOMAIN1d);
      call mpp_define_null_domain(NULL_DOMAIN2d);
      call mpp_define_null_UG_domain(NULL_DOMAINUG)

      if( domain_clocks_on )then
          pack_clock          = mpp_clock_id( 'Halo pack' )
          send_clock          = mpp_clock_id( 'Halo send' )
          recv_clock          = mpp_clock_id( 'Halo recv' )
          unpk_clock          = mpp_clock_id( 'Halo unpk' )
          wait_clock          = mpp_clock_id( 'Halo wait' )
          nest_pack_clock     = mpp_clock_id( 'nest pack' )
          nest_send_clock     = mpp_clock_id( 'nest send' )
          nest_recv_clock     = mpp_clock_id( 'nest recv' )
          nest_unpk_clock     = mpp_clock_id( 'nest unpk' )
          nest_wait_clock     = mpp_clock_id( 'nest wait' )
          group_pack_clock     = mpp_clock_id( 'group pack' )
          group_send_clock     = mpp_clock_id( 'group send' )
          group_recv_clock     = mpp_clock_id( 'group recv' )
          group_unpk_clock     = mpp_clock_id( 'group unpk' )
          group_wait_clock     = mpp_clock_id( 'group wait' )
      end if
      return
    end subroutine mpp_domains_init

!#####################################################################
! <SUBROUTINE NAME="mpp_domains_exit">
!  <OVERVIEW>
!    Exit <TT>mpp_domains_mod</TT>.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Serves no particular purpose, but is provided should you require to
!    re-initialize <TT>mpp_domains_mod</TT>, for some odd reason.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_domains_exit()
!  </TEMPLATE>
! </SUBROUTINE>
    subroutine mpp_domains_exit()
      integer :: unit 
      if( .NOT.module_is_initialized )return
      call mpp_max(mpp_domains_stack_hwm)
      unit = stdout()
      if( mpp_pe().EQ.mpp_root_pe() )write( unit,* )'MPP_DOMAINS_STACK high water mark=', mpp_domains_stack_hwm
      module_is_initialized = .FALSE.
      return
    end subroutine mpp_domains_exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_CHECK_FIELD: Check parallel                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <SUBROUTINE NAME="mpp_check_field_3D" INTERFACE="mpp_check_field">
!   <IN NAME="field_in" TYPE="real, dimension(:,:,:)" > </IN>
!   <IN NAME="pelist1, pelist2" TYPE="integer, dimension(:)" > </IN>
!   <IN NAME="domain" TYPE="type(domain2d)" > </IN>
!   <IN NAME="mesg" TYPE="character(len=*)" > </IN>
!   <IN NAME="w_halo, s_halo, e_halo, n_halo" TYPE="integer, optional" > </IN>
!   <IN NAME="force_abort" TYPE="logical,optional" > </IN>
! </SUBROUTINE>

  subroutine mpp_check_field_3D(field_in, pelist1, pelist2, domain, mesg, &
                                 w_halo, s_halo, e_halo, n_halo, force_abort, position  )
!  This routine is used to do parallel checking for 3d data between n and m pe. The comparison is
!  is done on pelist2. When size of pelist2 is 1, we can check the halo; otherwise,
!  halo can not be checked.

  real, dimension(:,:,:), intent(in)  :: field_in ! field to be checked
  integer, dimension(:), intent(in) :: pelist1, pelist2 ! pe list for the two groups
  type(domain2d),        intent(in) :: domain   ! domain for each pe
  character(len=*),     intent(in)  :: mesg     ! message to be printed out
                                               ! if differences found      
  integer, intent(in), optional     :: w_halo,  s_halo, e_halo, n_halo
                                       ! halo size for west, south, east and north
  logical, intent(in), optional     :: force_abort   ! when true, call mpp_error if any difference
                                                    ! found. default value is false.  
  integer, intent(in), optional     :: position    ! when domain is symmetry, only value = CENTER is 
                                                   ! implemented. 
            
    integer :: k
    character(len=256) :: temp_mesg


    do k = 1, size(field_in,3)
       write(temp_mesg, '(a, i3)') trim(mesg)//" at level " , k
       call mpp_check_field_2d(field_in(:,:,k), pelist1, pelist2, domain, temp_mesg, &
                                 w_halo, s_halo, e_halo, n_halo, force_abort, position )
    enddo
    
  end subroutine mpp_check_field_3D
  
  
!#####################################################################################
! <SUBROUTINE NAME="mpp_check_field_2D" INTERFACE="mpp_check_field">
!   <IN NAME="field_in" TYPE="real, dimension(:,:)" > </IN>
! </SUBROUTINE>
  subroutine mpp_check_field_2d(field_in, pelist1, pelist2, domainObj, mesg, &
                                 w_halo, s_halo, e_halo, n_halo,force_abort, position  )
!  This routine is used to do parallel checking for 2d data between n and m pe. The comparison is
!  is done on pelist2. When size of pelist2 is 1, we can check the halo; otherwise,
!  halo can not be checked.

  real, dimension(:,:), intent(in)  :: field_in ! field to be checked
  integer, dimension(:), intent(in) :: pelist1, pelist2 ! pe list for the two groups
  type(domain2d),        intent(in) :: domainObj   ! domain for each pe
  character(len=*),     intent(in)  :: mesg     ! message to be printed out
                                               ! if differences found
  integer, intent(in), optional     :: w_halo,  s_halo, e_halo, n_halo
                                       ! halo size for west, south, east and north
  logical, intent(in), optional     :: force_abort   ! when, call mpp_error if any difference
                                                    ! found. default value is false.
  integer, intent(in), optional     :: position    ! when domain is symmetry, only value = CENTER is 
                                                   ! implemented. 
  type(domain2d_private), pointer :: domain=>NULL()

  call check_domain2d_index(domainObj, "mpp_check_field_2D")
  domain => domainList(domainObj%index)
     if(present(position)) then
        if(position .NE. CENTER .AND. domain%symmetry) call mpp_error(FATAL,  &
           'mpp_check_field: when domain is symmetry, only value CENTER is implemented, contact author')
     endif

     if(size(pelist2(:)) == 1) then                 
        call check_field_2d_type1(field_in, pelist1, pelist2, domain, mesg, &
                                      w_halo, s_halo, e_halo, n_halo, force_abort )
     else if(size(pelist1(:)) == 1) then
        call check_field_2d_type1(field_in, pelist2, pelist1, domain, mesg, &
                                      w_halo, s_halo, e_halo, n_halo, force_abort )
     else if(size(pelist1(:)) .gt. 1 .and. size(pelist2(:)) .gt. 1) then
        call check_field_2d_type2(field_in, pelist1, pelist2, domain, mesg, force_abort )
     else
        call mpp_error(FATAL, 'mpp_check_field: size of both pelists should be greater than 0')
     endif

  end subroutine mpp_check_field_2D


!####################################################################################

  subroutine check_field_2d_type1(field_in, pelist1, pelist2, domain, mesg, &
                                 w_halo, s_halo, e_halo, n_halo,force_abort  )
!  This routine is used to check field between running on 1 pe (pelist2) and
!  n pe(pelist1). The need_to_be_checked data is sent to the pelist2 and All the
!  comparison is done on pelist2.

  real, dimension(:,:), intent(in)  :: field_in ! field to be checked
  integer, dimension(:), intent(in) :: pelist1, pelist2 ! pe list for the two groups
  type(domain2d_private),  intent(in) :: domain   ! domain for each pe
  character(len=*),     intent(in)  :: mesg     ! message to be printed out
                                               ! if differences found      
  integer, intent(in), optional     :: w_halo,  s_halo, e_halo, n_halo
                                       ! halo size for west, south, east and north
  logical, intent(in), optional     :: force_abort   ! when, call mpp_error if any difference
                                                    ! found. default value is false.         
! some local data

  integer                :: pe,npes, p
  integer                :: hwest, hsouth, heast, hnorth, isg, ieg, jsg, jeg, xhalo, yhalo
  integer                :: i,j,im,jm,l,is,ie,js,je,isc,iec,jsc,jec,isd,ied,jsd,jed
  real,dimension(:,:), allocatable :: field1,field2
  real,dimension(:),   allocatable :: send_buffer
  integer, dimension(4)  ::  ibounds
  logical                :: check_success,  error_exit

  check_success = .TRUE.
  error_exit    = .FALSE.
  if(present(force_abort)) error_exit = force_abort
  hwest  = 0; if(present(w_halo)) hwest  = w_halo
  heast  = 0; if(present(e_halo)) heast  = e_halo
  hsouth = 0; if(present(s_halo)) hsouth = s_halo
  hnorth = 0; if(present(n_halo)) hnorth = n_halo

  pe = mpp_pe ()
  npes = mpp_npes()

  call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
  xhalo = isc - isd
  yhalo = jsc - jsd
  !--- need to checked halo size should not be bigger than x_halo or y_halo
  if(hwest .gt. xhalo .or. heast .gt. xhalo .or. hsouth .gt. yhalo .or. hnorth .gt. yhalo) &
     call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//': The halo size is not correct')

  is = isc - hwest; ie = iec + heast; js = jsc - hsouth; je = jec + hnorth
  allocate(field2(is:ie,js:je))

  ! check if the field_in is on compute domain or data domain
  if((size(field_in,1) .eq. iec-isc+1) .and. (size(field_in,2) .eq. jec-jsc+1)) then
     !if field_in on compute domain, you can not check halo points
     if( hwest .ne. 0 .or. heast .ne. 0 .or. hsouth .ne. 0 .or. hnorth .ne. 0 ) &
         call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//': field is on compute domain, can not check halo')
     field2(:,:) = field_in(:,:)
  else if((size(field_in,1) .eq. ied-isd+1) .and. (size(field_in,2) .eq. jed-jsd+1)) then
     field2(is:ie,js:je) = field_in(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1)
  else if((size(field_in,1) .eq. ieg-isg+1) .and. (size(field_in,2) .eq. jeg-jsg+1)) then
     if( hwest .ne. 0 .or. heast .ne. 0 .or. hsouth .ne. 0 .or. hnorth .ne. 0 ) &
         call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//': field is on compute domain, can not check halo')
     field2(is:ie,js:je) = field_in(1:ie-is+1,1:je-js+1)
  else if((size(field_in,1) .eq. ieg-isg+1+2*xhalo) .and. (size(field_in,2) .eq. jeg-jsg+1+2*yhalo)) then
     field2(is:ie,js:je) = field_in(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1)
  else
     print*, 'on pe ', pe, 'domain: ', isc, iec, jsc, jec, isd, ied, jsd, jed, 'size of field: ', size(field_in,1), size(field_in,2)
     call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//':field is not on compute, data or global domain')
  endif
     
  call mpp_sync_self()
  
  if(any(pelist1 == pe)) then  ! send data to root pe
  
     im = ie-is+1; jm=je-js+1
     allocate(send_buffer(im*jm))
     
     ibounds(1) = is; ibounds(2) = ie; ibounds(3) = js; ibounds(4) = je
     l = 0
     do i = is,ie
     do j = js,je
        l = l+1
        send_buffer(l) = field2(i,j)
     enddo
     enddo
!  send the check bounds and data to the root pe
     ! Force use of "scalar", integer pointer mpp interface
     call mpp_send(ibounds(1), plen=4, to_pe=pelist2(1), tag=COMM_TAG_1)
     call mpp_send(send_buffer(1),plen=im*jm, to_pe=pelist2(1), tag=COMM_TAG_2)
     deallocate(send_buffer)
     
   else if(pelist2(1) == pe) then        ! receive data and compare
     do p = pelist1(1), pelist1(size(pelist1(:)))
     ! Force use of "scalar", integer pointer mpp interface
        call mpp_recv(ibounds(1), glen=4,from_pe=p, tag=COMM_TAG_1)
        is = ibounds(1); ie = ibounds(2); js=ibounds(3); je=ibounds(4)
        im = ie-is+1; jm=je-js+1
        if(allocated(field1)) deallocate(field1)
        if(allocated(send_buffer)) deallocate(send_buffer)
        allocate(field1(is:ie,js:je),send_buffer(im*jm))
     ! Force use of "scalar", integer pointer mpp interface
        call mpp_recv(send_buffer(1),glen=im*jm,from_pe=p, tag=COMM_TAG_2)
        l = 0
        
!  compare here, the comparison criteria can be changed according to need
        do i = is,ie
        do j = js,je
           l = l+1
           field1(i,j) = send_buffer(l)
           if(field1(i,j) .ne. field2(i,j)) then
         !   write to standard output
             print*,trim(mesg)//": ", i, j, field1(i,j), field2(i,j), field1(i,j) - field2(i,j)
!             write(stdout(),'(a,2i,2f)') trim(mesg), i, j, pass_field(i,j), field_check(i,j)
             check_success = .FALSE.
             if(error_exit) call mpp_error(FATAL,"mpp_check_field: can not reproduce at this point")
           endif
        enddo
        enddo
      enddo
        
      if(check_success) then
         print*, trim(mesg)//": ", 'comparison between 1 pe and ', npes-1, ' pes is ok'
      endif
  ! release memery
      deallocate(field1, send_buffer)
    endif
      
    deallocate(field2)
    
    call mpp_sync()
    
  end subroutine check_field_2d_type1
  
!####################################################################

  subroutine check_field_2d_type2(field_in, pelist1, pelist2, domain, mesg,force_abort)
!  This routine is used to check field between running on m pe (root pe) and
!  n pe. This routine can not check halo.
 
  real, dimension(:,:),  intent(in) :: field_in
  type(domain2d_private),intent(in) :: domain
  integer, dimension(:), intent(in) :: pelist1
  integer, dimension(:), intent(in) :: pelist2
  character(len=*),      intent(in) :: mesg
  logical, intent(in), optional     :: force_abort   ! when, call mpp_error if any difference
                                                    ! found. default value is false.
! some local variables                              
  logical                :: check_success, error_exit
  real, dimension(:,:), allocatable :: field1, field2
  integer :: i, j, pe, npes, isd,ied,jsd,jed, is, ie, js, je
  type(domain2d_private) :: domain1, domain2
  
    check_success = .TRUE.
    error_exit    = .FALSE.
    if(present(force_abort)) error_exit = force_abort
    pe = mpp_pe()
    npes = mpp_npes()
    call mpp_sync_self()
    if(any(pelist1 == pe)) domain1 = domain
    if(any(pelist2 == pe)) domain2 = domain
    
!  Comparison is made on pelist2.
    if(any(pelist2 == pe)) then
       call mpp_get_data_domain(domain2, isd, ied, jsd, jed)
       call mpp_get_compute_domain(domain2, is, ie, js, je)
       allocate(field1(isd:ied, jsd:jed),field2(isd:ied, jsd:jed))
      if((size(field_in,1) .ne. ied-isd+1) .or. (size(field_in,2) .ne. jed-jsd+1)) &
         call mpp_error(FATAL,'mpp_check_field: input field is not on the data domain')
      field2(isd:ied, jsd:jed) = field_in(:,:)
    endif
      
!  broadcast domain
    call broadcast_domain_1(domain1)
    call broadcast_domain_1(domain2)
    
    call redistribute_2d(domain1,field_in,domain2,field1)
    
    if(any(pelist2 == pe)) then
        do i =is,ie
        do j =js,je
          if(field1(i,j) .ne. field2(i,j)) then
             print*, trim(mesg)//": ", i, j, field1(i,j), field2(i,j), field1(i,j) - field2(i,j)
!             write(stdout(),'(a,2i,2f)') trim(mesg), i, j, field_check(i,j), field_out(i,j)
             check_success = .FALSE.
             if(error_exit) call mpp_error(FATAL,"mpp_check_field: can not reproduce at this point")
          endif
        enddo
        enddo
        if(check_success) &
             print*, trim(mesg)//": ", 'comparison between ', size(pelist1(:)), ' pes and ', &
                   size(pelist2(:)), ' pe on', pe, ' pes is ok'
    endif          
                   
    if(any(pelist2 == pe))    deallocate(field1, field2)
    
    call mpp_sync()
    call deallocate_domain2D_local(domain1)
    call deallocate_domain2D_local(domain2)    

    return
    
  end subroutine check_field_2d_type2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_BROADCAST_DOMAIN                                           !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_domain_1( domain )
!broadcast domain (useful only outside the context of its own pelist)
      type(domain2D), intent(inout) :: domain

      if(domain%index == DEFAULT_DOMAIN2D_INDEX ) then
         cur_domain_index = cur_domain_index+1
         domain%index = cur_domain_index
      endif
      call check_domain2d_index(domain, "mpp_broadcast_domain_1")
      call broadcast_domain_1(domainList(domain%index))

    end subroutine mpp_broadcast_domain_1

    subroutine broadcast_domain_1( domain )
!broadcast domain (useful only outside the context of its own pelist)
      type(domain2D_private), intent(inout) :: domain
      integer, allocatable :: pes(:)
      logical :: native         !true if I'm on the pelist of this domain
      integer :: listsize, listpos
      integer :: n
      integer, dimension(11) :: msg, info         !pe and compute domain of each item in list
      integer                :: errunit      

      errunit = stderr()
      if( .NOT.module_is_initialized ) &
                 call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN_1: You must first call mpp_domains_init.' ) 
      
!get the current pelist
      allocate( pes(0:mpp_npes()-1) )
      call mpp_get_current_pelist(pes)
      
!am I part of this domain?
      native = ASSOCIATED(domain%list)
      
!set local list size
      if( native )then
          listsize = size(domain%list(:))
      else
          listsize = 0
      end if
      call mpp_max(listsize)
      
      if( .NOT.native )then
!initialize domain%list and set null values in message
          allocate( domain%list(0:listsize-1) )
          domain%pe = NULL_PE
          domain%pos = -1
          allocate(domain%x(1), domain%y(1), domain%tile_id(1))
          do n = 0, listsize-1
             allocate(domain%list(n)%x(1), domain%list(n)%y(1), domain%list(n)%tile_id(1) )
          end do
          domain%x%compute%begin =  1
          domain%x%compute%end   = -1
          domain%y%compute%begin =  1
          domain%y%compute%end   = -1
          domain%x%global %begin = -1
          domain%x%global %end   = -1
          domain%y%global %begin = -1
          domain%y%global %end   = -1
          domain%tile_id         = -1
          domain%whalo           = -1
          domain%ehalo           = -1
          domain%shalo           = -1
          domain%nhalo           = -1
          domain%symmetry        = .false.
      end if
!initialize values in info
      info(1) = domain%pe
      call mpp_get_compute_domain( domain, info(2), info(3), info(4), info(5) )
      info(6) = domain%tile_id(1)
      info(7) = domain%whalo
      info(8) = domain%ehalo
      info(9) = domain%shalo
      info(10)= domain%nhalo
      if(domain%symmetry) then
         info(11) = 1
      else
         info(11) = 0
      endif
!broadcast your info across current pelist and unpack if needed
      listpos = 0
      do n = 0,mpp_npes()-1
         msg = info
         if( mpp_pe().EQ.pes(n) .AND. debug )write( errunit,* )'PE ', mpp_pe(), 'broadcasting msg ', msg
         call mpp_broadcast( msg, 11, pes(n) )
!no need to unpack message if native
!no need to unpack message from non-native PE
         if( .NOT.native .AND. msg(1).NE.NULL_PE )then
             domain%list(listpos)%pe = msg(1)
             domain%list(listpos)%x%compute%begin = msg(2)
             domain%list(listpos)%x%compute%end   = msg(3)
             domain%list(listpos)%y%compute%begin = msg(4)
             domain%list(listpos)%y%compute%end   = msg(5)
             domain%list(listpos)%tile_id(1)      = msg(6)
             if(domain%x(1)%global%begin < 0) then
                domain%x(1)%global%begin = msg(2)
                domain%x(1)%global%end   = msg(3)
                domain%y(1)%global%begin = msg(4)
                domain%y(1)%global%end   = msg(5)
                domain%whalo             = msg(7)
                domain%ehalo             = msg(8)
                domain%shalo             = msg(9)
                domain%nhalo             = msg(10)
                if(msg(11) == 1) then
                   domain%symmetry       = .true.
                else
                   domain%symmetry       = .false.
                endif
             else
                domain%x(1)%global%begin = min(domain%x(1)%global%begin, msg(2))
                domain%x(1)%global%end   = max(domain%x(1)%global%end,   msg(3))
                domain%y(1)%global%begin = min(domain%y(1)%global%begin, msg(4))
                domain%y(1)%global%end   = max(domain%y(1)%global%end,   msg(5)) 
             endif
             listpos = listpos + 1
             if( debug )write( errunit,* )'PE ', mpp_pe(), 'received domain from PE ', msg(1), 'is,ie,js,je=', msg(2:5)
         end if
      end do 

    end subroutine broadcast_domain_1


!##############################################################################
    subroutine mpp_broadcast_domain_2( domain_in, domain_out )
!broadcast domain (useful only outside the context of its own pelist)
      type(domain2D), intent(in)    :: domain_in
      type(domain2D), intent(inout) :: domain_out
  
      call check_domain2d_index(domain_in, "mpp_broadcast_domain_2")
      if(domain_in%index == DEFAULT_DOMAIN2D_INDEX)  &
            call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN_2: domain_in%index == DEFAULT_DOMAIN2D_INDEX')
      if(domain_out%index .NE. DEFAULT_DOMAIN2D_INDEX) &
          call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN_2: domain_out%index .NE. DEFAULT_DOMAIN2D_INDEX')
      cur_domain_index = cur_domain_index+1
      domain_out%index = cur_domain_index
      call check_domain2d_index(domain_out, "mpp_broadcast_domain_2")      

      call broadcast_domain_2(domainList(domain_in%index), domainList(domain_out%index))

    end subroutine mpp_broadcast_domain_2

    subroutine broadcast_domain_2( domain_in, domain_out )
!broadcast domain (useful only outside the context of its own pelist)
      type(domain2D_private), intent(in)    :: domain_in
      type(domain2D_private), intent(inout) :: domain_out
      integer, allocatable :: pes(:)
      integer :: listpos
      integer :: n
      integer, dimension(12) :: msg, info         !pe and compute domain of each item in list
      integer                :: errunit, npes_in, npes_out, pstart, pend

      errunit = stderr()
      if( .NOT.module_is_initialized ) &
                 call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN_2: You must first call mpp_domains_init.' )

!get the current pelist
      allocate( pes(0:mpp_npes()-1) )
      call mpp_get_current_pelist(pes)

! domain_in must be initialized
      if( .not. ASSOCIATED(domain_in%list) ) then
         call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN_2: domain_in is not initialized')
      endif
      if( ASSOCIATED(domain_out%list) ) then
         call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN_2: domain_out is already initialized')
      endif

      npes_in = size(domain_in%list(:))
      if( npes_in == mpp_npes() ) then
         call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN_2: size(domain_in%list(:)) == mpp_npes()')
      endif
      npes_out = mpp_npes() - npes_in

!initialize domain_out%list and set null values in message
      allocate( domain_out%list(0:npes_out-1) )
      domain_out%pe = NULL_PE
      domain_out%pos = -1
      allocate(domain_out%x(1), domain_out%y(1), domain_out%tile_id(1))
      do n = 0, npes_out-1
         allocate(domain_out%list(n)%x(1), domain_out%list(n)%y(1), domain_out%list(n)%tile_id(1) )
      end do
      domain_out%x%compute%begin =  1
      domain_out%x%compute%end   = -1
      domain_out%y%compute%begin =  1
      domain_out%y%compute%end   = -1
      domain_out%x%global %begin = -1
      domain_out%x%global %end   = -1
      domain_out%y%global %begin = -1
      domain_out%y%global %end   = -1
      domain_out%tile_id         = -1
      domain_out%whalo           = -1
      domain_out%ehalo           = -1
      domain_out%shalo           = -1
      domain_out%nhalo           = -1
      domain_out%symmetry        = .false.
!initialize values in info
      info(1) = domain_in%pe
      call mpp_get_compute_domain( domain_in, info(2), info(3), info(4), info(5) )
      info(6) = domain_in%tile_id(1)
      info(7) = domain_in%whalo
      info(8) = domain_in%ehalo
      info(9) = domain_in%shalo
      info(10)= domain_in%nhalo
      if(domain_in%symmetry) then
         info(11) = 1
      else
         info(11) = 0
      endif
      info(12) = domain_in%ntiles

!broadcast your info across current pelist and unpack if needed
      if( domain_in%list(0)%pe == mpp_root_pe() ) then
          pstart = npes_in
          pend   = mpp_npes()-1
      else
          pstart = 0
          pend   = npes_out-1
      endif
      do n = 0,mpp_npes()-1
         msg = info
         if( mpp_pe().EQ.pes(n) .AND. debug )write( errunit,* )'PE ', mpp_pe(), 'broadcasting msg ', msg
         call mpp_broadcast( msg, 12, pes(n) )
         !--- pack if from other domain
         if( n .GE. pstart .AND. n .LE. pend )then
             listpos = n - pstart
             domain_out%list(listpos)%pe = msg(1)
             domain_out%list(listpos)%x%compute%begin = msg(2)
             domain_out%list(listpos)%x%compute%end   = msg(3)
             domain_out%list(listpos)%y%compute%begin = msg(4)
             domain_out%list(listpos)%y%compute%end   = msg(5)
             domain_out%list(listpos)%tile_id(1)      = msg(6)
             if(domain_out%x(1)%global%begin < 0) then
                domain_out%x(1)%global%begin = msg(2)
                domain_out%x(1)%global%end   = msg(3)
                domain_out%y(1)%global%begin = msg(4)
                domain_out%y(1)%global%end   = msg(5)
                domain_out%whalo             = msg(7)
                domain_out%ehalo             = msg(8)
                domain_out%shalo             = msg(9)
                domain_out%nhalo             = msg(10)
                if(msg(11) == 1) then
                   domain_out%symmetry       = .true.
                else
                   domain_out%symmetry       = .false.
                endif
                domain_out%ntiles            = msg(12)
             else
                domain_out%x(1)%global%begin = min(domain_out%x(1)%global%begin, msg(2))
                domain_out%x(1)%global%end   = max(domain_out%x(1)%global%end,   msg(3))
                domain_out%y(1)%global%begin = min(domain_out%y(1)%global%begin, msg(4))
                domain_out%y(1)%global%end   = max(domain_out%y(1)%global%end,   msg(5))
             endif
             if( debug )write( errunit,* )'PE ', mpp_pe(), 'received domain from PE ', msg(1), 'is,ie,js,je=', msg(2:5)
         end if
      end do

    end subroutine broadcast_domain_2

