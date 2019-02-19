include(`macros.m4')dnl
!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
`subroutine domain_read_'NUM_DIMS`d(fileobj, &'
                                    variable_name, &
                                    vdata, &
                                    unlim_dim_level)

    !Inputs/outputs
    type(FmsNetcdfDomainFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dim_declare(NUM_DIMS) intent(inout) :: vdata !< Data that will
                                                          !! be written out
                                                          !! to the netcdf file.
    integer,intent(in),optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.

ifelse(NUM_DIMS,0,`define(`ATLEAST2D',0)',
       NUM_DIMS,1,`define(`ATLEAST2D',0)',
       `define(`ATLEAST2D',1)')dnl

    !Local variables.
    integer :: xdim_index
    integer :: ydim_index
ifelse(ATLEAST2D,1,
   `type(domain2d),pointer :: io_domain
    integer,dimension(:),allocatable :: dimension_sizes
    integer :: dpos
    integer :: i
    integer :: ds
    integer,dimension(:),allocatable :: isc
    integer,dimension(:),allocatable :: icsize
    integer,dimension(:),allocatable :: jsc
    integer,dimension(:),allocatable :: jcsize
    integer,dimension(NUM_DIMS) :: c
    integer,dimension(NUM_DIMS) :: e
    integer(kind=int32),dim_declare(NUM_DIMS) allocatable :: buf_int32
    integer(kind=int64),dim_declare(NUM_DIMS) allocatable :: buf_int64
    real(kind=real32),dim_declare(NUM_DIMS) allocatable :: buf_real32
    real(kind=real64),dim_declare(NUM_DIMS) allocatable :: buf_real64
    integer(kind=int32),dim_declare(NUM_DIMS) allocatable :: lbuf_int32
    integer(kind=int64),dim_declare(NUM_DIMS) allocatable :: lbuf_int64
    real(kind=real32),dim_declare(NUM_DIMS) allocatable :: lbuf_real32
    real(kind=real64),dim_declare(NUM_DIMS) allocatable :: lbuf_real64')

    if (.not. is_variable_domain_decomposed(fileobj,variable_name,.true., &
                                            xdim_index,ydim_index)) then
        call netcdf_read_data(fileobj, &
                              variable_name, &
                              vdata, &
                              unlim_dim_level=unlim_dim_level, &
                              broadcast=.true.)
        return
ifelse(ATLEAST2D,0,
   `else
        call error("this branch should never be reached.")')
    endif
ifelse(ATLEAST2D,1,
   `io_domain => mpp_get_io_domain(fileobj%domain)
    call get_variable_size(fileobj, &
                           variable_name, &
                           dimension_sizes, &
                           broadcast=.true.)
    dpos = get_domain_decomposed_variable_position(fileobj, &
                                                   variable_name)
    do i = 1,size(dimension_sizes)
        if (i .eq. xdim_index) then
            call mpp_get_compute_domain(io_domain, &
                                        xsize=ds, &
                                        position=dpos)
            if (size(vdata,i) .ne. ds) then
                call error("domain decomposed x dimension size mismatch.")
            endif
        elseif (i .eq. ydim_index) then
            call mpp_get_compute_domain(io_domain, &
                                        ysize=ds, &
                                        position=dpos)
            if (size(vdata,i) .ne. ds) then
                call error("domain decomposed y dimension size mismatch.")
            endif
        endif
    enddo

    !I/O root reads in the data and scatters it.
    if (fileobj%is_root) then
        allocate(isc(size(fileobj%pelist)))
        allocate(icsize(size(fileobj%pelist)))
        allocate(jsc(size(fileobj%pelist)))
        allocate(jcsize(size(fileobj%pelist)))
        call mpp_get_compute_domains(io_domain, &
                                     xbegin=isc, &
                                     xsize=icsize, &
                                     ybegin=jsc, &
                                     ysize=jcsize, &
                                     position=dpos)
        if (isc(1) .ne. 1) then
            isc(:) = isc(:) - isc(1) + 1
        endif
        if (jsc(1) .ne. 1) then
            jsc(:) = jsc(:) - jsc(1) + 1
        endif
        c(:) = 1
        c(xdim_index) = isc(1)
        c(ydim_index) = jsc(1)
        e(:) = shape(vdata)
        e(xdim_index) = icsize(1)
        e(ydim_index) = jcsize(1)
        select type(vdata)
            type is (integer(kind=int32))
                call allocate_array(buf_int32, &
                                    dimension_sizes)
                call netcdf_read_data(fileobj, &
                                      variable_name, &
                                      buf_int32, &
                                      unlim_dim_level=unlim_dim_level, &
                                      broadcast=.false.)
                call get_array_section(vdata, &
                                       buf_int32, &
                                       c, &
                                       e)
            type is (integer(kind=int64))
                call allocate_array(buf_int64, &
                                    dimension_sizes)
                call netcdf_read_data(fileobj, &
                                      variable_name, &
                                      buf_int64, &
                                      unlim_dim_level=unlim_dim_level, &
                                      broadcast=.false.)
                call get_array_section(vdata, &
                                       buf_int64, &
                                       c, &
                                       e)
            type is (real(kind=real32))
                call allocate_array(buf_real32, &
                                    dimension_sizes)
                call netcdf_read_data(fileobj, &
                                      variable_name, &
                                      buf_real32, &
                                      unlim_dim_level=unlim_dim_level, &
                                      broadcast=.false.)
                call get_array_section(vdata, &
                                       buf_real32, &
                                       c, &
                                       e)
            type is (real(kind=real64))
                call allocate_array(buf_real64, &
                                    dimension_sizes)
                call netcdf_read_data(fileobj, &
                                      variable_name, &
                                      buf_real64, &
                                      unlim_dim_level=unlim_dim_level, &
                                      broadcast=.false.)
                call get_array_section(vdata, &
                                       buf_real64, &
                                       c, &
                                       e)
            class default
                call error("unsupported type.")
        end select
        do i = 2,size(fileobj%pelist)
            c(xdim_index) = isc(i)
            c(ydim_index) = jsc(i)
            e(xdim_index) = icsize(i)
            e(ydim_index) = jcsize(i)
            select type(vdata)
                type is (integer(kind=int32))
                    call allocate_array(lbuf_int32, &
                                        e)
                    call get_array_section(lbuf_int32, &
                                           buf_int32, &
                                           c, &
                                           e)
                    call mpp_send(lbuf_int32, &
                                  size(lbuf_int32), &
                                  fileobj%pelist(i))
                    call mpp_sync_self(check=EVENT_SEND)
                    deallocate(lbuf_int32)
                type is (integer(kind=int64))
                    call allocate_array(lbuf_int64, &
                                        e)
                    call get_array_section(lbuf_int64, &
                                           buf_int64, &
                                           c, &
                                           e)
                    call mpp_send(lbuf_int64, &
                                  size(lbuf_int64), &
                                  fileobj%pelist(i))
                    call mpp_sync_self(check=EVENT_SEND)
                    deallocate(lbuf_int64)
                type is (real(kind=real32))
                    call allocate_array(lbuf_real32, &
                                        e)
                    call get_array_section(lbuf_real32, &
                                           buf_real32, &
                                           c, &
                                           e)
                    call mpp_send(lbuf_real32, &
                                  size(lbuf_real32), &
                                  fileobj%pelist(i))
                    call mpp_sync_self(check=EVENT_SEND)
                    deallocate(lbuf_real32)
                type is (real(kind=real64))
                    call allocate_array(lbuf_real64, &
                                        e)
                    call get_array_section(lbuf_real64, &
                                           buf_real64, &
                                           c, &
                                           e)
                    call mpp_send(lbuf_real64, &
                                  size(lbuf_real64), &
                                  fileobj%pelist(i))
                    call mpp_sync_self(check=EVENT_SEND)
                    deallocate(lbuf_real64)
                class default
                    call error("unsupported type.")
            end select
        enddo
        deallocate(isc)
        deallocate(icsize)
        deallocate(jsc)
        deallocate(jcsize)
        if (allocated(buf_int32)) then
            deallocate(buf_int32)
        endif
        if (allocated(buf_int64)) then
            deallocate(buf_int64)
        endif
        if (allocated(buf_real32)) then
            deallocate(buf_real32)
        endif
        if (allocated(buf_real64)) then
            deallocate(buf_real64)
        endif
    else
        select type(vdata)
            type is (integer(kind=int32))
                call mpp_recv(vdata, &
                              size(vdata), &
                              fileobj%io_root)
            type is (integer(kind=int64))
                call mpp_recv(vdata, &
                              size(vdata), &
                              fileobj%io_root)
            type is (real(kind=real32))
                call mpp_recv(vdata, &
                              size(vdata), &
                              fileobj%io_root)
            type is (real(kind=real64))
                call mpp_recv(vdata, &
                              size(vdata), &
                              fileobj%io_root)
            class default
                call error("unsupported type.")
        end select
        call mpp_sync_self(check=EVENT_RECV)
    endif
    deallocate(dimension_sizes)')
`end subroutine domain_read_'NUM_DIMS`d'
