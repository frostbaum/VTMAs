program calc_vtma
  use v3d_func_rep
  implicit none
  integer :: natm, nl, ierr
  integer :: i, j, k
  double precision, dimension(:,:,:), allocatable :: dipgrad, nmodes
  double precision, dimension(:,:), allocatable :: vibspec
  double precision, dimension(3) :: vtm, dipm
  character(len=99) :: line, dir
  character(len=1) :: dummy
  logical :: numread = .false., dipread = .false., rerr = .false.
  
  dir = ''
  call getarg(1,dir)
  
  !~ Reading file "control"
  !~ Getting: natm, dipm
  open(101,file=trim(join_fpath(trim(dir),'control')),status='old',iostat=ierr)
  
  if (ierr .ne. 0) then
    write(6,'(A)') 'Cannot open '//trim(join_fpath(trim(dir),'control'))
    stop
  end if
  
  do
    if (numread .and. dipread) exit
    read(101,'(A)',iostat=ierr) line
    
    if (ierr .ne. 0) then
      write(6,'(A)') trim(join_fpath(trim(dir),'control'))//' seems to be broken'
      rerr = .true.
      exit
    else
      
      if (line(4:9) .eq. 'natoms') then
        line(1:10) = ''
        read(line,*,iostat=ierr) natm
        
        if (ierr .ne. 0) then
          write(6,'(A)') 'Error reading "natoms"'
          rerr = .true.
          exit
        end if
        
        numread = .true.
        
      else if (line(1:7) .eq. '$dipole') then
        read(101,*,iostat=ierr) dummy, dipm(1), dummy, dipm(2), dummy, dipm(3)
        
        if (ierr .ne. 0) then
          write(6,'(A)') 'Error reading "$dipole"'
          rerr = .true.
          exit
        end if
        
        dipread = .true.
        
      end if
      
    end if
  end do
  
  close(101)
  
  if (rerr) stop
  
  allocate(nmodes(3,natm,3*natm),dipgrad(3,3,natm),vibspec(3,3*natm))
  
  nl = floor(3*natm/5.d0)
  
  !~ Reading file "vib_normal_modes"
  !~ Getting: nmodes
  open(102,file=trim(join_fpath(trim(dir),'vib_normal_modes')),status='old',iostat=ierr)
  read(102,*,iostat=ierr)
  
  if (ierr .ne. 0) then
    write(6,'(A)') 'Cannot open '//trim(join_fpath(trim(dir),'vib_normal_modes'))
    stop
  end if
  
  do i = 1, natm
    do j = 1, 3
      do k = 1, nl
        read(102,'(A5,A)',iostat=ierr) dummy, line
        
        if (ierr .ne. 0) then
          write(6,'(A)') trim(join_fpath(trim(dir),'vib_normal_modes'))//' seems to be broken'
          close(102)
          stop
        else
          read(line,*,iostat=ierr) nmodes(j,i,5*(k-1)+1:5*k)
          
          if (ierr .ne. 0) then
            write(6,'(A)') 'Error reading vib_normal_modes'
            close(102)
            stop
          end if
          
        end if
      end do
      
      if (modulo(natm*3,5) .ne. 0) then
        read(102,'(A5,A)',iostat=ierr) dummy, line
        
        if (ierr .ne. 0) then
          write(6,'(A)') trim(join_fpath(trim(dir),'vib_normal_modes'))//' seems to be broken'
          close(102)
          stop
        else
          read(line,*,iostat=ierr) nmodes(j,i,5*nl+1:natm*3)
          
          if (ierr .ne. 0) then
            write(6,'(A)') 'Error reading vib_normal_modes'
            close(102)
            stop
          end if
          
        end if
      end if
    end do
  end do
  
  close(102)
  
  !~ Reading file "dipgrad"
  !~ Getting: dipgrad
  open(103,file=trim(join_fpath(trim(dir),'dipgrad')),status='old',iostat=ierr)
  read(103,*,iostat=ierr)
  
  if (ierr .ne. 0) then
    write(6,'(A)') 'Cannot open '//trim(join_fpath(trim(dir),'dipgrad'))
    stop
  end if
  
  do i = 1, natm
    do j = 1, 3
      
      read(103,*,iostat=ierr) dipgrad(:,j,i)
      
      if (ierr .ne. 0) then
        write(6,'(A)') 'Error reading dipgrad'
        close(103)
        stop
      end if
      
    end do
  end do
  
  close(103)
  
  !~ Reading file "vibspectrum"
  !~ Getting: vibspec
  open(104,file=trim(join_fpath(trim(dir),'vibspectrum')),status='old',iostat=ierr)
  read(104,*,iostat=ierr)
  read(104,*,iostat=ierr)
  read(104,*,iostat=ierr)
  
  if (ierr .ne. 0) then
    write(6,'(A)') 'Cannot open '//trim(join_fpath(trim(dir),'vibspectrum'))
    stop
  end if
  
  do i = 1, 3*natm
    read(104,'(A20,A)',iostat=ierr) dummy,line
    
    if (ierr .ne. 0) then
      write(6,'(A)') 'Error reading vibspec'
      close(104)
      stop
    end if
    
    read(line,*,iostat=ierr) vibspec(1:2,i)
    
    if (ierr .ne. 0) then
      write(6,'(A)') 'Error reading vibspec'
      close(104)
      stop
    end if
  end do
  
  close(104)
  
  !~ Calculating VTMA
  do k = 1, 3*natm
    
    vtm = 0.d0
    
    do i=1,natm
      do j=1,3
        vtm(:) = vtm(:) + nmodes(j,i,k)*dipgrad(:,j,i)
      end do
    end do
    
    vibspec(3,k) = min(get_angle(dipm,-vtm),get_angle(dipm,vtm))*180/3.141592653589793d0
  end do
  
  !~ Printing frequencies, intensities and VTMAs
  write(6,'(A)') '# mode  wave number  IR intensity    VTMA'
  write(6,'(A)') '#         [cm^-1]      [km/mol]    [degree]'
  do i = 1, 3*natm
    write(6,'(I5,F12.2,F14.5,F10.1)') i, vibspec(:,i)
  end do
  
  contains
  
  function join_fpath(dir,fname) result(string)
    character(*) :: dir, fname
    character(len=99) :: string
    integer :: lpos
    
    if (dir .ne. '') then
      lpos = len_trim(dir)
      
      if (dir(lpos:lpos) .ne. '/') then
        string = trim(adjustl(dir))//'/'//trim(adjustl(fname))
      else
        string = trim(adjustl(dir))//trim(adjustl(fname))
      end if
    else
      string = trim(adjustl(fname))
    end if
  end function
end program
