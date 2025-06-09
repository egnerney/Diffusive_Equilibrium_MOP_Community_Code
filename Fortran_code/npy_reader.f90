!=====================================================================
!  npy_reader.f90  –  loader for *.npy input files
!=====================================================================
module npy_reader
   use, intrinsic :: iso_fortran_env, only : dp => real64, int8,int16,int32
   implicit none
contains
!---------------------------------------------------------------------
subroutine read_npy_1d(vec,fname)
   character(*),          intent(in)  :: fname
   real(dp), allocatable, intent(out) :: vec(:)

   integer :: fd, ios, hlen
   character(len=6) :: magic
   integer(int8)  :: ver(2)
   integer(int16) :: h16
   character(len=:), allocatable :: header
   integer, allocatable :: shape(:)
   logical :: f_order

   open(newunit=fd,file=fname,access='stream',status='old',action='read',iostat=ios)
   if (ios /= 0) then
      write(*,*) 'ERROR: cannot open ', trim(fname)
      stop 1
   end if
   read(fd) magic; read(fd) ver; read(fd) h16
   if (magic /= char(147)//'NUMPY') then
      write(*,*) 'ERROR: not a NumPy file ', trim(fname)
      stop 1
   end if
   hlen = int(h16,int32); if (hlen < 0) hlen = hlen + 65536

   allocate(character(len=hlen)::header); read(fd) header
   call parse_hdr(header, f_order, shape)

   select case (size(shape))
   case (1)
      allocate(vec(shape(1)))
   case (2)
      if (shape(1) == 1) then
         allocate(vec(shape(2)))
      else if (shape(2) == 1) then
         allocate(vec(shape(1)))
      else
         write(*,*) 'ERROR: expected 1-D in ', trim(fname)
         stop 1
      end if
   case default
      write(*,*) 'ERROR: expected 1-D in ', trim(fname)
      stop 1
   end select

   read(fd) vec
   close(fd)
end subroutine read_npy_1d
!---------------------------------------------------------------------
subroutine read_npy_2d(mat,fname)
   character(*),          intent(in)  :: fname
   real(dp), allocatable, intent(out) :: mat(:,:)

   integer :: fd, ios, hlen, nx, ny
   character(len=6) :: magic
   integer(int8)  :: ver(2)
   integer(int16) :: h16
   character(len=:), allocatable :: header
   integer, allocatable :: shape(:)
   logical :: f_order
   real(dp), allocatable :: tmp(:)

   open(newunit=fd,file=fname,access='stream',status='old',action='read',iostat=ios)
   if (ios /= 0) then
      write(*,*) 'ERROR: cannot open ', trim(fname)
      stop 1
   end if
   read(fd) magic; read(fd) ver; read(fd) h16
   if (magic /= char(147)//'NUMPY') then
      write(*,*) 'ERROR: not a NumPy file ', trim(fname)
      stop 1
   end if
   hlen = int(h16,int32); if (hlen < 0) hlen = hlen + 65536

   allocate(character(len=hlen)::header); read(fd) header
   call parse_hdr(header, f_order, shape)
   if (size(shape) /= 2) then
      write(*,*) 'ERROR: expected 2-D in ', trim(fname)
      stop 1
   end if

   nx = shape(1); ny = shape(2); allocate(mat(nx,ny))
   if (f_order) then
      read(fd) mat
   else
      allocate(tmp(nx*ny)); read(fd) tmp
      mat = transpose(reshape(tmp,[ny,nx])); deallocate(tmp)
   end if
   close(fd)
end subroutine read_npy_2d
!---------------------------------------------------------------------
subroutine parse_hdr(txt,f_order,shape)
   character(*), intent(in)  :: txt
   logical,      intent(out) :: f_order
   integer, allocatable, intent(out) :: shape(:)

   integer :: p1,p2,ios,n1,n2
   character(len=:),allocatable :: shp

   p1 = index(txt,'fortran_order') + 15
   f_order = (txt(p1:p1) == 'T')

   p1 = index(txt,'(') + 1
   p2 = index(txt,')') - 1
   shp = adjustl(txt(p1:p2))
   if (shp(len_trim(shp):len_trim(shp)) == ',') shp(len_trim(shp):) = ' '

   read(shp,*,iostat=ios) n1, n2
   if (ios == 0 .and. n2 /= 0) then
      allocate(shape(2)); shape = [n1, n2]
   else
      allocate(shape(1)); shape(1) = n1
   end if
end subroutine parse_hdr
!=====================================================================
end module npy_reader
