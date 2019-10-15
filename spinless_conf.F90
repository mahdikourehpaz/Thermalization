subroutine spinless_fermion_plus_special_particle

  implicit none

  integer :: i,j,N,p,kk
  integer*8 , allocatable , dimension(:) :: A
  integer*8  :: m , m_max , m_min 

  open(unit=1,file="initial.in",action="read")
  read(1,*)
  read(1,*) N , p
  close(1)
  print*, " Number of Sites         " , N 
  print*, " Number of Bath Particles" , p
  print*, " Test particle included"
  open(unit=1,file="configuration.dat")
  

  !All possible configuration of having "p" polarized fermion and the test particle in "N" sites. First p numbers are positions of
  !fermions and p+1 number corresponds to test particle position
  !Example: N=7 & p=3
  !        ( 1     2     4     6)           --->           -|-    -|-   --   -|-  --   -o-   -- 



  if(p .eq. 1) then
    do i=1,N
      do j=1,N
        write(1,333) i , j
      end do
    end do
  else
    m_min = (N+1)**(p-1)
    m_max = (N-p+2)*((N+1)**(p-1))

    allocate(A(1:p)) ; A = 0

    do m=m_min,m_max

      do j=1,p
        A(j) = mod(m/(N+1)**(p-j),N)
        if(A(j) .eq. 0) then
          A(j) = N
        end if
      end do


      do i=p,2,-1
        if(A(i) .le. A(i-1)) then
          goto 10
        else if(i .eq. 2) then

          do kk=1,N
            write(1,333) A(:) , kk
          end do
        else
          goto 20
        end if
        20 continue
      end do
      10 continue

    end do
  end if
  deallocate(A)
  333 format(100I4)
end subroutine spinless_fermion_plus_special_particle
