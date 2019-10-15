subroutine H_forming_and_Diag

  implicit none
  character ( len = 2 ) :: IC_type , IC_Dist , Onsite_P_Type , time_type
  character(len=80):: filename
  integer , allocatable , dimension(:,:) :: cf 
  integer :: L , Nb , info , wave_f , Itteration , init_sp_position
  integer :: i , j , k , whole_HSS , l1 ,p , m , kk , i_t , timee , iii_it , ie0
  integer*8 :: bn
  real*8 , allocatable , dimension(:,:)  ::  H
  real*8 , allocatable , dimension(:) :: E , work , onsiteP , lmda , ic , s ,inner , MF
  real*8 :: t1 , t2 , v1 , v2 , vr , ts , vs , a , t00 , tff , t_scale
  real*8 :: tf , t0 , dt , t , temp_randn , var=0.d0 , E_in , E_final , tmax
  complex*16 , allocatable , dimension(:,:) :: rho 
  complex*16 , allocatable , dimension(:)   :: wfn , rhowork 
  complex*16 , parameter :: ui = cmplx(0.0,1.0)


  open(1,file="initial.in",action="read")
  read(1,*)
  read(1,*) L , p
  read(1,*)
  read(1,*) t1 , v1 , t2 , v2 , vr
  read(1,*)
  read(1,*) ts , vs
  read(1,*)
  read(1,*) IC_type  ,  init_sp_position
  read(1,*) 
  read(1,*) wave_f
  read(1,*)
  read(1,*) timee
  read(1,*)
  read(1,*) t00 , tff , dt , time_type
  read(1,*) 
  read(1,*) IC_Dist , var
  read(1,*)
  read(1,*) Itteration
  read(1,*)
  read(1,*) E_in , E_final
  read(1,*)
  read(1,*) Onsite_P_Type
  close(1)
  E_in = E_in / 100.d0
  E_final = E_final / 100.0
  write(*,*)
  write(*,353) "       t1       =   " , -1.d0*t1
  write(*,353) "       u1       =   " , v1
  write(*,353) "       t2       =   " , -1.d0*t2
  write(*,353) "       u2       =   " , v2
  write(*,353) "       vr       =   ", vr
  write(*,353) "       ts       =   ", -1.d0*ts
  write(*,353) "       ws       =   ", vs
  write(*,*)
  t1 = -1.d0 * t1
  t2 = -1.d0 * t2
  ts = -1.d0 * ts
  if(timee .eq. 0 ) write(*,*) "  No time propagation"
  Nb = int(bn(L,p),kind=4)
  whole_HSS = L*Nb                 !!! whole Hilbert Space Size
  write(*,*) " The Hilbert Space Size ", whole_HSS
  write(*,*)
  write(*,*)

  allocate(cf(whole_HSS,p+1)) ; cf = 0
  open(unit=2,file="configuration.dat",action="read")

  do i=1,whole_HSS
    read(2,*) cf(i,:)
  end do
  close(2)



  !Random on-site Potential
  allocate(onsiteP(L)) ; onsiteP = 0.d0

  open(unit=11,file="Onsite_Random_Potential.dat")
  if (Onsite_P_Type .eq. "R") then
    call random_seed()
    call random_number(OnsiteP)
    OnsiteP = vr*(2.d0*OnsiteP-1.d0)
  else if (Onsite_P_Type .eq. "L") then
    do i=1,L
      onsiteP(i) = -0.5 + i*1.0/(L+1.0)
    end do
    onsiteP(:) = 0.5*vr*onsiteP(:)/onsiteP(L)
  else if (Onsite_P_Type .eq. "Q") then
    do i=1,L
      onsiteP(i) = (-0.5 + i*1.0/(L+1.0))**2
    end do
    onsiteP(:) = vr*onsiteP(:)/onsiteP(L)
  end if 
  do i=1,L
    write(11,313) OnsiteP(i)
  end do
  close(11)

  allocate(H(whole_HSS,whole_HSS)) ; H = 0.d0
  do i=1,whole_HSS
    !diagonal elements
    H(i,i) = H(i,i) + onsiteP(cf(i,p+1))
    do j=1,p
      H(i,i) = H(i,i) + onsiteP(cf(i,j))
      if (abs(cf(i,j)-cf(i,p+1)) .eq. 0) H(i,i) = H(i,i) + Vs
    end do
    do j=1,p-1
      do k=j,p
        if (abs(cf(i,j)-cf(i,k)) .eq. 1) H(i,i) = H(i,i) + V1
      end do
      do k=j,p
        if (abs(cf(i,j)-cf(i,k)) .eq. 2) H(i,i) = H(i,i) + V2
      end do
    end do

    ! in order to have a non zero off-diagonal Hamiltonian entity H(i,j), test particle position displacement must be  either zero
    ! (bath can jump) or one (bath must be freeze)

    do j=i+1,whole_HSS
      ! Bath can jump
      if (abs(cf(i,p+1) - cf(j,p+1)) .eq. 0) then
        l1=0
        do k=1,p
          do kk=1,p
            if(cf(i,k).eq.cf(j,kk)) l1=l1+1
          end do
        end do
        if (l1 .eq. p-1) then
          k  = 0
          kk = 0
          do m=1,p
            k  =  k + cf(i,m)
            kk = kk + cf(j,m)
          end do
          if (abs(k-kk) .eq. 1) H(i,j) = t1
          if (abs(k-kk) .eq. 2) H(i,j) = t2
        end if

        ! Test particle can jump
      else if (abs(cf(i,p+1) - cf(j,p+1)) .eq.1 ) then
        l1=0
        do k=1,p
          if(cf(i,k).ne.cf(j,k)) l1=l1+1
        end do
        if (l1 .eq. 0) H(i,j) = ts
      end if

      H(j,i) = H(i,j)
    end do

  end do
  allocate(E(whole_HSS)) ; E = 0.d0
  allocate(work(3*whole_HSS-1)) ; work = 0.d0

  call dsyev('V','U',whole_HSS,H,whole_HSS,E,work,3*whole_HSS-1,info)

  open(unit=3,file="eigenenergies.dat")
  do i=1,whole_HSS
    write(3,313) E(i)
  end do
  close(3)

  if (wave_f .eq. 1) then
    open(unit=4,file="eigenvectors.dat")
    do i=1,whole_HSS
      write(4,313) H(i,:)
    end do
  end if
  close(4)
  deallocate(work)
  deallocate(onsiteP)


  write(*,*)
  write(*,*) "       HUUURRRRRAAAAAAAHHHHH      "
  write(*,*) " HAMILTONIAN HAS BEEN DIAGONALIZED"
  write(*,*)

  write(*,*) '#############################################'
  write(*,*) '#############################################'
  ! Now we have eigenenergies and energy eigenstates and can do the time propagation
  ! But first we have to do sth about initial condition

  if (IC_type .ne. 'WP') itteration = 1

  if (timee .eq. 1) then
    do iii_it=0,Itteration-1

      if (iii_it .eq. 0) then
        print*, L
        allocate(ic(whole_HSS)) ; ic = 0.d0              !!!! Initial Condition in Site Basis
        allocate(wfn(whole_HSS)) ; wfn = (0.d0,0.d0)     !!!! Time Dependent Many-Body Wave Function
        allocate(inner(whole_HSS)) ; inner = 0.d0        !!!! inner(i)=<E(i)|wfn(0)>  initial state in energy eigen basis
        allocate(MF(L)) ; MF = 0.d0                      !!!! Averaged Density of Bath Particles in each Site
        allocate(rho(L,L)) ; rho = (0.d0,0.d0)           !!!! 1RDM 
        allocate(lmda(L)) ; lmda = 0.d0                  !!!! Occupation of Natural Orbitals
        !!!! LAPACK
        allocate(rhowork(2*L-1)) ; rhowork = 0.d0
        allocate(work(3*L-2)) ; work = 0.d0
      end if
      wfn(:) = (0.d0,0.d0)
      MF(:) = 0.d0
      rho(:,:) = (0.d0,0.d0)
      lmda(:) = 0.d0
      rhowork(:) = 0.d0
      work(:) = 0.d0
      inner(:) = 0.d0
      ic(:) = 0.d0


      if (IC_type .eq. 'WP') then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!! Wave Packet Initial Condition!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ie0 = int(e_in + (e_final-e_in)*iii_it/itteration)*whole_HSS  
        if (ie0 .eq. 0) then
          ie0 = 1
        end if
        temp_randn = 0.d0

        if (IC_Dist .eq. 'GS') then
          !Gaussian WP around ie0
          inner = 0.d0
          do i=1,whole_HSS
            temp_randn = 0.d0
            call random_number(temp_randn)
            temp_randn = 2.d0*temp_randn - 1.d0
            inner(i) = temp_randn * exp(-((E(i)-E(ie0))**2)*1.d0/(2.d0*var))
          end do


        else if (IC_Dist .eq. 'UN') then
          !Uniformly WP around ie0
          inner = 0.d0
          do i=1,whole_HSS
            if(E(i).lt.(E(ie0)+var*0.d5) .and. E(i).gt.(E(ie0)-var*0.d5)) then
              call random_number(temp_randn)
              inner(i) = 2.d0* temp_randn- 1.d0
            end if
          end do


        else if (IC_Dist .eq. 'CN') then
          !Constantly WP around ie0
          inner = 0.d0
          do i=1,whole_HSS
            if(E(i).lt.(E(ie0)+var*0.5) .and. E(i).gt.(E(ie0)-var*0.d5)) then
              inner(i) = (-1.d0)**i
            end if
          end do
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else if (IC_type .eq. 'AR') then
        !!!!!!!!A random initial condition in whole space
        call RANDOM_NUMBER(ic)
      else if (IC_type .eq. 'SR') then
        !!!!!!!!Non-Localized Test-Particle and Localized Bath Left-Aligned
        do j=1,L
          call random_number(temp_randn)
          ic(j) = 2.d0*temp_randn - 1.d0
        end do
      else if (IC_type .eq. 'BR') then
        !!!!!!!!Localized Test-Particle in site 1 and Non-Localized Bath
        do j=0,Nb-1
          call random_number(temp_randn)
          ic(j*L+1) = 2.d0*temp_randn - 1.d0
        end do

      else if (IC_type .eq. 'AL') then
        !!!!!!!!Localized Test-Particle in site 'IC_type' where Bath is Left-Aligned
        ic(init_sp_position) = 1.d0

      end if


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!IC Normalization and printing ir!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (IC_type .ne. 'WP') then
        do i=1,whole_HSS
          inner(i) = dot_product(H(:,i),ic(:))
        end do
      end if

      a = 0.d0
      do i=1,whole_HSS
        a = a + abs(inner(i))**2
      end do
      inner = inner / sqrt(a)

      write(filename,'(a,i0,a)')'Wave_Function_Initial_',iii_it,'.dat'
      open(unit=8,file=filename)
      do i=1,whole_HSS
        if (i .eq. 1) then
          write(8,*) "#Energy_Basis"
        end if
        write(8,*) inner(i)
      end do
      close(8)

      a = sum(abs(inner)**2)
      if (abs(a - 1.000000000) .gt. 1e-8) then
        write(*,*) 'Initial Condition is not Normalized'
        write(*,*) a
        goto 1010
      end if



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!Finding Time steps, Initial!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!and Final Time Using Smallest!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!Energy Difference!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (IC_type .eq. 'WP') then
        i = 1
        k = whole_HSS
        do
          if(E(i) .ge. E(ie0) - 0.5*var) exit
          i = i+1
        end do
        do
          if(E(k) .le. E(ie0) + 0.5*var) exit
          k = k-1
        end do

        if (k-i .ne. 0) then
          allocate(s(k-i)) ; s = 0.d0
          do j = 0,k-i-1
            s(j+1) = E(i+j+1) - E(i+j)
          end do
          t_scale = 1.0 / minval(s)
          deallocate(s)
        else if (k-i .eq. 0) then
          t_scale = 1.0 / abs(E(k))
        end if

      else
        a = 0.d0
        do j=1,whole_HSS-1
          a = a + (E(j+1) - E(j))
        end do
        a = a / float(whole_HSS-1)          
        t_scale = 500.d0 / a
      end if

      if (time_type .eq. "L") then
        t0 = t00*t_scale
        tf = tff*t_scale
        tmax = tf-t0
      else
        t0 = t00
        tf = tff
        tmax = tf-t0
      end if

      if (iii_it .eq. 0) then
        write(*,*) " Now Time Propagation"
        write(*,*) " Time Steps           " , (int(dt)) 
        write(*,*)
        write(*,*)
      end if
      write(*,*) '#############################################'
      write(*,*)   "    " , iii_it
      write(*,343) "  From      t   =   " , (t0)
      write(*,343) "  Until     t   =   " , (tf)
      if (IC_type .eq. 'WP') write(*,343) "   <H>          =   " , (E(ie0)) 
      write(*,*)
      write(*,*)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!! Time Propagation and !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Tracing!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      write(filename,'(a,i0,a)')'Occ_Nat_Orb_',iii_it,'.dat'
      open(unit=7,file=filename)
      write(filename,'(a,i0,a)')'Re_Natural_Orb_',iii_it,'.dat'
      open(unit=9,file=filename)
      write(filename,'(a,i0,a)')'Im_Natural_Orb_',iii_it,'.dat'
      open(unit=10,file=filename)

      do i_t=0,int(dt)
        t = t0*1.d0 + i_t*tmax/dt
        wfn = cmplx(0.0,0.0)
        rho = cmplx(0.0,0.0)
        !$omp parallel do private(i,j)
        do i=1,whole_HSS
          do j=1,whole_HSS
            wfn(i) = wfn(i) + inner(j)*H(i,j)*exp(-ui*1.d0*t*E(j))
          end do
        end do 
        !$omp end parallel do

        do i=1,whole_HSS
          do j=1,p
            MF(cf(i,j)) = MF(cf(i,j)) + vs*(abs(wfn(i)))**2
          end do
        end do

        do j=1,L
          do k=1,L
            do i=0,Nb-1
              rho(j,k) = rho(j,k) + conjg(wfn(i*L+k))*wfn(i*L+j)
            end do
          end do
        end do

        lmda(:) = 0.d0
        call zheev('V','U',L,rho,L,lmda,rhowork,2*L-1,work,info)


        do i=1,L
          write(9,313) real(rho(i,:))
          write(10,313) aimag(rho(i,:))
        end do


        write(7,313) t , lmda(:)
      end do

      MF = MF / dt
      write(filename,'(a,i0,a)')'Mean_Field_',iii_it,'.dat'
      open(unit=14,file=filename)
      write(14,313) MF(:)


      close(7)
      close(9)
      close(14)
      close(10)
    end do
  end if

  1010 continue
  313 format(10000000E22.13)
  343 format(A,E15.5)
  353 format(A,F9.3)
  deallocate(cf)
  deallocate(H)
  write(*,*) "   LEAVE ME ALONE"
end subroutine

integer*8 function bn(x,y)
  integer :: x , y
  integer*8 :: fac1 , fac2 , b , a , i_count
  a = max(y,x-y)
  b = min(y,x-y)
  i = 0
  fac1 = 1
  fac2 = 1
  do i_count=a+1,x
    fac1 = fac1 * i_count
  end do
  do i_count=1,b
    fac2 = fac2 * i_count
  end do
  bn = abs(fac1/fac2)
end function

