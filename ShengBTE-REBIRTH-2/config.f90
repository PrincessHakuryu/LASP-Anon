!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2023 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2012-2023 Jesús Carrete Montaña <jcarrete@gmail.com>
!  Copyright (C) 2012-2023 Nebil Ayape Katcho <nebil.ayapekatcho@cea.fr>
!  Copyright (C) 2012-2023 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
!  Copyright (C) 2021-2022 Fanchen Meng <fanchem@g.clemson.edu>
!  Copyright (C) 2022-2023 Ben Durham <bd740@york.ac.uk>
!  Copyright (C) 2023 Martí Raya Moreno <marti_94@hotmail.com>
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Variables and routines used to read and store the configuration.
module config
  use iso_fortran_env
  use misc
  use data
  use symmetry
  use mpi
  implicit none

  integer :: nelements,natoms,ngrid(3),norientations
  namelist /allocations/ nelements,natoms,ngrid,norientations
  real(kind=dp) :: lfactor,lattvec(3,3),epsilon(3,3),alpha,beta_pulay
  character(len=3),allocatable :: elements(:)
  integer,allocatable :: types(:),orientations(:,:)
  integer :: scell(3)
  real(kind=dp),allocatable :: positions(:,:),born(:,:,:)
  real(kind=dp),allocatable :: masses(:),gfactors(:)
  namelist /crystal/ lfactor,lattvec,elements,types,positions,masses,gfactors,&
       epsilon,born,scell,orientations
  integer       :: maxiter,nticks,mmax
  real :: T,scalebroad,rmin,rmax,dr,eps
  real :: T_min,T_max,T_step,omega_max
  namelist /parameters/ T,scalebroad,rmin,rmax,dr,maxiter,nticks,eps,&
           T_min,T_max,T_step,omega_max,alpha,beta_pulay,mmax
  logical :: nonanalytic,convergence,isotopes,autoisotopes,nanowires,onlyharmonic,espresso,itergpu,DIIS
  
  ! Number of threads requested by user
  ! Later overwritten to number of threads systems allows
  integer :: nthreads

  namelist /flags/ nonanalytic,convergence,isotopes,autoisotopes,&
       nanowires,onlyharmonic,espresso,nthreads,itergpu,DIIS

   !Type for isotopic scattering
   type isotopic_information
        ! Information about phonons
        integer, allocatable      :: Indof1stPhononIso(:)
        integer, allocatable      :: Indof2ndPhononIso(:)
        ! Matrix elements for isotopic processes in given proc
        real(kind=dp),allocatable :: Gamma_isotopic(:)
        ! Number of isotopic processes in given proc
        integer :: nisotopic = 0
   end type isotopic_information



  integer       :: nbands,nptk,nwires
  real(kind=dp) :: cgrid,V,rV,rlattvec(3,3),slattvec(3,3)
  real(kind=sp) :: rlattvec_sp(3,3)
  real(kind=dp),allocatable :: cartesian(:,:),uorientations(:,:)

  ! Factor common to all 3 phonon processes 1/sqrt(M_i*M_j*M_k)
  real(kind=dp),allocatable :: inv_root_mmm(:,:,:)

  ! symmetry
  integer                   :: nsymm,nsymm_rot
  integer,      allocatable :: rotations(:,:,:)
  integer,      allocatable :: rotations_orig(:,:,:)
  real(kind=dp),allocatable :: crotations(:,:,:),qrotations(:,:,:)
  real(kind=dp),allocatable :: crotations_orig(:,:,:)
  real(kind=dp),allocatable :: qrotations_orig(:,:,:)
  real(kind=dp),allocatable :: symmetrizers(:,:,:)
  character(len=10)         :: international
  ! Vp matrix
  real(kind=dp),allocatable :: Vp_plus_matrix(:), Vp_minus_matrix(:)
  ! Ind
  integer,      allocatable :: Indof2ndPhonon_plus(:),Indof3rdPhonon_plus(:)
  integer,      allocatable :: Indof2ndPhonon_minus(:),Indof3rdPhonon_minus(:)
  real(kind=dp),allocatable :: Gamma_plus(:),Gamma_minus(:)
  ! Index of Vp, Ind*, Gamma arrays for a particular state on this mpi process
  integer(kind=long),allocatable :: Naccum_plus_array(:),Naccum_minus_array(:)
  ! MPI variables, assigned in ShengBTE.f90.
  integer :: myid,numprocs,nstates
  ! Isotopic information
  type(isotopic_information) :: IsoInfo
!$acc declare create( rlattvec_sp,rlattvec,ngrid, nbands, nptk, myid, nstates, omega_max, scalebroad)
contains

  subroutine read_config()
    implicit none

    integer :: i,j,k,ii,jj,kk,ll,info,ierr
    integer :: P(3)
    integer,allocatable :: rtmp(:,:,:),ID_equi(:,:)
    logical,allocatable :: valid(:)
    real(kind=dp) :: dnrm2,tmp1(3,3),tmp2(3,3),tmp3(3,3)
    real(kind=dp),allocatable :: crtmp(:,:,:),qrtmp(:,:,:)
    real(kind=dp),allocatable :: translations(:,:),ctranslations(:,:)

    ! Set the defaults and read the namelist.
    nelements=0
    natoms=0
    ngrid=0
    norientations=0
    open(1,file="CONTROL",status="old")
    read(1,nml=allocations)
    if(norientations.lt.0) then
       if(myid.eq.0)write(error_unit,*) "Error: norientations must be >=0"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    if(nelements.lt.1.or.natoms.lt.1.or.natoms.lt.nelements) then
       if(myid.eq.0)write(error_unit,*) "Error: nelements,natoms must be >0, natoms must be >=nelements"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    if(any(ngrid.lt.1)) then
       if(myid.eq.0)write(error_unit,*) "Error: all components of ngrid must be must be >0"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    allocate(elements(nelements),types(natoms),positions(3,natoms),&
         masses(nelements),gfactors(nelements),born(3,3,natoms),&
         cartesian(3,natoms))
    if(norientations.ne.0) then
       allocate(orientations(3,norientations))
       allocate(uorientations(3,norientations))
       orientations=0
    end if
    lfactor=1.0_dp
    scell=-1
    epsilon=0.0_dp
    epsilon(1,1)=1.0_dp
    epsilon(2,2)=1.0_dp
    epsilon(3,3)=1.0_dp
    born=0.0_dp
    gfactors=0.0_dp
    types=0
    read(1,nml=crystal)
    if(.not.all(scell.gt.0)) then
       if(myid.eq.0)write(error_unit,*) "Error: all supercell sizes must be >0"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    if(.not.all(types.ne.0)) then
       if(myid.eq.0)write(error_unit,*) "Error: atom types must be initialized correctly"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    do i=1,norientations
       if(all(orientations(:,i).eq.0)) then
          if(myid.eq.0)write(error_unit,*) "Error: orientation uninitialized or zero"
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call MPI_FINALIZE(ierr)
       end if
    end do
    T=0.0_dp
    T_min=0.0_dp
    scalebroad=1.0_dp

    rmin=5.0_dp
    rmax=505.0_dp
    dr=100.0_dp
    maxiter=1000
    nticks=100
    eps=1e-5_dp
    omega_max=1.0e30_dp
    alpha=0.5
    mmax=6
    beta_pulay=0.2
    read(1,nml=parameters)
    if ((T.le.0.0_dp).and.(T_min.le.0.0_dp)) then
       if(myid.eq.0)write(error_unit,*) "Error: T must be >0 K"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    elseif (T.gt.0.0_dp) then 
      T_min=T
      T_max=T
      T_step=T
    end if
    if(rmin.le.0.0_dp.or.rmax.le.rmin.or.dr.le.0.0_dp) then
       if(myid.eq.0)write(error_unit,*) "Error: rmin and dr must be >0, and rmax must be > rmin"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    if(maxiter.le.0) then
       if(myid.eq.0)write(error_unit,*) "Error: maxiter must be >0"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    nonanalytic=.true.
    convergence=.true.
    isotopes=.true.
    autoisotopes=.true.
    nanowires=.false.
    onlyharmonic=.false.
    espresso=.false.
    itergpu=.true.
    DIIS=.true.
    nthreads=1
    read(1,nml=flags)
    if(nanowires.and.norientations.eq.0) then
       if(myid.eq.0)write(error_unit,*) "Error: nanowires=.TRUE. but norientations=0"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    close(1)

    nptk=product(ngrid) !网格总数
    cgrid=real(nptk,kind=dp)**(1.0_dp/3.0_dp) !每个方向的平均网格数
    nbands=3*natoms !一共三N个简正频率 3N个自由度
    nwires=ceiling((rmax-rmin)/dr)!纳米线 暂时不管

    lattvec=lfactor*lattvec !单位转换系数

    do i=1,3
       j=mod(i,3)+1
       k=mod(j,3)+1
       call cross_product(lattvec(:,j),lattvec(:,k),rlattvec(:,i))  !构造了倒晶格方向矢量
    end do

    V=abs(dot_product(lattvec(:,1),rlattvec(:,1))) !归一化 b1 = 2π (a2 × a3) / (a1 · (a2 × a3))
    rV=2.0_dp*pi/V
    rlattvec=rV*rlattvec
    rlattvec_sp = real(rlattvec,sp)
    cartesian=matmul(lattvec,positions) !输入的是分数坐标 这是变成笛卡尔坐标

    if(nanowires) then !纳米线不管
       uorientations=matmul(lattvec,real(orientations,kind=dp))
       do i=1,norientations
          uorientations(:,i)=uorientations(:,i)/&
               dnrm2(3,uorientations(:,i),1)
       end do
    end if

    ! Compute the average masses and g-factors automatically, if
    ! requested.
    if(autoisotopes) then !计算同位素效应 获取了masses（每种元素对应的质量数组）和gfactors
       call data_fill_isotopes()
       do i=1,nelements
          call data_calc_mandg(elements(i),masses(i),gfactors(i))
       end do
       call data_free_isotopes()
    end if

    allocate(inv_root_mmm(nelements,nelements,nelements)) !预计算好不同元素之间的质量约化矩阵 用于后续的三阶力计算
    do i = 1, nelements
      do j = 1, nelements
        do k = 1, nelements
          inv_root_mmm(k,j,i) = 1.0_dp/sqrt(masses(k)*masses(j)*masses(i))
        end do
      end do
    end do

    ! Find out the symmetries of the system.
    nsymm=get_num_operations(lattvec,natoms,types,positions) !types是每个原子对应的元素标号 1 1 1 2 2 2 类似 
    ! We need do double that number to take time reversal symetry into
    ! account.
    nsymm_rot=2*nsymm
    allocate(translations(3,nsymm),rotations(3,3,nsymm_rot),&
         ctranslations(3,nsymm),crotations(3,3,nsymm_rot),&
         qrotations(3,3,nsymm_rot))
    allocate(rotations_orig(3,3,nsymm),&
         crotations_orig(3,3,nsymm),&
         qrotations_orig(3,3,nsymm),symmetrizers(3,3,nptk))
    call get_operations(lattvec,natoms,types,positions,nsymm,&
         rotations_orig,translations,international)
    rotations(:,:,1:nsymm)=rotations_orig
    if(myid.eq.0)write(*,*) "Info: symmetry group ",trim(international)," detected"
    if(myid.eq.0)write(*,*) "Info: ",nsymm," symmetry operations"
    call get_cartesian_operations(lattvec,nsymm,rotations_orig,translations,&
         crotations_orig,ctranslations)
    crotations(:,:,1:nsymm)=crotations_orig
    deallocate(translations,ctranslations)

    ! Transform the rotation matrices to the reciprocal-space basis.
    do i=1,nsymm
       tmp1=matmul(transpose(lattvec),lattvec)
       tmp2=transpose(rotations_orig(:,:,i))
       tmp3=tmp1
       call dgesv(3,3,tmp1,3,P,tmp2,3,info)
       qrotations_orig(:,:,i)=transpose(matmul(tmp2,tmp3))
    end do
    qrotations(:,:,1:nsymm)=qrotations_orig

    ! Fill the second half of the rotation matrix list
    ! using time reversal symmetry.
    rotations(:,:,nsymm+1:2*nsymm)=-rotations_orig(:,:,1:nsymm)
    qrotations(:,:,nsymm+1:2*nsymm)=-qrotations_orig(:,:,1:nsymm)
    crotations(:,:,nsymm+1:2*nsymm)=-crotations_orig(:,:,1:nsymm)

    allocate(ID_Equi(nsymm_rot,nptk))
    call symmetry_map_notransl(ID_equi)
    ! Create the "symmetrizers", linear operators that extract the
    ! component of a vector compatible with the symmetries at each q point.
    symmetrizers=0.0_dp
    do ii=1,nptk
       kk=0
       do jj=1,nsymm_rot
          if(ID_equi(jj,ii).eq.ii) then
             symmetrizers(:,:,ii)=symmetrizers(:,:,ii)+&
                  crotations(:,:,jj)
             kk=kk+1
          end if
       end do
       if(kk.gt.1) then
          symmetrizers(:,:,ii)=symmetrizers(:,:,ii)/real(kk,kind=dp)
       end if
    end do
    ! Find rotations that are either duplicated or incompatible with 有些时候 原子坐标得出的对称操作并不能完美的映射到q网络上  因此要舍弃一定的对称性
    ! the q-point grid.    
    call symmetry_map(ID_equi)
    allocate(valid(nsymm_rot))
    valid=.TRUE.
    jj=0
    do ii=1,nsymm_rot
       if(valid(ii).and.any(ID_equi(ii,:).eq.-1)) then
          valid(ii)=.FALSE.
          jj=jj+1
       end if
    end do
    if(myid.eq.0.and.jj.ne.0) then
       write(*,*) "Info:",jj,&
            "rotations are incompatible with the q-point grid and will be discarded"
    end if
    ll=0
    do ii=2,nsymm_rot
       do i=1,ii-1
          if(.not.valid(i))cycle
          if(all(rotations(:,:,ii).eq.rotations(:,:,i))) then
             valid(ii)=.FALSE.
             ll=ll+1
             exit
          end if
       end do
    end do
    if(myid.eq.0.and.ll.ne.0) then
       write(*,*) "Info:",ll,"duplicated rotations will be discarded"
    end if
    ! Filter out those rotations through a series of move_alloc calls.
    ! Arrays to take into account: rotations,crotations,qrotations.
    if(ll+jj.ne.0) then
       allocate(rtmp(3,3,nsymm_rot-ll-jj))
       allocate(crtmp(3,3,nsymm_rot-ll-jj))
       allocate(qrtmp(3,3,nsymm_rot-ll-jj))
       kk=0
       do ii=1,nsymm_rot
          if(valid(ii)) then
             kk=kk+1
             rtmp(:,:,kk)=rotations(:,:,ii)
             crtmp(:,:,kk)=crotations(:,:,ii)
             qrtmp(:,:,kk)=qrotations(:,:,ii)
          end if
       end do
       nsymm_rot=nsymm_rot-ll-jj
       call move_alloc(rtmp,rotations)
       call move_alloc(crtmp,crotations)
       call move_alloc(qrtmp,qrotations)
    end if
    deallocate(ID_Equi,valid)
  end subroutine read_config

  ! Free the memory used by all config structures.
  subroutine free_config()
    deallocate(elements,types,positions,masses,gfactors,born,cartesian,&
         rotations,crotations,qrotations,symmetrizers,inv_root_mmm)
    if(nanowires) then
       deallocate(orientations,uorientations)
    end if
  end subroutine free_config

  ! Compute all images of a reciprocal-space point using the
  ! rotational part of the symmetry operations. Everything is
  ! performed in lattice coordinates.
  subroutine symm(r_in,r_out)
    implicit none
    integer,intent(in) :: r_in(3)
    real(kind=dp),intent(out) :: r_out(3,nsymm_rot)

    integer :: ii

    do ii=1,nsymm_rot
       r_out(:,ii)=real(ngrid,kind=dp)*matmul(qrotations(:,:,ii),real(r_in,kind=dp)/real(ngrid,kind=dp))
    end do
  end subroutine symm

  ! Find the equivalences among points.
  subroutine symmetry_map(ID_equi)
    implicit none
    integer,intent(out) :: ID_equi(nsymm_rot,nptk)

    integer :: Ind_cell(3,nptk)
    integer :: i,isym,ivec(3)
    real(kind=dp) :: vec(3),vec_symm(3,nsymm_rot),dnrm2

    call Id2Ind(Ind_cell)
    do i=1,nptk
       call symm(Ind_cell(:,i),vec_symm)
       do isym=1,nsymm_rot
          vec=vec_symm(:,isym)
          ivec=nint(vec)
          if(dnrm2(3,abs(vec-real(ivec,kind=dp)),1).gt.1.0e-2_dp) then
             ID_equi(isym,i)=-1
          else
             ID_equi(isym,i)=Ind2Id(modulo(ivec,ngrid))
          end if
       end do
    end do
  end subroutine symmetry_map

  ! Find symmetry operations excluding translations that can bring q to itself.
  subroutine symmetry_map_notransl(ID_equi)
    implicit none
    integer,intent(out) :: ID_equi(nsymm_rot,nptk)

    integer :: Ind_cell(3,nptk)
    integer :: i,isym,ivec(3)
    real(kind=dp) :: vec(3),vec_symm(3,nsymm_rot),dnrm2

    call Id2Ind(Ind_cell)
    do i=1,nptk
       call symm(Ind_cell(:,i),vec_symm)
       do isym=1,nsymm_rot
          vec=vec_symm(:,isym)
          ivec=nint(vec)
          if(dnrm2(3,abs(real(Ind_cell(:,i),kind=dp)-vec),1).le.1.0e-5_dp) then
             ID_equi(isym,i)=i
          else
             ID_equi(isym,i)=-1
          endif
       end do
    end do
  end subroutine symmetry_map_notransl

  ! Create a table that can be used to demultiplex cell indices.
  subroutine Id2ind(Ind_cell)
    implicit none
    integer,intent(out) :: Ind_cell(3,nptk)

    integer :: ii,tmp1

    do ii=1,nptk
       call divmod(ii-1,Ngrid(1),tmp1,Ind_cell(1,ii))
       call divmod(tmp1,Ngrid(2),Ind_cell(3,ii),Ind_cell(2,ii))
    end do
  end subroutine Id2ind

  ! Multiplex three cell indices into one.
  function Ind2Id(Ind_cell)
    implicit none
    integer,intent(in) :: Ind_cell(3)

    integer :: Ind2Id

    Ind2Id=1+Ind_cell(1)+(Ind_cell(2)+Ind_cell(3)*Ngrid(2))*Ngrid(1)
  end function Ind2Id

  ! Return the base broadening (without prefactor) for a mode.

  function base_sigma(v)
    implicit none
    real(kind=dp),intent(in) :: v(3)

    real(kind=dp) :: base_sigma

    integer :: nu

    base_sigma=0.0_dp
    do nu=1,3
       base_sigma=base_sigma+(dot_product(rlattvec(:,nu),v)/real(ngrid(nu),kind=dp))**2
    end do

    base_sigma=sqrt(base_sigma/6.0_dp)
  end function base_sigma
  

  ! Force a real 3x3 Cartesian tensor to fulfill all the symmetries.
  subroutine symmetrize_tensor(tensor)
    implicit none

    real(kind=dp),intent(inout) :: tensor(3,3)

    integer :: isym
    real(kind=dp) :: tmp(3,3)

    tmp = 0.0_dp
    do isym=1,nsymm_rot
       tmp = tmp + matmul(crotations(:,:,isym),&
            matmul(tensor,transpose(crotations(:,:,isym))))
    end do

    tensor=tmp/real(nsymm_rot,kind=dp)
  end subroutine symmetrize_tensor
end module config
