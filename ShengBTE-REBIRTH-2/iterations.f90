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

! Code used to initialize F_n for the zeroth-order iteration of the
! BTE and to perform successive iterations.
module iterations
  use data
  use config
  use wedgetc
  use mpi
  use, intrinsic :: ieee_arithmetic
  implicit none
!================ 持久化（仅主机端） ================
  logical :: host_inited    = .false.
  integer :: saved_nbands   = 0
  integer :: saved_nptk     = 0
  integer :: saved_nlist    = 0
  integer :: saved_nsymmrot = 0
  integer :: IsoN_saved     = 0
  integer :: sz_Gp_saved    = 0
  integer :: sz_Gm_saved    = 0
  logical :: itergpu_debug_header_printed = .false.

  

  ! 主机端持久化影子与辅助表（均在主机侧，仅在需要时更新）
  integer, allocatable :: ID_equi(:,:)         ! (Nsymm_rot, nptk)

  real(sp), allocatable :: omega_sp(:,:)       ! (Nbands, nptk)
  real(sp), allocatable :: velocity_sp(:,:,:)  ! (3, Nbands, nptk)
  real(sp), allocatable :: tau_zero_sp(:,:)    ! (Nbands, Nlist)

  real(sp), allocatable :: F_sp(:,:,:)         ! (3, Nbands, nptk)
  real(sp), allocatable :: DeltaF_sp(:,:,:)    ! (3, Nbands, nptk)

  real(sp), allocatable :: Gamma_plus_sp(:), Gamma_minus_sp(:), Gamma_iso_sp(:)
  integer,  allocatable :: Iso_Ind1(:), Iso_Ind2(:)

contains
  subroutine iterationgpu_host_finalize()
    implicit none
    if (.not. host_inited) return

    if (allocated(ID_equi))        deallocate(ID_equi)
    if (allocated(omega_sp))       deallocate(omega_sp)
    if (allocated(velocity_sp))    deallocate(velocity_sp)
    if (allocated(tau_zero_sp))    deallocate(tau_zero_sp)
    if (allocated(F_sp))           deallocate(F_sp)
    if (allocated(DeltaF_sp))      deallocate(DeltaF_sp)
    if (allocated(Gamma_plus_sp))  deallocate(Gamma_plus_sp)
    if (allocated(Gamma_minus_sp)) deallocate(Gamma_minus_sp)
    if (allocated(Gamma_iso_sp))   deallocate(Gamma_iso_sp)
    if (allocated(Iso_Ind1))       deallocate(Iso_Ind1)
    if (allocated(Iso_Ind2))       deallocate(Iso_Ind2)

    host_inited    = .false.
    saved_nbands   = 0
    saved_nptk     = 0
    saved_nlist    = 0
    saved_nsymmrot = 0
    IsoN_saved     = 0
    sz_Gp_saved    = 0
    sz_Gm_saved    = 0
  end subroutine iterationgpu_host_finalize
  ! Fill F_n with its initial values for the iteration, computed from the
  ! RTA approximation to tau.
  subroutine iteration0(Nlist,Nequi,ALLEquiList,omega,velocity,tau_zero,F_n)
    implicit none

    integer,intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm_rot,nptk)
    real(kind=dp),intent(in) :: omega(Nbands,nptk),velocity(3,Nbands,nptk),tau_zero(Nbands,Nlist)
    real(kind=dp),intent(out) :: F_n(3,Nbands,nptk)

    integer :: ii,kk,ll

    !$OMP PARALLEL DO default(none) schedule(static) shared(tau_zero,omega,velocity) &
    !$OMP & shared(F_n,Nlist,nbands,Nequi,ALLEquiList) private(ll,ii,kk)
    do ll=1,Nlist
       do kk=1,Nequi(ll)
          do ii=1,Nbands
             F_n(:,ii,ALLEquiList(kk,ll))=tau_zero(ii,ll)*velocity(:,ii,ALLEquiList(kk,ll))*&
                  omega(ii,ALLEquiList(kk,ll))
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine iteration0

  !Advance the algorithm one iteration. F_n is used both as the input
  !and as the output.
  subroutine iteration(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,&
                      & omega,velocity,tau_zero,F_n)
    implicit none
    integer,intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm_rot,nptk)
    integer,intent(in) :: TypeofSymmetry(Nsymm_rot,nptk)
    integer,intent(in) :: N_plus(Nlist*Nbands),N_minus(Nlist*Nbands)
    real(kind=dp),intent(in) :: omega(nbands,nptk),velocity(3,nbands,nptk)
    real(kind=dp),intent(in) :: tau_zero(nbands,nlist)
    real(kind=dp),intent(inout) :: F_n(3,Nbands,nptk)

    integer :: ID_equi(Nsymm_rot,nptk)
    integer (kind=long) :: Naccum_plus,Naccum_minus
    integer :: i,j,k,jj,kk,ll,mm,nn,mmm,nnn
    real(kind=dp) :: DeltaF(3,Nbands,nptk)

    call symmetry_map(ID_equi)
    DeltaF=0.0_dp
    !$OMP PARALLEL DO default(none) schedule(dynamic,1) shared(nstates,myid,nbands,nlist) &
    !$OMP & shared(Nequi,F_n,ID_equi,TypeofSymmetry,ALLEquiList,DeltaF) &
    !$OMP & shared(N_plus,Naccum_plus_array,Indof2ndPhonon_plus,Indof3rdPhonon_plus,Gamma_plus) &
    !$OMP & shared(N_minus,Naccum_minus_array,Indof2ndPhonon_minus,Indof3rdPhonon_minus,Gamma_minus) &
    !$OMP & private(i,j,k,kk,ll,mm,nn,mmm,nnn,Naccum_plus,Naccum_minus)
    do nnn=1,nstates
        mmm=myid*nstates+nnn
        i=modulo(mmm-1,Nbands)+1
        ll=int((mmm-1)/Nbands)+1

        if (mmm.gt.nlist*nbands) cycle

        Naccum_plus=Naccum_plus_array(nnn)
        Naccum_minus=Naccum_minus_array(nnn)

        do kk=1,Nequi(ll)
           if ((N_plus(mmm).ne.0)) then
              do jj=1,N_plus(mmm) 
                 j=modulo(Indof2ndPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                 mm=(Indof2ndPhonon_plus(Naccum_plus+jj)-1)/Nbands + 1
                 k=modulo(Indof3rdPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                 nn=(Indof3rdPhonon_plus(Naccum_plus+jj)-1)/Nbands + 1
                 DeltaF(:,i,ALLEquiList(kk,ll))=DeltaF(:,i,ALLEquiList(kk,ll))+&               !这里是要求所有Q点 不只是不可约Q点了
                      Gamma_plus(Naccum_plus+jj)*(F_n(:,k,ID_equi(TypeofSymmetry(kk,ll),nn))-& !TypeofSymmetry保存的是旋转操作的对象ll到kk是什么旋转过去的
                      F_n(:,j,ID_equi(TypeofSymmetry(kk,ll),mm)))                              !然后ID_equi就可以把对应的lamda‘’也旋转回去 
              end do !jj 相当于蹭了不可约q点的旋转方式回到了不可约q点转回的q点对应的q''点上
           end if
           if ((N_minus(mmm).ne.0)) then
              do jj=1,N_minus(mmm)
                 j=modulo(Indof2ndPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                 mm=(Indof2ndPhonon_minus(Naccum_minus+jj)-1)/Nbands + 1
                 k=modulo(Indof3rdPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                 nn=(Indof3rdPhonon_minus(Naccum_minus+jj)-1)/Nbands + 1
                 DeltaF(:,i,ALLEquiList(kk,ll))=DeltaF(:,i,ALLEquiList(kk,ll))+&
                      Gamma_minus(Naccum_minus+jj)*(F_n(:,k,ID_equi(TypeofSymmetry(kk,ll),nn))+&
                      F_n(:,j,ID_equi(TypeofSymmetry(kk,ll),mm)))*0.5_dp
              end do !jj
           end if
        end do !kk
    end do !nnn
    !$OMP END PARALLEL DO

    ! Isotopic contribution to iterative value
    if ((IsoInfo%nisotopic .ne. 0)) then
        !$OMP PARALLEL DO default(none) schedule(static) shared(nbands,IsoInfo) &
        !$OMP & shared(Nequi,F_n,ID_equi,TypeofSymmetry,ALLEquiList,DeltaF) &
        !$OMP & private(i,j,kk,ll,nn,jj)
        do jj=1,IsoInfo%nisotopic
            ! q and ib of input phonon
            i  = modulo(IsoInfo%Indof1stPhononIso(jj)-1,Nbands)+1
            ll = int((IsoInfo%Indof1stPhononIso(jj)-1)/Nbands)+1
            ! q and ib of outcome phonon
            j  = modulo(IsoInfo%Indof2ndPhononIso(jj)-1,Nbands)+1
            kk = int((IsoInfo%Indof2ndPhononIso(jj)-1)/Nbands)+1
            !Iterate over equivalence class of first phonon
            do nn=1,Nequi(ll)
               DeltaF(:,i,ALLEquiList(nn,ll))=DeltaF(:,i,ALLEquiList(nn,ll))+&
                  IsoInfo%gamma_isotopic(jj) * F_n(:,j,ID_equi(TypeofSymmetry(nn,ll),kk))
            end do !nn
        end do !jj
        !$OMP END PARALLEL DO
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,mmm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,DeltaF,nptk*nbands*3,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,mm)

         

    !$OMP PARALLEL DO default(none) schedule(static) shared(tau_zero,omega,velocity,ALLEquiList) &
    !$OMP & shared(F_n,DeltaF,Nlist,nbands,Nequi) private(ll,i,kk)
    do ll=1,Nlist
       do kk=1,Nequi(ll)
          do i=1,Nbands
             F_n(:,i,ALLEquiList(kk,ll))=tau_zero(i,ll)*velocity(:,i,ALLEquiList(kk,ll))*&
                   omega(i,ALLEquiList(kk,ll))+tau_zero(i,ll)*DeltaF(:,i,ALLEquiList(kk,ll))
          enddo
       enddo
    enddo  
    !$OMP END PARALLEL DO       

  end subroutine iteration


subroutine iterationgpu(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus, &
                       & omega,velocity,tau_zero,F_n,start_time)
    implicit none
    ! 形参（保持双精度接口）
    integer,              intent(in)    :: Nlist
    integer,              intent(in)    :: Nequi(Nlist), ALLEquiList(Nsymm_rot,nptk)
    integer,              intent(in)    :: TypeofSymmetry(Nsymm_rot,nptk)
    integer,              intent(in)    :: N_plus(Nlist*Nbands), N_minus(Nlist*Nbands)
    real(kind=dp),        intent(in)    :: omega(Nbands,nptk), velocity(3,Nbands,nptk)
    real(kind=dp),        intent(in)    :: tau_zero(Nbands,Nlist)
    real(kind=dp),        intent(in)    :: start_time
    real(kind=dp),        intent(inout) :: F_n(3,Nbands,nptk)

    ! 局部量
    integer(kind=long) :: Naccum_plus, Naccum_minus
    integer :: i,j,k,jj,kk,ll,mm,nn,mmm,nnn,c,q_idx,idx_sym_j,idx_sym_k,sym_idx
    integer :: ierr, Nequi_max, IsoN, N_plus_max, N_minus_max
    integer :: bad_map_count, bad_sym_count, bad_q_count
    integer, allocatable :: TypeofSymmetry_safe(:,:), ALLEquiList_safe(:,:)
    integer :: finite_local, finite_global
    real(kind=dp) :: deltaf_absmax_local, deltaf_absmax_global
    integer :: c_bad, i_bad, q_bad, dbg_state, dbg_kk, dbg_jj
    real(kind=sp) :: dbg_term
    logical :: dbg_found

    !---------------- 主机端：一次性分配/构建 或 当规模变化时重建 ----------------
    if (.not. host_inited .or. saved_nbands/=Nbands .or. saved_nptk/=nptk .or.  &
        saved_nlist/=Nlist .or. saved_nsymmrot/=Nsymm_rot .or.                 &
        sz_Gp_saved/=size(Gamma_plus) .or. sz_Gm_saved/=size(Gamma_minus) .or. &
        IsoN_saved /= IsoInfo%nisotopic) then

      call iterationgpu_host_finalize()

      saved_nbands   = Nbands
      saved_nptk     = nptk
      saved_nlist    = Nlist
      saved_nsymmrot = Nsymm_rot
      sz_Gp_saved    = size(Gamma_plus)
      sz_Gm_saved    = size(Gamma_minus)
      IsoN_saved     = IsoInfo%nisotopic

      ! 对称映射（只做一次或规模变化时重建）
      allocate(ID_equi(Nsymm_rot, nptk))
      call symmetry_map(ID_equi)

      ! 常量影子（仅构建一次）
      allocate(omega_sp(Nbands,nptk))
      allocate(velocity_sp(3,Nbands,nptk))
      allocate(tau_zero_sp(Nbands,Nlist))
      omega_sp    = real(omega,    kind=sp)
      velocity_sp = real(velocity, kind=sp)
      tau_zero_sp = real(tau_zero, kind=sp)

      ! 逐步更新的影子（每步前会轻量刷新 F_sp、清零 DeltaF_sp）
      allocate(F_sp(3,Nbands,nptk))
      allocate(DeltaF_sp(3,Nbands,nptk))

      ! 散射表影子（一次性）
      allocate(Gamma_plus_sp(sz_Gp_saved))
      allocate(Gamma_minus_sp(sz_Gm_saved))
      Gamma_plus_sp  = real(Gamma_plus,  kind=sp)
      Gamma_minus_sp = real(Gamma_minus, kind=sp)

      if (IsoN_saved > 0) then
        allocate(Iso_Ind1(IsoN_saved), Iso_Ind2(IsoN_saved), Gamma_iso_sp(IsoN_saved))
        Iso_Ind1     = IsoInfo%Indof1stPhononIso
        Iso_Ind2     = IsoInfo%Indof2ndPhononIso
        Gamma_iso_sp = real(IsoInfo%Gamma_isotopic, kind=sp)
      end if

      host_inited = .true.
    end if

    

    !---------------- 每次调用：轻量刷新/清零（主机端） ----------------
    ! if (myid.eq.0) write(*,'(A,F0.3,A)') " Info: host-prep begin at ", MPI_Wtime()-start_time, " seconds"
    F_sp      = real(F_n, kind=sp)   ! 用上一步主机态的 F_n 更新本步初值
    DeltaF_sp = 0.0_sp
    IsoN      = IsoN_saved
    ! if (myid.eq.0) write(*,'(A,F0.3,A)') " Info: host-prep end at   ", MPI_Wtime()-start_time, " seconds"

    Nequi_max   = maxval(Nequi)
    N_plus_max  = maxval(N_plus)
    N_minus_max = maxval(N_minus)
    bad_map_count = 0
    bad_sym_count = 0
    bad_q_count   = 0
    allocate(TypeofSymmetry_safe(size(TypeofSymmetry,1),size(TypeofSymmetry,2)))
    allocate(ALLEquiList_safe(size(ALLEquiList,1),size(ALLEquiList,2)))
    TypeofSymmetry_safe = TypeofSymmetry
    ALLEquiList_safe    = ALLEquiList
    do ll=1,size(TypeofSymmetry,2)
      do kk=1,size(TypeofSymmetry,1)
        if (TypeofSymmetry_safe(kk,ll) < 1 .or. TypeofSymmetry_safe(kk,ll) > Nsymm_rot) then
          TypeofSymmetry_safe(kk,ll) = 1
          bad_sym_count = bad_sym_count + 1
        end if
        if (ALLEquiList_safe(kk,ll) < 1 .or. ALLEquiList_safe(kk,ll) > nptk) then
          ALLEquiList_safe(kk,ll) = 1
          bad_q_count = bad_q_count + 1
        end if
      end do
    end do
    bad_map_count = bad_sym_count + bad_q_count
    call MPI_ALLREDUCE(MPI_IN_PLACE,bad_map_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    if (myid.eq.0 .and. bad_map_count.gt.0) then
      write(*,*) "WARN_ITERGPU: invalid mapping entries detected(sym+q) =",bad_map_count
    end if

    if (myid.eq.0 .and. .not. itergpu_debug_header_printed) then
      write(*,*) "DEBUG_ITERGPU_INPUT: min/max(TypeofSymmetry)=",minval(TypeofSymmetry),maxval(TypeofSymmetry)
      write(*,*) "DEBUG_ITERGPU_INPUT: min/max(ID_equi)=",minval(ID_equi),maxval(ID_equi)
      write(*,*) "DEBUG_ITERGPU_INPUT: min/max(ALLEquiList)=",minval(ALLEquiList),maxval(ALLEquiList)
      write(*,*) "DEBUG_ITERGPU_INPUT: min/max(N_plus,N_minus)=",minval(N_plus),maxval(N_plus),minval(N_minus),maxval(N_minus)
      write(*,*) "DEBUG_ITERGPU_INPUT: finite(Gamma_plus_sp,Gamma_minus_sp,tau_zero_sp,F_sp)=", &
           ieee_is_finite(sum(Gamma_plus_sp)), ieee_is_finite(sum(Gamma_minus_sp)), &
           ieee_is_finite(sum(tau_zero_sp)), ieee_is_finite(sum(F_sp))
      itergpu_debug_header_printed = .true.
    end if

    ! if (myid.eq.0) write(*,'(A,4(I0,1X))') " Info: dims Nequi_max nstates N+max N-max = ", &
    !                                        Nequi_max, nstates, N_plus_max, N_minus_max

    
    !$acc data                                                                 &
    !$acc& copyin(Nequi, ALLEquiList_safe, TypeofSymmetry_safe, ID_equi)                 &
    !$acc& copyin(N_plus, N_minus,                                            &
    !$acc&       Naccum_plus_array, Indof2ndPhonon_plus,                      &
    !$acc&       Indof3rdPhonon_plus, Gamma_plus_sp,                          &
    !$acc&       Naccum_minus_array, Indof2ndPhonon_minus,                    &
    !$acc&       Indof3rdPhonon_minus, Gamma_minus_sp)                        &
    !$acc& copyin(omega_sp, velocity_sp, tau_zero_sp)                          &
    !$acc& copy(DeltaF_sp, F_sp)

    if (IsoN > 0) then
      !$acc enter data copyin(Iso_Ind1, Iso_Ind2, Gamma_iso_sp)
    end if


    !$acc parallel loop gang collapse(2)                                      &
    !$acc& present(DeltaF_sp,Nequi,ALLEquiList_safe,ID_equi,F_sp,                   &
    !$acc&         N_plus,N_minus,                                            &
    !$acc&         Naccum_plus_array,Indof2ndPhonon_plus,                      &
    !$acc&         Indof3rdPhonon_plus,Gamma_plus_sp,                          &
    !$acc&         Naccum_minus_array,Indof2ndPhonon_minus,                    &
    !$acc&         Indof3rdPhonon_minus,Gamma_minus_sp, TypeofSymmetry_safe)
    do nnn = 1, nstates
      do kk  = 1, Nequi_max

        mmm = myid*nstates + nnn
        if (mmm <= Nlist*Nbands) then
          i  = modulo(mmm-1,Nbands) + 1
          ll = int((mmm-1)/Nbands) + 1

          if (kk <= Nequi(ll)) then
            sym_idx = TypeofSymmetry_safe(kk,ll)
            q_idx = ALLEquiList_safe(kk,ll)
            if (q_idx < 1 .or. q_idx > nptk .or. sym_idx < 1 .or. sym_idx > Nsymm_rot) then
              bad_map_count = bad_map_count + 1
              cycle
            end if
            Naccum_plus  = Naccum_plus_array(nnn)
            Naccum_minus = Naccum_minus_array(nnn)

            ! plus 通道
            if (N_plus(mmm) .ne. 0) then
              !$acc loop seq
              do jj = 1, N_plus(mmm)
                j  = modulo(Indof2ndPhonon_plus(Naccum_plus+jj)-1, Nbands) + 1
                mm = int((Indof2ndPhonon_plus(Naccum_plus+jj)-1)/Nbands) + 1
                k  = modulo(Indof3rdPhonon_plus(Naccum_plus+jj)-1, Nbands) + 1
                nn = int((Indof3rdPhonon_plus(Naccum_plus+jj)-1)/Nbands) + 1

                idx_sym_k = ID_equi(sym_idx, nn)
                idx_sym_j = ID_equi(sym_idx, mm)
                if (idx_sym_k < 1 .or. idx_sym_k > nptk .or. idx_sym_j < 1 .or. idx_sym_j > nptk) cycle
                !$acc loop seq
                do c = 1, 3
                  !$acc atomic update
                  DeltaF_sp(c,i,q_idx) = DeltaF_sp(c,i,q_idx) + &
                       Gamma_plus_sp(Naccum_plus+jj) * ( F_sp(c,k,idx_sym_k) - F_sp(c,j,idx_sym_j) )
                end do
              end do
            end if

            ! minus 通道
            if (N_minus(mmm) .ne. 0) then
              !$acc loop seq
              do jj = 1, N_minus(mmm)
                j  = modulo(Indof2ndPhonon_minus(Naccum_minus+jj)-1, Nbands) + 1
                mm = int((Indof2ndPhonon_minus(Naccum_minus+jj)-1)/Nbands) + 1
                k  = modulo(Indof3rdPhonon_minus(Naccum_minus+jj)-1, Nbands) + 1
                nn = int((Indof3rdPhonon_minus(Naccum_minus+jj)-1)/Nbands) + 1

                idx_sym_k = ID_equi(sym_idx, nn)
                idx_sym_j = ID_equi(sym_idx, mm)
                if (idx_sym_k < 1 .or. idx_sym_k > nptk .or. idx_sym_j < 1 .or. idx_sym_j > nptk) cycle
                !$acc loop seq
                do c = 1, 3
                  !$acc atomic update
                  DeltaF_sp(c,i,q_idx) = DeltaF_sp(c,i,q_idx) + &
                       0.5_sp*Gamma_minus_sp(Naccum_minus+jj) * ( F_sp(c,k,idx_sym_k) + F_sp(c,j,idx_sym_j) )
                end do
              end do
            end if

          end if
        end if

      end do
    end do
    !$acc end parallel

    if (IsoN .ne. 0) then
      !$acc parallel loop gang collapse(2)                                   &
      !$acc& present(DeltaF_sp, F_sp, ALLEquiList_safe, ID_equi,                  &
      !$acc&         TypeofSymmetry_safe, Nequi, Iso_Ind1, Iso_Ind2, Gamma_iso_sp)
      do jj = 1, IsoN
        do nn = 1, Nequi_max

          i  = modulo(Iso_Ind1(jj)-1, Nbands) + 1
          ll = int((Iso_Ind1(jj)-1)/Nbands) + 1

          if (nn <= Nequi(ll)) then
            j  = modulo(Iso_Ind2(jj)-1, Nbands) + 1
            kk = int((Iso_Ind2(jj)-1)/Nbands) + 1

            q_idx     = ALLEquiList_safe(nn, ll)
            sym_idx = TypeofSymmetry_safe(nn,ll)
            if (q_idx < 1 .or. q_idx > nptk .or. sym_idx < 1 .or. sym_idx > Nsymm_rot) cycle
            idx_sym_k = ID_equi(sym_idx, kk)
            if (idx_sym_k < 1 .or. idx_sym_k > nptk) cycle

            do c = 1, 3
              !$acc atomic update
              DeltaF_sp(c,i,q_idx) = DeltaF_sp(c,i,q_idx) + Gamma_iso_sp(jj) * F_sp(c,j,idx_sym_k)
            end do
          end if

        end do
      end do
      !$acc end parallel
    end if

    !$acc update self(DeltaF_sp)
    finite_local = merge(1,0,ieee_is_finite(sum(DeltaF_sp)))
    deltaf_absmax_local = maxval(abs(real(DeltaF_sp,kind=dp)))
    call MPI_ALLREDUCE(finite_local,finite_global,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(deltaf_absmax_local,deltaf_absmax_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (myid.eq.0) then
      if (finite_global.eq.0) then
        write(*,*) "ERROR_ITERGPU: at least one rank has non-finite DeltaF_sp before MPI_ALLREDUCE"
      end if
      write(*,*) "DEBUG_ITERGPU: global max|DeltaF_sp| before reduce =",deltaf_absmax_global
    end if
    if (finite_local.eq.0) then
      c_bad = -1; i_bad = -1; q_bad = -1
      do q_idx=1,nptk
        do i=1,Nbands
          do c=1,3
            if (.not. ieee_is_finite(DeltaF_sp(c,i,q_idx))) then
              c_bad=c; i_bad=i; q_bad=q_idx
              exit
            end if
          end do
          if (c_bad.gt.0) exit
        end do
        if (c_bad.gt.0) exit
      end do
      write(*,*) "ERROR_ITERGPU: rank",myid," local DeltaF_sp non-finite at (c,i,q)=",&
                 c_bad,i_bad,q_bad," value=",DeltaF_sp(c_bad,i_bad,q_bad)
      dbg_found = .false.
      do dbg_state = 1, nstates
        mmm = myid*nstates + dbg_state
        if (mmm <= Nlist*Nbands) then
          i  = modulo(mmm-1,Nbands) + 1
          ll = int((mmm-1)/Nbands) + 1
          if (i /= i_bad) cycle
          do dbg_kk = 1, Nequi(ll)
            if (ALLEquiList_safe(dbg_kk,ll) /= q_bad) cycle
            sym_idx = TypeofSymmetry_safe(dbg_kk,ll)
            Naccum_plus  = Naccum_plus_array(dbg_state)
            Naccum_minus = Naccum_minus_array(dbg_state)
            do dbg_jj = 1, N_plus(mmm)
              j  = modulo(Indof2ndPhonon_plus(Naccum_plus+dbg_jj)-1, Nbands) + 1
              mm = int((Indof2ndPhonon_plus(Naccum_plus+dbg_jj)-1)/Nbands) + 1
              k  = modulo(Indof3rdPhonon_plus(Naccum_plus+dbg_jj)-1, Nbands) + 1
              nn = int((Indof3rdPhonon_plus(Naccum_plus+dbg_jj)-1)/Nbands) + 1
              if (mm<1 .or. mm>nptk .or. nn<1 .or. nn>nptk) then
                write(*,*) "DEBUG_TERM: invalid mm/nn in plus",mm,nn," at state/jj=",dbg_state,dbg_jj
                dbg_found = .true.; exit
              end if
              idx_sym_k = ID_equi(sym_idx, nn); idx_sym_j = ID_equi(sym_idx, mm)
              dbg_term = Gamma_plus_sp(Naccum_plus+dbg_jj) * ( F_sp(c_bad,k,idx_sym_k) - F_sp(c_bad,j,idx_sym_j) )
              if (.not. ieee_is_finite(dbg_term)) then
                write(*,*) "DEBUG_TERM: non-finite PLUS term at state/kk/jj=",dbg_state,dbg_kk,dbg_jj
                write(*,*) "DEBUG_TERM: sym,mm,nn,idxj,idxk=",sym_idx,mm,nn,idx_sym_j,idx_sym_k
                write(*,*) "DEBUG_TERM: raw Ind2/Ind3=",Indof2ndPhonon_plus(Naccum_plus+dbg_jj),&
                           Indof3rdPhonon_plus(Naccum_plus+dbg_jj)
                write(*,*) "DEBUG_TERM: gamma_plus=",Gamma_plus_sp(Naccum_plus+dbg_jj),&
                           " Fk=",F_sp(c_bad,k,idx_sym_k)," Fj=",F_sp(c_bad,j,idx_sym_j)
                write(*,*) "DEBUG_TERM: finite(gamma,Fk,Fj,diff)=", &
                           ieee_is_finite(Gamma_plus_sp(Naccum_plus+dbg_jj)), &
                           ieee_is_finite(F_sp(c_bad,k,idx_sym_k)), &
                           ieee_is_finite(F_sp(c_bad,j,idx_sym_j)), &
                           ieee_is_finite(F_sp(c_bad,k,idx_sym_k) - F_sp(c_bad,j,idx_sym_j))
                dbg_found = .true.; exit
              end if
            end do
            if (dbg_found) exit
            do dbg_jj = 1, N_minus(mmm)
              j  = modulo(Indof2ndPhonon_minus(Naccum_minus+dbg_jj)-1, Nbands) + 1
              mm = int((Indof2ndPhonon_minus(Naccum_minus+dbg_jj)-1)/Nbands) + 1
              k  = modulo(Indof3rdPhonon_minus(Naccum_minus+dbg_jj)-1, Nbands) + 1
              nn = int((Indof3rdPhonon_minus(Naccum_minus+dbg_jj)-1)/Nbands) + 1
              if (mm<1 .or. mm>nptk .or. nn<1 .or. nn>nptk) then
                write(*,*) "DEBUG_TERM: invalid mm/nn in minus",mm,nn," at state/jj=",dbg_state,dbg_jj
                dbg_found = .true.; exit
              end if
              idx_sym_k = ID_equi(sym_idx, nn); idx_sym_j = ID_equi(sym_idx, mm)
              dbg_term = 0.5_sp*Gamma_minus_sp(Naccum_minus+dbg_jj) * ( F_sp(c_bad,k,idx_sym_k) + F_sp(c_bad,j,idx_sym_j) )
              if (.not. ieee_is_finite(dbg_term)) then
                write(*,*) "DEBUG_TERM: non-finite MINUS term at state/kk/jj=",dbg_state,dbg_kk,dbg_jj
                write(*,*) "DEBUG_TERM: sym,mm,nn,idxj,idxk=",sym_idx,mm,nn,idx_sym_j,idx_sym_k
                write(*,*) "DEBUG_TERM: raw Ind2/Ind3=",Indof2ndPhonon_minus(Naccum_minus+dbg_jj),&
                           Indof3rdPhonon_minus(Naccum_minus+dbg_jj)
                write(*,*) "DEBUG_TERM: gamma_minus=",Gamma_minus_sp(Naccum_minus+dbg_jj),&
                           " Fk=",F_sp(c_bad,k,idx_sym_k)," Fj=",F_sp(c_bad,j,idx_sym_j)
                write(*,*) "DEBUG_TERM: finite(gamma,Fk,Fj,sum)=", &
                           ieee_is_finite(Gamma_minus_sp(Naccum_minus+dbg_jj)), &
                           ieee_is_finite(F_sp(c_bad,k,idx_sym_k)), &
                           ieee_is_finite(F_sp(c_bad,j,idx_sym_j)), &
                           ieee_is_finite(F_sp(c_bad,k,idx_sym_k) + F_sp(c_bad,j,idx_sym_j))
                dbg_found = .true.; exit
              end if
            end do
            if (dbg_found) exit
          end do
        end if
        if (dbg_found) exit
      end do
    end if
    if (myid.eq.0) then
      if (.not. ieee_is_finite(sum(DeltaF_sp))) then
        write(*,*) "ERROR_ITERGPU: DeltaF_sp becomes non-finite before MPI_ALLREDUCE"
        write(*,*) "DEBUG_ITERGPU: max|DeltaF_sp|=", maxval(abs(DeltaF_sp))
      end if
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, DeltaF_sp, nptk*Nbands*3, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (myid.eq.0) then
      if (.not. ieee_is_finite(sum(DeltaF_sp))) then
        write(*,*) "ERROR_ITERGPU: DeltaF_sp non-finite after MPI_ALLREDUCE"
      end if
    end if
    !$acc update device(DeltaF_sp)

    !$acc parallel loop gang collapse(3)                                      &
    !$acc& present(F_sp, tau_zero_sp, velocity_sp, omega_sp, DeltaF_sp,        &
    !$acc&         ALLEquiList_safe, Nequi)
    do ll = 1, Nlist
      do kk = 1, Nequi_max
        do i = 1, Nbands
          if (kk <= Nequi(ll)) then
            q_idx = ALLEquiList_safe(kk, ll)
            do c = 1, 3
              F_sp(c,i,q_idx) = tau_zero_sp(i,ll) * ( velocity_sp(c,i,q_idx) * omega_sp(i,q_idx) + &
                                   DeltaF_sp(c,i,q_idx) )
            end do
          end if
        end do
      end do
    end do
    !$acc end parallel
    !$acc update self(F_sp)
    if (myid.eq.0) then
      if (.not. ieee_is_finite(sum(F_sp))) then
        write(*,*) "ERROR_ITERGPU: F_sp non-finite right after F update kernel"
        write(*,*) "DEBUG_ITERGPU: max|tau_zero_sp|=", maxval(abs(tau_zero_sp)), " max|DeltaF_sp|=", maxval(abs(DeltaF_sp))
      end if
    end if

    if (IsoN > 0) then
      !$acc exit data delete(Iso_Ind1, Iso_Ind2, Gamma_iso_sp)
    end if

    !$acc end data
    if (allocated(TypeofSymmetry_safe)) deallocate(TypeofSymmetry_safe)
    if (allocated(ALLEquiList_safe)) deallocate(ALLEquiList_safe)

    ! 回写双精度接口
    F_n = real(F_sp,kind=dp)

  end subroutine iterationgpu


  ! Restricted variation of the above subroutine, limited to cases where kappa
  ! is an scalar. Used in the nanowire calculation.
  subroutine iteration_scalar(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,&
      & Ntotal_plus,Ntotal_minus,&
      & omega,velocity,Gamma_plus,Gamma_minus,tau_zero,F_n)
    implicit none

    integer,intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm_rot,nptk)
    integer,intent(in) :: TypeofSymmetry(Nsymm_rot,nptk)
    integer,intent(in) :: N_plus(Nlist*Nbands),N_minus(Nlist*Nbands)
    integer(kind=long),intent(in) :: Ntotal_plus,Ntotal_minus
    real(kind=dp),intent(in) :: omega(nbands,nptk),velocity(nbands,nptk)
    real(kind=dp),intent(in) :: Gamma_plus(Ntotal_plus),Gamma_minus(Ntotal_minus),tau_zero(nbands,nlist)
    real(kind=dp),intent(inout) :: F_n(Nbands,nptk)

    integer :: ID_equi(Nsymm_Rot,nptk)
    integer (kind=long) :: Naccum_plus,Naccum_minus
    integer :: i,j,k,jj,kk,ll,mm,nn,mmm,nnn
    real(kind=dp) :: DeltaF(Nbands,nptk)

    DeltaF=0.0_dp
    call symmetry_map(ID_equi)
    !$OMP PARALLEL DO default(none) schedule(dynamic,1) shared(nstates,myid,nbands,nlist) &
    !$OMP & shared(Nequi,F_n,ID_equi,TypeofSymmetry,ALLEquiList,DeltaF) &
    !$OMP & shared(N_plus,Naccum_plus_array,Indof2ndPhonon_plus,Indof3rdPhonon_plus,Gamma_plus) &
    !$OMP & shared(N_minus,Naccum_minus_array,Indof2ndPhonon_minus,Indof3rdPhonon_minus,Gamma_minus) &
    !$OMP & private(i,j,k,kk,ll,mm,nn,mmm,nnn,Naccum_plus,Naccum_minus)
    do nnn=1,nstates
        mmm=myid*nstates+nnn
        i=modulo(mmm-1,Nbands)+1
        ll=int((mmm-1)/Nbands)+1

        if (mmm.gt.nlist*nbands) cycle
        Naccum_plus=Naccum_plus_array(nnn)
        Naccum_minus=Naccum_minus_array(nnn)

        do kk=1,Nequi(ll)
           if ((N_plus(mmm).ne.0)) then
              do jj=1,N_plus(mmm) 
                  j=modulo(Indof2ndPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                  mm=(Indof2ndPhonon_plus(Naccum_plus+jj)-1)/Nbands + 1
                  k=modulo(Indof3rdPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                  nn=(Indof3rdPhonon_plus(Naccum_plus+jj)-1)/Nbands + 1
                  DeltaF(i,ALLEquiList(kk,ll))=DeltaF(i,ALLEquiList(kk,ll))+&
                        Gamma_plus(Naccum_plus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn))-&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm)))
              end do !jj
           end if
           if ((N_minus(mmm).ne.0)) then
               do jj=1,N_minus(mmm)
                  j=modulo(Indof2ndPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                  mm=(Indof2ndPhonon_minus(Naccum_minus+jj)-1)/Nbands + 1
                  k=modulo(Indof3rdPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                  nn=(Indof3rdPhonon_minus(Naccum_minus+jj)-1)/Nbands + 1
                  DeltaF(i,ALLEquiList(kk,ll))=DeltaF(i,ALLEquiList(kk,ll))+&
                        Gamma_minus(Naccum_minus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn))+&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm)))*0.5_dp
               end do !jj
           end if
        end do !kk
    end do ! nnn
    !$OMP END PARALLEL DO

    if ((IsoInfo%nisotopic .ne. 0)) then
        !$OMP PARALLEL DO default(none) schedule(static) shared(nbands,IsoInfo) &
        !$OMP & shared(Nequi,F_n,ID_equi,TypeofSymmetry,ALLEquiList,DeltaF) &
        !$OMP & private(i,j,kk,ll,nn,jj)
        do jj=1,IsoInfo%nisotopic
            ! q and ib of input phonon
            i  = modulo(IsoInfo%Indof1stPhononIso(jj)-1,Nbands)+1
            ll = int((IsoInfo%Indof1stPhononIso(jj)-1)/Nbands)+1
            ! q and ib of outcome phonon
            j  = modulo(IsoInfo%Indof2ndPhononIso(jj)-1,Nbands)+1
            kk = int((IsoInfo%Indof2ndPhononIso(jj)-1)/Nbands)+1
            !Iterate over equivalence class of first phonon
            do nn=1,Nequi(ll)
               DeltaF(i,ALLEquiList(nn,ll))=DeltaF(i,ALLEquiList(nn,ll))+&
                  IsoInfo%gamma_isotopic(jj) * F_n(j,ID_equi(TypeofSymmetry(nn,ll),kk))
            end do !nn
        end do !jj
        !$OMP END PARALLEL DO
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,mmm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,DeltaF,nptk*nbands,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,mm)

    !$OMP PARALLEL DO default(none) schedule(static) shared(tau_zero,omega,velocity,ALLEquiList) &
    !$OMP & shared(F_n,DeltaF,Nlist,nbands,Nequi) private(ll,i,kk)
    do ll=1,Nlist
       do kk=1,Nequi(ll)
          do i=1,Nbands
             F_n(i,ALLEquiList(kk,ll))=tau_zero(i,ll)*velocity(i,ALLEquiList(kk,ll))*&
                  omega(i,ALLEquiList(kk,ll))+tau_zero(i,ll)*DeltaF(i,ALLEquiList(kk,ll))
          enddo
       enddo
    enddo   
    !$OMP END PARALLEL DO 
  end subroutine iteration_scalar
end module iterations
