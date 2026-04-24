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

! Compute the number of allowed three-phonon processes, their
! scattering amplitudes and their phase-space volume.

module processes
  use iso_fortran_env, only : output_unit
  use iso_c_binding,  only : c_backspace
  use omp_lib,        only : omp_get_thread_num
  use misc
  use data
  use config
  use mpi
  implicit none

  real(kind=dp), parameter :: hbarp = hbar*1.0e22_dp

  !==================== 跨迭代复用的全局事件列表缓存 ====================

  integer(kind=long),save :: plus_count   = 0_long
  integer(kind=long),save :: minus_count  = 0_long

  ! 事件列表（与 Vp_*_matrix 的顺序一一对应；不在 Ind_driver 中释放）
  integer, allocatable, save :: plus_i(:),  plus_kp(:),  plus_ii(:),  plus_j(:),  plus_k(:),  plus_ss(:)
  integer, allocatable, save :: minus_i(:), minus_kp(:), minus_ii(:), minus_j(:), minus_k(:), minus_ss(:)

contains

!======================================================================
! 构建事件列表 + 计算 Vp（GPU）；事件列表缓存为模块级且不释放
!======================================================================
  subroutine calculate_Vp(energy, velocity, eigenvect, Nlist, List, &
                          Ntri, Phi, phase_Rj, phase_Rk, Index_i, Index_j, Index_k, IJK, &
                          plus_size, minus_size, start_time)
    use mpi, only : MPI_BARRIER, MPI_COMM_WORLD, MPI_Wtime
    implicit none

    !---------------- 实参（接口保持双精度） ----------------
    integer,            intent(in) :: NList, List(Nlist), IJK(3,nptk), Ntri
    integer(kind=long), intent(in) :: plus_size, minus_size
    integer,            intent(in) :: Index_i(Ntri), Index_j(Ntri), Index_k(Ntri)
    real(kind=dp),      intent(in) :: energy(Nbands,nptk), velocity(3,Nbands,nptk)
    real(kind=dp),      intent(in) :: Phi(3,3,3,Ntri)
    complex(kind=dp),   intent(in) :: phase_Rj(Ntri,nptk), phase_Rk(Ntri,nptk)
    complex(kind=dp),   intent(in) :: eigenvect(Nbands,Nbands,nptk)
    real(kind=dp),      intent(in) :: start_time

    !---------------- 局部变量 ----------------
    integer :: nn, ii, j, k, ll, i, kp_ll
    integer :: q1, q2, q3, qp1, qp2, qp3
    integer :: qdp1, qdp2, qdp3, qdm1, qdm2, qdm3, ss_p, ss_m
    integer :: mm, ierr
    integer(kind=long) :: xx
    integer(kind=long) :: plus_write, minus_write
    integer(kind=long) :: pos, hcode, plus_counter, minus_counter
    integer(kind=long), allocatable :: plus_hash(:), minus_hash(:)

    ! Step1 判据（sp）
    real(kind=sp) :: omega, omegap, omegadp
    real(kind=sp) :: vx, vy, vz, sigma, two_sigma, dvx,dvy,dvz
    real(kind=sp) :: scalebroad_sp, omega_max_sp
    real(kind=sp), allocatable :: energy_sp(:,:), velocity_sp(:,:,:)

    ! Step2/3 GPU 计算（sp）
    real(kind=sp),    allocatable :: Phi_sp(:,:,:,:)
    complex(kind=sp), allocatable :: phase_Rj_sp(:,:), phase_Rk_sp(:,:)
    complex(kind=sp), allocatable :: eigenvect_sp(:,:,:)
    real(kind=sp),    allocatable :: Vp_plus_matrix_sp(:), Vp_minus_matrix_sp(:)

    !---------------- 输出矩阵（双精度） ----------------
   
    allocate(Vp_plus_matrix(plus_size))
    allocate(Vp_minus_matrix(minus_size))
    Vp_plus_matrix  = 0.0_dp
    Vp_minus_matrix = 0.0_dp
    allocate(plus_i(plus_size))
    allocate(plus_ii(plus_size))
    allocate(plus_k(plus_size))
    allocate(plus_kp(plus_size))
    allocate(plus_j(plus_size))
    allocate(plus_ss(plus_size))
    allocate(minus_i(minus_size))
    allocate(minus_ii(minus_size))
    allocate(minus_k(minus_size))
    allocate(minus_kp(minus_size))
    allocate(minus_j(minus_size))
    allocate(minus_ss(minus_size))

    !---------------- Step1：准备单精度镜像 & hash ----------------
    allocate(energy_sp(Nbands,nptk));       energy_sp   = real(energy,   kind=sp)
    allocate(velocity_sp(3,Nbands,nptk));   velocity_sp = real(velocity, kind=sp)
    scalebroad_sp = real(scalebroad, kind=sp)
    omega_max_sp  = real(omega_max,  kind=sp)

    allocate(plus_hash(plus_size), minus_hash(minus_size))
    plus_hash  = 0_long
    minus_hash = 0_long
    plus_counter  = 0_long
    minus_counter = 0_long

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (myid == 0) then
      write(*,'(A)') " Info (Step1 build list, OpenACC, sp): start."
      flush(output_unit)
    end if

    !$acc data copyin(energy_sp, velocity_sp, IJK, Ngrid, List, Nlist, nptk, Nbands, nstates, myid) &
    !$acc      copyin(omega_max_sp, scalebroad_sp, plus_size, minus_size)                             &
    !$acc      copy  (plus_i, plus_kp, plus_ii, plus_j, plus_k, plus_ss, plus_hash, plus_counter,     &
    !$acc            minus_i, minus_kp, minus_ii, minus_j, minus_k, minus_ss, minus_hash, minus_counter)

      !$acc parallel default(present)
        !$acc loop gang collapse(3) private(nn,ii,j,k,mm,i,ll,kp_ll,q1,q2,q3,qp1,qp2,qp3, &
        !$acc                                qdp1,qdp2,qdp3,qdm1,qdm2,qdm3,ss_p,ss_m,omega,omegap,omegadp, &
        !$acc                                vx,vy,vz,dvx,dvy,dvz,sigma,two_sigma,pos,hcode)
        do nn = 1, nstates   !nbands*等价点的个数（Nlist)除以mpi进程数的向上取整=nstates
          do ii = 1, nptk
            do j = 1, Nbands
            mm = myid*nstates + nn
            if (mm > Nlist*Nbands) cycle !由于向上取整了 mm会大于 Nlist*Nbands会有轮空 要去掉

            i     = modulo(mm-1, Nbands) + 1
            ll    = int((mm-1)/Nbands) + 1
            kp_ll = List(ll)

            omega = energy_sp(i, kp_ll)
            if (omega == 0.0_sp .or. omega > omega_max_sp) cycle

            q1 = IJK(1, kp_ll);  q2 = IJK(2, kp_ll);  q3 = IJK(3, kp_ll)
            qp1 = IJK(1, ii);    qp2 = IJK(2, ii);    qp3 = IJK(3, ii)

            ! q + q'
            qdp1 = modulo(q1 + qp1, Ngrid(1))
            qdp2 = modulo(q2 + qp2, Ngrid(2))
            qdp3 = modulo(q3 + qp3, Ngrid(3))
            ss_p = (qdp3*Ngrid(2) + qdp2)*Ngrid(1) + qdp1 + 1

            ! q - q'
            qdm1 = modulo(q1 - qp1, Ngrid(1))
            qdm2 = modulo(q2 - qp2, Ngrid(2))
            qdm3 = modulo(q3 - qp3, Ngrid(3))
            ss_m = (qdm3*Ngrid(2) + qdm2)*Ngrid(1) + qdm1 + 1

            omegap = energy_sp(j, ii)
            if (omegap == 0.0_sp) cycle

            vx = velocity_sp(1, j, ii)
            vy = velocity_sp(2, j, ii)
            vz = velocity_sp(3, j, ii)

              !$acc loop seq
              do k = 1, Nbands

                ! ---- plus ----
                omegadp = energy_sp(k, ss_p)
                if (omegadp /= 0.0_sp) then
                  dvx = vx - velocity_sp(1, k, ss_p)
                  dvy = vy - velocity_sp(2, k, ss_p)
                  dvz = vz - velocity_sp(3, k, ss_p)
                  sigma     = scalebroad_sp * base_sigma_sp(dvx,dvy,dvz)
                  two_sigma = 2.0_sp * sigma

                  if (abs(omega + omegap - omegadp) <= two_sigma) then
                    !$acc atomic capture
                      pos = plus_counter; plus_counter = plus_counter + 1_long
                    !$acc end atomic
                    if (pos < plus_size) then
                      hcode = (((mm-1_long)*nptk + (ii-1_long)) * int(Nbands,kind=long) + (j-1_long)) * &
                              int(Nbands,kind=long) + (k-1_long)
                      plus_hash(pos+1_long) = hcode
                      plus_i   (pos+1_long) = i
                      plus_kp  (pos+1_long) = kp_ll
                      plus_ii  (pos+1_long) = ii
                      plus_j   (pos+1_long) = j
                      plus_k   (pos+1_long) = k
                      plus_ss  (pos+1_long) = ss_p
                    end if
                  end if
                end if

                ! ---- minus ----
                omegadp = energy_sp(k, ss_m)
                if (omegadp /= 0.0_sp) then
                  dvx = vx - velocity_sp(1, k, ss_m)
                  dvy = vy - velocity_sp(2, k, ss_m)
                  dvz = vz - velocity_sp(3, k, ss_m)
                  sigma     = scalebroad_sp * base_sigma_sp(dvx,dvy,dvz)
                  two_sigma = 2.0_sp * sigma

                  if (abs(omega - omegap - omegadp) <= two_sigma) then
                    !$acc atomic capture
                      pos = minus_counter; minus_counter = minus_counter + 1_long
                    !$acc end atomic
                    if (pos < minus_size) then
                      hcode = (((mm-1_long)*int(nptk,kind=long) + (ii-1_long)) * int(Nbands,kind=long) + (j-1_long)) * &
                              int(Nbands,kind=long) + (k-1_long)
                      minus_hash(pos+1_long) = hcode
                      minus_i   (pos+1_long) = i
                      minus_kp  (pos+1_long) = kp_ll
                      minus_ii  (pos+1_long) = ii
                      minus_j   (pos+1_long) = j
                      minus_k   (pos+1_long) = k
                      minus_ss  (pos+1_long) = ss_m
                    end if
                  end if
                end if

              end do
            end do

          end do
        end do
      !$acc end parallel
    !$acc end data

    plus_write  = min(plus_counter,  plus_size)
    minus_write = min(minus_counter, minus_size)
    if (myid == 0) write(*,'(A,F0.3,A)') " Info (Step1 build list, OpenACC, sp): kernel done; start CPU  radix sort by hash at", MPI_Wtime()-start_time,"s"
    if (plus_counter  /= plus_size)  write(*,*) "Warning: plus_counter=",  plus_counter,  " != plus_size=",  plus_size
    if (minus_counter /= minus_size) write(*,*) "Warning: minus_counter=", minus_counter, " != minus_size=", minus_size

    if (plus_write  > 1_long) call sort_pack7_long_int(plus_hash,  plus_i,  plus_kp,  plus_ii,  plus_j,  plus_k,  plus_ss,  plus_write)
    if (minus_write > 1_long) call sort_pack7_long_int(minus_hash, minus_i, minus_kp, minus_ii, minus_j, minus_k, minus_ss, minus_write)

    plus_count  = plus_write
    minus_count = minus_write

    deallocate(energy_sp, velocity_sp, plus_hash, minus_hash)

    if (myid == 0) write(*,'(A,F0.3,A)') " Info: matrices start filling with GPU at ", MPI_Wtime()-start_time, " s"

    !==================== Step 2：GPU 计算 Vp_plus(xx) 和 Vp_minus(xx)  ====================
    if (plus_write > 0_long .and. minus_write > 0_long) then
      allocate(Phi_sp(3,3,3,Ntri));               Phi_sp       = real(Phi,       kind=sp)
      allocate(phase_Rj_sp(Ntri,nptk));           phase_Rj_sp  = cmplx(phase_Rj, kind=sp)
      allocate(phase_Rk_sp(Ntri,nptk));           phase_Rk_sp  = cmplx(phase_Rk, kind=sp)
      allocate(eigenvect_sp(Nbands,Nbands,nptk)); eigenvect_sp = cmplx(eigenvect,kind=sp)
      allocate(Vp_plus_matrix_sp(plus_write));    Vp_plus_matrix_sp = 0.0_sp
      allocate(Vp_minus_matrix_sp(minus_write));  Vp_minus_matrix_sp = 0.0_sp
      !$acc data copyin(Phi_sp, phase_Rj_sp, phase_Rk_sp, eigenvect_sp, Index_i, Index_j, Index_k,plus_i, plus_kp, plus_ii, plus_j, plus_k, plus_ss, minus_i, minus_kp, minus_ii, minus_j, minus_k, minus_ss) &
      !$acc      copy(Vp_plus_matrix_sp,Vp_minus_matrix_sp)
        !$acc parallel default(present)
          !$acc loop gang vector private(i,kp_ll,ii,j,k,ss_p,xx)
          do xx = 1_long, plus_write
            i     = plus_i(xx)
            kp_ll = plus_kp(xx)
            ii    = plus_ii(xx)
            j     = plus_j(xx)
            k     = plus_k(xx)
            ss_p  = plus_ss(xx)

            Vp_plus_matrix_sp(xx) = Vp_plus_sp( &
                 eigenvect_sp(:, i,  kp_ll), &
                 eigenvect_sp(:, j,  ii),    &
                 eigenvect_sp(:, k,  ss_p),  &
                 Ntri, Phi_sp,               &
                 phase_Rj_sp(:, ii),         &
                 phase_Rk_sp(:, ss_p),       &
                 Index_i, Index_j, Index_k )
          end do
        !$acc end parallel

        !$acc parallel default(present)
          !$acc loop gang vector private(i,kp_ll,ii,j,k,ss_m,xx)
          do xx = 1_long, minus_write
            i     = minus_i(xx)
            kp_ll = minus_kp(xx)
            ii    = minus_ii(xx)
            j     = minus_j(xx)
            k     = minus_k(xx)
            ss_m  = minus_ss(xx)

            Vp_minus_matrix_sp(xx) = Vp_minus_sp( &
                 eigenvect_sp(:, i,  kp_ll), &
                 eigenvect_sp(:, j,  ii),    &
                 eigenvect_sp(:, k,  ss_m),  &
                 Ntri, Phi_sp,               &
                 phase_Rj_sp(:, ii),         &
                 phase_Rk_sp(:, ss_m),       &
                 Index_i, Index_j, Index_k )
          end do
        !$acc end parallel
      !$acc end data
      Vp_plus_matrix = real(Vp_plus_matrix_sp, dp)
      Vp_minus_matrix = real(Vp_minus_matrix_sp, dp)
      deallocate(Vp_plus_matrix_sp, Vp_minus_matrix_sp, Phi_sp, phase_Rj_sp, phase_Rk_sp, eigenvect_sp)
    end if

    if (myid == 0) write(*,'(A,F0.3,A)') " Info: matrices filled with GPU at ", MPI_Wtime()-start_time, " s"

  end subroutine calculate_Vp


  !==================== CPU：按 hash 基数排序，同步重排 6 个 int 数组 ====================




!!!!并行化基数排序
!!!! int64拆成8个8bit 从最末端得8bit开始 一共排序八次
!!!! 构造256个桶（因为8bit 最多256种情况）
!!!! 将hash值根据并行核数分割成几段 每一段都有自己得桶 桶的编号是(a,b) a是余数类似，b是线程编号
!!!! 将这些值的最后八个bit分到所有桶里 然后对同一个a的桶里的元素总量计数 得到start(1)=xx start(2)=xx start1是0 start2是线程上所有1号桶内元素的总和 标记了写入位置 然后每个桶都有自己的数量 根据start得到pos pos(2,0)=6 依次写入元素 然后结束了以后下一个8bit再开始.
subroutine sort_pack7_long_int(keys, a1, a2, a3, a4, a5, a6, n)
    use omp_lib, only : omp_get_thread_num, omp_get_num_threads
    implicit none

    integer(kind=long), intent(inout), target :: keys(:)
    integer,            intent(inout), target :: a1(:), a2(:), a3(:), a4(:), a5(:), a6(:)
    integer(kind=long), intent(in)            :: n

    integer :: nt
    integer :: pass, shift, npass
    integer :: t, tid
    integer :: bb, b

    integer(kind=long), parameter :: mask8 = 255_long

    integer(kind=long), allocatable, target :: keys_w(:)
    integer,            allocatable, target :: a1_w(:), a2_w(:), a3_w(:), a4_w(:), a5_w(:), a6_w(:)

    integer(kind=long), allocatable :: hist(:, :)     ! (0:255, 0:nt-1)
    integer(kind=long), allocatable :: tstart(:, :)   ! (0:255, 0:nt-1)
    integer(kind=long), allocatable :: total(:)       ! (0:255)
    integer(kind=long), allocatable :: start(:)       ! (0:255)

    integer(kind=long) :: i, j, istart, iend
    integer(kind=long) :: run
    integer(kind=long) :: pos(0:255)

    integer(kind=long), pointer :: k_in(:), k_out(:), k_tmp(:)
    integer,            pointer :: p1(:), p2(:), p3(:), p4(:), p5(:), p6(:)
    integer,            pointer :: q1(:), q2(:), q3(:), q4(:), q5(:), q6(:)
    integer,            pointer :: tpi(:)

    if (n <= 1_long) return

    ! radix 趟数：按 long 的位宽决定（一般 64bit -> 8 趟；32bit -> 4 趟）
    npass = (bit_size(0_long) + 7) / 8

    ! 获取实际线程数（若未开启 OpenMP，nt=1）
    nt = 1
    !$omp parallel
    !$omp single
      nt = omp_get_num_threads()
    !$omp end single
    !$omp end parallel

    allocate(keys_w(n), a1_w(n), a2_w(n), a3_w(n), a4_w(n), a5_w(n), a6_w(n))
    allocate(hist(0:255, 0:nt-1), tstart(0:255, 0:nt-1))
    allocate(total(0:255), start(0:255))

    ! 指针指向“当前输入/输出缓冲”，每趟结束交换指针，避免来回拷贝
    k_in  => keys(1:n);   k_out => keys_w
    p1    => a1(1:n);     q1    => a1_w
    p2    => a2(1:n);     q2    => a2_w
    p3    => a3(1:n);     q3    => a3_w
    p4    => a4(1:n);     q4    => a4_w
    p5    => a5(1:n);     q5    => a5_w
    p6    => a6(1:n);     q6    => a6_w

    do pass = 0, npass - 1
      shift = 8 * pass
      hist  = 0_long

      ! 1) 每线程局部直方图（256 桶）
      !$omp parallel default(none) num_threads(nt) &
      !$omp shared(hist, k_in, n, shift, mask8, nt) &
      !$omp private(tid, i, istart, iend, bb)
        tid = omp_get_thread_num()
        call chunk_bounds_1based(n, nt, tid, istart, iend)

        do i = istart, iend
          bb = int(iand(ishft(k_in(i), -shift), mask8))   ! 0..255
          hist(bb, tid) = hist(bb, tid) + 1_long
        end do
      !$omp end parallel

      ! 2) 汇总每个 bucket 的总数，并计算全局起始位置 start(b)
      do b = 0, 255
        total(b) = sum(hist(b, 0:nt-1))
      end do

      start(0) = 1_long
      do b = 1, 255
        start(b) = start(b-1) + total(b-1)
      end do

      ! 3) 计算每线程每 bucket 的写入起点（按 tid 顺序拼接 -> 稳定）
      do b = 0, 255
        run = start(b)
        do t = 0, nt-1
          tstart(b, t) = run
          run = run + hist(b, t)
        end do
      end do

      ! 4) stable scatter：每线程按自己的 tstart 写入，不冲突
      !$omp parallel default(none) num_threads(nt) &
      !$omp shared(tstart, k_in, p1,p2,p3,p4,p5,p6, k_out, q1,q2,q3,q4,q5,q6, n, shift, mask8, nt) &
      !$omp private(tid, i, j, istart, iend, bb, pos)
        tid = omp_get_thread_num()
        call chunk_bounds_1based(n, nt, tid, istart, iend)

        pos = tstart(:, tid)

        do i = istart, iend
          bb = int(iand(ishft(k_in(i), -shift), mask8))
          j  = pos(bb)
          pos(bb) = j + 1_long

          k_out(j) = k_in(i)
          q1(j)    = p1(i)
          q2(j)    = p2(i)
          q3(j)    = p3(i)
          q4(j)    = p4(i)
          q5(j)    = p5(i)
          q6(j)    = p6(i)
        end do
      !$omp end parallel

      ! 5) 交换输入/输出缓冲（下一趟继续）
      k_tmp => k_in;  k_in  => k_out;  k_out => k_tmp

      tpi => p1; p1 => q1; q1 => tpi
      tpi => p2; p2 => q2; q2 => tpi
      tpi => p3; p3 => q3; q3 => tpi
      tpi => p4; p4 => q4; q4 => tpi
      tpi => p5; p5 => q5; q5 => tpi
      tpi => p6; p6 => q6; q6 => tpi
    end do

    ! 如果趟数是奇数，结果在工作数组里，需要拷回原数组
    if (mod(npass, 2) == 1) then
      keys(1:n) = k_in
      a1(1:n)   = p1
      a2(1:n)   = p2
      a3(1:n)   = p3
      a4(1:n)   = p4
      a5(1:n)   = p5
      a6(1:n)   = p6
    end if

    deallocate(keys_w, a1_w, a2_w, a3_w, a4_w, a5_w, a6_w)
    deallocate(hist, tstart, total, start)

  contains

    subroutine chunk_bounds_1based(nloc, ntloc, tidloc, istartloc, iendloc)
      implicit none
      integer(kind=long), intent(in)  :: nloc
      integer,            intent(in)  :: ntloc, tidloc
      integer(kind=long), intent(out) :: istartloc, iendloc

      integer(kind=long) :: base, rem, t64

      t64  = int(tidloc, kind=long)
      base = nloc / int(ntloc, kind=long)
      rem  = mod(nloc, int(ntloc, kind=long))

      istartloc = 1_long + t64 * base + min(t64, rem)
      iendloc   = istartloc + base - 1_long
      if (t64 < rem) iendloc = iendloc + 1_long

      if (iendloc < istartloc) then
        istartloc = 1_long
        iendloc   = 0_long
      end if
    end subroutine chunk_bounds_1based

  end subroutine sort_pack7_long_int

  
  

!======================================================================
! 复用事件列表计算 Ind（GPU 计算 WP3/Gamma，CPU 分段累加）
!======================================================================
subroutine Ind_driver( energy,velocity,Nlist,List,IJK,N_plus,N_minus, Naccum_plus_array, &
                       Naccum_minus_array, rate_scatt_plus,rate_scatt_minus,              &
                       Pspace_plus_total,Pspace_minus_total,start_time)
  use iso_fortran_env, only : output_unit
  use omp_lib
  use config
  use misc
  use data
  implicit none

  !---------------- arguments ----------------
  real(kind=dp),intent(in)    :: energy(nbands,nptk)
  real(kind=dp),intent(in)    :: velocity(3,nbands,nptk)
  real(kind=dp),intent(in)    :: start_time
  integer,      intent(in)    :: NList
  integer,      intent(in)    :: List(Nlist)
  integer,      intent(in)    :: IJK(3,nptk)
  integer,      intent(in)    :: N_plus(Nlist*Nbands), N_minus(Nlist*Nbands)
  integer(kind=long),intent(in) :: Naccum_plus_array(nstates), Naccum_minus_array(nstates)
  real(kind=dp),intent(inout) :: rate_scatt_plus(nbands,nlist), rate_scatt_minus(nbands,nlist)
  real(kind=dp),intent(inout) :: Pspace_plus_total(nbands,nlist), Pspace_minus_total(nbands,nlist)

  !---------------- locals ----------------
  real(kind=sp), allocatable :: energy_sp(:,:), velocity_sp(:,:,:)
  integer :: nn, mm, i, ll
  integer :: ii, j, k, ss
  integer :: kp, ib, kp_ev
  integer(kind=long) :: base_plus, base_minus
  integer(kind=long) :: nloc_plus, nloc_minus
  integer(kind=long) :: e, idx
  integer(kind=long) :: totalp, totalm

  real(kind=sp) :: omega_sp, omegap_sp, omegadp_sp
  real(kind=sp) :: fBEprime_sp, fBEdprime_sp
  real(kind=sp) :: vel_diff_sp(3)
  real(kind=sp) :: sigma_sp, delta_sp
  real(kind=sp) :: WP3_sp
  real(kind=sp) :: Vp_sp, gamma_plus_sp, gamma_minus_sp

  real(kind=sp) :: nptk_sp, scalebroad_sp, hbar_sp, Kb_sp, T_sp, hbarp_sp, pi_sp, sqrt_pi_sp, omega_max_sp
  real(kind=sp), parameter :: zero_sp = 0.0_sp, one_sp = 1.0_sp, two_sp = 2.0_sp
  real(kind=sp), parameter :: half_sp = 0.5_sp, four_sp = 4.0_sp
  real(kind=sp), parameter :: five_point_six_zero_sp = 5.60626442_sp
  real(kind=sp), parameter :: one_eight_sp = 1.0e8_sp

  ! 逐事件缓存（全局）
  real(kind=sp), allocatable :: WP3_plus_sp_total(:), WP3_minus_sp_total(:)

  ! 每个本地态的段信息与分段累加结果
  integer :: countm, g
  integer, allocatable :: band_of_state(:), list_of_state(:)     ! 映射到 (i,ll)
  integer(kind=long), allocatable :: base_plus_list(:), base_minus_list(:)
  integer(kind=long), allocatable :: nloc_plus_total(:), nloc_minus_total(:)
  real(kind=sp),    allocatable :: WP3_plus_sp_partition(:), WP3_minus_sp_partition(:)
  real(kind=sp),    allocatable :: gamma_plus_sp_partition(:), gamma_minus_sp_partition(:)

  ! 分段累加的局部标量（在并行区内设为 private）
  real(kind=sp) :: wsum_p, gsum_p, wsum_m, gsum_m

  !---------------- prepare SP mirrors (serial) ----------------
  allocate(energy_sp(nbands,nptk))
  allocate(velocity_sp(3,nbands,nptk))
  energy_sp   = real(energy,   kind=sp)
  velocity_sp = real(velocity, kind=sp)

  nptk_sp       = real(nptk,       kind=sp)
  scalebroad_sp = real(scalebroad, kind=sp)
  hbar_sp       = real(hbar,       kind=sp)
  Kb_sp         = real(Kb,         kind=sp)
  T_sp          = real(T,          kind=sp)
  hbarp_sp      = real(hbarp,      kind=sp)
  pi_sp         = real(pi,         kind=sp)
  sqrt_pi_sp    = sqrt(pi_sp)
  omega_max_sp  = real(omega_max,  kind=sp)

  !---------------- sanity check: totals must match cached list sizes ----------------
  totalp = 0_long
  totalm = 0_long
  do nn = 1, nstates
    mm = myid*nstates + nn
    if (mm > Nlist*Nbands) cycle
    totalp = totalp + int(N_plus(mm),  kind=long)
    totalm = totalm + int(N_minus(mm), kind=long)
  end do

  if (totalp /= plus_count) then
    write(output_unit,*) "ERROR: sum(N_plus(local))=", totalp, " but plus_count=", plus_count
    error stop "Ind_driver mapping mismatch: plus"
  end if
  if (totalm /= minus_count) then
    write(output_unit,*) "ERROR: sum(N_minus(local))=", totalm, " but minus_count=", minus_count
    error stop "Ind_driver mapping mismatch: minus"
  end if

  !---------------- phase 1: compute per-event WP3 and Gamma (global) ----------------
  allocate(WP3_plus_sp_total(plus_count), WP3_minus_sp_total(minus_count))
  WP3_plus_sp_total  = 0.0_sp
  WP3_minus_sp_total = 0.0_sp

  if (plus_count  > 0_long) Gamma_plus (:) = 0.0_dp
  if (minus_count > 0_long) Gamma_minus(:) = 0.0_dp

  !================ OpenACC: 仅并行化两个大事件循环 ================
  !$acc data                                                                                                        &
  !$acc& copyin(energy_sp(1:nbands,1:nptk), velocity_sp(1:3,1:nbands,1:nptk))                                       &
  !$acc& copyin(plus_i(1:plus_count), plus_kp(1:plus_count), plus_ii(1:plus_count), plus_j(1:plus_count),           &
  !$acc&        plus_k(1:plus_count), plus_ss(1:plus_count), Vp_plus_matrix(1:plus_count))                          &
  !$acc& copyin(minus_i(1:minus_count), minus_kp(1:minus_count), minus_ii(1:minus_count), minus_j(1:minus_count),   &
  !$acc&        minus_k(1:minus_count), minus_ss(1:minus_count), Vp_minus_matrix(1:minus_count))                    &
  !$acc& copyin(nptk_sp, scalebroad_sp, hbar_sp, Kb_sp, T_sp, hbarp_sp, pi_sp, sqrt_pi_sp, omega_max_sp)            &
  !$acc& copyout(WP3_plus_sp_total(1:plus_count), WP3_minus_sp_total(1:minus_count),                                &
  !$acc&         Gamma_plus(1:plus_count), Gamma_minus(1:minus_count),                                              &
  !$acc&         Indof2ndPhonon_plus(1:plus_count), Indof3rdPhonon_plus(1:plus_count),                              &
  !$acc&         Indof2ndPhonon_minus(1:minus_count), Indof3rdPhonon_minus(1:minus_count))

    ! plus events
    !$acc parallel loop gang vector present(energy_sp,velocity_sp,plus_i,plus_kp,plus_ii,plus_j,plus_k,plus_ss,      &
    !$acc&                                  Vp_plus_matrix,WP3_plus_sp_total,Gamma_plus,                              &
    !$acc&                                  Indof2ndPhonon_plus,Indof3rdPhonon_plus)                                  &
    !$acc& private(ib,kp_ev,ii,j,k,ss,omega_sp,omegap_sp,omegadp_sp,fBEprime_sp,fBEdprime_sp,                         &
    !$acc&         vel_diff_sp,sigma_sp,delta_sp,WP3_sp,Vp_sp,gamma_plus_sp)
    do idx = 1_long, plus_count
      ib    = plus_i(idx)
      kp_ev = plus_kp(idx)
      ii    = plus_ii(idx)
      j     = plus_j(idx)
      k     = plus_k(idx)
      ss    = plus_ss(idx)

      omega_sp   = energy_sp(ib, kp_ev)
      omegap_sp  = energy_sp(j,  ii)
      omegadp_sp = energy_sp(k,  ss)

      fBEprime_sp  = one_sp/(exp(hbar_sp*omegap_sp /(Kb_sp*T_sp)) - one_sp)
      fBEdprime_sp = one_sp/(exp(hbar_sp*omegadp_sp/(Kb_sp*T_sp)) - one_sp)

      vel_diff_sp(:) = velocity_sp(:, j, ii) - velocity_sp(:, k, ss)
      sigma_sp = scalebroad_sp * base_sigma_sp(vel_diff_sp(1), vel_diff_sp(2), vel_diff_sp(3))

      delta_sp = omega_sp + omegap_sp - omegadp_sp
      Vp_sp = real(Vp_plus_matrix(idx), kind=sp)

      if (sigma_sp <= 0.0_sp) then
         ! Keep indexing intact, but kill the singular event contribution
         ! Be careful!! FP32 would cause sigma_sp=0 more often in extreme situations than FP64.... 
         ! GUGUGAGA!!
         WP3_sp        = 0.0_sp
         gamma_plus_sp = 0.0_sp
      else
         WP3_sp = (fBEprime_sp - fBEdprime_sp) * &
                  exp(-(delta_sp*delta_sp)/(sigma_sp*sigma_sp)) / &
                  (sigma_sp*sqrt_pi_sp*omega_sp*omegap_sp*omegadp_sp)

         gamma_plus_sp = hbarp_sp*pi_sp/four_sp * WP3_sp * Vp_sp*Vp_sp
         gamma_plus_sp = gamma_plus_sp * five_point_six_zero_sp * one_eight_sp / nptk_sp
      end if
      
      WP3_plus_sp_total(idx) = WP3_sp

      Indof2ndPhonon_plus(idx) = (ii-1)*Nbands + j
      Indof3rdPhonon_plus(idx) = (ss-1)*Nbands + k
      Gamma_plus(idx)          = real(gamma_plus_sp, kind=dp)
    end do

    ! minus events
    !$acc parallel loop gang vector present(energy_sp,velocity_sp,minus_i,minus_kp,minus_ii,minus_j,minus_k,minus_ss, &
    !$acc&                                  Vp_minus_matrix,WP3_minus_sp_total,Gamma_minus,                           &
    !$acc&                                  Indof2ndPhonon_minus,Indof3rdPhonon_minus)                                 &
    !$acc& private(ib,kp_ev,ii,j,k,ss,omega_sp,omegap_sp,omegadp_sp,fBEprime_sp,fBEdprime_sp,                          &
    !$acc&         vel_diff_sp,sigma_sp,delta_sp,WP3_sp,Vp_sp,gamma_minus_sp)
    do idx = 1_long, minus_count
      ib    = minus_i(idx)
      kp_ev = minus_kp(idx)
      ii    = minus_ii(idx)
      j     = minus_j(idx)
      k     = minus_k(idx)
      ss    = minus_ss(idx)

      omega_sp   = energy_sp(ib, kp_ev)
      omegap_sp  = energy_sp(j,  ii)
      omegadp_sp = energy_sp(k,  ss)

      fBEprime_sp  = one_sp/(exp(hbar_sp*omegap_sp /(Kb_sp*T_sp)) - one_sp)
      fBEdprime_sp = one_sp/(exp(hbar_sp*omegadp_sp/(Kb_sp*T_sp)) - one_sp)

      vel_diff_sp(:) = velocity_sp(:, j, ii) - velocity_sp(:, k, ss)
      sigma_sp = scalebroad_sp * base_sigma_sp(vel_diff_sp(1), vel_diff_sp(2), vel_diff_sp(3))

      delta_sp = omega_sp - omegap_sp - omegadp_sp
      Vp_sp = real(Vp_minus_matrix(idx), kind=sp)

      if (sigma_sp <= 0.0_sp) then
         WP3_sp         = 0.0_sp
         gamma_minus_sp = 0.0_sp
      else
         WP3_sp = (fBEprime_sp + fBEdprime_sp + one_sp) * &
                  exp(-(delta_sp*delta_sp)/(sigma_sp*sigma_sp)) / &
                  (sigma_sp*sqrt_pi_sp*omega_sp*omegap_sp*omegadp_sp)

         gamma_minus_sp = hbarp_sp*pi_sp/four_sp * WP3_sp * Vp_sp*Vp_sp
         gamma_minus_sp = gamma_minus_sp * five_point_six_zero_sp * one_eight_sp / nptk_sp
      end if

     
      WP3_minus_sp_total(idx) = WP3_sp
      


      Indof2ndPhonon_minus(idx) = (ii-1)*Nbands + j
      Indof3rdPhonon_minus(idx) = (ss-1)*Nbands + k
      Gamma_minus(idx)          = real(gamma_minus_sp, kind=dp)
    end do
  !$acc end data
  !================ OpenACC: 结束并行区 ===============================
  if (myid == 0) write(*,'(A,F0.3,A)') " Info: GPU done,starting scattering results on CPU at ", MPI_Wtime()-start_time, " s"
 
  !---------------- phase 2: scatter (segment sum by counts/offsets) ----------------
  countm = 0
  do nn = 1, nstates
    mm = myid*nstates + nn
    if (mm > Nlist*Nbands) cycle
    countm = countm + 1
  end do

  allocate(band_of_state(countm), list_of_state(countm))
  allocate(base_plus_list(countm), base_minus_list(countm))
  allocate(nloc_plus_total(countm), nloc_minus_total(countm))

  g = 0
  do nn = 1, nstates
    mm = myid*nstates + nn
    if (mm > Nlist*Nbands) cycle
    g  = g + 1

    i  = modulo(mm-1, Nbands) + 1
    ll = int((mm-1)/Nbands) + 1

    band_of_state(g) = i
    list_of_state(g) = ll

    base_plus_list(g)  = Naccum_plus_array(nn)
    base_minus_list(g) = Naccum_minus_array(nn)

    nloc_plus_total(g)  = int(N_plus(mm),  kind=long)
    nloc_minus_total(g) = int(N_minus(mm), kind=long)

    if (nloc_plus_total(g)  > 0_long) then
      if (base_plus_list(g)  + nloc_plus_total(g)  > plus_count)  error stop "plus segment out of range"
    end if
    if (nloc_minus_total(g) > 0_long) then
      if (base_minus_list(g) + nloc_minus_total(g) > minus_count) error stop "minus segment out of range"
    end if
  end do

  allocate(WP3_plus_sp_partition(countm),  WP3_minus_sp_partition(countm))
  allocate(gamma_plus_sp_partition(countm),gamma_minus_sp_partition(countm))
  WP3_plus_sp_partition   = 0.0_sp
  WP3_minus_sp_partition  = 0.0_sp
  gamma_plus_sp_partition = 0.0_sp
  gamma_minus_sp_partition= 0.0_sp

  !================ OpenMP: g 维度并行的分段累加（plus，64 位 idx） ================
  !$omp parallel do default(shared) schedule(static) private(g, base_plus, nloc_plus, idx, wsum_p, gsum_p)
  do g = 1, countm
    base_plus = base_plus_list(g)
    nloc_plus = nloc_plus_total(g)

    wsum_p = 0.0_sp
    gsum_p = 0.0_sp

    if (nloc_plus > 0_long) then
      do idx = base_plus + 1_long, base_plus + nloc_plus
        wsum_p = wsum_p + WP3_plus_sp_total(idx)
        gsum_p = gsum_p + real(Gamma_plus(idx), kind=sp)
      end do
    end if

    WP3_plus_sp_partition(g)   = wsum_p
    gamma_plus_sp_partition(g) = gsum_p
  end do
  !$omp end parallel do

  !================ OpenMP: g 维度并行的分段累加（minus，64 位 idx） ================
  !$omp parallel do default(shared) schedule(static) private(g, base_minus, nloc_minus, idx, wsum_m, gsum_m)
  do g = 1, countm
    base_minus = base_minus_list(g)
    nloc_minus = nloc_minus_total(g)

    wsum_m = 0.0_sp
    gsum_m = 0.0_sp

    if (nloc_minus > 0_long) then
      do idx = base_minus + 1_long, base_minus + nloc_minus
        wsum_m = wsum_m + WP3_minus_sp_total(idx)
        gsum_m = gsum_m + real(Gamma_minus(idx), kind=sp) * half_sp
      end do
    end if

    WP3_minus_sp_partition(g)   = wsum_m
    gamma_minus_sp_partition(g) = gsum_m
  end do
  !$omp end parallel do

  !================ OpenMP: g 维度并行写回 (i,ll) =======================
  !$omp parallel do default(shared) schedule(static) private(g,i,ll)
  do g = 1, countm
    i  = band_of_state(g)
    ll = list_of_state(g)
    Pspace_plus_total(i,ll)  = real(WP3_plus_sp_partition(g)              / nptk_sp, kind=dp)
    Pspace_minus_total(i,ll) = real((WP3_minus_sp_partition(g) * half_sp) / nptk_sp, kind=dp)
    rate_scatt_plus(i,ll)    = real(gamma_plus_sp_partition(g),  kind=dp)
    rate_scatt_minus(i,ll)   = real(gamma_minus_sp_partition(g), kind=dp)
  end do
  !$omp end parallel do

  !---------------- cleanup ----------------
  deallocate(energy_sp, velocity_sp)
  deallocate(WP3_plus_sp_total, WP3_minus_sp_total)
  deallocate(band_of_state, list_of_state, base_plus_list, base_minus_list)
  deallocate(nloc_plus_total, nloc_minus_total)
  deallocate(WP3_plus_sp_partition, WP3_minus_sp_partition)
  deallocate(gamma_plus_sp_partition, gamma_minus_sp_partition)
end subroutine Ind_driver





subroutine NP_driver(energy,velocity,Nlist,List,IJK, &
     N_plus,Pspace_plus_total,N_minus,Pspace_minus_total)
  use config
  use misc            ! 确保 base_sigma_sp 有 !$acc routine seq
  implicit none
  !------------------------ 形参（接口双精度） ------------------------
  real(kind=dp),intent(in) :: energy(nbands,nptk)
  real(kind=dp),intent(in) :: velocity(3,nbands,nptk)
  integer,     intent(in)  :: NList
  integer,     intent(in)  :: List(Nlist)
  integer,     intent(in)  :: IJK(3,nptk)
  integer,     intent(out) :: N_plus(Nlist*Nbands)
  integer,     intent(out) :: N_minus(Nlist*Nbands)
  real(kind=dp),intent(out):: Pspace_plus_total(Nbands,Nlist)
  real(kind=dp),intent(out):: Pspace_minus_total(Nbands,Nlist)

  !------------------------ 局部变量 ------------------------
  integer :: mm                        ! 也被用作 MPI ierr
  integer :: i, ll
  integer :: ii, j, k, ss_p,ss_m
  integer :: q(3), qprime(3), qdprime(3)
  integer :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
  integer :: iii, jjj, kkk

  ! 单精度工作变量
  real(kind=sp) :: omega_sp, omegap_sp, omegadp_sp
  real(kind=sp) :: sigma_sp, delta_sp, contrib_sp
  real(kind=sp) :: dvx,dvy,dvz
  real(kind=sp) :: inv_norm_sp, omega_max_sp, scalebroad_sp, Pi_sp

  ! 单精度镜像与临时累加器
  real(kind=sp), allocatable :: energy_sp(:,:), velocity_sp(:,:,:)
  real(kind=sp), allocatable :: Pspace_plus_total_sp(:,:), Pspace_minus_total_sp(:,:)

  !------------------------ 主机侧初始化 ------------------------
  Pspace_plus_total  = 0.0_dp
  Pspace_minus_total = 0.0_dp
  N_plus             = 0
  N_minus            = 0

  ! 0 基 -> 1 基的三维到线性索引映射
  do iii=0, Ngrid(1)-1
     do jjj=0, Ngrid(2)-1
        do kkk=0, Ngrid(3)-1
           Index_N(iii,jjj,kkk) = (kkk*Ngrid(2) + jjj)*Ngrid(1) + iii + 1
        end do
     end do
  end do

  ! 单精度镜像
  allocate(energy_sp(nbands,nptk));       energy_sp   = real(energy,   kind=sp)
  allocate(velocity_sp(3,nbands,nptk));   velocity_sp = real(velocity, kind=sp)
  omega_max_sp  = real(omega_max,  kind=sp)
  scalebroad_sp = real(scalebroad, kind=sp)
  Pi_sp         = real(Pi,         kind=sp)
  inv_norm_sp   = 1.0_sp / ( sqrt(Pi_sp) * real(nptk,sp)**2 * real(nbands,sp)**3 )

  allocate(Pspace_plus_total_sp(Nbands,Nlist));  Pspace_plus_total_sp  = 0.0_sp
  allocate(Pspace_minus_total_sp(Nbands,Nlist)); Pspace_minus_total_sp = 0.0_sp

  !------------------------ OpenACC 数据区（单精度） ------------------------
  !$acc data copyin(energy_sp,velocity_sp,List,IJK,Index_N)                                &
  !$acc&     copyin(Ngrid,nptk,nbands,nlist,omega_max_sp,scalebroad_sp,myid,numprocs)      &
  !$acc&     copyin(inv_norm_sp)                                                           &
  !$acc&     copy(N_plus,N_minus,Pspace_plus_total_sp,Pspace_minus_total_sp)

  !$acc parallel default(present)
    !$acc loop gang collapse(3) private(i,ll,q,omega_sp,ii,j,k,qprime,qdprime,ss_p,ss_m,omegap_sp,omegadp_sp,dvx,dvy,dvz,sigma_sp,delta_sp,contrib_sp)
    do mm = myid+1, Nbands*Nlist, numprocs
      do ii = 1, nptk
        do j = 1, Nbands
          if (energy_sp(modulo(mm-1,Nbands)+1, List(int((mm-1)/Nbands)+1)) .le. omega_max_sp) then

          i  = modulo(mm-1, Nbands) + 1
          ll = int((mm-1)/Nbands) + 1
          q  = IJK(:, List(ll))
          omega_sp = energy_sp(i, List(ll))

          if (omega_sp .ne. 0.0_sp) then
            qprime  = IJK(:, ii)
            qdprime = q + qprime
            qdprime = modulo(qdprime, Ngrid)
            ss_p = Index_N(qdprime(1), qdprime(2), qdprime(3))  
            qdprime = q - qprime
            qdprime = modulo(qdprime, Ngrid)
            ss_m = Index_N(qdprime(1), qdprime(2), qdprime(3))
            omegap_sp   = energy_sp(j, ii)
              !$acc loop seq
              do k = 1, Nbands
                omegadp_sp = energy_sp(k, ss_p)
                if ((omegap_sp .ne. 0.0_sp) .and. (omegadp_sp .ne. 0.0_sp)) then
                  dvx = velocity_sp(1, j, ii) - velocity_sp(1, k, ss_p)
                  dvy = velocity_sp(2, j, ii) - velocity_sp(2, k, ss_p)
                  dvz = velocity_sp(3, j, ii) - velocity_sp(3, k, ss_p)
                  sigma_sp = scalebroad_sp * base_sigma_sp(dvx,dvy,dvz)
                  delta_sp = omega_sp + omegap_sp - omegadp_sp
                  if (abs(delta_sp) .le. (2.0_sp*sigma_sp)) then
                    !$acc atomic update
                    N_plus(mm) = N_plus(mm) + 1
                    contrib_sp = exp( - (delta_sp*delta_sp) / (sigma_sp*sigma_sp) ) * (inv_norm_sp / sigma_sp)
                    !$acc atomic update
                    Pspace_plus_total_sp(i, ll) = Pspace_plus_total_sp(i, ll) + contrib_sp
                  end if
                end if 
                omegadp_sp = energy_sp(k, ss_m)
                if ((omegap_sp .ne. 0.0_sp) .and. (omegadp_sp .ne. 0.0_sp)) then
                  dvx = velocity_sp(1, j, ii) - velocity_sp(1, k, ss_m)
                  dvy = velocity_sp(2, j, ii) - velocity_sp(2, k, ss_m)
                  dvz = velocity_sp(3, j, ii) - velocity_sp(3, k, ss_m)
                  sigma_sp = scalebroad_sp * base_sigma_sp(dvx,dvy,dvz)
                  delta_sp = omega_sp - omegap_sp - omegadp_sp
                  if (abs(delta_sp) .le. (2.0_sp*sigma_sp)) then
                      !$acc atomic update
                      N_minus(mm) = N_minus(mm) + 1
                      contrib_sp = exp( - (delta_sp*delta_sp) / (sigma_sp*sigma_sp) ) * (inv_norm_sp / sigma_sp)
                      !$acc atomic update
                      Pspace_minus_total_sp(i, ll) = Pspace_minus_total_sp(i, ll) + contrib_sp
                  end if
                end if 
              end do
            end if  ! omega_sp /= 0
          end if    ! energy_sp(...) <= omega_max_sp
        end do
      end do
    end do       ! mm
  !$acc end parallel

  !$acc end data

  !------------------------ 回写到双精度并 MPI 归并 ------------------------
  Pspace_plus_total  = real(Pspace_plus_total_sp,  kind=dp)
  Pspace_minus_total = real(Pspace_minus_total_sp, kind=dp)

  call MPI_ALLREDUCE(MPI_IN_PLACE,N_plus,            Nbands*Nlist,MPI_INTEGER,          MPI_SUM,MPI_COMM_WORLD,mm)
  call MPI_ALLREDUCE(MPI_IN_PLACE,N_minus,           Nbands*Nlist,MPI_INTEGER,          MPI_SUM,MPI_COMM_WORLD,mm)
  call MPI_ALLREDUCE(MPI_IN_PLACE,Pspace_plus_total, Nbands*Nlist,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,mm)
  call MPI_ALLREDUCE(MPI_IN_PLACE,Pspace_minus_total,Nbands*Nlist,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,mm)

  if (allocated(energy_sp))              deallocate(energy_sp)
  if (allocated(velocity_sp))            deallocate(velocity_sp)
  if (allocated(Pspace_plus_total_sp))   deallocate(Pspace_plus_total_sp)
  if (allocated(Pspace_minus_total_sp))  deallocate(Pspace_minus_total_sp)
end subroutine NP_driver

  ! RTA-only version of Ind_plus.
  subroutine RTA_plus(nn,q,omega,energy,velocity,IJK,&
       Gamma_plus,WP3_plus)
    implicit none

    integer,intent(in) :: nn,IJK(3,nptk), q(3)
    real(kind=dp),intent(in) :: energy(Nbands,nptk),velocity(3,Nbands,nptk), omega
    real(kind=dp),intent(out) :: Gamma_plus,WP3_plus

    integer :: qprime(3),qdprime(3),i,j,k
    integer :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer :: ii,jj,kk,ll,ss,mm
    integer(kind=long) :: N_plus_count
    real(kind=dp) :: sigma, vel_j(3)
    real(kind=dp) :: fBEprime,fBEdprime
    real(kind=dp) :: omegap,omegadp
    real(kind=dp) :: realqprime(3),realqdprime(3)
    real(kind=dp) :: WP3
    real(kind=dp) :: Vp

    Gamma_plus=0.0_dp
    WP3_plus=0.0_dp
    N_plus_count = Naccum_plus_array(nn)
    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    mm=myid*nstates+nn
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    if(omega .ne. 0.0_dp) then
       do ii=1,nptk
          qprime=IJK(:,ii)
          realqprime=matmul(rlattvec,real(qprime,kind=dp)/real(ngrid,kind=dp))
          qdprime=q+qprime
          qdprime=modulo(qdprime,Ngrid)
          realqdprime=matmul(rlattvec,real(qdprime,kind=dp)/real(ngrid,kind=dp))
          ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
          do j=1,Nbands
             omegap=energy(j,ii)
             vel_j(:) = velocity(:,j,ii)
             if (omegap .eq. 0.0_dp) cycle
             fBEprime=1.0_dp/(exp(hbar*omegap/Kb/T)-1.0_dp)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                omegadp=energy(k,ss)
                if (omegadp .eq. 0.0_dp) cycle
                sigma=scalebroad*base_sigma(vel_j(:)-velocity(:,k,ss))
                if(abs(omega+omegap-omegadp).le.(2.0_dp*sigma)) then
                   N_plus_count = N_plus_count + 1
                   fBEdprime=1.0_dp/(exp(hbar*omegadp/Kb/T)-1.0_dp)
                   WP3=(fBEprime-fBEdprime)*&
                      exp(-(omega+omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                      (omega*omegap*omegadp)
                   WP3_plus=WP3_plus+WP3
                   if (.not.onlyharmonic) then
                   Vp=Vp_plus_matrix(N_plus_count)
                   Gamma_plus=Gamma_plus+hbarp*pi/4.0_dp*WP3*Vp**2
                   endif
                end if
             end do ! k
             !--------END absorption process-------------!
          end do ! ii
       end do  ! j
       WP3_plus=WP3_plus/real(nptk,kind=dp)
    end if
    Gamma_plus=Gamma_plus*5.60626442_dp*1.0e8_dp/real(nptk,kind=dp) ! THz
  end subroutine RTA_plus

  ! RTA-only version of Ind_minus.
  subroutine RTA_minus(nn,q,omega,energy,velocity,&
       IJK,Gamma_minus,WP3_minus)
    implicit none

    integer,intent(in) :: nn,IJK(3,nptk),q(3)
    real(kind=dp),intent(in) :: energy(Nbands,nptk),velocity(3,Nbands,nptk), omega
    real(kind=dp),intent(out) :: Gamma_minus,WP3_minus

    integer :: qprime(3),qdprime(3),i,j,k
    integer(kind=long) :: N_minus_count
    integer :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer :: ii,jj,kk,ll,ss,mm
    real(kind=dp) :: sigma, vel_j(3)
    real(kind=dp) :: fBEprime,fBEdprime
    real(kind=dp) :: omegap,omegadp
    real(kind=dp) :: realqprime(3),realqdprime(3)
    real(kind=dp) :: WP3
    real(kind=dp) :: Vp

    Gamma_minus=0.0_dp
    WP3_minus=0.0_dp
    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_minus_count=Naccum_minus_array(nn)
    mm=myid*nstates+nn
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    if(omega .ne. 0.0_dp) then
       do ii=1,nptk
          qprime=IJK(:,ii)
          realqprime=matmul(rlattvec,real(qprime,kind=dp)/real(ngrid,kind=dp))
          qdprime=q-qprime
          qdprime=modulo(qdprime,Ngrid)
          realqdprime=matmul(rlattvec,real(qdprime,kind=dp)/real(ngrid,kind=dp))
          ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
          do j=1,Nbands
             omegap=energy(j,ii)
             vel_j(:) = velocity(:,j,ii)
             if (omegap .eq. 0.0_dp) cycle
             fBEprime=1.0_dp/(exp(hbar*omegap/Kb/T)-1.0_dp)
             !--------BEGIN emission process-----------
             do k=1,Nbands
                omegadp=energy(k,ss)
                if (omegadp .eq. 0.0_dp) cycle
                sigma=scalebroad*base_sigma(vel_j(:)-velocity(:,k,ss))
                if (abs(omega-omegap-omegadp).le.(2.0_dp*sigma)) then
                   N_minus_count = N_minus_count + 1
                   fBEdprime=1.0_dp/(exp(hbar*omegadp/Kb/T)-1.0_dp)
                   WP3=(fBEprime+fBEdprime+1)*&
                      exp(-(omega-omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                      (omega*omegap*omegadp)
                   WP3_minus=WP3_minus+WP3
                   if (.not.onlyharmonic) then
                   Vp=Vp_minus_matrix(N_minus_count)
                   Gamma_minus=Gamma_minus+hbarp*pi/4.0_dp*WP3*Vp**2
                   endif
                end if
             end do ! k
             !--------END emission process-------------
          end do ! ii
       end do  ! j
       WP3_minus=WP3_minus*0.5_dp/real(nptk,kind=dp)
    end if
    Gamma_minus=Gamma_minus*5.60626442_dp*1.0e8_dp/real(nptk,kind=dp)
  end subroutine RTA_minus

  ! Wrapper around RTA_plus and RTA_minus that splits the work among processors.
  subroutine RTA_driver(energy,velocity,Nlist,List,IJK,&
       rate_scatt,rate_scatt_plus,rate_scatt_minus,WP3_plus,WP3_minus)
    implicit none

    real(kind=dp),intent(in) :: energy(nbands,nptk)
    real(kind=dp),intent(in) :: velocity(3,nbands,nptk)
    integer,intent(in) :: NList
    integer,intent(in) :: List(Nlist)
    integer,intent(in) :: IJK(3,nptk)
    real(kind=dp),intent(out) :: rate_scatt(Nbands,Nlist),rate_scatt_plus(Nbands,Nlist),rate_scatt_minus(Nbands,Nlist)
    real(kind=dp),intent(out) :: WP3_plus(Nbands,Nlist)
    real(kind=dp),intent(out) :: WP3_minus(Nbands,Nlist)

    integer :: i
    integer :: ll, mm, nn
    integer :: ii, j, k, ss
    integer :: q(3), qprime(3), qdprime(3)
    integer :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer :: iii, jjj, kkk
    integer(kind=long) :: plus_offset, minus_offset

    real(kind=sp), allocatable :: energy_sp(:,:), velocity_sp(:,:,:)
    real(kind=sp), allocatable :: WP3_plus_sp(:,:), WP3_minus_sp(:,:)
    real(kind=sp), allocatable :: rate_scatt_plus_sp(:,:), rate_scatt_minus_sp(:,:)
    real(kind=sp) :: omega_sp, omegap_sp, omegadp_sp
    real(kind=sp) :: sigma_sp, delta_sp
    real(kind=sp) :: fBEprime_sp, fBEdprime_sp
    real(kind=sp) :: vel_diff_sp(3)
    real(kind=sp) :: WP3_local_sp, Gamma_plus_sp, Gamma_minus_sp, Vp_sp
    real(kind=sp) :: hbar_sp, Kb_sp, T_sp, hbarp_sp, pi_sp, sqrt_pi_sp, nptk_sp
    real(kind=sp), parameter :: one_sp = 1.0_sp, two_sp = 2.0_sp
    real(kind=sp), parameter :: half_sp = 0.5_sp, four_sp = 4.0_sp
    real(kind=sp), parameter :: five_point_six_zero_sp = 5.60626442_sp
    real(kind=sp), parameter :: one_eight_sp = 1.0e8_sp

    rate_scatt=0.0_dp
    WP3_plus=0.0_dp
    WP3_minus=0.0_dp
    rate_scatt_plus=0.0_dp
    rate_scatt_minus=0.0_dp

    allocate(energy_sp(nbands,nptk))
    allocate(velocity_sp(3,nbands,nptk))
    allocate(WP3_plus_sp(nbands,nlist), WP3_minus_sp(nbands,nlist))
    allocate(rate_scatt_plus_sp(nbands,nlist), rate_scatt_minus_sp(nbands,nlist))

    energy_sp   = real(energy, kind=sp)
    velocity_sp = real(velocity, kind=sp)
    WP3_plus_sp = 0.0_sp
    WP3_minus_sp = 0.0_sp
    rate_scatt_plus_sp = 0.0_sp
    rate_scatt_minus_sp = 0.0_sp

    hbar_sp    = real(hbar,  kind=sp)
    Kb_sp      = real(Kb,    kind=sp)
    T_sp       = real(T,     kind=sp)
    hbarp_sp   = real(hbarp, kind=sp)
    pi_sp      = real(pi,    kind=sp)
    sqrt_pi_sp = sqrt(pi_sp)
    nptk_sp    = real(nptk,  kind=sp)

    do iii=0, Ngrid(1)-1
       do jjj=0, Ngrid(2)-1
          do kkk=0, Ngrid(3)-1
             Index_N(iii,jjj,kkk) = (kkk*Ngrid(2)+jjj)*Ngrid(1)+iii+1
          end do
       end do
    end do

    !$OMP PARALLEL DO default(shared) schedule(dynamic,1) shared(Nbands,NList,myid,nstates) &
    !$OMP & shared(energy_sp,velocity_sp,List,IJK,Naccum_plus_array,Naccum_minus_array,omega_max) &
    !$OMP & shared(WP3_plus_sp,WP3_minus_sp,rate_scatt_plus_sp,rate_scatt_minus_sp,Index_N) &
    !$OMP & shared(Vp_plus_matrix,Vp_minus_matrix,scalebroad,onlyharmonic,Ngrid) &
    !$OMP & private(nn,mm,i,ll,q,plus_offset,minus_offset,omega_sp,ii,j,k,qprime,qdprime,ss,omegap_sp,omegadp_sp) &
    !$OMP & private(fBEprime_sp,fBEdprime_sp,vel_diff_sp,sigma_sp,delta_sp,WP3_local_sp,Gamma_plus_sp,Gamma_minus_sp,Vp_sp)
    do nn=1,nstates
       mm=myid*nstates+nn
       if (mm.gt.nlist*nbands) cycle
       i=modulo(mm-1,Nbands)+1
       ll=int((mm-1)/Nbands)+1
       q=IJK(:,list(ll))

       omega_sp = energy_sp(i,List(ll))
       if ((omega_sp .ne. 0.0_sp) .and. (real(omega_sp,kind=dp).le.omega_max)) then
          Gamma_plus_sp = 0.0_sp
          Gamma_minus_sp = 0.0_sp
          WP3_plus_sp(i,ll) = 0.0_sp
          WP3_minus_sp(i,ll) = 0.0_sp
          plus_offset = Naccum_plus_array(nn)
          minus_offset = Naccum_minus_array(nn)

          do ii=1,nptk
             qprime=IJK(:,ii)

             qdprime=q+qprime
             qdprime=modulo(qdprime,Ngrid)
             ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
             do j=1,Nbands
                omegap_sp=energy_sp(j,ii)
                if (omegap_sp .eq. 0.0_sp) cycle
                fBEprime_sp=one_sp/(exp(hbar_sp*omegap_sp/(Kb_sp*T_sp))-one_sp)

                do k=1,Nbands
                   omegadp_sp=energy_sp(k,ss)
                   if (omegadp_sp .eq. 0.0_sp) cycle
                   vel_diff_sp(:)=velocity_sp(:,j,ii)-velocity_sp(:,k,ss)
                   sigma_sp=real(scalebroad,kind=sp)*base_sigma_sp(vel_diff_sp(1),vel_diff_sp(2),vel_diff_sp(3))
                   if (sigma_sp .le. 0.0_sp) cycle
                   delta_sp=omega_sp+omegap_sp-omegadp_sp
                   if(abs(delta_sp).le.(two_sp*sigma_sp)) then
                      plus_offset = plus_offset + 1_long
                      fBEdprime_sp=one_sp/(exp(hbar_sp*omegadp_sp/(Kb_sp*T_sp))-one_sp)
                      WP3_local_sp=(fBEprime_sp-fBEdprime_sp)* &
                           exp(-(delta_sp*delta_sp)/(sigma_sp*sigma_sp))/ &
                           (sigma_sp*sqrt_pi_sp*omega_sp*omegap_sp*omegadp_sp)
                      WP3_plus_sp(i,ll)=WP3_plus_sp(i,ll)+WP3_local_sp
                      if (.not.onlyharmonic) then
                         Vp_sp=real(Vp_plus_matrix(plus_offset),kind=sp)
                         Gamma_plus_sp=Gamma_plus_sp+ hbarp_sp*pi_sp/four_sp*WP3_local_sp*Vp_sp*Vp_sp
                      end if
                   end if
                end do
             end do

             qdprime=q-qprime
             qdprime=modulo(qdprime,Ngrid)
             ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
             do j=1,Nbands
                omegap_sp=energy_sp(j,ii)
                if (omegap_sp .eq. 0.0_sp) cycle
                fBEprime_sp=one_sp/(exp(hbar_sp*omegap_sp/(Kb_sp*T_sp))-one_sp)

                do k=1,Nbands
                   omegadp_sp=energy_sp(k,ss)
                   if (omegadp_sp .eq. 0.0_sp) cycle
                   vel_diff_sp(:)=velocity_sp(:,j,ii)-velocity_sp(:,k,ss)
                   sigma_sp=real(scalebroad,kind=sp)*base_sigma_sp(vel_diff_sp(1),vel_diff_sp(2),vel_diff_sp(3))
                   if (sigma_sp .le. 0.0_sp) cycle
                   delta_sp=omega_sp-omegap_sp-omegadp_sp
                   if(abs(delta_sp).le.(two_sp*sigma_sp)) then
                      minus_offset = minus_offset + 1_long
                      fBEdprime_sp=one_sp/(exp(hbar_sp*omegadp_sp/(Kb_sp*T_sp))-one_sp)
                      WP3_local_sp=(fBEprime_sp+fBEdprime_sp+one_sp)* &
                           exp(-(delta_sp*delta_sp)/(sigma_sp*sigma_sp))/ &
                           (sigma_sp*sqrt_pi_sp*omega_sp*omegap_sp*omegadp_sp)
                      WP3_minus_sp(i,ll)=WP3_minus_sp(i,ll)+WP3_local_sp
                      if (.not.onlyharmonic) then
                         Vp_sp=real(Vp_minus_matrix(minus_offset),kind=sp)
                         Gamma_minus_sp=Gamma_minus_sp+hbarp_sp*pi_sp/four_sp*WP3_local_sp*Vp_sp*Vp_sp
                      end if
                   end if
                end do
             end do
          end do
          rate_scatt_plus_sp(i,ll)=Gamma_plus_sp*five_point_six_zero_sp*one_eight_sp/nptk_sp
          rate_scatt_minus_sp(i,ll)=Gamma_minus_sp*five_point_six_zero_sp*one_eight_sp/nptk_sp*half_sp
          WP3_plus_sp(i,ll)=WP3_plus_sp(i,ll)/nptk_sp
          WP3_minus_sp(i,ll)=WP3_minus_sp(i,ll)*half_sp/nptk_sp
       endif
    end do
    !$OMP END PARALLEL DO

    rate_scatt_plus = real(rate_scatt_plus_sp, kind=dp)
    rate_scatt_minus = real(rate_scatt_minus_sp, kind=dp)
    WP3_plus = real(WP3_plus_sp, kind=dp)
    WP3_minus = real(WP3_minus_sp, kind=dp)

    call MPI_ALLREDUCE(MPI_IN_PLACE,rate_scatt_plus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,rate_scatt_minus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,WP3_plus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(MPI_IN_PLACE,WP3_minus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    rate_scatt=rate_scatt_plus+rate_scatt_minus

    deallocate(energy_sp,velocity_sp,WP3_plus_sp,WP3_minus_sp)
    deallocate(rate_scatt_plus_sp,rate_scatt_minus_sp)
  end subroutine RTA_driver

  subroutine free_cached_event_lists()
    if (allocated(plus_i))  deallocate(plus_i,plus_kp,plus_ii,plus_j,plus_k,plus_ss)
    if (allocated(minus_i)) deallocate(minus_i,minus_kp,minus_ii,minus_j,minus_k,minus_ss)
    plus_count   = 0_long
    minus_count  = 0_long
  end subroutine free_cached_event_lists
!-------------------------In one fortran file to avoid rdc and improve the chance of inline---------------------------------------
  !$acc routine(base_sigma_sp) seq
  function base_sigma_sp(vx,vy,vz)
    implicit none
    real(kind=sp),intent(in) :: vx,vy,vz

    real(kind=sp) :: base_sigma_sp

    integer :: nu
    
    base_sigma_sp=0.0_sp
    do nu=1,3
       base_sigma_sp=base_sigma_sp+((rlattvec_sp(1,nu)*vx + rlattvec_sp(2,nu)*vy + rlattvec_sp(3,nu)*vz)/real(ngrid(nu),kind=sp))**2
    end do

    base_sigma_sp=sqrt(base_sigma_sp/6.0)
  end function base_sigma_sp
!$acc routine(Vp_plus_sp) seq
real(kind=sp) function Vp_plus_sp(eigenvect_i,eigenvect_j,eigenvect_k, &
     Ntri,Phi,phase_Rj,phase_Rk,Index_i,Index_j,Index_k)
  implicit none
  
  complex(kind=sp),intent(in) :: eigenvect_i(Nbands),eigenvect_j(Nbands),eigenvect_k(Nbands)
  integer,intent(in) :: Ntri
  real(kind=sp),intent(in) :: Phi(3,3,3,Ntri)
  complex(kind=sp),intent(in) :: phase_Rj(Ntri),phase_Rk(Ntri)
  integer,intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)

  integer :: ll, rr, bi, bj, bk
  complex(kind=sp) :: xi(3), xj(3), xk(3)
  complex(kind=sp) :: t1, f2, term, acc
  complex(kind=sp) :: phase_prod

  acc = (0.0_sp, 0.0_sp)

  do ll=1,Ntri
     bi = 3*(Index_i(ll)-1); bj = 3*(Index_j(ll)-1); bk = 3*(Index_k(ll)-1)
     xi(1)=eigenvect_i(bi+1); xi(2)=eigenvect_i(bi+2); xi(3)=eigenvect_i(bi+3)
     xj(1)=eigenvect_j(bj+1); xj(2)=eigenvect_j(bj+2); xj(3)=eigenvect_j(bj+3)
     xk(1)=conjg(eigenvect_k(bk+1)); xk(2)=conjg(eigenvect_k(bk+2)); xk(3)=conjg(eigenvect_k(bk+3))

     phase_prod = phase_Rj(ll) * conjg(phase_Rk(ll))

     term = (0.0_sp, 0.0_sp)
     do rr=1,3
        t1 = Phi(1,1,rr,ll)*xi(1) + Phi(2,1,rr,ll)*xi(2) + Phi(3,1,rr,ll)*xi(3)
        f2 = t1*xj(1)
        t1 = Phi(1,2,rr,ll)*xi(1) + Phi(2,2,rr,ll)*xi(2) + Phi(3,2,rr,ll)*xi(3)
        f2 = f2 + t1*xj(2)
        t1 = Phi(1,3,rr,ll)*xi(1) + Phi(2,3,rr,ll)*xi(2) + Phi(3,3,rr,ll)*xi(3)
        f2 = f2 + t1*xj(3)
        term = term + f2*xk(rr)
     end do

     acc = acc + term*phase_prod
  end do

  Vp_plus_sp = abs(acc)
end function Vp_plus_sp


!$acc routine(Vp_minus_sp) seq
real(kind=sp) function Vp_minus_sp(eigenvect_i,eigenvect_j,eigenvect_k, &
     Ntri,Phi,phase_Rj,phase_Rk,Index_i,Index_j,Index_k)
  implicit none

  complex(kind=sp),intent(in) :: eigenvect_i(Nbands),eigenvect_j(Nbands),eigenvect_k(Nbands)
  integer,intent(in) :: Ntri
  real(kind=sp),intent(in) :: Phi(3,3,3,Ntri)
  complex(kind=sp),intent(in) :: phase_Rj(Ntri),phase_Rk(Ntri)
  integer,intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)

  integer :: ll, rr, bi, bj, bk
  complex(kind=sp) :: xi(3), xj(3), xk(3)
  complex(kind=sp) :: t1, f2, term, acc
  complex(kind=sp) :: phase_prod

  acc = (0.0_sp, 0.0_sp)

  do ll=1,Ntri
     bi = 3*(Index_i(ll)-1); bj = 3*(Index_j(ll)-1); bk = 3*(Index_k(ll)-1)
     xi(1)=eigenvect_i(bi+1); xi(2)=eigenvect_i(bi+2); xi(3)=eigenvect_i(bi+3)
     xj(1)=conjg(eigenvect_j(bj+1)); xj(2)=conjg(eigenvect_j(bj+2)); xj(3)=conjg(eigenvect_j(bj+3))
     xk(1)=conjg(eigenvect_k(bk+1)); xk(2)=conjg(eigenvect_k(bk+2)); xk(3)=conjg(eigenvect_k(bk+3))

     ! conjg(phase_Rj)*conjg(phase_Rk) = conjg(phase_Rj*phase_Rk)
     phase_prod = conjg( phase_Rj(ll) * phase_Rk(ll) )

     term = (0.0_sp, 0.0_sp)
     do rr=1,3
        t1 = Phi(1,1,rr,ll)*xi(1) + Phi(2,1,rr,ll)*xi(2) + Phi(3,1,rr,ll)*xi(3)
        f2 = t1*xj(1)
        t1 = Phi(1,2,rr,ll)*xi(1) + Phi(2,2,rr,ll)*xi(2) + Phi(3,2,rr,ll)*xi(3)
        f2 = f2 + t1*xj(2)
        t1 = Phi(1,3,rr,ll)*xi(1) + Phi(2,3,rr,ll)*xi(2) + Phi(3,3,rr,ll)*xi(3)
        f2 = f2 + t1*xj(3)
        term = term + f2*xk(rr)
     end do

     acc = acc + term*phase_prod
  end do

  Vp_minus_sp = abs(acc)
end function Vp_minus_sp
!--------------------------------------------------------------------------------------------

end module processes
