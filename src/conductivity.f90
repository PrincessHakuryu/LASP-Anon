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

! Compute the thermal conductivity and its cumulative version.
module conductivity
  use config
  use data, only : dp, hbar, Kb
  implicit none

contains

  ! Straightforward implementation of the thermal conductivity as an integral
  ! over the whole Brillouin zone in terms of frequencies, velocities and F_n.
  subroutine TConduct(omega,velocity,F_n,ThConductivity,ThConductivityMode)
    implicit none

    real(kind=dp),intent(in) :: omega(Nbands,nptk),velocity(3,Nbands,nptk),F_n(3,Nbands,nptk)
    real(kind=dp),intent(out) :: ThConductivity(3,3,nbands)
    real(kind=dp),intent(out) :: ThConductivityMode(3,3,nbands,nptk)

    real(kind=dp) :: fBE,tmp(3,3)
    integer :: ii,jj,dir1,dir2

    ThConductivity=0.0_dp
    ThConductivityMode=0.0_dp
    !$OMP PARALLEL DO default(none) collapse(2) schedule(static) &
    !$OMP & shared(nptk,nbands,omega,velocity,F_n,ThConductivityMode,T) &
    !$OMP & private(jj,ii,dir1,dir2,tmp,fBE) reduction(+:ThConductivity)
    do ii=2,nptk
       do jj=1,Nbands
          do dir1=1,3
             do dir2=1,3
                tmp(dir1,dir2)=velocity(dir1,jj,ii)*F_n(dir2,jj,ii)
             end do
          end do
          fBE=1.0_dp/(exp(hbar*omega(jj,ii)/Kb/T)-1.0_dp)
          ThConductivityMode(:,:,jj,ii)=fBE*(fBE+1.0_dp)*omega(jj,ii)*tmp
          ThConductivity(:,:,jj)=ThConductivity(:,:,jj)+ThConductivityMode(:,:,jj,ii)
       end do
    end do
    !$OMP END PARALLEL DO
    ThConductivity=1.0e21_dp*hbar**2*ThConductivity/(kB*T*T*V*nptk)
    ThConductivityMode=1.0e21_dp*hbar**2*ThConductivityMode/(kB*T*T*V*nptk)
  end subroutine TConduct

  ! Specialized version of the above subroutine for those cases where kappa
  ! can be treated as a scalar.
  subroutine TConductScalar(omega,velocity,F_n,ThConductivity)
    implicit none

    real(kind=dp),intent(in) :: omega(Nbands,nptk),velocity(Nbands,nptk),F_n(Nbands,nptk)
    real(kind=dp),intent(out) :: ThConductivity(Nbands)

    real(kind=dp) :: fBE
    integer :: ii,jj

    ThConductivity=0.0_dp
    !$OMP PARALLEL DO default(none) collapse(2) schedule(static) &
    !$OMP & shared(nptk,nbands,omega,velocity,F_n,T) &
    !$OMP & private(jj,ii,fBE) reduction(+:ThConductivity)
    do ii=2,nptk
       do jj=1,Nbands
          fBE=1.0_dp/(exp(hbar*omega(jj,ii)/Kb/T)-1.0_dp)
          ThConductivity(jj)=ThConductivity(jj)+fBE*(fBE+1.0_dp)*omega(jj,ii)*&
               velocity(jj,ii)*F_n(jj,ii)
       end do
    end do
    !$OMP END PARALLEL DO
    ThConductivity=1.0e21_dp*hbar**2*ThConductivity/(kB*T*T*V*nptk)
  end subroutine TConductScalar

  ! "Cumulative thermal conductivity": value of kappa obtained when
  ! only phonon with mean free paths below a threshold are considered.
  ! ticks is a list of the thresholds to be employed.
  subroutine CumulativeTConduct(omega,velocity,F_n,ticks,results)
    implicit none

    real(kind=dp),intent(in) :: omega(Nbands,nptk),velocity(3,Nbands,nptk),F_n(3,Nbands,nptk)
    real(kind=dp),intent(out) :: ticks(nticks),results(3,3,nbands,Nticks)

    real(kind=dp) :: fBE,tmp(3,3),lambda
    integer :: ii,jj,kk,dir1,dir2

    real(kind=dp) :: dnrm2

    do ii=1,nticks
       ticks(ii)=10.0_dp**(8.0_dp*(real(ii,kind=dp)-1.0_dp)/(real(nticks,kind=dp)-1.0_dp)-2.0_dp)
    end do

    results=0.0_dp
    !$OMP PARALLEL DO default(none) collapse(2) schedule(static) &
    !$OMP & shared(nptk,nbands,nticks,omega,velocity,F_n,T,ticks) &
    !$OMP & private(jj,ii,kk,dir1,dir2,tmp,fBE,lambda) reduction(+:results)
    do ii=2,nptk
       do jj=1,Nbands
          do dir1=1,3
             do dir2=1,3
                tmp(dir1,dir2)=velocity(dir1,jj,ii)*F_n(dir2,jj,ii)
             end do
          end do
          fBE=1.0_dp/(exp(hbar*omega(jj,ii)/Kb/T)-1.0_dp)
          tmp=fBE*(fBE+1.0_dp)*omega(jj,ii)*tmp
          lambda=dot_product(F_n(:,jj,ii),velocity(:,jj,ii))/(&
               (omega(jj,ii)+1.0e-12_dp)*dnrm2(3,velocity(:,jj,ii),1))
          do kk=1,nticks
             if(ticks(kk).gt.lambda) then
                results(:,:,jj,kk)=results(:,:,jj,kk)+tmp
             end if
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    results=1.0e21_dp*hbar**2*results/(kB*T*T*V*nptk)
  end subroutine CumulativeTConduct

  ! "Cumulative thermal conductivity vs angular frequency": value of kappa obtained when
  ! only frequencies below a threshold are considered.
  ! ticks is a list of the thresholds to be employed.
  subroutine CumulativeTConductOmega(omega,velocity,F_n,ticks,results)
    implicit none

    real(kind=dp),intent(in) :: omega(Nbands,nptk),velocity(3,Nbands,nptk),F_n(3,Nbands,nptk)
    real(kind=dp),intent(out) :: ticks(nticks),results(3,3,nbands,Nticks)

    real(kind=dp) :: fBE,tmp(3,3),lambda
    integer :: ii,jj,kk,dir1,dir2

    REAL(kind=dp)  EMIN,EMAX

    EMIN=1.0e10_dp
    EMAX=-1.0e10_dp
    DO ii=1,NPTK
    DO jj=1,Nbands
       emin = min(emin, omega(jj,ii))
       emax = max(emax, omega(jj,ii))
    ENDDO
    ENDDO
    emax=emax*1.1_dp

    do ii=1,nticks
       ticks(ii)=emin+(emax-emin)*ii/nticks
    end do

    results=0.0_dp
    !$OMP PARALLEL DO default(none) collapse(2) schedule(static) &
    !$OMP & shared(nptk,nbands,nticks,omega,velocity,F_n,T,ticks) &
    !$OMP & private(jj,ii,kk,dir1,dir2,tmp,fBE,lambda) reduction(+:results)
    do ii=2,nptk
       do jj=1,Nbands
          do dir1=1,3
             do dir2=1,3
                tmp(dir1,dir2)=velocity(dir1,jj,ii)*F_n(dir2,jj,ii)
             end do
          end do
          fBE=1.0_dp/(exp(hbar*omega(jj,ii)/Kb/T)-1.0_dp)
          tmp=fBE*(fBE+1.0_dp)*omega(jj,ii)*tmp
          lambda=omega(jj,ii)
          do kk=1,nticks
             if(ticks(kk).gt.lambda) then
                results(:,:,jj,kk)=results(:,:,jj,kk)+tmp
             end if
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    results=1e21_dp*hbar**2*results/(kB*T*T*V*nptk)
  end subroutine CumulativeTConductOmega
end module conductivity
