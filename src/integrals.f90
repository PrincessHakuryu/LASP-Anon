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

! Lattice specific heat and small-grain-limit reduced thermal
! conductivity. Both are calculated as integrals over the BZ and contain
! no anharmonic information.
module integrals
  use config
  use data
  implicit none

contains

  ! Lattice specific heat per unit volume.
  function cv(omega)
    implicit none
    real(kind=dp),intent(in) :: omega(nbands,nptk)

    integer :: ii,jj
    real(kind=dp) :: cv,dBE,x

    cv=0.0_dp
    !$OMP PARALLEL DO default(none) schedule(static) collapse(2) &
    !$OMP & shared(omega,T,nbands,nptk) private(x,dBE,ii,jj) reduction(+:cv)
    do jj=1,nptk
       do ii=1,nbands
          x=hbar*omega(ii,jj)/(2.0_dp*kB*T)
          if(x.eq.0.0_dp) then
             dBE=1.0_dp
          else
             dBE=(x/sinh(x))**2
          end if
          cv=cv+dBE
       end do
    end do
    !$OMP END PARALLEL DO
    cv=kB*cv/(1.0e-27_dp*V*real(nptk,kind=dp)) ! J/(K m^3)
  end function cv

  ! Small-grain-limit reduced thermal conductivity. In contrast with cv,
  ! group velocities play a role in this integral.
  subroutine kappasg(omega,velocity,nruter)
    implicit none
    real(kind=dp),intent(in) :: omega(nbands,nptk),velocity(3,nbands,nptk)
    real(kind=dp),intent(out) :: nruter(3,3)

    integer :: ii,jj,dir1,dir2
    real(kind=dp) :: x,dBE,tmp

    real(kind=dp) :: dnrm2

    nruter=0.0_dp
    !$OMP PARALLEL DO default(none) schedule(static) collapse(2)  &
    !$OMP & shared(omega,velocity,T,nbands,nptk) private(x,dBE,tmp,ii,jj,dir1,dir2) &
    !$OMP & reduction(+:nruter)
    do jj=2,nptk
       do ii=1,nbands
          x=hbar*omega(ii,jj)/(2.0_dp*kB*T)
          dBE=(x/sinh(x))**2
          tmp=dnrm2(3,velocity(:,ii,jj),1)
          if(tmp.lt.1.0e-12_dp) then
             cycle
          else
             dBE=dBE/tmp
          end if
          do dir1=1,3
             do dir2=1,3
                nruter(dir1,dir2)=nruter(dir1,dir2)+dBE*&
                     velocity(dir1,ii,jj)*velocity(dir2,ii,jj)
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    nruter=1.0e21_dp*kB*nruter/(real(nptk,kind=dp)*V)
  end subroutine kappasg
end module integrals
