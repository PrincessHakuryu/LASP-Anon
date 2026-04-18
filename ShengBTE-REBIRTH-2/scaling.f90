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

! Correction to the bulk tau required for nanowires.

module scaling
  use config
  use data
  implicit none

contains
  ! Correction to tau for nanowires, based on a radial average. Note
  ! that a direction needs to be selected before calling this
  ! subroutine: velocity_z is the projection of velocity over that
  ! axis.
  subroutine ScalingOfTau(Nlist,Nequi,ALLEquiList,velocity_z,&
       velocity,tauzero_wedge,radnw,ffunc)
    implicit none
    integer,intent(in) :: Nlist,Nequi(nptk),ALLEquiList(Nsymm_rot,nptk)
    real(kind=dp),intent(in) :: velocity(3,Nbands,nptk),velocity_z(Nbands,nptk)
    real(kind=dp),intent(in) :: tauzero_wedge(Nbands,Nlist),radnw
    real(kind=dp),intent(out) :: ffunc(Nbands,nptk)

    real(kind=dp) :: vrho,xsum,tauint,vmod,vaxial,tauzero(Nbands,nptk)
    integer,parameter :: nsum=500
    integer :: iisum,ii,jj

    real(kind=dp) :: dnrm2

    do ii=1,Nlist
      do jj=1,Nequi(ii)
          tauzero(:,ALLEquiList(jj,ii))=tauzero_wedge(:,ii)
       end do
    end do

    ffunc=0.0_dp
    !$OMP PARALLEL DO default(none) collapse(2) schedule(static) &
    !$OMP & shared(nptk,nbands,T,tauzero,ffunc,radnw,cgrid,velocity,velocity_z) &
    !$OMP & private(jj,ii,vrho,vmod,tauint,vaxial,xsum)
    do ii=1,nptk
       do jj=1,Nbands
          vmod=dnrm2(3,velocity(:,jj,ii),1)
          tauint=abs(tauzero(jj,ii))
          vaxial=velocity_z(jj,ii)
          if ((vmod**2-vaxial**2).le.vaxial**2*1.0e-8_dp) then
             vaxial=vaxial*2.0_dp*cgrid/sqrt(4.0_dp*(cgrid**2)+1.0_dp)
          end if
          vrho=tauint*sqrt(vmod**2-vaxial**2)
          if ((vrho /= 0.0_dp)) then
             ffunc(jj,ii)=0.0_dp
             do iisum=1,nsum
                xsum=real(iisum,kind=dp)*radnw/real(nsum,kind=dp)
                ffunc(jj,ii)=ffunc(jj,ii)+2.0_dp*sqrt(abs(radnw**2-xsum**2))/vrho-&
                     1.0_dp+exp(-2.0_dp*sqrt(abs(radnw**2-xsum**2))/vrho)
             end do
             ffunc(jj,ii)=(vrho*2.0_dp/(pi*radnw*real(nsum,kind=dp)))*ffunc(jj,ii)
          end if
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine ScalingOfTau
end module scaling
