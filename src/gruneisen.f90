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

! Subroutines and functions related to Grüneisen parameters.
module gruneisen
  use misc
  use data, only : dp ,hbar, kb
  use config
  use input
  implicit none

contains

  ! Subroutine to compute the mode Grüneisen parameters.
  subroutine mode_grun(omega,eigenvect,Ntri,Phi,&
       R_j,R_k,Index_i,Index_j,Index_k,grun)
    implicit none

    integer,intent(in) :: Ntri,Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=dp),intent(in) :: omega(nbands,nptk)
    real(kind=dp),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=dp),intent(in) :: eigenvect(nbands,nbands,nptk)
    real(kind=dp),intent(out) :: grun(nbands,nptk)

    real(kind=dp),parameter :: unitfactor=9.6472e4_dp ! From nm*eV/(amu*A^3*THz^2) to 1.

    integer :: ik,ii,jj,kk,iband,itri,ialpha,ibeta
    real(kind=dp) :: kspace(nptk,3)
    complex(kind=dp) :: factor1,factor2,factor3,g

    do ii=1,Ngrid(1)
       do jj=1,Ngrid(2)
          do kk=1,Ngrid(3)
             ik=((kk-1)*Ngrid(2)+(jj-1))*Ngrid(1)+ii
             kspace(ik,:)=rlattvec(:,1)*(ii-1.0_dp)/real(ngrid(1),kind=dp)+&
                  rlattvec(:,2)*(jj-1.0_dp)/real(ngrid(2),kind=dp)+&
                  rlattvec(:,3)*(kk-1.0_dp)/real(ngrid(3),kind=dp)
          end do
       end do
    end do    
    !$OMP PARALLEL DO default(none) collapse(2) schedule(static) &
    !$OMP & shared(nptk,nbands,Ntri,kspace,grun,omega,Phi,types,cartesian) &
    !$OMP & shared(Index_i,Index_j,Index_k,eigenvect,masses,R_k,R_j) &
    !$OMP & private(ik,iband,itri,ialpha,ibeta,factor1,factor2,factor3,g)
    do ik=1,nptk
       do iband=1,nbands
          g=0.0_dp
          do itri=1,Ntri
             factor1=phexp(dot_product(kspace(ik,:),R_j(:,itri)))/&
                  sqrt(masses(types(Index_i(itri)))*&
                  masses(types(Index_j(itri))))
             do ialpha=1,3
                factor2=factor1*conjg(eigenvect(3*(Index_i(itri)-1)+ialpha,&
                                                iband,ik))
                do ibeta=1,3
                   factor3=factor2*eigenvect(3*(Index_j(itri)-1)+ibeta,& 
                                             iband,ik)
                   g=g+factor3*dot_product(&
                        Phi(ialpha,ibeta,:,itri),&
                        (cartesian(:,Index_k(itri))+R_k(:,itri)))
                end do
             end do
          end do
         if (omega(iband,ik).eq.0.0_dp) then
            grun(iband,ik)=0.0_dp
         else
            g=-unitfactor*g/6.0_dp/omega(iband,ik)**2
            grun(iband,ik)=real(g)
         endif
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine mode_grun

  ! Obtain the total Grüneisen parameter as a weighted sum over modes.
  function total_grun(omega,grun)
    implicit none
    real(kind=dp),intent(in) :: omega(nbands,nptk),grun(nbands,nptk)

    integer :: ii,jj
    real(kind=dp) :: total_grun,weight,dBE,x

    total_grun=0.0_dp
    weight=0.0_dp
    !$OMP PARALLEL DO default(none) collapse(2) schedule(static) &
    !$OMP & shared(nptk,nbands,omega,grun,T) private(x,dBE,jj,ii) &
    !$OMP & reduction(+:total_grun) reduction(+:weight)
    do ii=1,nptk
       do jj=1,nbands
          x=hbar*omega(jj,ii)/(2.0_dp*kB*T)
          if(x.gt.1e-6_dp) then
             dBE=(x/sinh(x))**2
             weight=weight+dBE
             total_grun=total_grun+dBE*grun(jj,ii)
          end if
       end do
    end do
    !$OMP END PARALLEL DO
    total_grun=total_grun/weight
  end function total_grun
end module gruneisen
