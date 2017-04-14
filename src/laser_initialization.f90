!*****************************************************************************************************!
!             Copyright 2014-2016 Alberto Marocchino, Francesco Massimo                               !
!*****************************************************************************************************!

!*****************************************************************************************************!
!  This file is part of architect.                                                                    !
!                                                                                                     !
!  Architect is free software: you can redistribute it and/or modify                                  !
!  it under the terms of the GNU General Public License as published by                               !
!  the Free Software Foundation, either version 3 of the License, or                                  !
!  (at your option) any later version.                                                                !
!                                                                                                     !
!  Architect is distributed in the hope that it will be useful,                                       !
!  but WITHOUT ANY WARRANTY; without even the implied warranty of                                     !
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                      !
!  GNU General Public License for more details.                                                       !
!                                                                                                     !
!  You should have received a copy of the GNU General Public License                                  !
!  along with architect.  If not, see <http://www.gnu.org/licenses/>.                                 !
!*****************************************************************************************************!

 module initialize_laser

	USE my_types
	USE use_my_types


 implicit none

 !--- --- ---!
 contains
 !--- --- ---!
 subroutine init_gaussian_laser_envelope

 	REAL(8)    :: sigma_r_laser,sigma_z_laser,a0_Architect_units,conversion_factor
 	INTEGER :: i,j

  write(*,'(A)') ' --- Laser initialisation'

 	! adimensional laser pulse dimensions
 	sigma_r_laser 			= laser_initialization%w0*plasma%k_p
 	sigma_z_laser 			= laser_initialization%FWHM_fs*0.3/2.355*plasma%k_p
 	conversion_factor 	= 1. !(2.*pi)/LaserPulse%lambda0/plasma%k_p


  write(*,'(A,f11.3)')'Laser wavelength (um) =',laser_initialization%lambda0
  write(*,'(A,f11.3)')'Laser a0 =',laser_initialization%a0                  ! normalized peak vector potential
  write(*,'(A,f11.3)')'Laser w0 (um) =',laser_initialization%w0             ! 2*sigma_r_laser
  write(*,'(A,f11.3)')'Laser duration (fs) =',laser_initialization%FWHM_fs  ! FWHM in intensity

 	!write(*,*) LaserPulse%lambda0,plasma%k_p,conversion_factor,LaserPulse%a0
 	! laser parameter a0, with proper units
 	! remember that usually a0 is defined as a0 = me * omega0 * c /e * E0,
 	! but Architect uses omegap instead of omega0 in the normalization, thus conversion is necessary
 	a0_Architect_units 	= laser_initialization%a0*conversion_factor
 	write(*,'(A,f11.3)')'Laser a0 (architect units)=',laser_initialization%a0 ! 2*sigma_r_laser

 	! Laser envelope, centered in z=0
 	do i= Node_min_z,Node_max_z
     	do j= Node_min_r,Node_max_r
       	mesh(i,j)%a    = a0_Architect_units*sqrt(exp(-x_mesh(j)**2./2./sigma_r_laser**2)*exp(-z_mesh(i)**2./2./sigma_z_laser**2))
     	enddo
 	enddo

 	! BC on lower boundary
 	mesh(:,1)%a = mesh(:,2)%a

 	return

 end subroutine init_gaussian_laser_envelope


 !--- --- ---!
 end module initialize_laser
 !--- --- ---!
