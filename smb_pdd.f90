module smb_pdd
    ! This module contains all subroutines related to the calculation of
    ! surface mass balance using the PDD method.

    use smbpal_precision 

    implicit none

    real(prec), parameter :: sec_day = 86400.0   ! [sec]
    real(prec), parameter :: rho_w   = 1.d3      ! Density of pure water [kg/m3]
    real(prec), parameter :: L_m     = 3.35e5    ! Latent heat of melting [J/kg]

contains 


    elemental subroutine calc_ablation_pdd(abl,sif,pdds,acc,csnow,cice,csi)
        ! Determine total ablation based on input pdds
        
        implicit none 

        real(prec), intent(INOUT) :: abl, sif 
        real(prec), intent(IN)    :: pdds, acc
        real(prec), intent(IN)    :: csnow, csi, cice 

        real(prec) :: simax
        real(prec) :: abl_snow_pot, abl_ice_pot 

        ! (* maximum amount of super. ice that can be formed *)
        simax = acc*csi
       
        ! Determine the potential snowmelt (m) from the 
        ! available pdds, then determine the potential ice melt (m)
        ! from the left over pdds, if there are any.
        abl_snow_pot = pdds*csnow 
        abl_ice_pot  = max(0.0, (abl_snow_pot-acc)*cice/csnow)

        ! Get the mass balance
        if (abl_snow_pot .le. simax) then 
            abl = 0.0 
        else if (abl_snow_pot .gt. simax .and. abl_snow_pot .le. acc) then 
            abl = abl_snow_pot - simax 
        else if (abl_snow_pot .gt. acc .and. abl_ice_pot .le. simax) then 
            abl = abl_snow_pot + (simax-abl_ice_pot)
        else 
            abl = abl_snow_pot + (abl_ice_pot-simax) 
        end if  


        ! Calculate the actual superimposed ice from refreezing 
        if (abl_snow_pot .le. simax) then   
            sif=pdds*csnow
        else 
            sif=simax 
        endif

        return 

    end subroutine calc_ablation_pdd


    elemental function calc_temp_effective(temp, sigma) result(teff)
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Subroutine : e f f e c t i v e T
        ! Author     : Reinhard Calov
        ! Purpose    : Computation of the positive degree days (PDD) with
        !              statistical temperature fluctuations;
        !              based on semi-analytical solution by Reinhard Calov.
        !              This subroutine uses days as time unit, each day
        !              is added individually
        !              (the same sigma as for pdd monthly can be used)
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! temp = [degrees Celcius] !!! 

        implicit none

        real(prec), intent(IN) :: temp, sigma
        real(prec) :: teff
        
        real(prec) :: temp_c, inv_sigma
        real(prec), parameter :: inv_sqrt2   = 1.0/sqrt(2.0)
        real(prec), parameter :: inv_sqrt2pi = 1.0/sqrt(2.0*pi)

        inv_sigma   = 1.0/sigma

        teff = sigma*inv_sqrt2pi*exp(-0.5*(temp*inv_sigma)**2)  &
                  + temp*0.5*erfcc(-temp*inv_sigma*inv_sqrt2)

        ! Result is the assumed/felt/effective positive degrees, 
        ! given the actual temperature (accounting for fluctuations in day/month/etc, 
        ! based on the sigma chosen)
                     
        return

    end function calc_temp_effective 

        
    elemental function erfcc(x)
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Function :  e r f c c
        ! Author   :  Reinhard Calov and Ralf Greve
        ! Purpose  :  Returns the complementary error function erfc(x) with 
        !             fractional error everywhere less than 1.2 x 10^(-7).
        !             Credit: Press et al., 'Numerical recipes in Fortran 77'.
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        implicit none
            
        real(prec), intent(IN) :: x
        real(prec) :: erfcc

        real(prec) :: t, z

        z = abs(x)
        t = 1.0/(1.0+0.5*z)

        erfcc = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+  &
        t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*  &
        (-1.13520398+t*(1.48851587+t*(-0.82215223+  &
        t*0.17087277)))))))))

        if (x .lt. 0.0) erfcc = 2.0-erfcc

        return
      
    end function erfcc

end module smb_pdd 
