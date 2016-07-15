module smb_itm
    ! This module contains all subroutines related to the calculation of
    ! surface mass balance. The external interface to this module is via
    ! the melt_budget() subroutine, which takes arrays from the main
    ! program as input and outputs arrays of the calculated variables.

    use smbpal_precision 

    implicit none

    real(prec), parameter :: sec_day = 86400.0   ! [sec]
    real(prec), parameter :: rho_w   = 1.d3      ! Density of pure water [kg/m3]
    real(prec), parameter :: L_m     = 3.35e5    ! Latent heat of melting [J/kg]

    type(itm_par)

    end type 

contains

  elemental subroutine snowpack_budget(par,H_snow,alb_s,smbi,smb,melt,runoff,refrz, &
                                        z_srf,H_ice,S,t2m,PDDs, pr, sf )
    ! Determine the total melt, accumulation and surface mass balance at a given point
    !  * Modified from rembo subroutine `melt_budget`
    !  * input in mm water equivalent
    ! Note: Definitions as in Ettema et al (2009) supplementary information,
    !     SMB  = snow + rain - runoff     [kg m2 / d ] == [mm / d]
    !   runoff = rain + melt - refrozen   [kg m2 / d ] == [mm / d]

    implicit none
    
    type(itm_par), intent(IN)    :: par 
    real(prec),    intent(INOUT) :: H_snow, alb_s 
    real(prec),    intent(INOUT) :: smbi, smb, melt, runoff, refrz 
    real(prec),    intent(IN)    :: z_srf, H_ice, S, t2m, PDDs 
    real(prec),    intent(IN)    :: pr, sf
    real(prec) :: refreezing_max, refrozen_snow_max
    real(prec) :: smb, smbi
    real(prec) :: melt_pot
    
    ! Local variables
    real(prec) :: rf, atrans 
    real(prec) :: rfac  
    real(prec) ::              refrozen_rain, runoff_rain
    real(prec) :: melted_snow, refrozen_snow, runoff_snow
    real(prec) :: melted_firn, refrozen_firn, runoff_firn
    real(prec) :: melted_ice, snow_to_ice, new_ice
    real(prec) :: melt, refrozen, runoff 

    ! Determine rainfall from precip and snowfall 
    rf = pr - sf

    ! Determine preliminary surface and planetary albedo
    alb_s = calc_albedo_surface(par,z_srf,H_ice,H_snow,pdds)

    ! Add additional snowfall 
    H_snow  = H_snow + sf
    
    ! Get amount of potential melt from ITM scheme
    atrans   = calc_atmos_transmissivity(z_srf,par%trans_a,par%trans_b)
    melt_pot = calc_itm(S,t2m,alb_s,atrans,par%itm_c,par%itm_t)

    ! Determine how much snow and ice would be melted today
    if (melt_pot .gt. H_snow)
      
      ! All snow is melted, the rest of energy converted to melt some ice
      ! The rest of energy will go into melting ice
      melted_snow = H_snow
      melted_ice  = melt_pot - H_snow

    else
      
      ! Snow melt will use all energy, none left for ice melt
      melted_snow = melt_pot
      melted_ice  = 0.d0
      
    end if    
    
    ! Now calculate the actual melt (total ablation)
    melt   = melted_snow + melted_ice
    
    ! Remove melted snow, if any, from the snow height budget
    H_snow = H_snow - melted_snow
    
    ! Adjust the albedo (accounting for actual amount of melt)
    alb_s = calc_albedo_surface(par,z_srf,H_ice,H_snow,pdds,melt=melt)

    ! Determine how much new ice is made from compression of remaining snow
    snow_to_ice = 0.d0
    if (H_snow .gt. par%H_snow_max)
      ! Assume excess contributes to new ice
      snow_to_ice = (H_snow - par%H_snow_max)
      
      ! Reset snow height down to maximum height
      H_snow = par%H_snow_max
    end if
    
    ! Determine what fraction of the melted snow and rain will refreeze, 
    ! (Note: rf is zero if not on ice sheet or there is no snow cover)
    if ( H_ice .gt. 0.d0 .and. H_snow .gt. 0.d0 )
    
        rfac = par%Pmaxfrac * sf / max(1d-5,pr)    ! max() here ensures no division by zero, if (snow+rain)==0, then rfac is zero anyway

        ! Modify refreezing factor based on height of snow
        if ( H_snow .gt. 2e3 )
            rfac = 1.d0                                         ! refreezing factor is 1 for large snow heights.      
        else if ( H_snow .gt. 1e3 )
            rfac = rfac + ((H_snow-1e3)/(2e3-1e3) ) * (1.d0 - rfac) ! linear function increasing to rf=1 as H_snow increases to 2m.
        end if
    
    else 
        rfac = 0.0 

    end if
   
    ! Determine the actual maximum amount of refreezing
    refreezing_max    = H_snow                                ! Total refreezing depends on amount of snow left!
    refrozen_rain     = min(rf*rfac,refreezing_max)           ! First rain takes up refreezing capacity
    refrozen_snow_max = refreezing_max - refrozen_rain        ! Subtract rain from refreezing capacity to determine what's left for snow
    refrozen_snow     = min(melted_snow*rfac,refrozen_snow_max) ! melted_snow uses remaining capacity if it needs it
    refrozen          = refrozen_snow + refrozen_rain         ! Amount of ice created from refreezing

    ! Determine how much water will runoff for each component
    runoff_snow = melted_snow - refrozen_snow                 ! Net snow melt
    runoff_rain = rf - refrozen_rain                          ! Net rainfall
    runoff      = runoff_snow + runoff_rain + melted_ice      ! Total runoff
    
    ! Get the total ice accumulated for the day
    new_ice = snow_to_ice + refrozen
    
    ! Get the global surface mass balance
    smb = snow + rain - runoff
    
    ! Get the internal surface mass balance (what ice sheet model needs)
    ! These should apply at separate places, so
    ! smbi = -melted_ice, for negative mass balance
    ! smbi =     new_ice, for positive mass balance
    smbi = new_ice - melted_ice
    
    ! Determine the actual melted amount of snow or snow+ice,
    ! depending on where melt is taking place (for energy-balance)
    ! Note: no refreezing off of ice sheet, so there total melt can
    ! only equal the amount of snow melted
    ! Note: this value can be negative too (negative contribution to energy flux)
    if ( mask .eq. 0 ) 
      dh = melt-refrozen       ! Technically refrozen=0 here, but include it anyway
    else
      dh = melted_snow-refrozen
    end if
    
    return
  
  end subroutine snowpack_budget
    
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine : g e t _ a l b e d o
    ! Author     : Alex Robinson
    ! Purpose    : Determine the current surface and planetary albedos
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    elemental function calc_albedo_surface(par,z_srf,H_ice,H_snow,PDDs,melt) result(alb)
    
        implicit none
        
        type(itm_par), intent(IN) :: par 
        real(prec),   intent(IN) :: z_srf, H_ice, H_snow, PDDs 
        real(prec),   intent(IN), optional :: melt
        real(prec) :: alb 

        ! Local variables
        real(prec) :: H_snow_crit, depth, as_snow, alb_bg, melt_now
        integer :: n, k
        
        ! Determine the critical snow height based on number of PDDs,
        ! which correspond to vegetation growth
        ! 0-100 PDDs    = desert, 10mm
        ! 100-1000 PDDs = tundra (grass), linearly increasing between 10mm => 100mm
        ! 1000 PDDs +   = forest, 100mm
        if ( PDDs .le. 100.0 ) 
            H_snow_crit = par%H_snow_crit_desert ! 10 mm
        else if (PDDs .le. 1000.0) 
            H_snow_crit = par%H_snow_crit_desert +   &
                (par%H_snow_crit_forest-par%H_snow_crit_desert) * (PDDs-100.0)/(1000.0-100.0)
        else
            H_snow_crit = par%H_snow_crit_forest ! 100 mm
        end if
        
        ! Determine the scaled snowdepth
        depth = min( H_snow / H_snow_crit, 1.0 ) 
        
        ! Determine the amount of melt to affect abledo calculation.
        ! Default is high melt everywhere, so that jump in albedo is avoided
        ! at onset of melt season. After calculating actual melt, this is included 
        ! as an argument and albedo is recalculated.
        melt_now = par%melt_crit+1.0
        if ( present(melt) ) melt_now = melt
    
        ! Figure out what the surface albedo would be with no snow, 
        ! based on the type of ground underneath ( ice or land )
        if (z_srf .le. 0.0) then 
            alb_bg = par%alb_ocean

        else if (z_srf .gt. 0.0 .and. H_ice .eq. 0.0) then 
             alb_bg = par%alb_land*(1d3-min(PDDs,1d3))/(1d3-0d0) &
                                        + par%alb_forest*(min(PDDs,1d3))/(1d3-0d0)
        else 
            alb_bg = par%alb_ice 

        end if

        ! Determine current snow albedo: if melt gt eg, 1mm/day, then change albedo to melting snow!
        as_snow = par%alb_snow_dry
        if (melt_now .gt. par%melt_crit) as_snow = par%alb_snow_wet  
        
        ! Get current surface albedo to be used for ITM melt scheme
        ! It will either be the maximum albedo (that of dry snow: as_snow)
        ! or the wet snow albedo plus a fraction depending on height of snow
        ! minimum albedo now should be that of ice / wet snow minimum
        alb = alb_bg + depth*(as_snow-alb_bg)
    
        return
    
    end function calc_albedo_surface
    
    elemental function calc_albedo_planet(alb_s,a,b) result(alb_p)

        implicit none 

        real(prec), intent(IN) :: alb_s, a, b 
        real(prec) :: alb_p 

        alb_p = a + b*alb_s
    
        return 

    end function calc_albedo_planet

    elemental function calc_atmos_transmissivity(z_srf,a,b) result(at)

        implicit none 

        real(prec), intent(IN) :: z_srf, a, b 
        real(prec) :: at 

        at = a + b*max(z_srf,0.d0)
    
        return 

    end function calc_atmos_transmissivity

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine : i t m (insolation temperature melt)
    ! Author     : Alex Robinson
    ! Purpose    : Determine the total melt, accumulation and 
    !              surface mass balance at a given point
    !              * output in mm water equivalent per day 
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    elemental function calc_itm(S,t2m,alb_s,atrans,c,t) result(melt)
        ! Determine the potential melt rate
        implicit none

        real(prec), intent(IN) :: S, t2m, alb_s, atrans, c, t
        real(prec) :: melt

        ! Calculate potential melt [m/s]
        melt = (atrans*(1.d0 - alb_s)*S + c + t*t2m) / (rho_w*L_m)

        ! Convert: [m/s] => [mm/day], only positive melt
        melt = max( melt, 0.d0 ) * sec_day * 1d3  

        return

    end function calc_itm

end module smb_itm 

