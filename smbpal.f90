 module smbpal

    use smbpal_precision
    use insolation
    use ncio
    use smb_itm 

    implicit none 

    real(4), parameter :: pi = 3.14159265359

    type smbpal_param_class
        type(itm_par_class) :: itm
        character(len=16)   :: abl_method 
        real(prec) :: Teff_sigma, sf_a, sf_b  

        real(prec), allocatable :: x(:), y(:)
        real(prec), allocatable :: lats(:,:)           ! Latitude of domain [deg N]
        real(prec) :: rho_sw
        real(prec) :: rho_ice

    end type 

    type smbpal_state_class 
        real(prec), allocatable   :: t2m(:,:)            ! Surface temperature [K]
        real(prec), allocatable   :: pr(:,:), sf(:,:)    ! Precip, snowfall [mm/a or mm/d]
        real(prec), allocatable   :: S(:,:)              ! Insolation [W/m2]
        real(prec), allocatable   :: teff(:,:)           ! Effective temp. (ie, PDDs) [num. of days]
        
        ! Prognostic variables
        real(prec), allocatable   :: H_snow(:,:)         ! Snow thickness [mm]
        real(prec), allocatable   :: alb_s(:,:)          ! Surface albedo 
        real(prec), allocatable   :: smbi(:,:), smb(:,:) ! Surface mass balance [mm/a or mm/d]
        real(prec), allocatable   :: melt(:,:), runoff(:,:), refrz(:,:)   ! smb components
    end type 

    type smbpal_class
        type(smbpal_param_class) :: par 
        type(smbpal_state_class) :: now, mon(12), ann
    end type

    interface smbpal_update 
        module procedure smbpal_update_2temp
        module procedure smbpal_update_monthly
        module procedure smbpal_update_daily
    end interface 

    private
    public :: smbpal_class
    public :: smbpal_init 
    public :: smbpal_update
    public :: smbpal_end 
    public :: smbpal_write_init, smbpal_write

contains 

    subroutine smbpal_init(smb,filename,x,y,lats)

        implicit none 

        type(smbpal_class) :: smb
        character(len=*), intent(IN)  :: filename  ! Parameter file 
        real(prec) :: x(:), y(:), lats(:,:)

        ! Local variables
        integer :: nx, ny, m  
        
        nx = size(x,1)
        ny = size(y,1)

        ! Load smbpal parameters
        call smbpal_par_load(smb%par,filename)

        ! Additionally define dimension info 
        if (allocated(smb%par%x)) deallocate(smb%par%x)
        if (allocated(smb%par%y)) deallocate(smb%par%y)
        if (allocated(smb%par%lats)) deallocate(smb%par%lats)
        allocate(smb%par%x(nx),smb%par%y(ny),smb%par%lats(nx,ny))

        smb%par%x    = x 
        smb%par%y    = y 
        smb%par%lats = lats 

        ! Allocate the smbpal object 
        call smbpal_allocate(smb%now,nx,ny)
        call smbpal_allocate(smb%ann,nx,ny)
    
        do m = 1, 12 
            call smbpal_allocate(smb%mon(m),nx,ny)
        end do 

        return 

    end subroutine smbpal_init

    subroutine smbpal_update_2temp(smb,t2m_ann,t2m_sum,pr_ann, &
                                   z_srf,H_ice,time_bp,sf_ann,file_out,file_out_mon,file_out_day)
        ! Generate climate using two points in year (Tsum,Tann)

        implicit none 
        
        type(smbpal_class), intent(INOUT) :: smb
        real(prec), intent(IN) :: t2m_ann(:,:), t2m_sum(:,:)
        real(prec), intent(IN) ::  pr_ann(:,:), z_srf(:,:), H_ice(:,:)
        real(prec), intent(IN) :: time_bp       ! years BP 
        real(prec), intent(IN), optional :: sf_ann(:,:)
        character(len=*), intent(IN), optional :: file_out      ! Annual output
        character(len=*), intent(IN), optional :: file_out_mon  ! Monthly output
        character(len=*), intent(IN), optional :: file_out_day  ! Daily output 

        ! Local variables
        integer, parameter :: ndays = 360   ! 360-day year
        integer, parameter :: ndays_mon = 30   ! 30 days per month  
        integer :: day, m, nx, ny, mnow, mday  
        real(prec) :: PDDs(size(t2m_ann,1),size(t2m_ann,2))

        type(smbpal_param_class) :: par
        type(smbpal_state_class) :: now

        ! Fill in local versions for easier access 
        par = smb%par 
        now = smb%now 

        if (present(file_out_day)) then
            call smbpal_write_init(par,file_out_day,z_srf,H_ice)
        end if 

        ! First calculate PDDs for the whole year (input to itm)
        PDDs = 0.0 
        do day = 1, ndays
            now%t2m = t2m_ann-(t2m_sum-t2m_ann)*cos(2.0*pi*real(day-15)/real(ndays))
            PDDs    = PDDs + calc_temp_effective(now%t2m,par%Teff_sigma)
        end do 

        ! Initialize averaging 
        call smbpal_average(smb%ann,now,step="init")
        do m = 1, 12 
            call smbpal_average(smb%mon(m),now,step="init")
        end do 

        mnow = 1 
        mday = 1 

        do day = 1, ndays

            ! Determine t2m, teff, pr, sf and S today 
            now%t2m  = t2m_ann-(t2m_sum-t2m_ann)*cos(2.0*pi*real(day-15)/real(ndays))
            now%teff = calc_temp_effective(now%t2m,par%Teff_sigma)
            now%pr = pr_ann / real(ndays)

            if (present(sf_ann)) then 
                now%sf = sf_ann / real(ndays)
            else 
                now%sf = calc_snowfrac(now%t2m,par%sf_a,par%sf_b)
            end if 

            now%S = calc_insol_day(day,dble(par%lats),dble(time_bp),fldr="libs/insol/input")

            ! Call mass budget for today
            call calc_snowpack_budget_day(par%itm,z_srf,H_ice,now%S,now%t2m,PDDs, &
                                          now%pr,now%sf,now%H_snow,now%alb_s,now%smbi, &
                                          now%smb,now%melt,now%runoff,now%refrz)
        

            ! Get averages 
            call smbpal_average(smb%ann,now,step="step")
            call smbpal_average(smb%mon(mnow),now,step="step")

            mday = mday + 1 
            if (mday .eq. ndays_mon) then 
                call smbpal_average(smb%mon(mnow),now,step="end",nt=real(ndays_mon))
                mnow = mnow + 1
                mday = 1 
            end if 

            if (present(file_out_day)) then 
                ! Write daily output for this year 
                call smbpal_write(now,file_out,ndat=day,time_bp=time_bp,step="day")
            end if 
    
        end do 

        ! Finalize annual average 
        call smbpal_average(smb%ann,now,step="end",nt=real(ndays))

        ! Repopulate global now variable (in case it is needed)
        smb%now = now 

        ! Annual I/O 
        if (present(file_out)) then
            call smbpal_write_init(par,file_out,z_srf,H_ice)
            call smbpal_write(smb%ann,file_out,ndat=1,time_bp=time_bp,step="ann")

        end if 


        return 

    end subroutine smbpal_update_2temp

    subroutine smbpal_update_monthly(par,now,t2m,pr,time_bp)
        ! Generate climate using monthly input data [nx,ny,nmon]
        
        implicit none 
        
        type(smbpal_param_class), intent(IN)    :: par
        type(smbpal_state_class), intent(INOUT) :: now
        real(prec), intent(IN) :: t2m(:,:,:), pr(:,:,:)
        real(prec), intent(IN) :: time_bp       ! years BP 

        ! Local variables
        real(prec), allocatable :: t2m_now(:,:), pr_now(:,:)

        allocate(t2m_now(size(t2m,1),size(t2m,2)))
        allocate(pr_now(size(pr,1),size(pr,2)))


        return 

    end subroutine smbpal_update_monthly



    subroutine smbpal_update_daily(par,now,t2m,pr,time_bp)
        ! Generate climate using input data for each day

        implicit none 
        
        type(smbpal_param_class), intent(IN)    :: par
        type(smbpal_state_class), intent(INOUT) :: now
        real(prec), intent(IN) :: t2m(:,:), pr(:,:)
        real(prec), intent(IN) :: time_bp       ! years BP  

        ! ** TO DO ** 
        
        return 

    end subroutine smbpal_update_daily


    subroutine smbpal_end(smbpal)

        implicit none 

        type(smbpal_class) :: smbpal 

        ! Deallocate smbpal state object
        call smbpal_deallocate(smbpal%now)
	
        return 

    end subroutine smbpal_end


    subroutine smbpal_par_load(par,filename)

        type(smbpal_param_class)     :: par
        character(len=*), intent(IN) :: filename 

        ! Local parameter definitions (identical to object)
        character(len=16) :: abl_method
        real(prec)        :: Teff_sigma, sf_a, sf_b

        namelist /smbpal_par/ abl_method, Teff_sigma, sf_a, sf_b
                
        ! Store initial values in local parameter values 
        abl_method = par%abl_method
        Teff_sigma = par%Teff_sigma 
        sf_a       = par%sf_a 
        sf_b       = par%sf_b 

        ! Read parameters from input namelist file
        open(7,file=trim(filename))
        read(7,nml=smbpal_par)
        close(7)

        ! Store local parameter values in output object
        par%abl_method = abl_method
        par%Teff_sigma = Teff_sigma 
        par%sf_a       = sf_a 
        par%sf_b       = sf_b 

        ! Also load itm parameters
        call itm_par_load(par%itm,filename)

        return

    end subroutine smbpal_par_load

   
    ! =======================================================
    !
    ! smb physics
    !
    ! =======================================================

    elemental subroutine calc_ablation_pdd(abl,sif,pdds,acc,csnow,csi,cice)
        ! Determine total ablation based on input pdds
        
        implicit none 

        real(4), intent(INOUT) :: abl, sif 
        real(4), intent(IN)    :: pdds, acc
        real(4), intent(IN)    :: csnow, csi, cice 

        real(4) :: simax
        real(4) :: abl_snow_pot, abl_ice_pot 

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

    elemental function calc_temp_surf(tann,sif) result(ts)
        ! Surface temperature is equal to the annual mean
        ! near-surface temperature + warming due to 
        ! freezing of superimposed ice
        implicit none 

        real(4), intent(IN) :: tann, sif 
        real(4) :: ts 

        ts = (tann+26.6*sif)
        ts = min(0.0,ts)

        return 

    end function calc_temp_surf
    
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
    elemental function calc_temp_effective(temp, sigma) result(teff)

        implicit none

        real(4), intent(IN) :: temp, sigma
        real(4) :: teff
        
        real(4) :: inv_sigma
        real(4), parameter :: inv_sqrt2   = 1.0/sqrt(2.0)
        real(4), parameter :: inv_sqrt2pi = 1.0/sqrt(2.0*pi)

        inv_sigma   = 1.0/sigma

        teff = sigma*inv_sqrt2pi*exp(-0.5*(temp*inv_sigma)**2)  &
                  + temp*0.5*erfcc(-temp*inv_sigma*inv_sqrt2)

        ! Result is the assumed/felt/effective positive degrees, 
        ! given the actual temperature (accounting for fluctuations in day/month/etc, 
        ! based on the sigma chosen)
                     
        return

    end function calc_temp_effective 

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Function :  e r f c c
    ! Author   :  Reinhard Calov and Ralf Greve
    ! Purpose  :  Returns the complementary error function erfc(x) with 
    !             fractional error everywhere less than 1.2 x 10^(-7).
    !             Credit: Press et al., 'Numerical recipes in Fortran 77'.
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    elemental function erfcc(x)

        implicit none
            
        real(4), intent(IN) :: x
        real(4) :: erfcc

        real(4) :: t, z

        z = abs(x)
        t = 1.0/(1.0+0.5*z)

        erfcc = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+  &
        t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*  &
        (-1.13520398+t*(1.48851587+t*(-0.82215223+  &
        t*0.17087277)))))))))

        if (x .lt. 0.0) erfcc = 2.0-erfcc

        return
      
    end function erfcc

    elemental function calc_snowfrac(t2m,a,b) result(f)
        ! Return the fraction of snow from total precipitation
        ! expected for a given temperature
        
        implicit none 

        real(4), intent(IN) :: t2m, a, b 
        real(4)             :: f 

        f = -0.5*tanh(a*(t2m-b))+0.5 

        return 

    end function calc_snowfrac 

    subroutine smbpal_average(ave,now,step,nt)
        implicit none 

        type(smbpal_state_class), intent(INOUT) :: ave
        type(smbpal_state_class), intent(IN)    :: now 
        character(len=*)  :: step
        real(prec), optional :: nt 
        
        call field_average(ave%t2m,    now%t2m,    step,nt)
        call field_average(ave%teff,   now%teff,   step,nt)
        call field_average(ave%pr,     now%pr,     step,nt)
        call field_average(ave%sf,     now%sf,     step,nt)
        call field_average(ave%S,      now%S,      step,nt)
        
        call field_average(ave%H_snow, now%H_snow, step,nt)
        call field_average(ave%alb_s,  now%alb_s,  step,nt)
        call field_average(ave%smbi,   now%smbi,   step,nt)
        call field_average(ave%smb,    now%smb,    step,nt)
        call field_average(ave%melt,   now%melt,   step,nt)
        call field_average(ave%runoff, now%runoff, step,nt)
        call field_average(ave%refrz,  now%refrz,  step,nt)
        
        return

    end subroutine smbpal_average

    subroutine field_average(ave,now,step,nt)
        ! Generic routine to average a field through time 

        implicit none 
        real(prec), intent(INOUT)    :: ave(:,:)
        real(prec), intent(IN)       :: now(:,:)
        character(len=*), intent(IN) :: step
        real(prec), intent(IN), optional :: nt 

        if (trim(step) .eq. "init") then
            ! Initialize field to zero  
            ave = 0.0 
        else if (trim(step) .eq. "step") then 
            ! Sum intermediate steps
            ave = ave + now 
        else if (trim(step) .eq. "end") then
            if (.not.  present(nt)) then 
                write(*,*) "Averaging step total not provided."
                stop 
            end if 
            ! Divide by total steps
            ave = ave / nt 
        else
            write(*,*) "Step not recognized: ",trim(step)
            stop 
        end if 

        return 

    end subroutine field_average 

    ! =======================================================
    !
    ! smbpal I/O
    !
    ! =======================================================

    subroutine smbpal_write_init(par,filename,z_srf,H_ice)

        implicit none 

        type(smbpal_param_class), intent(IN) :: par 
        character(len=*),         intent(IN) :: filename 
        real(prec), intent(IN), optional :: z_srf(:,:), H_ice(:,:) 

        call nc_create(filename)
        call nc_write_dim(filename,"xc",x=par%x)
        call nc_write_dim(filename,"yc",x=par%y)
        call nc_write_dim(filename,"day",  x=1,nx=360,dx=1)
        call nc_write_dim(filename,"month",x=1,nx=12,dx=1)
        call nc_write_dim(filename,"time",x=0.0,units="ka BP",unlimited=.TRUE.)
        
        ! Write the 2D latitude field to file
        call nc_write(filename,"lat2D",par%lats,dim1="xc",dim2="yc")

        if (present(z_srf)) call nc_write(filename,"z_srf",z_srf,dim1="xc",dim2="yc")
        if (present(H_ice)) call nc_write(filename,"H_ice",H_ice,dim1="xc",dim2="yc")

        return 

    end subroutine smbpal_write_init 

    subroutine smbpal_write(now,filename,ndat,time_bp,step)

        implicit none 

        type(smbpal_state_class), intent(IN) :: now 
        character(len=*),         intent(IN) :: filename 
        integer,                  intent(IN) :: ndat 
        real(prec),               intent(IN) :: time_bp  
        character(len=*),         intent(IN) :: step  

        ! Local variables 
        real(prec) :: ka_bp 
        integer :: dim3_val, nx, ny  

        ka_bp = time_bp * 1e-3 

        if (trim(step) .ne. "ann" .and. trim(step) .ne. "mon" .and. trim(step) .ne. "day") then 
            write(*,*) "smbpal_write:: error: step should be one of: ann, mon or day."
            stop 
        end if 

        nx = size(now%t2m,1)
        ny = size(now%t2m,2)

            
        ! Write the variables
        if (trim(step) .eq. "ann") then 
            ! Write the annual mean with time as 3rd dimension 

            ! Update the timestep 
            call nc_write(filename,"time",ka_bp,dim1="time",start=[ndat],count=[1])

            call nc_write(filename,"t2m",now%t2m,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1])
            call nc_write(filename,"S",now%S,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1])

        else 
            ! Write the step along the year (mon or day)

            ! Update the timestep 
            call nc_write(filename,"time",ka_bp,dim1="time",start=[1],count=[1])

            call nc_write(filename,"t2m",now%t2m,dim1="xc",dim2="yc",dim3=trim(step), &
                          start=[1,1,ndat],count=[nx,ny,1])
            call nc_write(filename,"S",now%S,dim1="xc",dim2="yc",dim3=trim(step), &
                          start=[1,1,ndat],count=[nx,ny,1])

        end if 

        return 

    end subroutine smbpal_write 


    ! =======================================================
    !
    ! smbpal memory management
    !
    ! =======================================================

    subroutine smbpal_allocate(now,nx,ny)

        implicit none 

        type(smbpal_state_class) :: now 
        integer :: nx, ny 

        ! Make object is deallocated
        call smbpal_deallocate(now)

        ! Allocate variables
        allocate(now%t2m(nx,ny))
        allocate(now%pr(nx,ny))
        allocate(now%sf(nx,ny))
        allocate(now%S(nx,ny))
        allocate(now%teff(nx,ny))
        allocate(now%H_snow(nx,ny))
        allocate(now%alb_s(nx,ny))
        allocate(now%smbi(nx,ny))
        allocate(now%smb(nx,ny))
        allocate(now%melt(nx,ny))
        allocate(now%runoff(nx,ny))
        allocate(now%refrz(nx,ny))

        return

    end subroutine smbpal_allocate

    subroutine smbpal_deallocate(now)

        implicit none 

        type(smbpal_state_class) :: now 

        ! Allocate state objects
        if (allocated(now%t2m))      deallocate(now%t2m)
        if (allocated(now%pr))       deallocate(now%pr)
        if (allocated(now%sf))       deallocate(now%sf)
        if (allocated(now%S))        deallocate(now%S)
        if (allocated(now%teff))     deallocate(now%teff)
        if (allocated(now%H_snow))   deallocate(now%H_snow)
        if (allocated(now%alb_s))    deallocate(now%alb_s)
        if (allocated(now%smbi))     deallocate(now%smbi)
        if (allocated(now%smb))      deallocate(now%smb)
        if (allocated(now%melt))     deallocate(now%melt)
        if (allocated(now%runoff))   deallocate(now%runoff)
        if (allocated(now%refrz))    deallocate(now%refrz)
        
        return

    end subroutine smbpal_deallocate

end module smbpal


