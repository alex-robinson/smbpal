 module smbpal

    use smbpal_precision
    use insolation
    use smb_itm 

    implicit none 

    type smbpal_param_class
        type(itm_par) :: itm
        logical :: use_subgrid 
        real(prec) :: par1 

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
        type(smbpal_state_class) :: now, mon, ann
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

contains 

    subroutine smbpal_init(smb,filename,nx,ny)

        implicit none 

        type(smbpal_class) :: smb
        character(len=*), intent(IN)  :: filename  ! Parameter file 

        ! Local variables
        integer :: nx, ny 
        
        ! Load smbpal parameters
        call smbpal_par_load(smb%par,filename)

        ! Allocate the smbpal object 
        call smbpal_allocate(smb%now,nx,ny)
 
        return 

    end subroutine smbpal_init

    subroutine smbpal_update_2temp(par,now,t2m_ann,t2m_sum,pr_ann, &
                                   lats,z_srf,H_ice,time_bp)
        ! Generate climate using two points in year (Tsum,Tann)

        implicit none 
        
        type(smbpal_param_class), intent(IN)    :: par
        type(smbpal_state_class), intent(INOUT) :: now
        real(prec), intent(IN) :: t2m_ann(:,:), t2m_sum(:,:)
        real(prec), intent(IN) ::  pr_ann(:,:), lats(:,:), z_srf(:,:), H_ice(:,:)
        real(prec), intent(IN) :: time_bp       ! years BP 

        ! Local variables
        integer, parameter :: ndays = 360   ! 360-day year 
        integer :: day 

        do day = 1, ndays 

            ! Determine t2m, pr, S and PDDs of today 
!             now%t2m = cos(2*pi*day/real(ndays))
            
            now%pr = pr_ann / real(ndays)   ! Evenly divided precip over the year (account for temp?)

            now%S = calc_insol_day(day,dble(lats),dble(time_bp),fldr="libs/insol/input")

            ! Call mass budget for today
            call snowpack_budget(par%itm,z_srf,H_ice,now%S,now%t2m,now%teff,now%pr,now%sf, &
                                 now%H_snow,now%alb_s,now%smbi,now%smb, &
                                 now%melt,now%runoff,now%refrz)
        end do 


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
        logical :: use_subgrid
        real(prec) :: par1

        namelist /smbpal_par/ use_subgrid, par1
                
        ! Store initial values in local parameter values 
        use_subgrid = par%use_subgrid
        par1        = par%par1 

        ! Read parameters from input namelist file
        open(7,file=trim(filename))
        read(7,nml=smbpal_par)
        close(7)

        ! Store local parameter values in output object
        par%use_subgrid = use_subgrid
        par%par1        = par1 

        return

    end subroutine smbpal_par_load

   
    ! =======================================================
    !
    ! smb physics
    !
    ! =======================================================



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


