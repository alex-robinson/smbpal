 module smbpal

    use smbpal_precision
    
    implicit none 

    type smbpal_param_class
        logical :: use_subgrid 
        real (prec) :: par1 

        real (prec) :: rho_sw
        real (prec) :: rho_ice

    end type 

    type smbpal_state_class 
        integer, allocatable   :: mask(:,:)           ! Ocean-land-ice mask
        real (prec), allocatable   :: mask_ice(:,:)       ! Ice cover fraction (0.0-1.0)

        real (prec), allocatable   :: z_srf(:,:)          ! Surface elevation [m]
        real (prec), allocatable   :: t2m(:,:)            ! Surface temperature [K]
        real (prec), allocatable   :: S(:,:)              ! Insolation [W/m2]
        
        ! Prognostic variables
        real (prec), allocatable   :: H(:,:)              ! Snow thickness
        real (prec), allocatable   :: smb(:,:)           ! Surface mass balance (mm/a)

    end type 

    type smbpal_class
        type(smbpal_param_class) :: par 
        type(smbpal_state_class) :: now 
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

    subroutine smbpal_update_2temp(par,now)
        ! Generate climate using two points in year (Tsum,Twin)

        implicit none 
        
        type(smbpal_param_class), intent(IN)    :: par
        type(smbpal_state_class), intent(INOUT) :: now

        call smbpal_update_daily(par,now)
        
        return 

    end subroutine smbpal_update_2temp

    subroutine smbpal_update_monthly(par,now)
        ! Generate climate using monthly input data
        
        implicit none 
        
        type(smbpal_param_class), intent(IN)    :: par
        type(smbpal_state_class), intent(INOUT) :: now

        call smbpal_update_daily(par,now)

        return 

    end subroutine smbpal_update_monthly



    subroutine smbpal_update_daily(par,now)
        ! Generate climate using daily input data (default)

        implicit none 
        
        type(smbpal_param_class), intent(IN)    :: par
        type(smbpal_state_class), intent(INOUT) :: now

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
        real (prec) :: par1

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
        allocate(now%mask(nx,ny))
        allocate(now%mask_ice(nx,ny))
        allocate(now%z_srf(nx,ny))
        allocate(now%t2m(nx,ny))
        allocate(now%S(nx,ny))
        allocate(now%H(nx,ny))
        allocate(now%smb(nx,ny))

        return

    end subroutine smbpal_allocate

    subroutine smbpal_deallocate(now)

        implicit none 

        type(smbpal_state_class) :: now 

        ! Allocate state objects
        if (allocated(now%mask))     deallocate(now%mask)
        if (allocated(now%mask_ice)) deallocate(now%mask_ice)
        if (allocated(now%z_srf))    deallocate(now%z_srf)
        if (allocated(now%S))        deallocate(now%S)
        if (allocated(now%H))        deallocate(now%H)
        if (allocated(now%smb))      deallocate(now%smb)

        return

    end subroutine smbpal_deallocate

end module smbpal


