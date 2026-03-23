module smbpal_precision

    implicit none 

    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)
    
    ! Precision used here
    integer,  parameter :: wp = sp 

    ! Missing value aliases 
    real(wp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 
    real(wp), parameter :: MV     = MISSING_VALUE_DEFAULT
    real(wp), parameter :: MV_INT = int(MISSING_VALUE_DEFAULT)

    ! Error values
    real(wp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 
        
    ! Constants
    real(wp), parameter :: pi = 3.14159265359

end module smbpal_precision 