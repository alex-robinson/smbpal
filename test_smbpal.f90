program test 
    
    use smbpal 
    use insolation 

    implicit none 


    type(smbpal_class) :: smb1 
    real(4) :: insol_pt 

    ! Initialize the smbpal object 
    call smbpal_init(smb1,"Greenland.nml",nx=10,ny=10)

    ! Get the present-day insolation for the domain 
    insol_pt = calc_insol_day(170,65.d0,0.d0,fldr="libs/insol/input")
    write(*,*) "insol = ", insol_pt 



    ! Finalize the smbpal object 
    call smbpal_end(smb1)


end program test 

