program test 
    
    use smbpal 
    use insolation 
    use ncio 

    implicit none 


    type(smbpal_class) :: smb1 
    real(4) :: insol_pt 

    ! Input data information
    character(len=256) :: file_topo, file_clim, file_out  
    integer :: nx, ny 
    real(4), allocatable :: x(:), y(:)
    real(4), allocatable :: lats(:,:), z_srf(:,:), H_ice(:,:)
    real(4), allocatable :: t2m(:,:,:), pr(:,:,:), sf(:,:,:)
    real(4), allocatable :: t2m_ann(:,:), t2m_sum(:,:), pr_ann(:,:), sf_ann(:,:)

    ! Output file 
    file_out = "output/smbpal_eraint_pdd.nc"

    ! Load domain test data (Greenland)
    file_topo = "data/GRL-20KM_TOPO-B13_gl0.05.nc"
    nx = nc_size(file_topo,"xc")
    ny = nc_size(file_topo,"yc")
    allocate(x(nx),y(ny),lats(nx,ny),z_srf(nx,ny),H_ice(nx,ny))
    call nc_read(file_topo,"xc",x)
    call nc_read(file_topo,"yc",y)
    call nc_read(file_topo,"lat2D",lats)
    call nc_read(file_topo,"zs",   z_srf)
    call nc_read(file_topo,"H",    H_ice)
    where(z_srf .lt. 0.0) z_srf = 0.0 

    file_clim = "data/GRL-20KM_MARv3.5-ERA-30km-monthly_1981-2010.nc"
    allocate(t2m(nx,ny,12),pr(nx,ny,12),sf(nx,ny,12))
    call nc_read(file_clim,"rf",t2m) ! Hold rainfall here temporarily
    call nc_read(file_clim,"sf",sf)
    pr = t2m + sf 
    call nc_read(file_clim,"t3m",t2m)
    t2m = t2m + 273.15                  ! [C] => [K]

    ! Derived climate variables for testing
    allocate(t2m_ann(nx,ny),t2m_sum(nx,ny),pr_ann(nx,ny),sf_ann(nx,ny))
    t2m_ann = sum(t2m,dim=3) / 12.0
    t2m_sum = sum(t2m(:,:,6:8),dim=3) / 3.0
    pr_ann  = sum(pr,dim=3) / 12.0
    sf_ann  = sum(sf,dim=3) / 12.0

    ! Check input data 
    write(*,"(5x,a)")      "======== Input data ========"
    write(*,"(a15,2f10.3)") "lats: ",  minval(lats),  maxval(lats)
    write(*,"(a15,2f10.3)") "z_srf: ", minval(z_srf), maxval(z_srf)
    write(*,"(a15,2f10.3)") "H_ice: ", minval(H_ice), maxval(H_ice)
    
    write(*,"(a15,2f10.3)") "t2m_ann: ", minval(t2m_ann), maxval(t2m_ann)
    write(*,"(a15,2f10.3)") "t2m_sum: ", minval(t2m_sum), maxval(t2m_sum)
    write(*,"(a15,2f10.3)") "pr_ann: ",  minval(pr_ann),  maxval(pr_ann)
    write(*,"(a15,2f10.3)") "sf_ann: ",  minval(sf_ann),  maxval(sf_ann)
    
    ! Initialize the smbpal object and output file
    call smbpal_init(smb1,"Greenland.nml",x,y,lats)
    
    call smbpal_update(smb1%par,smb1%now,t2m_ann,t2m_sum,pr_ann, &
                            z_srf,H_ice,time_bp=0.0,file_out=file_out)
    
    ! Finalize the smbpal object 
    call smbpal_end(smb1)


end program test 

