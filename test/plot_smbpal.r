library(RNetCDF)

rmse = function (x, ii = c(1:length(x)), round = NA) 
{
    x <- x[ii]
    ii1 <- which(!is.na(x))
    rmse <- sqrt(sum(x[ii1]^2)/length(x[ii1]))
    if (!is.na(round)) 
        rmse <- round(rmse, round)
    return(rmse)
}

# Load data
if (TRUE) {

    a = open.nc("data/GRL-20KM_TOPO-B13_gl0.05.nc")
    topo = read.nc(a)
    close.nc(a)

    a = open.nc("data/GRL-20KM_REGIONS.nc")
    regs = read.nc(a)
    close.nc(a)

    topo$zs[which(regs$mask < 3.09 | regs$mask > 3.21)] = 0.0 


    a = open.nc("data/GRL-20KM_MARv3.5-ERA-30km-monthly_1981-2010.nc")
    mar = read.nc(a)
    close.nc(a)
    mar$t3m     = mar$t3m     + 273.15 
    mar$t3m_min = mar$t3m_min + 273.15 
    mar$t3m_max = mar$t3m_max + 273.15 

    # Calculate some quantities 
    marann = mar 
    for (q in 1:length(marann)) {
        if (length(dim(mar[[q]]))==3) {
            marann[[q]] = apply(mar[[q]],c(1,2),mean)
        }
    }
    marann$pdd = marann$pdd*365 
}

a = open.nc("output/smbpal_eraint_pdd.nc")
pdd = read.nc(a)
close.nc(a)

a = open.nc("output/smbpal_eraint_itm.nc")
itm = read.nc(a)
close.nc(a)

# Get errors of interest 
tbl = data.frame(name=c("t2m","tsrf","PDDs","smb","melt","refrz","runoff"),pdd=NA,itm=NA,mar=NA,
                 err_pdd=NA,err_itm=NA,
                 nm_mar=c("t3m","ts","pdd","smb","me","rz","ru")) 
for (q in 1:dim(tbl)[1]) {
    tbl$err_pdd[q] = rmse(marann[[tbl$nm_mar[q]]]-pdd[[tbl$name[q]]],ii=which(topo$zs>0),round=3)
    tbl$err_itm[q] = rmse(marann[[tbl$nm_mar[q]]]-itm[[tbl$name[q]]],ii=which(topo$zs>0),round=3)
}

# image.plot(smb$smb-marann$smb,zlim=c(-1,1))
# image.plot(smb$melt-marann$me,zlim=c(-5,5))
# image.plot(smb$alb_s-marann$al,zlim=c(-0.2,0.2))


