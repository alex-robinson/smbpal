library(RNetCDF)
source("functions.r")

# Load data
if (FALSE) {

    a = open.nc("../data/GRL-20KM_TOPO-B13_gl0.05.nc")
    topo = read.nc(a)
    close.nc(a)

    a = open.nc("../data/GRL-20KM_REGIONS.nc")
    regs = read.nc(a)
    close.nc(a)

    topo$zs[which(regs$mask < 3.09 | regs$mask > 3.21)] = 0.0 

    dx = diff(topo$xc)[1]
    dy = diff(topo$yc)[1]

    a = open.nc("../data/GRL-20KM_MARv3.5-ERA-30km-monthly_1981-2010.nc")
    mar = read.nc(a)
    close.nc(a)
    mar$t3m     = mar$t3m     + 273.15 
    mar$t3m_min = mar$t3m_min + 273.15 
    mar$t3m_max = mar$t3m_max + 273.15 
    mar$ts      = mar$ts      + 273.15 

    # Calculate some quantities 
    marann = mar 
    for (q in 1:length(marann)) {
        if (length(dim(mar[[q]]))==3) {
            marann[[q]] = apply(mar[[q]],c(1,2),mean)
        }
    }
    marann$pdd = marann$pdd*365 
}

a = open.nc("../output/smbpal_eraint_pdd.nc")
pdd = read.nc(a)
close.nc(a)

a = open.nc("../output/smbpal_eraint_itm.nc")
itm = read.nc(a)
close.nc(a)

# Get errors of interest 
mask = topo$zs*NA 
mask[topo$zs>0] = 1.0 
mask[topo$zs>0 & topo$H >0] = 1.0 
nmask = sum(mask,na.rm=TRUE)

tbl = data.frame(name=c("t2m","tsrf","PDDs","smb","melt","refrz","runoff"),pdd=NA,itm=NA,mar=NA,
                 err_pdd=NA,err_itm=NA,
                 nm_mar=c("t3m","ts","pdd","smb","me","rz","ru")) 
for (q in 1:dim(tbl)[1]) {

    nm     = tbl$name[q]
    nm_mar = tbl$nm_mar[q]
        
    if (nm %in% c("smb","melt","refrz","runoff")) {
        tbl$pdd[q] = sum(pdd[[nm]]*360*topo$area*mask,na.rm=TRUE)*1e-12
        tbl$itm[q] = sum(itm[[nm]]*360*topo$area*mask,na.rm=TRUE)*1e-12
        tbl$mar[q] = sum(marann[[nm_mar]]*365*topo$area*mask,na.rm=TRUE)*1e-12
    } else {
        tbl$pdd[q] = sum(pdd[[nm]]*mask,na.rm=TRUE)/nmask
        tbl$itm[q] = sum(itm[[nm]]*mask,na.rm=TRUE)/nmask
        tbl$mar[q] = sum(marann[[nm_mar]]*mask,na.rm=TRUE)/nmask
    }

    tbl$err_pdd[q] = rmse(marann[[tbl$nm_mar[q]]]-pdd[[tbl$name[q]]],ii=which(mask>0))
    tbl$err_itm[q] = rmse(marann[[tbl$nm_mar[q]]]-itm[[tbl$name[q]]],ii=which(mask>0))
}

for (k in 1:dim(tbl)[2]) {
    if (names(tbl)[k] %in% c("pdd","itm","mar","err_pdd","err_itm")) tbl[[k]] = round(tbl[[k]],3)
}

ptype = "png"

mask = topo$zs*NA 
mask[topo$zs>0] = 1.0 

myfigure("plots","compare_pdd",asp=1.0,pointsize=12,type=ptype)
plot_comparison(topo,mask,marann$pdd,pdd$PDDs,itm$PDDs,zlim=c(0,500),dzlim=c(-100,100))
graphics.off()

# myfigure("plots","compare_smb",asp=1.0,pointsize=12,type=ptype)
# plot_comparison(topo,mask,marann$smb,pdd$smb,itm$smb,zlim=c(-20,10),dzlim=c(-5,5))
# graphics.off()

# myfigure("plots","compare_melt",asp=1.0,pointsize=12,type=ptype)
# plot_comparison(topo,mask,marann$me,pdd$melt,itm$melt,zlim=c(0,30),dzlim=c(-5,5))
# graphics.off()

# myfigure("plots","compare_alb",asp=1.0,pointsize=12,type=ptype)
# plot_comparison(topo,mask,marann$al,pdd$alb_s,itm$alb_s,zlim=c(0,1),dzlim=c(-0.2,0.2))
# graphics.off()

# myfigure("plots","compare_refrz",asp=1.0,pointsize=12,type=ptype)
# plot_comparison(topo,mask,marann$rz,pdd$refrz,itm$refrz,zlim=c(0,5),dzlim=c(-1,1))
# graphics.off()

# myfigure("plots","compare_tsrf",asp=1.0,pointsize=12,type=ptype)
# plot_comparison(topo,mask,marann$ts,pdd$tsrf,itm$tsrf,zlim=c(240,280),dzlim=c(-5,5))
# graphics.off()


