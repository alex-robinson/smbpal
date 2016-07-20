library(RNetCDF)
library(myr)

rmse = function (x, ii = c(1:length(x)), round = NA) 
{
    x <- x[ii]
    ii1 <- which(!is.na(x))
    rmse <- sqrt(sum(x[ii1]^2)/length(x[ii1]))
    if (!is.na(round)) 
        rmse <- round(rmse, round)
    return(rmse)
}

plot_comparison = function(topo,mask,rcm,pdd,itm,zlim=NULL,dzlim=NULL,
                           xlim=range(topo$xc),ylim=range(topo$yc))
{   
    colpalette  = rev(c('#d7191c','#fdae61','#ffffbf','#abdda4','#2b83ba'))
    dcolpalette = rev(c('#ca0020','#f4a582','#f7f7f7','#92c5de','#0571b0'))

    var0 = rcm
    var1 = pdd
    var2 = itm

    # Limit to mask 
    var0[which(mask!=1 | is.na(mask))] = NA 
    var1[which(mask!=1 | is.na(mask))] = NA 
    var2[which(mask!=1 | is.na(mask))] = NA 

    dvar1 = var1-var0 
    dvar2 = var2-var0 

    # Determine zlims
    if (is.null(zlim))  zlim  = range(var0,var1,var2,na.rm=TRUE)
    if (is.null(dzlim)) dzlim = range(dvar1,dvar2,na.rm=TRUE)

    # Generate breaks based on desired zlims
    breaks  = pretty(zlim,11)
    dbreaks = pretty(dzlim,11)

    # Make sure zlims are consistent with actual breaks
    zlim  = range(breaks)
    dzlim = range(dbreaks)

    col  = colorRampPalette(colpalette)(length(breaks)-1)
    dcol = colorRampPalette(dcolpalette)(length(dbreaks)-1)
    
    var0[var0<zlim[1]] = zlim[1]
    var0[var0>zlim[2]] = zlim[2]
    var1[var1<zlim[1]] = zlim[1]
    var1[var1>zlim[2]] = zlim[2]
    var2[var2<zlim[1]] = zlim[1]
    var2[var2>zlim[2]] = zlim[2]

    par(plt=c(0.02,0.3,0.52,0.98))
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    image(topo$xc,topo$yc,var0,add=TRUE,breaks=breaks,col=col)
    par(plt=c(0.32,0.6,0.52,0.98),new=TRUE)
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    image(topo$xc,topo$yc,var1,add=TRUE,breaks=breaks,col=col)
    par(plt=c(0.62,0.9,0.52,0.98),new=TRUE)
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    image(topo$xc,topo$yc,var2,add=TRUE,breaks=breaks,col=col)
    par(plt=c(0.92,0.94,0.65,0.85),new=TRUE)
    mylegend(breaks=breaks,col=col)

    par(plt=c(0.02,0.3,0.02,0.48),new=TRUE)
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    par(plt=c(0.32,0.6,0.02,0.48),new=TRUE)
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    image(topo$xc,topo$yc,dvar1,add=TRUE,breaks=dbreaks,col=dcol)
    par(plt=c(0.62,0.9,0.02,0.48),new=TRUE)
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    image(topo$xc,topo$yc,dvar2,add=TRUE,breaks=dbreaks,col=dcol)
    par(plt=c(0.92,0.94,0.15,0.35),new=TRUE)
    mylegend(breaks=dbreaks,col=dcol)

}

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
# mask[topo$zs>0 & topo$H >0] = 1.0 
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

myfigure("plots","compare_smb",asp=1.0,pointsize=12,type=ptype)
plot_comparison(topo,mask,marann$smb,pdd$smb,itm$smb,zlim=c(-20,10),dzlim=c(-5,5))
graphics.off()

myfigure("plots","compare_melt",asp=1.0,pointsize=12,type=ptype)
plot_comparison(topo,mask,marann$me,pdd$melt,itm$melt,zlim=c(0,30),dzlim=c(-5,5))
graphics.off()

myfigure("plots","compare_alb",asp=1.0,pointsize=12,type=ptype)
plot_comparison(topo,mask,marann$al,pdd$alb_s,itm$alb_s,zlim=c(0,1),dzlim=c(-0.5,0.5))
graphics.off()


