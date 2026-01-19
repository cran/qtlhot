## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 5,
                      collapse = TRUE, comment = "#>")

## ----loadlib------------------------------------------------------------------
library(qtlhot)

## -----------------------------------------------------------------------------
ncross1 <- sim.null.cross(chr.len = rep(100, 4),
                          n.mar = 51,
                          n.ind = 100,
                          type = "bc",
                          n.phe = 1000,
                          latent.eff = 3,
                          res.var = 1,
                          init.seed = 123457)

## -----------------------------------------------------------------------------
cross1 <- include.hotspots(cross = ncross1,
                           hchr = c(2, 3, 4),
                           hpos = c(25, 75, 50),
                           hsize = c(100, 50, 20),
                           Q.eff = 2,
                           latent.eff = 3,
                           lod.range.1 = c(2.5, 2.5),
                           lod.range.2 = c(5, 8),
                           lod.range.3 = c(10, 15),
                           res.var = 1,
                           n.phe = 1000,
                           init.seed = 12345)

## -----------------------------------------------------------------------------
ncor1 <- cor(cross1$pheno)
summary(ncor1[lower.tri(ncor1)])
rm(ncor1)

## ----echo=FALSE---------------------------------------------------------------
if(file.exists("savedperms.RData")) load("savedperms.RData")

## ----echo=FALSE---------------------------------------------------------------
if(!exists("pt.scanone")) {
  set.seed(123)
  pt.scanone <- scanone(ncross1, method = "hk", n.perm = 1000)
}

## -----------------------------------------------------------------------------
alphas <- seq(0.01, 0.10, by=0.01)
# Following works but breaks CRAN build:
lod.thrs <- summary(pt.scanone, alphas)
#lod.thrs <- qtl:::summary.scanoneperm(pt.scanone, alphas)
lod.thrs
lod.thr <- lod.thrs[5]

## -----------------------------------------------------------------------------
scan1 <- qtl::scanone(cross1, pheno.col = 1:1000, method = "hk")

## -----------------------------------------------------------------------------
high1 <- highlod(scan1, lod.thr = min(lod.thrs), drop.lod = 1.5)
max(high1, lod.thr = lod.thrs)

## -----------------------------------------------------------------------------
hots1 <- hotsize(high1, lod.thr = lod.thr)
summary(hots1)

## ----label=plotex1counts,width=6,height=3.5,fig=TRUE,echo=FALSE---------------
par(mar=c(4.1,4.1,0.5,0.1))
plot(hots1, cex.lab = 1.5, cex.axis = 1.5)

## ----echo=FALSE---------------------------------------------------------------
if(!exists("hotperm1")) {
  set.seed(12345)
  hotperm1 <- hotperm(cross = cross1,
                  n.quant = 300,
                  n.perm = 100,
                  lod.thrs = lod.thrs,
                  alpha.levels = alphas,
                  drop.lod = 1.5,
                  verbose = FALSE)
}

## -----------------------------------------------------------------------------
names(hotperm1)
summary(hotperm1)

## ----label=figex1slidingbar,include=FALSE,echo=FALSE--------------------------
par(mar=c(4.1,4.1,0.5,4.1))
quant1 <- quantile_hotperm(hotperm1, 0.05, lod.thr = lod.thr)
plot(high1, quant.level = quant1, sliding = TRUE)

## ----label=plotex1slidingbar,fig=TRUE,width=6,height=3.5,echo=FALSE-----------
par(mar=c(4.1,4.1,0.5,4.1))
quant1 <- quantile_hotperm(hotperm1, 0.05, lod.thr = lod.thr)
plot(high1, quant.level = quant1, sliding = TRUE)

## ----echo=FALSE---------------------------------------------------------------
hotsq1 <- hotsize(high1, lod = lod.thr, window = 5, quant.level = quant1)
quant.axis <- pmax(1, pretty(c(0, qtl:::summary.scanone(hotsq1)[3,"quant"])))
quant.level <- round(attr(hotsq1, "quant.level")[quant.axis], 2)

## ----label=figex1elias,include=FALSE,echo=FALSE-------------------------------
par(mar=c(4.1,4.1,0.5,4.1))
plot(hotsq1)

## ----label=plotex1elias,fig=TRUE,width=6,height=3.5,echo=FALSE----------------
par(mar=c(4.1,4.1,0.5,4.1))
plot(hotsq1)

## -----------------------------------------------------------------------------
summary(hotsq1)

## ----echo=FALSE---------------------------------------------------------------
ncross2 <- sim.null.cross(chr.len = rep(100,4), 
                          n.mar = 51, 
                          n.ind = 100,
                          type = "bc", 
                          n.phe = 1000, 
                          latent.eff = 0, 
                          res.var = 1, 
                          init.seed = 123457)
cross2 <- include.hotspots(cross = ncross2,
                           hchr = c(2, 3, 4),
                           hpos = c(25, 75, 50),
                           hsize = c(100, 50, 20),
                           Q.eff = 2,
                           latent.eff = 0,
                           lod.range.1 = c(2.5, 2.5),
                           lod.range.2 = c(5, 8),
                           lod.range.3 = c(10, 15),
                           res.var = 1,
                           n.phe = 1000,
                           init.seed = 12345)

## ----label=figex2counts,echo=FALSE,include=FALSE------------------------------
scan2 <- scanone(cross2, pheno.col = 1:1000, method = "hk")
high2 <- highlod(scan2, lod.thr = lod.thr, drop.lod = 1.5)
hots2 <- hotsize(high2)
par(mar=c(4.1,4.1,0.1,0.1))
plot(hots2, cex.lab = 1.5, cex.axis = 1.5)

## ----label=plotex2counts,fig=TRUE,width=6,height=3.5,echo=FALSE---------------
scan2 <- scanone(cross2, pheno.col = 1:1000, method = "hk")
high2 <- highlod(scan2, lod.thr = lod.thr, drop.lod = 1.5)
hots2 <- hotsize(high2)
par(mar=c(4.1,4.1,0.1,0.1))
plot(hots2, cex.lab = 1.5, cex.axis = 1.5)

## -----------------------------------------------------------------------------
ncor2 <- cor(cross2$pheno)
summary(ncor2[lower.tri(ncor2)])
rm(ncor2)

## ----echo=FALSE---------------------------------------------------------------
if(!exists("hotperm2")) {
  set.seed(12345)
  hotperm2 <- hotperm(cross = cross2, 
                  n.quant = 300, 
                  n.perm = 100, 
                  lod.thrs = lod.thrs, 
                  alpha.levels = alphas,
                  drop.lod = 1.5, 
                  verbose = FALSE) 
}

## ----echo=FALSE---------------------------------------------------------------
quant2 <- quantile_hotperm(hotperm2, 0.05, lod.thr = lod.thr)

## ----label=plotex2slidingbar,width=6,height=3.5,fig=TRUE,echo=FALSE,include=FALSE----
par(mar=c(4.1,4.1,0.5,0.1))
plot(high2, lod.thr = lod.thr, quant.level = quant2, sliding = TRUE)

## ----label=plotex2sigct,width=6,height=3.5,fig=TRUE,echo=FALSE,include=FALSE----
par(mar=c(4.1,4.1,0.5,0.1))
plot(high2, quant.level = quant2)

## ----echo=FALSE---------------------------------------------------------------
if(!file.exists("savedperms.RData"))
  save(pt.scanone, hotperm1, hotperm2, file = "savedperms.RData", compress = TRUE)

