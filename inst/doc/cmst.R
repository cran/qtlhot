## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 5,
                      collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(qtlhot)

## -----------------------------------------------------------------------------
set.seed(987654321)
CMSTCross <- SimCrossCausal(n.ind = 100,
                        len = rep(100, 3),
                        n.mar = 101,
                        beta = rep(0.5, 2),
                        add.eff = 1,
                        dom.eff = 0,
                        sig2.1 = 0.4,
                        sig2.2 = 0.1,
                        eq.spacing = FALSE,
                        cross.type = "bc",
                        normalize = TRUE)

## -----------------------------------------------------------------------------
CMSTCross <- calc.genoprob(CMSTCross, step = 1)

## -----------------------------------------------------------------------------
Scan <- scanone(CMSTCross, pheno.col = 1 : 3, method = "hk")
summary(Scan[, c(1, 2, 3)], thr = 3)
summary(Scan[, c(1, 2, 4)], thr = 3)
summary(Scan[, c(1, 2, 5)], thr = 3)

## ----label=lodprofiles,width=6,height=6,fig=TRUE------------------------------
plot(Scan, lodcolumn = 1 : 3, ylab = "LOD")

## -----------------------------------------------------------------------------
commqtls <- GetCommonQtls(CMSTCross,
                          pheno1 = "y1",
                          pheno2 = "y3",
                          thr = 3,
                          peak.dist = 5,
                          addcov1 = NULL,
                          addcov2 = NULL,
                          intcov1 = NULL,
                          intcov2 = NULL)
commqtls

## -----------------------------------------------------------------------------
nms <- names(CMSTCross$pheno)
out1 <- CMSTtests(CMSTCross,
                  pheno1 = nms[1],
                  pheno2 = nms[2],
                  Q.chr = 1,
                  Q.pos = 55,
                  addcov1 = NULL,
                  addcov2 = NULL,
                  intcov1 = NULL,
                  intcov2 = NULL,
                  method = "all",
                  penalty = "both")

## -----------------------------------------------------------------------------
out1[1:3]

## -----------------------------------------------------------------------------
out1[4]

## -----------------------------------------------------------------------------
out1[5]

## -----------------------------------------------------------------------------
out1[6]

## -----------------------------------------------------------------------------
out1[7]

## -----------------------------------------------------------------------------
out1[8]

## -----------------------------------------------------------------------------
out1[9]

## -----------------------------------------------------------------------------
out1[10:12]

## -----------------------------------------------------------------------------
out1[13:17]

## -----------------------------------------------------------------------------
out2 <- CMSTtests(CMSTCross,
                  pheno1 = nms[1],
                  pheno2 = nms[-1],
                  Q.chr = 1,
                  Q.pos = 55.5,
                  addcov1 = NULL,
                  addcov2 = NULL,
                  intcov1 = NULL,
                  intcov2 = NULL,
                  method = "all",
                  penalty = "both")
out2

## -----------------------------------------------------------------------------
CMSTscan <- scanone(CMSTCross, pheno.col = 1:3, method = "hk")
CMSThigh <- highlod(CMSTscan)

## -----------------------------------------------------------------------------
traits <- names(CMSTCross$pheno)
annot <- data.frame(name = traits, traits = traits, chr = rep(1, 3),
 Mb.pos = c(55,10,100))
annot$cM.pos <- annot$Mb.pos
annot
targets <- list(y1 = c("y2","y3"))

## -----------------------------------------------------------------------------
cand.reg <- GetCandReg(CMSThigh, annot, traits)
cand.reg
cis.cand.reg <- GetCisCandReg(CMSThigh, cand.reg)
cis.cand.reg
comap.targets <- GetCoMappingTraits(CMSThigh, cand.reg)
comap.targets

## -----------------------------------------------------------------------------
tests <- list()
for(k in seq(names(comap.targets))) {
  tests[[k]] <- FitAllTests(CMSTCross, pheno1 = names(comap.targets)[k],
                      pheno2 = comap.targets[[k]],
                      Q.chr = cand.reg[k, 4],
                      Q.pos = cand.reg[k, 5])
}
names(tests) <- names(comap.targets)
tests <- JoinTestOutputs(comap.targets, tests)
tests

## -----------------------------------------------------------------------------
PrecTpFpMatrix(alpha = seq(0.01, 0.10, by = 0.01),
  val.targets = targets, all.orfs = CMSThigh$names, tests = tests,
  cand.reg = cand.reg, cis.cand.reg = cis.cand.reg)

