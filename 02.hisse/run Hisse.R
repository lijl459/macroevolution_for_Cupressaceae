
library(ape)
library(geiger)
library(diversitree)
library(hisse)
library(parallel)
library(deSolve)
library(GenSA)
library(subplex)
library(nloptr)

#download HiSSE.null4.9rate.R from Dryad repository https://datadryad.org/stash/dataset/doi:10.5061/dryad.c494p
#load("HiSSE.null4.9rate.R")???????Ի???source??source֮??û??TransMatMaker???????ĳ?TransMatMaker.old



# 设置工作目录
setwd("F:\\实验室\\project\\Cupressaceae\\taxa204\\04.hisse\\all_trait\\24hisse/")
getwd()


source("F:\\实验室\\project\\Cupressaceae\\script\\HiSSE.null4.9rate.R")

# 读取数据

TopA <- read.tree(file="../Cup171.noOUT.regular.Ma.tre")
#orthop.data <- TopA

dat <- read.csv("../code.bin.csv",header = T, row.names = 2)
head(dat)

traits <- colnames(dat)[-1]
traits

trait <- traits[1]
trait

Sting<-cbind(rownames(dat), dat[, trait])
head(Sting)


sampfrac_s <-c(1, 0.95) # sampling frequency for stinger


##Set up all relevant transition matrices##

trans.rates.hisse <- TransMatMaker.old(hidden.states=TRUE)
#exclude dual transitions between both the observed trait and the hidden trait (e.g., q0A<->q1B)
trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3,5,8,10))
trans.rates.hisse.test <- TransMatMaker.old(hidden.states=TRUE)
trans.rates.hisse.test <- ParDrop(trans.rates.hisse.test, c(3,5,8,9,10,12))
trans.rates.hisse.red <- trans.rates.hisse
trans.rates.hisse.red.test <- trans.rates.hisse.test
trans.rates.hisse.red[!is.na(trans.rates.hisse.red) & !trans.rates.hisse.red == 0] = 1
trans.rates.hisse.red.test[!is.na(trans.rates.hisse.red.test) & !trans.rates.hisse.red.test == 0] = 1
##bisse model without hidden states
trans.rates.bisse <- TransMatMaker.old(hidden.states=FALSE)
trans.rates.bisse.red <- trans.rates.bisse
trans.rates.bisse.red[!is.na(trans.rates.bisse.red)] = 1
###hisse model with irreversible states
trans.rates.hisse.irrev <- TransMatMaker.old(hidden.states=TRUE)
trans.rates.hisse.irrev <- ParDrop(trans.rates.hisse.irrev, c(1,3,5,8,9,10))
###hisse model with three trans rates, one for 0 -> 1, one for 1 -> 0, and one for transitions among hidden states####
### the following is straight from the HiSSE Vignette:
trans.rates.nodual.threerates <- trans.rates.hisse.red
# Set all transitions from 0->1 to be governed by a single rate:
to.change <- cbind(c(1,3), c(2,4))
trans.rates.nodual.threerates[to.change] = 1
# Now set all transitions from 1->0 to be governed by a single rate:
to.change <- cbind(c(2,4), c(1,3))
trans.rates.nodual.threerates[to.change] = 2
# Finally, set all transitions between the hidden state to be a single rate (essentially giving 
# you an estimate of the rate by which shifts in diversification occur):
to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.nodual.threerates[to.change] = 3


##Fit the 24 models described in Beaulieu and O'Meara 2016, plus 6 additional models used by Harrington & Reeder (2016)
##The example we use here is topology A and trait=Sting; for running the remaining analyses, change the phy, data and f arguments

hisse.fit1 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.bisse, sann=FALSE)	
hisse.fit2 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse, sann=FALSE)
hisse.fit3 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.bisse.red, sann=FALSE)	
hisse.fit4 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse.red, sann=FALSE)
hisse.fit5 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit6 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit7 <- hisse.null4.old(TopA, Sting, f=sampfrac_s, sann=FALSE)
hisse.fit8 <- hisse.null4.old(TopA, Sting, f=sampfrac_s, eps.anc=rep(1,8), sann=FALSE)
hisse.fit9 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit10 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit11 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,2), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit12 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit13 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse.red.test, sann=FALSE)
hisse.fit14 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test, sann=FALSE)
hisse.fit15 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,2), trans.rate=trans.rates.hisse.red.test, sann=FALSE)
hisse.fit16 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test, sann=FALSE)	
hisse.fit17 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,2,1,3), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit18 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red, sann=FALSE)	
hisse.fit19 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,2,1,3), trans.rate=trans.rates.hisse.red.test, sann=FALSE)
hisse.fit20 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test, sann=FALSE)	
hisse.fit21 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,2,3), trans.rate=trans.rates.hisse.red, sann=FALSE)
hisse.fit22 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red, sann=FALSE)	
hisse.fit23 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,2,3), trans.rate=trans.rates.hisse.red.test, sann=FALSE)
hisse.fit24 <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test, sann=FALSE)	
#hisse_full <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse, sann=FALSE)
#hisse_full_irrev <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse.irrev, sann=FALSE)
#CID2_3rate <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual.threerates, sann=FALSE)
#CID2_8rate <- hisse.old(TopA, Sting, f=sampfrac_s, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.hisse, sann=FALSE)
#CID4_3rate <- hisse.null4.old(TopA, Sting, f=sampfrac_s, trans.type="three.rate", sann=FALSE)
#CID4_9rate <- hisse.null4.mod.9.rates(TopA, Sting, f=sampfrac_s, trans.type="All.no.dual", sann=FALSE)


## =========================
## 3. 汇总模型并保存 AIC 表
## =========================
fits <- list(
  hisse.fit1, hisse.fit2, hisse.fit3, hisse.fit4, hisse.fit5, hisse.fit6,
  hisse.fit7, hisse.fit8, hisse.fit9, hisse.fit10, hisse.fit11, hisse.fit12,
  hisse.fit13, hisse.fit14, hisse.fit15, hisse.fit16, hisse.fit17, hisse.fit18,
  hisse.fit19, hisse.fit20, hisse.fit21, hisse.fit22, hisse.fit23, hisse.fit24
)

modelnames <- c(
  "hisse.fit1","hisse.fit2","hisse.fit3","hisse.fit4","hisse.fit5","hisse.fit6",
  "hisse.fit7","hisse.fit8",
  "hisse.fit9","hisse.fit10","hisse.fit11","hisse.fit12",
  "hisse.fit13","hisse.fit14","hisse.fit15","hisse.fit16",
  "hisse.fit17","hisse.fit18","hisse.fit19","hisse.fit20",
  "hisse.fit21","hisse.fit22","hisse.fit23","hisse.fit24"
)


results <- lapply(fits, function(fit) fit[c("loglik","AIC","AICc")])
best_hebing <- do.call(rbind, results)
rownames(best_hebing) <- modelnames

write.csv(best_hebing, paste0("HiSSE_AIC_info_", trait, ".csv"))

best_index <- which.min(best_hebing[, "AICc"])
best_model_name <- rownames(best_hebing)[best_index]
hisse.best <- fits[[best_index]]

cat("最优模型是：", best_model_name, "\n")
saveRDS(hisse.best, paste0(trait, "_", best_model_name, "_hisse_best_model.rds"))
capture.output(hisse.best, file=paste0(trait, "_", best_model_name, "_hisse_best_summary.txt"))






traits <- colnames(dat)[-1]
traits

trait <- traits[8]
trait

Sting<-cbind(rownames(dat), dat[, trait])
head(Sting)


hisse.best <- hisse.fit22

if (is.null(hisse.best$hidden.states)) {
  hs <- TRUE  # 对 null4 / CID4 就默认 TRUE
} else {
  hs <- hisse.best$hidden.states
}

hisse_best.recon <- MarginRecon.old(TopA, Sting, f = sampfrac_s, 
                                    pars=hisse.best$solution, 
                                    hidden.states=hs, 
                                    AIC=hisse.best$AIC, n.cores=1)


#dir.create("pdf")

pdf(paste0("./pdf/", trait, "_hisse_speciation_recon.pdf"), width = 10, height = 12)
recon_plot <- plot.hisse.states(hisse_best.recon, rate.param = "speciation", 
                                #type = "phylogram",
                                show.tip.label = FALSE, 
                                label.offset = 2, fsize = 0.2, 
                                legend = "tips", legend.cex = 0.6, 
                                edge.width = 8, width.factor=0.4, 
                                rate.colors = c("darkcyan", "goldenrod"), 
                                state.colors = c("grey", "black"))
dev.off()

pdf(paste0("./pdf/", trait, "_hisse_netdiv_recon.pdf"), width = 10, height = 12)
recon_plot <- plot.hisse.states(hisse_best.recon, rate.param = "net.div", 
                                #type = "phylogram",
                                show.tip.label = FALSE, 
                                label.offset = 2, fsize = 0.2, 
                                legend = "tips", legend.cex = 0.6, 
                                edge.width = 8, width.factor=0.4, 
                                rate.colors = c("darkcyan", "goldenrod"), 
                                state.colors = c("grey", "black"))
dev.off()

pdf(paste0("./pdf/", trait, "_hisse_turnover_recon.pdf"), width = 10, height = 12)
recon_plot <- plot.hisse.states(hisse_best.recon, rate.param = "turnover", 
                                #type = "phylogram",
                                show.tip.label = FALSE, 
                                label.offset = 2, fsize = 0.2, 
                                legend = "tips", legend.cex = 0.6, 
                                edge.width = 8, width.factor=0.4, 
                                rate.colors = c("darkcyan", "goldenrod"), 
                                state.colors = c("grey", "black"))

dev.off()

pdf(paste0("./pdf/", trait, "_hisse_speication_boxplot.pdf"), width = 6, height = 6)

model.ave.rates <- GetModelAveRates(hisse_best.recon, AIC.weights = NULL, type="tips",
                                    bound.par.matrix=cbind(c(-2000,-2000,-2000,-2000,-2000),c(2000,2000,2000,130000,2000)) )
boxplot(model.ave.rates$speciation~model.ave.rates$state,boxwex=0.5, 
        notch = TRUE,main="speciation", 
        col = c("steelblue","tomato"),
        xlab=" ", ylab="Speciation rate")
dev.off()

write.csv(model.ave.rates, paste0(trait, "_average_rate.csv"))

