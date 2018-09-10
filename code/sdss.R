rm(list=ls())
set.seed(1234)
library('parallel')
source('funcs.R')

folder.data <- "data/"
tmss <- get_tmss(folder.data)
length(tmss)

## scatterplot of g band verus \sigma_g
tmss.g <- lapply(tmss,function(x){x[['g']]})
ms <- unlist(lapply(tmss.g,function(x){x[,2]}))
ts <- unlist(lapply(tmss.g,function(x){x[,3]}))
pdf("figs/mag_error.pdf",width=8,height=6)
par(mar=c(4.5,4.5,1,1))
plot(ms,ts,col="#00000040",ylim=c(0,.1),
     xlab="Magnitude",ylab="Magnitude Error Standard Deviation",
     cex.lab=1.3)
abline(v=18,lwd=2)
dev.off()

## select well sampled curves
tmss <- get_well_sampled(tmss,min_point=40,min_bands=5)
bright <- vapply(tmss,function(x){max(x[["g"]][,2]) < 18},c(TRUE))
tmss <- tmss[bright]
length(tmss)

## get true frequencies for tmss
stars <- read.table(paste(folder.data,"apj326724t2_mrt.txt",sep="/"),skip=42)
stars[,1] <- paste("LC_",stars[,1],".dat",sep="")
stars <- stars[,c(1,3)]
stars[,2] <- (2*pi)/stars[,2]
names(stars) <- c("ID","omega")
stars <- stars[stars[,1] %in% names(tmss),]
tmss <- tmss[names(tmss) %in% stars[,1]]
stars <- stars[order(stars[,1]),]
tmss <- tmss[order(names(tmss))]
identical(names(tmss),stars[,1])

## plot one lightcurve, both folded and unfolded
ii <- 1
jj <- 'g'
pdf("figs/unfolded_single.pdf",height=4,width=10)
w <- stars[ii,2]
par(mar=c(4,4,1,.5))
ran <- range(tmss[[ii]][[jj]][,2])
plot(tmss[[ii]][[jj]][,1],tmss[[ii]][[jj]][,2],ylim=c(ran[2],ran[1]),
     main="",
     xlab="Time (Days)",ylab="Magnitudes",cex.lab=1.5)
segments(tmss[[ii]][[jj]][,1],tmss[[ii]][[jj]][,2] - 2*tmss[[ii]][[jj]][,3],
         tmss[[ii]][[jj]][,1],tmss[[ii]][[jj]][,2] + 2*tmss[[ii]][[jj]][,3])
dev.off()




greater2 <- c(0,0)
pdf("figs/folded_single.pdf",height=4,width=10)
par(mar=c(4,4,1,.5))
preds1 <- predict_lc(w,1,tmss[[ii]][[jj]])
preds2 <- predict_lc(w,2,tmss[[ii]][[jj]])
ran <- range(c(tmss[[ii]][[jj]][,2],preds1[,2],preds2[,2]))
plot(tmss[[ii]][[jj]][,1] %% ((2*pi)/w),tmss[[ii]][[jj]][,2],ylim=c(ran[2],ran[1]),
     main="",
     xlab="Phase (Days)",
     ylab="Magnitudes",
     xaxs='i',
     cex.lab=1.5)
segments(tmss[[ii]][[jj]][,1] %% ((2*pi)/w),tmss[[ii]][[jj]][,2] - 2*tmss[[ii]][[jj]][,3],
         tmss[[ii]][[jj]][,1] %% ((2*pi)/w),tmss[[ii]][[jj]][,2] + 2*tmss[[ii]][[jj]][,3])
points(preds1[,1],preds1[,2],type='l',col='orange',lwd=2,lty=1)
points(preds2[,1],preds2[,2],type='l',col='blue',lwd=2,lty=2)
dev.off()
t <- tmss[[ii]][[jj]][,1] %% ((2*pi) / w)
m <- tmss[[ii]][[jj]][,2]
greater2[1] <- sum(abs(m - predict_lc(w,1,tmss[[ii]][[jj]],t)[,2]) > 2*tmss[[ii]][[jj]][,3])
greater2[2] <- sum(abs(m - predict_lc(w,8,tmss[[ii]][[jj]],t)[,2]) > 2*tmss[[ii]][[jj]][,3])
1 - (greater2 / nrow(tmss[[ii]][[jj]]))



## downsampled to n.points
n.points <- (1:4)*10
tmss.down <- list()
for(ii in 1:length(n.points)){
    tmss.down[[ii]] <- list()
    for(jj in 1:length(tmss)){
        tmss.down[[ii]][[jj]] <- list()
        for(kk in 1:length(tmss[[jj]])){
            tmss.down[[ii]][[jj]][[kk]] <- tmss[[jj]][[kk]][sample(1:nrow(tmss[[jj]][[kk]]),n.points[ii]),]
        }
        names(tmss.down[[ii]][[jj]]) <- names(tmss[[jj]]) 
    }
   names(tmss.down[[ii]]) <- names(tmss)
}
names(tmss.down) <- as.character(n.points)


period_min <- 0.2
period_max <- 1


## run all light curves through pipeline
mc.cores <- 7

## tmss.down <- tmss.down[1:2]
## for(ii in 1:length(tmss.down)){
##     tmss.down[[ii]] <- tmss.down[[ii]][1:2]
## }

out <- list()
for(ii in 1:length(tmss.down)){
    print(n.points[ii])
    out[[ii]] <- mclapply(tmss.down[[ii]],
                          all_freq_estimate,
                          period_min,period_max,
                          mc.cores=mc.cores)
}


names(out) <- n.points
save(tmss.down,tmss,out,stars,file="sdss_out.RData")
