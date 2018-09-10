get_well_sampled <- function(tmss,min_points=50,min_bands=2){
    to.keep <- vapply(tmss,function(x){length(x)>=min_bands},c(TRUE))
    tmss <- tmss[to.keep]
    to.keep <- vapply(tmss,function(x){all(vapply(x,nrow,c(0))>min_points)},
                      c(TRUE))
    tmss <- tmss[to.keep]
    return(tmss)
}

get_tmss <- function(folder.data){
    folders <- list.dirs(folder.data)
    folders <- folders[folders != folder.data]
    folders <- paste(folders,"/",sep="")
    folders.names <- vapply(strsplit(folders,"/"),
                            function(x){x[length(x)]},c("hello"))
    ## get all light curves in tms form
    fnames <- list()
    for(ii in 1:length(folders)){
        fnames[[ii]] <- list.files(folders[ii])
    }
    fnames <- unique(unlist(fnames))
    tmss <- list()
    for(ii in 1:length(fnames)){
        tmss[[ii]] <- list()
        for(jj in 1:length(folders)){
            fname <- paste(folders[jj],fnames[ii],sep="")
            if(file.exists(fname))  tmss[[ii]][[folders.names[jj]]] <- read.table(fname)
        }
    }
    names(tmss) <- fnames
    return(tmss)
}

# makes design matrix for weighted least squares
construct_design <- function(w,K,t){
    predesign <- w*outer(t,1:K)
    return(cbind(1,cos(predesign),sin(predesign)))
}

compute_cs <- function(s,Xs,resids){
    p <- ncol(Xs[[1]])
    B <- matrix(0,nrow=p,ncol=p)
    A <- matrix(0,nrow=p,ncol=p)
    for(ii in 1:length(Xs)){
        n <- nrow(Xs[[ii]])
        B.temp <- try(solve(t(Xs[[ii]])%*%Xs[[ii]]/n),silent=TRUE)
        if (!inherits(B.temp, "try-error")) {
            B <- B + B.temp
            Y <- ((t(Xs[[ii]])%*%diag(s[[ii]]^{-4}*(resids[[ii]]^2 - s[[ii]]^2))%*%Xs[[ii]])
                  / sum(s[[ii]]^{-4}))
            A <- A + B.temp%*%Y%*%B.temp
        }
    }
    if(identical(B,matrix(0,nrow=p,ncol=p)) | sum(diag(A)) < 0){
        out <- 0
    }
    else {
        out <- sum(diag(A))/sum(diag(B))
    }
    return(out)
}

compute_residuals <- function(w,K,mag,weights,X){
    B <- t(X) %*% (X * weights)
    d <- t(X) %*% (mag * weights)
    tmp <- try(solve(B, d),silent=TRUE)
    if (!inherits(tmp, "try-error")) {
        r <- mag - X %*% tmp
    }
    else {
        r <- rep(0,length(weights))
    }
    return(r)
}   

reweight_band <- function(lc,w,K){
    Xs <- list()
    resids <- list()
    for(ii in 1:length(lc)){
        weights <- rep(1,nrow(lc[[ii]]))
        Xs[[ii]] <- construct_design(w,K,lc[[ii]][,1])
        resids[[ii]] <- as.numeric(compute_residuals(w,K,lc[[ii]][,2],weights,Xs[[ii]]))
    }
    return(list(Xs,resids))
}

reweight_lightcurve <- function(w,K,lc){
    Xsresids <- reweight_band(lc,w,K)
    s <- lapply(lc,function(x){x[,3]})
    return(compute_cs(s,Xsresids[[1]],Xsresids[[2]]))
}

update_sds <- function(lc,cs){
    lc.new <- lc
    for(ii in 1:length(lc.new)){
        lc.new[[ii]][,3] <- sqrt(lc.new[[ii]][,3]^2 + cs)
    }
    return(lc.new)
}

compute_rss <- function(w,K,lc,use_weights=TRUE){
    if(use_weights){
        weights <- 1 / lc[,3]^2
    }
    else {
        weights <- rep(1,nrow(lc))
    }
    X <- construct_design(w,K,lc[,1])
    r <- compute_residuals(w,K,lc[,2],weights,X)
    return(sum(weights * (r^2)))
}

## used for plotting
predict_lc <- function(w,K,lc,t=seq(0,(2*pi)/w,length=100)){
    weights <- 1 / lc[,3]^2
    X <- construct_design(w,K,lc[,1])
    B <- t(X) %*% (X * weights)
    d <- t(X) %*% (lc[,2] * weights)
    z <- solve(B, d)
    X <- construct_design(w,K,t)
    preds <- X %*% z
    return(cbind(t,preds))
}

get_freqs <-function(period_min,period_max,freq_del = 0.1/4000){
    freq_max <- 1/period_min
    freq_min <- 1/period_max
    return(2 * pi * seq(freq_min, freq_max, freq_del))
}

get_time_range <- function(tms){
    time_range <- 0
    for (tm in tms) {
        time_range <- max(max(tm[, 1]) - min(tm[, 1]), time_range)
    }
    return(time_range)
}

freq_estimate <- function(lc,K,period_min,period_max,use_weights=TRUE){
    omegas <- get_freqs(period_min, period_max,.1/get_time_range(lc))
    B <- length(lc)
    rss <- matrix(0,nrow=length(omegas),ncol=B)
    for(bb in 1:B){
        rss[,bb] <- vapply(omegas,compute_rss,c(0),K,lc[[bb]],use_weights)
    }
    rss <- rowSums(rss)
    return(omegas[which.min(rss)])
}

## find freq estimate:  K=1,2,3 / error=no,yes,c
all_freq_estimate <- function(lc,period_min,period_max){
    print("running lc")
    K <- 3 ## how many harmonics to try
    out <- matrix(0,nrow=4,ncol=K)
    colnames(out) <- 1:K
    rownames(out) <- c("error","no error","c est","c value")
    for(kk in 1:K){
        out[1,kk] <- freq_estimate(lc['g'],kk,period_min,period_max,use_weights=TRUE)
        out[2,kk] <- freq_estimate(lc['g'],kk,period_min,period_max,use_weights=FALSE)
        out[4,kk] <- reweight_lightcurve(out[2,kk],kk,lc['g'])
        lc.new <- update_sds(lc['g'],out[4,kk])
        out[3,kk] <- freq_estimate(lc.new,kk,period_min,period_max,use_weights=TRUE)
    }
    return(out)
}
