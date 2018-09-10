rm(list=ls())
options(width=100)
source('funcs.R')
library(xtable)
load("sdss_out.RData")

## get stars and tmss in same order
stars <- stars[order(stars[,1]),]
for(ii in 1:length(out)){
    out[[ii]] <- out[[ii]][order(names(out[[ii]]))]
}

## is the estimate within 1\% of truth
is_correct <- function(truth,estimate){
    return((estimate < 1.01*truth) & (estimate > .99*truth))
}
    

results <- matrix(0,nrow=4,ncol=9)
rownames(results) <- names(out)
for(ii in 1:nrow(results)){
    for(jj in 1:3){
        for(kk in 1:3){
        ests <- vapply(out[[ii]],function(x){x[jj,kk]},c(0))
        results[ii,jj+3*(kk-1)] <- mean(is_correct(stars[,2],ests))
    }
    }
}

construct_table <- function(results){
    rnames <- rownames(results)
    results <- matrix(substr(paste(as.character(round(results,2)),"00",sep=""),0,4),nrow=nrow(results))
    rows <- rep("a",nrow(results))
    for(ii in 1:nrow(results)){
        rows[[ii]] <- paste(rnames[ii],paste(results[ii,],collapse="&"),sep="&")
    }
    rows <- paste(paste(rows,collapse="\\\\"),"\\\\",sep="")
    pre <- "\\begin{table}[ht]
\\centering
\\begin{tabular}{c|ccc|ccc|ccc}
 &  & $K=1$ & &   & $K=2$ &  &  & $K=3$ &  \\\\ 
  \\hline
 & $\\Sigma^{-1}$ &  $I$ & $\\Delta$ & $\\Sigma^{-1}$ &  $I$ & $\\Delta$ & $\\Sigma^{-1}$ &  $I$ & $\\Delta$\\\\
  \\hline"
    post <- "\\hline
\\end{tabular}
\\caption{Fraction of periods estimated correctly using different weightings for models with $K=1,2,3$ harmonics. Ignoring the observation uncertainties ($I$) in the fitting is superior to using them ($\\Sigma^{-1}$). The strategy for determining an optimal weight function ($\\Delta$) does not provide much improvement over ignoring the weights. More complex models ($K=3$) perform worse than simple models ($K=1$) when there is limited data ($n=10$), but better when the functions are better sampled ($n=40$). The standard deviation on these accuracies is no larger than $\\sqrt{0.5(1-0.5)/238} \\approx 0.032$ .}
\\label{tab:period_est_results}
\\end{table}"
    rows <- paste(pre,rows,post,sep="")
    return(rows)
}


construct_table_simple <- function(results){
    rnames <- rownames(results)
    results <- matrix(substr(paste(as.character(round(results,2)),"00",sep=""),0,4),nrow=nrow(results))
    rows <- rep("a",nrow(results))
    for(ii in 1:nrow(results)){
        rows[[ii]] <- paste(rnames[ii],paste(results[ii,],collapse="&"),sep="&")
    }
    rows <- paste(paste(rows,collapse="\\\\"),"\\\\",sep="")
    pre <- "\\begin{table}[ht]
\\centering
\\begin{tabular}{c|cc|cc|cc}
 &   $K=1$ &    & $K=2$ &    & $K=3$ &  \\\\ 
  \\hline
n & $\\Sigma^{-1}$ &  $I$ & $\\Sigma^{-1}$ &  $I$ & $\\Sigma^{-1}$ &  $I$ \\\\
  \\hline"
    post <- "\\hline
\\end{tabular}
\\end{table}"
    rows <- paste(pre,rows,post,sep="")
    return(rows)
}




construct_table_very_simple <- function(results){
    rnames <- rownames(results)
    results <- matrix(substr(paste(as.character(round(results,2)),"00",sep=""),0,4),nrow=nrow(results))
    rows <- rep("a",nrow(results))
    for(ii in 1:nrow(results)){
        rows[[ii]] <- paste(rnames[ii],paste(results[ii,],collapse="&"),sep="&")
    }
    rows <- paste(paste(rows,collapse="\\\\"),"\\\\",sep="")
    pre <- "\\begin{table}[ht]
\\centering
\\begin{tabular}{c|cc}
n & $\\widehat{\\omega}(\\Sigma^{-1})$ &  $\\widehat{\\omega}(I)$ \\\\
  \\hline"
    post <- "\\hline
\\end{tabular}
\\caption{Fraction of periods estimated correctly (within 1\\% of truth). The standard deviation on these accuracies is no larger than $\\sqrt{0.5(1-0.5)/238} \\approx 0.032$ .}
\\label{tab:period_est_results}
\\end{table}"
    rows <- paste(pre,rows,post,sep="")
    return(rows)
}




tab <- construct_table(results)
writeLines(tab,"figs/period_est_results.tex")


tab2 <- construct_table_simple(results[,c(1,2,4,5,7,8)])
writeLines(tab2,"figs/period_est_results_simple.tex")



tab3 <- construct_table_very_simple(results[,c(1,2)])
writeLines(tab3,"figs/period_est_results_very_simple.tex")
