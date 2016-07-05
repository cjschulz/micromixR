rank_abund_sim <- function(richness, samplesize, power_exp){
 rcurve <- c(1:richness)
 rcurve <- rcurve^power_exp
 rcurve <- rcurve/(sum(rcurve))
 rcurve <- round(rcurve*samplesize)
}

#try1 <- rank_abund_sim(100, 10000, -0.8)
#sum(try1)
#try1
 