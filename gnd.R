#######################################################################
# R FUNCTION TO CALCULATE GREENWOOD-NAM-D'AGOSTINO CALIBRATION TEST FOR SURVIVAL MODEL
# Most up-to date version of this code is available at http://ncook.bwh.harvard.edu/r-code.html
# Author: Olga Demler, BWH, HMS
# Version 2 - Updated 8/4/2015
# FOR MORE DETAILS SEE Demler, Paynter, Cook "Tests of Calibration and Goodness of Fit 
# in the Survival Setting" Stat Med 2015; 34(10):1659-80. PMID: 25684707
# TO RUN:
# GND.calib(pred,tvar,out,cens.t, groups, adm.cens)
# PARAMETERS:
# pred - PREDICTED PROBABILITIES OF AN EVENT CALCULATED FOR THE FIXED TIME WHICH IS THE SAME FOR ALL OBSERVATIONS (=adm.cens)
# out  - OUTCOME 0/1 1=EVENT
# cens.t - CENSORED/NOT CENSORED INDICATOR 1=CENSORED
# groups - GROUPING ASSIGNMENT FOR EACH OBSERVATION
# adm.cens - END OF STUDY TIME 
# REQUIRES AT LEAST 2 EVENTS PER GROUP, AT LEAST 5 EVENTS PER GROUP IS RECOMMENDED
# IF <2 EVENTS PER GROUP THEN QUITS
#######################################################################
kmdec=function(dec.num,dec.name, datain, adm.cens){
  stopped=0
  data.sub=datain[datain[,dec.name]==dec.num,]
  if (sum(data.sub$out)>1){
    avsurv=survfit(Surv(tvar,out) ~ 1, data=datain[datain[,dec.name]==dec.num,], error="g")
    avsurv.est=ifelse(min(avsurv$time)<=adm.cens,avsurv$surv[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],1)
    
    avsurv.stderr=ifelse(min(avsurv$time)<=adm.cens,avsurv$std.err[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],0)
    avsurv.stderr=avsurv.stderr*avsurv.est
    
    avsurv.num=ifelse(min(avsurv$time)<=adm.cens,avsurv$n.risk[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],0)
    
  } else {
    return(c(0,0,0,0,stopped=-1))
  }
  
  if (sum(data.sub$out)<5) stopped=1
  c(avsurv.est, avsurv.stderr, avsurv.num, dec.num, stopped) 
}#kmdec

GND.calib = function(pred, tvar, out, cens.t, groups, adm.cens){
  
  tvar.t=ifelse(tvar>adm.cens, adm.cens, tvar)
  out.t=ifelse(tvar>adm.cens, 0, out)
  
  datause=data.frame(pred=pred, tvar=tvar.t, out=out.t, count=1, cens.t=cens.t, dec=groups)
  numcat=length(unique(datause$dec))
  groups=sort(unique(datause$dec))
  
  kmtab=matrix(unlist(lapply(groups,kmdec,"dec",datain=datause, adm.cens)),ncol=5, byrow=TRUE)
  
  if (any(kmtab[,5] == -1)) stop("Stopped because at least one of the groups contains <2 events. Consider collapsing some groups.")
  else if (any(kmtab[,5] == 1)) warning("At least one of the groups contains < 5 events. GND can become unstable.\ 
(see Demler, Paynter, Cook 'Tests of Calibration and Goodness of Fit in the Survival Setting' DOI: 10.1002/sim.6428) \
Consider collapsing some groups to avoid this problem.")
  
  hltab=data.frame(group=kmtab[,4],
                   totaln=tapply(datause$count,datause$dec,sum),
                   censn=tapply(datause$cens.t,datause$dec,sum),
                   numevents=tapply(datause$out,datause$dec,sum),
                   expected=tapply(datause$pred,datause$dec,sum),
                   kmperc=1-kmtab[,1], 
                   kmvar=kmtab[,2]^2, 
                   kmnrisk=kmtab[,3],
                   expectedperc=tapply(datause$pred,datause$dec,mean))
  
  hltab$kmnum=hltab$kmperc*hltab$totaln
  hltab$GND_component=ifelse(hltab$kmvar==0, 0,(hltab$kmperc-hltab$expectedperc)^2/(hltab$kmvar))
  
  res_table = print(hltab[c(1,2,3,4,10,5,6,9,7,11)], digits=4)
  
  
  res_values = c(df=numcat-1, chi2gw=sum(hltab$GND_component),pvalgw=1-pchisq(sum(hltab$GND_component),numcat-1))
  
  return(list = c(res_table,res_values))
}#GND.calib