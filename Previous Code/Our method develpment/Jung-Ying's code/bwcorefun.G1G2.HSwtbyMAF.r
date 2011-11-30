########################################################################################
## (2)  FUNCTIONS
########################################################################################
library(MASS)
library(haplo.stats)

##=========================================================================================
##  for generate QT using random sampling in Simu1
##  (i.e., generate Y from N(X_Gamma, phi*diag(XX)%*%S%*%diag(XX) + tau*S + simga*I)
##         where XX exclude intercept)
##=========================================================================================

  
  ## n.adjloci=2;  win=655;  Dlocus    = (win-1)*10+5;  totloci   = n.adjloci*2+1;  pos       = (Dlocus-n.adjloci):(Dlocus+n.adjloci);  hapbank   = t(hapmat[pos,])
  ## xx          = rnorm(ngeno, 0,1)
  ## hap=hapbank;Gamma0.QT=0.5; Gamma1=2; pphi=comb01[1,"phi"]; ttau = comb01[1,"tau"]; ssigma2e = comb01[1,"sigma"];  XX=xx
  getYY.qt.simu1.fun<-function(hap=hapbank, Gamma0.QT=0.5, Gamma1=2, pphi=phi, ttau=tau, ssigma2e = sigma, XX=xx)
  {
    ## XX must be a vector excluding intercept
    NNhap         = nrow(hap)
    draw          = t(sapply(rep(1, 2*ngeno), sample, x=(1:NNhap), replace=T))
    hapmat.sample = hap[draw, ] ## 2*ngeno by nloci;
    ##    V31 V32 V33 V34 V35
    ## 8    0   0   1   0   1
    ## 67   1   0   0   0   0
    ## 90   1   1   0   0   0
    ## 86   1   1   0   0   0
    ## .....
    nloci       = ncol(hapmat.sample)
    hap.sample  = apply(hapmat.sample, 1, paste, collapse="")  ## vector with len=2*ngeno
    ##       8      67      90      86 ...
    ## "00101" "10000" "11000" "11000" ...
    uniqhap     = sort(unique(hap.sample))
    ## [1] "00100" "00101" "10000" "10010" "11000"
    tmp         = matrix(hap.sample, ncol=2, byrow=T);
    hap1        = tmp[,1]; ## [1] "00101" "11000" "10010"
    hap1.design = model.matrix(~factor(hap1, levels=uniqhap)); hap1.design[,1]=0; hap1.design[,1]=1-rowSums(hap1.design)
    colnames(hap1.design) = uniqhap
    hap2        = tmp[,2]  ## [1] "10000" "11000" "11000"
    hap2.design = model.matrix(~factor(hap2, levels=uniqhap)); hap2.design[,1]=0; hap2.design[,1]=1-rowSums(hap2.design)
    colnames(hap2.design) = uniqhap
    ##----------------------------------------
    ## get RR
    ##----------------------------------------
    nhap        = length(uniqhap)
    uniqhap.mat = matrix(unlist(strsplit(uniqhap, "")), ncol=nloci, byrow=T)
    ij.lst      = cbind(rep(1:nhap, each=nhap), rep(1:nhap, nhap)); ij.lst = ij.lst[ij.lst[,1]<ij.lst[,2],]
    RR          = matrix(0, nhap, nhap); colnames(RR) = rownames(RR) = uniqhap
    for(kkk in 1:nrow(ij.lst))
      {
        ij=ij.lst[kkk,]; iii=ij[1]; jjj = ij[2]
        RR[iii,jjj] = sum(uniqhap.mat[iii,]==uniqhap.mat[jjj,])/nloci
      }
    tRR=t(RR)
    RR[lower.tri(RR)]= tRR[lower.tri(tRR)]
    diag(RR) = 1
    ##----------------------------------------
    ## get SSS when phase known
    ##----------------------------------------
    hap.design = hap1.design+hap2.design
    SSS        = hap.design  %*% RR %*% t(hap.design)
    ##----------------------------------------
    ## get yy.qt
    ##----------------------------------------
    SIGMA    = diag(XX) %*% SSS %*% diag(XX)
    yy.qt    = mvrnorm(n = 1, mu=Gamma0.QT+Gamma1*XX, Sigma=ssigma2e * diag(rep(1, ngeno))+ ttau * SSS + pphi * SIGMA, tol = 1e-6, empirical = FALSE)
    ##----------------------------------------
    ## get geno from hap
    ##----------------------------------------
    hap1.mat = hapmat.sample[seq(1, 2*ngeno, 2),]
    hap2.mat = hapmat.sample[seq(1, 2*ngeno, 2)+1,]
    geno = NULL;for(cc in 1:nloci){  geno = cbind(geno, hap1.mat[,cc], hap2.mat[,cc] )        }
    return(list("yy.QT"=yy.qt, "geno"=geno, "HH"=hap.design, "RR"=RR, "SSS"=SSS))
  }
  

  ## hap=hapbank;Gamma0.QT=0.5; Gamma1=2; ssigma2e = sigma2.e; XX=xx; GammaG=1; GammaGE = 0; Gmode="dom"
  getYY.qt.simu2.fun<-function(hap=hapbank, Gamma0.QT=0.5, Gamma1=2, ssigma2e = sigma2.e, XX=xx, GammaG=1, GammaGE = 0, Gmode="add")
  {
    nloci    = ncol(hap)
    ## --- gen draw ----
    draw     = sample(1:NNhap, 2*ngeno, replace=T)
    draw.mat = matrix(draw, ncol=2, byrow=T)
    ## --- gen yy ----
    preG = matrix((hap[draw, (n.adjloci+1)]), ncol=2, byrow=T)
    if(Gmode=="add")
      {
        GG   = rowSums(  preG )
      }else if (effect.int=="dom")
        {
          GG = rowSums(  preG )>0
        }else if (effect.int=="rec")
          {
            GG = (preG[,1]*preG[,2]) > 0
          }
    Ey   = Gamma0.QT + Gamma1*XX + GammaG*GG +GammaGE*GG*XX
    yy   = sapply(Ey, rnorm, n=1, sd=sqrt(ssigma2e))
    ## --- get geno ----
    hap1.mat     = hapbank[draw.mat[,1],]
    hap2.mat     = hapbank[draw.mat[,2],]
    geno.hapstat = hap1.mat+ hap2.mat
    geno = NULL;for(cc in 1:nloci){  geno = cbind(geno, hap1.mat[,cc], hap2.mat[,cc] )        }
    return(list("yy"=yy, "geno"=geno, "geno.hapstat"=geno.hapstat))
  }
  




  ##hap=hapbank; Gamma0.QT=0.5; Gamma1=2; ssigma2e = sigma2.e; XX=xx; GammaG1=gaG1; GammaG2=gaG2; GammaG1E = gaG1xE; GammaG2E=gaG2xE; Gmode="add"
  getYY.qt.simu3.G1G2.fun<-function(hap=hapbank, Gamma0.QT=0.5, Gamma1=2, ssigma2e = sigma2.e, XX=xx, GammaG1=gaG1, GammaG2=gaG2, GammaG1E = gaG1E, GammaG2E=gaG2E, Gmode="add", n.DL = nDL)
  {
    nloci    = ncol(hap)
    ## --- gen draw ----
    draw     = sample(1:NNhap, 2*ngeno, replace=T)
    draw.mat = matrix(draw, ncol=2, byrow=T)
    ## --- gen yy ----
    preG1 = matrix((hap[draw, (n.adjloci+1)]),     ncol=2, byrow=T)
    preG2 = matrix((hap[draw, (n.adjloci*2+1+1)]), ncol=2, byrow=T)
    if(Gmode=="add")
      {
        GG1   = rowSums(  preG1 )
        GG2   = rowSums(  preG2 )
      }else if (effect.int=="dom")
        {
          GG1 = rowSums(  preG1 )>0
          GG2 = rowSums(  preG2 )>0
        }else if (effect.int=="rec")
          {
            GG1 = (preG1[,1]*preG1[,2]) > 0
            GG2 = (preG2[,1]*preG2[,2]) > 0
          }
    Ey   = Gamma0.QT + Gamma1*XX + GammaG1*GG1 + GammaG2*GG2 +GammaG1E*GG1*XX+GammaG2E*GG2*XX
    yy   = sapply(Ey, rnorm, n=1, sd=sqrt(ssigma2e))
    ## --- get geno ----
    hap1.mat     = hapbank[draw.mat[,1],]
    hap2.mat     = hapbank[draw.mat[,2],]
    geno.hapstat = hap1.mat+ hap2.mat
    geno = NULL;for(cc in 1:nloci){  geno = cbind(geno, hap1.mat[,cc], hap2.mat[,cc] )        }
    return(list("yy"=yy, "geno"=geno, "geno.hapstat"=geno.hapstat))
  }
  

##==========================================================
##  EM algorithm for REML
##==========================================================

  ##    sigma2.new = 1; tau.new=1; dif.tol=DIF.TOL; RRR=RRRR1; HHH=HHHH1; XXX=newx.adj; Sss=SS1; nsubj=n.subj; n.uniqhap=ll1; y.vec=y
  ##    sigma2.new = 1; tau.new=1; dif.tol=DIF.TOL; RRR=RRRR0; HHH=HHHH0; XXX=newx.adj; Sss=SS; nsubj=n.subj; n.uniqhap=ll0; y.vec=y

  getEME.fun<-function(sigma2.new=0.5, tau.new=0.5, dif.tol =DIF.TOL, RRR = RRRR, HHH=HHHH, XXX = newx.adj, Sss=SSS, nsubj=n.subj, n.uniqhap=ll, y.vec=y)
       {
# START CSCS
# CSCS Move unchanging values out of loop.
# CSCS Pre-calculate some unchanging values. 
         HR = HHH %*% RRR
         tHR = t(HR)
         HHR     = crossprod(HHH, HR)
         invRRR = solve(RRR)
         tXXX = t(XXX)
         I.nsubj = diag(rep(1,nsubj)) # CSCS Now a global
         AAA        = I.nsubj - XXX %*% solve(crossprod(XXX)) %*% tXXX
         HAH        = crossprod(HHH, AAA) %*% HHH
         I.uniqhap = diag(rep(1,n.uniqhap))
# END CSCS
         
         dif        = 1
         j.EM       = 1
         while(dif > dif.tol)
           {
             ## cat(j.EM, "tau=", tau.new, "sigma2=",sigma2.new,"\n")
             sigma2.old = sigma2.new
             tau.old    = tau.new
             #VVV        = tau.old *Sss + sigma2.old*diag(rep(1,nsubj)) # CSCS
             VVV        = tau.old * Sss + sigma2.old * I.nsubj          # CSCS
             ##------------------
             ## get V^{-1}
             ##------------------

             # CSCS HR      = HHH %*% RRR
             # CSCS HHR     = t(HHH) %*% HR
             # IHHR    = sigma2.old *diag(rep(1, n.uniqhap)) +tau.old * HHR

             IHHR    = sigma2.old * I.uniqhap + tau.old * HHR
             invIHHR = solve(IHHR)

             # CSCS invVVV  = diag(rep(1, nsubj))/sigma2.old - tau.old/sigma2.old *  HR %*% invIHHR %*% t(HHH)
             invVVV  = I.nsubj / sigma2.old - tau.old / sigma2.old * HR %*% tcrossprod(invIHHR, HHH)

             ##------------------
             ## get PPP
             ##------------------
             invVX      = invVVV %*% XXX
             XinvVX     = tXXX %*% invVX
             inv.XinvVX = solve(XinvVX)
             inv.XinvVX.tinvVX = tcrossprod(inv.XinvVX, invVX)
             # CSCS PPP        = invVVV - invVX %*% inv.XinvVX %*% t(invVX)
             PPP        = invVVV - invVX %*% inv.XinvVX.tinvVX
             RHPHR      = tHR %*% PPP %*% HR # CSCS Moved up.

             ##-----------------------
             ## get Varbeta 
             ##-----------------------
             Varbeta  = tau.old * RRR - tau.old^2 * RHPHR

             ##-----------------------
             ## get Ebeta = beta.hat
             ##-----------------------
             gammahat = inv.XinvVX.tinvVX %*% y.vec           # CSCS
             residual = y.vec - XXX %*% gammahat
             Ebeta    = tau.old * tHR %*% invVVV %*% residual # CSCS
             Ebeta2    = tau.old * tHR %*% PPP %*% y.vec
             cat(Ebeta, Ebeta2,"\n")
             ##-----------------------
             ## get tau.new
             ##-----------------------
             # CSCS invRRR   = solve(RRR)
             # CSCS BRB      = t(Ebeta) %*% invRRR %*% Ebeta
             # CSCS tmp      = invRRR %*% Varbeta
             # CSCS tau.new  = c( (BRB + sum(diag(tmp)))/n.uniqhap )
             BRB      = crossprod(Ebeta, invRRR) %*% Ebeta
             tmp      = get.sum.diag(invRRR, Varbeta) # invRRR %*% Varbeta
             tau.new  = c( (BRB + tmp) / n.uniqhap )

             ##-----------------------
             ## get sigma2.new
             ##-----------------------
             ddd        = ncol(XXX)
             r1         = y.vec - HHH %*% Ebeta
             # CSCS AAA        = diag(rep(1, nsubj)) - XXX %*% solve(t(XXX)%*% XXX) %*% t(XXX)
             # CSCS r1Ar1      = t(r1)%*% AAA %*% r1
             # CSCS HAH        = t(HHH) %*% AAA %*% HHH
             # CSCS RHPHR      = t(HR)%*% PPP %*% HR
             # CSCS tmp        = HAH %*% (tau.old*RRR - tau.old^2 * RHPHR) 
             # CSCS sigma2.new = c(1/(nsubj-ddd)*( r1Ar1 + sum(diag(tmp))   ))
             r1Ar1      = crossprod(r1, AAA) %*% r1
             tmp        = get.sum.diag(HAH, (tau.old*RRR - tau.old^2 * RHPHR)) # HAH %*% (tau.old*RRR - tau.old^2 * RHPHR) 
             sigma2.new = c(1/(nsubj-ddd)*( r1Ar1 + tmp ))
             ## sigma2.new = sigma2.e
             ##-----------------------
             ## dif1        = abs(sigma2.new-sigma2.old)
             ## dif2        = abs(tau.new   -tau.old)
             dif1       = abs(sigma2.new-sigma2.old)/(sigma2.old+1e-10)
             dif2       = abs(tau.new   -tau.old   )/(tau.old   +1e-10)
             dif        = max(dif1, dif2)
             j.EM       = j.EM+1

             if ((j.EM %% 10000) == 0) cat("\n(EME loop ", j.EM, ", dif = ", dif, ")\n") 
           }
         ## cat(j.EM, "tau=", tau.new, "sigma2=",sigma2.new,"\n")
         ## cat("------ TRUE: tau=", tau, "sigma2=",sigma2.e,"-----\n")
         ## CSCS Already calculated above: gammahat = inv.XinvVX %*% t(invVX) %*% y.vec; ## cat("gamma.hat=", gammahat)
         return(list("tau.EME"=tau.new,
                     "sigma2.EME"=sigma2.new,
                     "P" = PPP,
                     "V" = VVV,
                     "invV" = invVVV,
                     ## "XinvVX"= XinvVX,
                     "gammahat" = gammahat
                     ))
       }

##=============================
## single snp analysis.
##=============================

## CSCS Put this function back.
        snp.pval.xxGE.xxE.fun<-function(vec=Geno.hapstat[,1],Y=yy, xxxGE=xx, xxxE=xxE, family="gaussian")
          {
            ngeno = length(vec) # CSCS
 
            ## xxxGE is the variable for GxE interaction
            ## xxxE is the additional covariates that need to be adjusted in residuals. (Although
            ##      xxxGE will also be included in the residual adjustment, no need to include xxGE again
            ##      because xxxGE will be automatically in the model
            ## ps. the "xxx" in hsreg should be = cbind(xxxGE, xxxE)
            if(is.null(xxxE))
            {
                ##-----------------
                ## 1df, G effect
                ##-----------------
                fit1.G      = glm(Y~xxxGE+vec,      family=family)
                pval1.G     = anova(fit1.G,test="Chisq")["vec","P(>|Chi|)"]; names(pval1.G) = "pval1.G"
                ## pval1.G  = summary(fit1.G)$coef["vec", "Pr(>|t|)"]; names(pval1.G) = "pval.G"  ## wald's test

                ##-----------------------
                ## 1df, G and GxE effect 
                ##-----------------------
                fit1.GxE    = glm(Y~xxxGE+vec+xxxGE*vec,family=family)
                anova1      = anova(fit1.GxE,test="Chisq")
                pval1.GxE   = anova1["xxxGE:vec","P(>|Chi|)"]; names(pval1.GxE) = "pval1.GxE"
                ## pval1.GxE= anova1[c("vec","xxxGE:vec"),"P(>|Chi|)"]; names(pval1.GxE) = c("pval1.G.inGxE", "pval1.GxE.inGxE")
                ## pval1.GxE= summary(fit1.GxE)$coef[c("vec","xxxGE:vec"), "Pr(>|t|)"]; names(pval1.GxE) = c("pval.G", "pval.GxE")
                mse2.snp1   = fit1.GxE$deviance/(ngeno-length(fit1.GxE$coef))
                ## sqrt(mse2.snp1) should be equal to "summary(lm(Y~xxxGE+vec+xxxGE*vec))"
                pval1.joint = pchisq(abs(diff( anova1[c("xxxGE","xxxGE:vec"),"Resid. Dev"])) / mse2.snp1,
                  abs(diff(anova1[c("xxxGE","xxxGE:vec"),"Resid. Df"])), lower=F)
                names(pval1.joint)="pval1.joint"
                ## ---- for debug: the following should match with the answer given in anova1
                ## pchisq(abs(diff( anova1[c("xxxGE",    "vec"),"Resid. Dev"])) / mse2.snp1, abs(diff(anova1[c("xxxGE",    "vec"),"Resid. Df"])), lower=F)
                ## pchisq(abs(diff( anova1[c("xxxGE:vec","vec"),"Resid. Dev"])) / mse2.snp1, abs(diff(anova1[c("xxxGE:vec","vec"),"Resid. Df"])), lower=F)
                
                ##-----------------------
                ## 2df, G effect 
                ##-----------------------
                fit2.G      = glm(Y~xxxGE+factor(vec), family=family);
                pval2.G     = anova(fit2.G,test="Chisq")["factor(vec)","P(>|Chi|)"]; names(pval2.G) = "pval2.G"
                ##-----------------------
                ## 2df, G and GxE effect 
                ##-----------------------
                fit2.GxE    = glm(Y~xxxGE+factor(vec)+xxxGE*factor(vec), family=family);
                anova2      = anova(fit2.GxE,test="Chisq")
                pval2.GxE   = anova2["xxxGE:factor(vec)","P(>|Chi|)"]; names(pval2.GxE) = c("pval2.GxE")
                mse2.snp2   = fit2.GxE$deviance/(ngeno-length(fit2.GxE$coef))
                ## sqrt(mse2.snp2) should be equal to "summary(lm(Y~xxxGE+factor(vec)+xxxGE*factor(vec)))"
                pval2.joint = pchisq(abs(diff( anova2[c("xxxGE","xxxGE:factor(vec)"),"Resid. Dev"]))/ mse2.snp2, abs(diff(anova2[c("xxxGE","xxxGE:factor(vec)"),"Resid. Df"])), lower=F)
                names(pval2.joint)="pval2.joint"
                ## ---- for debug: the following should match with the answer given in anova2
                ## pchisq(abs(diff( anova2[c("xxxGE",            "factor(vec)"), "Resid. Dev"])) / mse2.snp2, abs(diff(anova2[c("xxxGE",            "factor(vec)"),"Resid. Df"])), lower=F)
                ## pchisq(abs(diff( anova2[c("xxxGE:factor(vec)","factor(vec)"), "Resid. Dev"])) / mse2.snp2, abs(diff(anova2[c("xxxGE:factor(vec)","factor(vec)"),"Resid. Df"])), lower=F)
                
            } else
            {## these are with additional xxxE
              ##-----------------
              ## 1df, G effect
              ##-----------------
              fit1.G      = glm(Y~xxxE+xxxGE+vec,      family=family)
              ## here it is importnat to put "xxxE+xxxGE" instead of "xxxGE+xxxE", as
              ## in anova analysis, we need the deviance that include all xx (and only xx)
              pval1.G     = anova(fit1.G,test="Chisq")["vec","P(>|Chi|)"]; names(pval1.G) = "pval1.G"
              ## pval1.G  = summary(fit1.G)$coef["vec", "Pr(>|t|)"]; names(pval1.G) = "pval.G"  ## wald's test
              ##-----------------------
              ## 1df, G and GxE effect 
              ##-----------------------
              fit1.GxE    = glm(Y~xxxE+xxxGE+vec+xxxGE*vec,family=family)
              anova1      = anova(fit1.GxE,test="Chisq")
              pval1.GxE   = anova1["xxxGE:vec","P(>|Chi|)"]; names(pval1.GxE) = "pval1.GxE"
              ## pval1.GxE= anova1[c("vec","xxxGE:vec"),"P(>|Chi|)"]; names(pval1.GxE) = c("pval1.G.inGxE", "pval1.GxE.inGxE")
              ## pval1.GxE= summary(fit1.GxE)$coef[c("vec","xxxGE:vec"), "Pr(>|t|)"]; names(pval1.GxE) = c("pval.G", "pval.GxE")
              mse2.snp1   = fit1.GxE$deviance/(ngeno-length(fit1.GxE$coef))
              ## sqrt(mse2.snp1) should be equal to "summary(lm(Y~xxxGE+xxxE+vec+xxxGE*vec))"
              pval1.joint = pchisq(abs(diff( anova1[c("xxxGE","xxxGE:vec"),"Resid. Dev"])) / mse2.snp1, abs(diff(anova1[c("xxxGE","xxxGE:vec"),"Resid. Df"])), lower=F)
              names(pval1.joint)="pval1.joint"
              ## ---- for debug: the following should match with the answer given in anova1
              ## pchisq(abs(diff( anova1[c("xxxGE",    "vec"),"Resid. Dev"])) / mse2.snp1, abs(diff(anova1[c("xxxGE",    "vec"),"Resid. Df"])), lower=F)
              ## pchisq(abs(diff( anova1[c("xxxGE:vec","vec"),"Resid. Dev"])) / mse2.snp1, abs(diff(anova1[c("xxxGE:vec","vec"),"Resid. Df"])), lower=F)
              
              ##-----------------------
              ## 2df, G effect 
              ##-----------------------
              fit2.G      = glm(Y~xxxE+xxxGE+factor(vec), family=family);
              pval2.G     = anova(fit2.G,test="Chisq")["factor(vec)","P(>|Chi|)"]; names(pval2.G) = "pval2.G"
              ##-----------------------
              ## 2df, G and GxE effect 
              ##-----------------------
              fit2.GxE    = glm(Y~xxxE+xxxGE+factor(vec)+xxxGE*factor(vec), family=family);
              anova2      = anova(fit2.GxE,test="Chisq")
              pval2.GxE   = anova2["xxxGE:factor(vec)","P(>|Chi|)"]; names(pval2.GxE) = c("pval2.GxE")
              mse2.snp2   = fit2.GxE$deviance/(ngeno-length(fit2.GxE$coef))
              ## sqrt(mse2.snp2) should be equal to "summary(lm(Y~xxxE+xxxGE+factor(vec)+xxxGE*factor(vec)))"
              pval2.joint = pchisq(abs(diff( anova2[c("xxxGE","xxxGE:factor(vec)"),"Resid. Dev"]))/ mse2.snp2, abs(diff(anova2[c("xxxGE","xxxGE:factor(vec)"),"Resid. Df"])), lower=F)
              names(pval2.joint)="pval2.joint"
              ## ---- for debug: the following should match with the answer given in anova2
              ## pchisq(abs(diff( anova2[c("xxxGE",            "factor(vec)"), "Resid. Dev"])) / mse2.snp2, abs(diff(anova2[c("xxxGE",            "factor(vec)"),"Resid. Df"])), lower=F)
              ## pchisq(abs(diff( anova2[c("xxxGE:factor(vec)","factor(vec)"), "Resid. Dev"])) / mse2.snp2, abs(diff(anova2[c("xxxGE:factor(vec)","factor(vec)"),"Resid. Df"])), lower=F)
            }
            return(c(pval1.G,pval1.GxE,pval1.joint, pval2.G,pval2.GxE, pval2.joint))
          }

        ## vec=Geno.hapstat[,3];Y=yy; xxx=xx;family="gaussian"
        snp.pval.fun<-function(vec=Geno.hapstat[,1],Y=yy, xxx=xx,family="gaussian")
          {
            if(length(unique(vec))<2)
              {
                pval1.G=pval1.GxE=pval1.joint=pval2.G=pval2.GxE= pval2.joint=1
                names(pval1.G) =     "pval1.G"
                names(pval1.GxE) =   "pval1.GxE"
                names(pval1.joint) = "pval1.joint"
                names(pval2.G) =     "pval2.G"
                names(pval2.GxE) =   "pval2.GxE"
                names(pval2.joint) = "pval2.joint"
              }else
            {
              ##-----------------
              ## 1df, G effect
              ##-----------------
              fit1.G      = glm(Y~xxx+vec,      family=family)
              pval1.G     = anova(fit1.G,test="Chisq")["vec","P(>|Chi|)"]; names(pval1.G) = "pval1.G"
              ## pval1.G  = summary(fit1.G)$coef["vec", "Pr(>|t|)"]; names(pval1.G) = "pval.G"  ## wald's test
              
              ##-----------------------
              ## 1df, G and GxE effect 
              ##-----------------------
              fit1.GxE    = glm(Y~xxx+vec+xxx*vec,family=family)
              anova1      = anova(fit1.GxE,test="Chisq")
              pval1.GxE   = anova1["xxx:vec","P(>|Chi|)"]; names(pval1.GxE) = "pval1.GxE"
              ## pval1.GxE= anova1[c("vec","xxx:vec"),"P(>|Chi|)"]; names(pval1.GxE) = c("pval1.G.inGxE", "pval1.GxE.inGxE")
              ## pval1.GxE= summary(fit1.GxE)$coef[c("vec","xxx:vec"), "Pr(>|t|)"]; names(pval1.GxE) = c("pval.G", "pval.GxE")
              mse2.snp1   = fit1.GxE$deviance/(ngeno-length(fit1.GxE$coef))
              ## sqrt(mse2.snp1) should be equal to "summary(lm(Y~xxx+vec+xxx*vec))"
              pval1.joint = pchisq(abs(diff( anova1[c("xxx","xxx:vec"),"Resid. Dev"])) / mse2.snp1, abs(diff(anova1[c("xxx","xxx:vec"),"Resid. Df"])), lower=F)
              names(pval1.joint)="pval1.joint"
              ## ---- for debug: the following should match with the answer given in anova1
              ## pchisq(abs(diff( anova1[c("xxx",    "vec"),"Resid. Dev"])) / mse2.snp1, abs(diff(anova1[c("xxx",    "vec"),"Resid. Df"])), lower=F)
              ## pchisq(abs(diff( anova1[c("xxx:vec","vec"),"Resid. Dev"])) / mse2.snp1, abs(diff(anova1[c("xxx:vec","vec"),"Resid. Df"])), lower=F)
              
              ##-----------------------
              ## 2df, G effect 
              ##-----------------------
              fit2.G      = glm(Y~xxx+factor(vec), family=family);
              pval2.G     = anova(fit2.G,test="Chisq")["factor(vec)","P(>|Chi|)"]; names(pval2.G) = "pval2.G"
              ##-----------------------
              ## 2df, G and GxE effect 
              ##-----------------------
              fit2.GxE    = glm(Y~xxx+factor(vec)+xxx*factor(vec), family=family);
              anova2      = anova(fit2.GxE,test="Chisq")
              pval2.GxE   = anova2["xxx:factor(vec)","P(>|Chi|)"]; names(pval2.GxE) = c("pval2.GxE")
              mse2.snp2   = fit2.GxE$deviance/(ngeno-length(fit2.GxE$coef))
              ## sqrt(mse2.snp2) should be equal to "summary(lm(Y~xxx+factor(vec)+xxx*factor(vec)))"
              pval2.joint = pchisq(abs(diff( anova2[c("xxx","xxx:factor(vec)"),"Resid. Dev"]))/ mse2.snp2, abs(diff(anova2[c("xxx","xxx:factor(vec)"),"Resid. Df"])), lower=F)
              names(pval2.joint)="pval2.joint"
              ## ---- for debug: the following should match with the answer given in anova2
              ## pchisq(abs(diff( anova2[c("xxx",            "factor(vec)"), "Resid. Dev"])) / mse2.snp2, abs(diff(anova2[c("xxx",            "factor(vec)"),"Resid. Df"])), lower=F)
              ## pchisq(abs(diff( anova2[c("xxx:factor(vec)","factor(vec)"), "Resid. Dev"])) / mse2.snp2, abs(diff(anova2[c("xxx:factor(vec)","factor(vec)"),"Resid. Df"])), lower=F)
            }
            return(c(pval1.G,pval1.GxE,pval1.joint, pval2.G,pval2.GxE, pval2.joint))
          }

##============================================
##    VC test
##============================================
# START CSCS
get.smat.fun = function(geno)
{
  n.cols = ncol(geno)
  odds = seq(1,n.cols,2)
  evens = odds+1
  geno.add = geno[,odds] + geno[,evens]

  geno.add = 1 - as.matrix(geno.add)

  n.loci = ncol(geno.add)
  SSS = ((2/n.loci) * tcrossprod(geno.add)) + 2

  return (SSS)
}

## geno=Geno
get.smat.wtbyMAF.34.fun = function(geno)
{
  n.cols = ncol(geno)
  odds = seq(1,n.cols,2)
  evens = odds+1

  geno.A =  geno[,odds] + geno[,evens]
  geno.a =  2-geno.A
  n.loci = ncol(geno.A)

  AFA = colSums(geno.A)/2/ngeno ## Allle Freq of "A"
  pseudo1 = 0.99999
  AFA[AFA==1]=pseudo1
  AFA[AFA==0]=(1-pseudo1)
  AFa = 1-AFA ## Allle Freq of "a"
  ## AFa = colSums(geno.a)/2/ngeno ## Allle Freq of "a"

  pre.A = 1/AFA ;   pre.a = 1/AFa;
  wt3.A = (pre.A)
  wt3.a = (pre.a)
  
  pre.A = 1/sqrt(AFA);   pre.a = 1/sqrt(AFa)
  wt4.A = (pre.A)
  wt4.a = (pre.a)
  ## cbind(AFA, wt3.A, wt4.A);   cbind(AFa, wt3.a, wt4.a)

  wtgeno.A = t( t(geno.A) * sqrt(wt3.A) )
  wtgeno.a = t( t(geno.a) * sqrt(wt3.a) )
  SSS.wt3  = (  tcrossprod(wtgeno.A) +tcrossprod(wtgeno.a)  ) /n.loci
  ## range((geno.A %*% diag(wt3.A) %*% t(geno.A) + geno.a %*% diag(wt3.a) %*% t(geno.a))/n.loci - SSS.wt3)

  wtgeno.A = t( t(geno.A) * sqrt(wt4.A) )
  wtgeno.a = t( t(geno.a) * sqrt(wt4.a) )
  SSS.wt4  = (  tcrossprod(wtgeno.A) +tcrossprod(wtgeno.a)  ) /n.loci
  ## range((geno.A %*% diag(wt4.A) %*% t(geno.A) + geno.a %*% diag(wt4.a) %*% t(geno.a))/n.loci - SSS.wt4)
  return (list("SSS.wt3"=SSS.wt3, "SSS.wt4"=SSS.wt4))
  ##pre.A = 1     ;   pre.a = 1 ;   wt0.A = (pre.A);  wt0.a = (pre.a)
  ## wtgeno.A = t( t(geno.A) * sqrt(wt0.A) );   wtgeno.a = t( t(geno.a) * sqrt(wt0.a) )
  ## SSS.wt0  = (  tcrossprod(wtgeno.A) +tcrossprod(wtgeno.a)  ) /n.loci
  ## range(SSS.wt0-get.smat.fun(geno))

  ## return (list("SSS.wt0"=SSS.wt0,"SSS.wt3"=SSS.wt3, "SSS.wt4"=SSS.wt4))
}


## geno=Geno
get.smat.wtbyMAF.fun = function(geno)
{
  n.cols = ncol(geno)
  odds = seq(1,n.cols,2)
  evens = odds+1

  geno.A =  geno[,odds] + geno[,evens]
  geno.a =  2-geno.A
  n.loci = ncol(geno.A)

  AFA = colSums(geno.A)/2/ngeno ## Allle Freq of "A"
  pseudo1 = 0.99999
  AFA[AFA==1]=pseudo1
  AFA[AFA==0]=(1-pseudo1)
  AFa = 1-AFA ## Allle Freq of "a"
  ## AFa = colSums(geno.a)/2/ngeno ## Allle Freq of "a"

  pre.A = 1/AFA ;   pre.a = 1/AFa;
  wt1.A = (pre.A)/(pre.A+pre.a) *2
  wt1.a = (pre.a)/(pre.A+pre.a) *2
  
  pre.A = 1/sqrt(AFA);   pre.a = 1/sqrt(AFa)
  wt2.A = (pre.A)/(pre.A+pre.a) *2
  wt2.a = (pre.a)/(pre.A+pre.a) *2
  ## cbind(AFA, wt1.A, wt2.A);   cbind(AFa, wt1.a, wt2.a)

  wtgeno.A = t( t(geno.A) * sqrt(wt1.A) )
  wtgeno.a = t( t(geno.a) * sqrt(wt1.a) )
  SSS.wt1  = (  tcrossprod(wtgeno.A) +tcrossprod(wtgeno.a)  ) /n.loci
  ## range((geno.A %*% diag(wt1.A) %*% t(geno.A) + geno.a %*% diag(wt1.a) %*% t(geno.a))/n.loci - SSS.wt1)

  wtgeno.A = t( t(geno.A) * sqrt(wt2.A) )
  wtgeno.a = t( t(geno.a) * sqrt(wt2.a) )
  SSS.wt2  = (  tcrossprod(wtgeno.A) +tcrossprod(wtgeno.a)  ) /n.loci
  ## range((geno.A %*% diag(wt2.A) %*% t(geno.A) + geno.a %*% diag(wt2.a) %*% t(geno.a))/n.loci - SSS.wt2)
  return (list("SSS.wt1"=SSS.wt1, "SSS.wt2"=SSS.wt2))
  ##pre.A = 1     ;   pre.a = 1 ;   wt0.A = (pre.A)/(pre.A+pre.a) *2;  wt0.a = (pre.a)/(pre.A+pre.a) *2
  ##wtgeno.A = t( t(geno.A) * sqrt(wt0.A) );   wtgeno.a = t( t(geno.a) * sqrt(wt0.a) )
  ## SSS.wt0  = (  tcrossprod(wtgeno.A) +tcrossprod(wtgeno.a)  ) /n.loci
  ## range(SSS.wt0-get.smat.fun(geno))

  ## return (list("SSS.wt0"=SSS.wt0,"SSS.wt1"=SSS.wt1, "SSS.wt2"=SSS.wt2))
}

## This function retruns the sum of the diagonal elements in the square of a matrix.
## For a large matrix this is faster than calculating the entire square of the matrix
## and then summing the diagonals.
get.sum.diag <- function(m,n=NULL)
{
  if (length(n) == 0)
  {
    r = nrow(m)
    c = ncol(m)
    if ((r<10) && (c<10)) return (sum(diag(m%*%m)))
    
    s = 0
    for (i in 1:r)
    {
      s = s + drop(m[i,] %*% m[,i])
    }
  } else {
    r = nrow(m)
    c = ncol(n)
    if ((r<10) && (c<10)) return (sum(diag(m%*%n)))
    s = 0
    for (i in 1:r)
    {
      s = s + drop(m[i,] %*% n[,i])
    }
  }
  return(s)
}
# END CSCS

##----------------------
## library(haplo.stats,lib.loc=lib.dir)
##=============================
## library("foreign")

    ## y=yy; geno=Geno; x.adj=xx.center; trait.type="gaussian"; SSS=SS
    ## y=yy; geno=Geno; x.adj=xx.std; trait.type="gaussian"; SSS=SS1
    vctest.qtG.GxE.jointwt.fun<-function(y=yy, geno=Geno, x.adj=xx, trait.type="gaussian", SSS=SS)
     {
      ## This function take the estimation error of sigma2 and tau into account
      trait.int <- charmatch(trait.type, c("gaussian", "binomial"))
      if (is.na(trait.int))        stop("Invalid trait type")
      if (trait.int == 0)          stop("Ambiguous trait type")
#CSCS if (length(y) != nrow(geno)) stop("Dims of y and geno are not compatible")
      if (length(y) != ncol(geno)) stop(paste("Dims of y and geno are not compatible", length(y), ncol(geno)) )

## CSCS Changed for genotype data in row-per-snp format.      
#CSCS n.subj <- length(y)
#CSCS n.ij   <- n.subj *(n.subj-1)/2
#CSCS n.loci <- ncol(geno)/2
#CSCS if (n.loci != (floor(ncol(geno)/2)))  stop("Odd number of cols of geno")

      n.subj <- length(y)
      n.loci <- nrow(geno)

      adjusted <- TRUE ;        if (all(is.na(x.adj)))       adjusted <- FALSE
      if (adjusted)
        {
          ## ## standardized X
          ## x.adj = (x.adj-mean(x.adj))/sd(x.adj)
          x.adj <- as.matrix(x.adj)
          if (nrow(x.adj) != length(y)) stop("Dims of y and x.adj are not compatible")
        }
       miss <- is.na(y)
       if (adjusted)       miss <- miss | apply(is.na(x.adj), 1, any)
       y    <- as.numeric(y[!miss])
       geno <- geno[,!miss] # CSCS For different genotype data format.
       if (adjusted)       x.adj <- x.adj[!miss, , drop = FALSE]
      ##-----------------------------
      ##-----------------------------
      ##-----------------------------
      if (!adjusted){ newx.adj = as.matrix(rep(1, n.subj))      }
      if (adjusted){  newx.adj = as.matrix( cbind(rep(1, n.subj), x.adj) )}
      ##########################################
      ## (I) marginal GxE test
      ##########################################
      ##------------------
      ## obtain SSS's FDF decomposition or eigen decomposition
      ##------------------
      ##--FDF decomposition--------?????????--------
      ## out = qr(SSS)
      ## qqq=qr.Q(out)
      ## rrr=qr.R(out)
      ##------------------------
      ##  eigen decomposition
      ##------------------------
      ##eigen.stm = proc.time()
      eg     = eigen(SSS, symmetric=T);
      tmpkey = round(eg$values, 10)>0;
      ll     = sum(tmpkey) ## rank of SSS
      RRRR   = diag(eg$values[tmpkey])
      HHHH   = eg$vectors[,tmpkey]
      ##show.time("Eigen time = ", eigen.stm)
      ## RRR=gendata$RR; HHH = gendata$HH
      ##------------------------
      ##  EM algorithm
      ##------------------------
      eme        = getEME.fun(sigma2.new = 1, tau.new=1, dif.tol=DIF.TOL, RRR=RRRR, HHH=HHHH, XXX=newx.adj, Sss=SSS, nsubj=n.subj, n.uniqhap=ll, y.vec=y)
      tau.eme    = eme$tau.EME
      sigma2.eme = eme$sigma2.EME
      ##-----------------------
      ## get T.GxE
      ##-----------------------
      P.mat.GxEm        = eme$P
      ## debug: range( P.mat.GxEm -      eme$invV +      eme$invV %*% (newx.adj) %*% solve(t(newx.adj)%*% eme$invV%*% newx.adj) %*% t(newx.adj)%*% eme$invV)
      xxGE             = c(x.adj[,1])                       ## CSCS 
      Sigma.mat         = diag(xxGE) %*% SSS %*% diag(xxGE) ## CSCS
      PSigma.GxEm       = P.mat.GxEm %*% Sigma.mat
      T.GxEm            = 1/2 * t(y) %*% PSigma.GxEm %*% P.mat.GxEm %*% y
      E.T.GxEm          = 1/2 * sum(diag(PSigma.GxEm))
      PS.GxEm           = P.mat.GxEm %*% SSS
      I.GxEm.phiphi     = 1/2 * get.sum.diag(PSigma.GxEm)
      I.GxEm.phitau     = 1/2 * get.sum.diag(PSigma.GxEm, PS.GxEm)
      I.GxEm.phisigma   = 1/2 * get.sum.diag(PSigma.GxEm, P.mat.GxEm)
      I.GxEm.tautau     = 1/2 * get.sum.diag(PS.GxEm)
      I.GxEm.tausigma   = 1/2 * get.sum.diag(PS.GxEm, P.mat.GxEm)
      I.GxEm.sigmasigma = 1/2 * get.sum.diag(P.mat.GxEm)
      I.GxEm.phitheta   = as.matrix(c(I.GxEm.phitau, I.GxEm.phisigma))
      I.GxEm.thetatheta = matrix(c(I.GxEm.tautau, I.GxEm.tausigma, I.GxEm.tausigma, I.GxEm.sigmasigma), 2,2)
      V.T.GxEm     = I.GxEm.phiphi - t(I.GxEm.phitheta) %*% solve(I.GxEm.thetatheta) %*% (I.GxEm.phitheta)
      ## ---- Gamma approximation ----
      aa.GxEm      = E.T.GxEm^2/V.T.GxEm
      bb.GxEm      = E.T.GxEm  /V.T.GxEm
      pval.GxEm    = 1-pgamma(T.GxEm, aa.GxEm, bb.GxEm)
      rm(P.mat.GxEm, PSigma.GxEm,PS.GxEm, I.GxEm.phiphi, I.GxEm.phitau,I.GxEm.phisigma, I.GxEm.tautau,I.GxEm.tausigma, I.GxEm.sigmasigma, I.GxEm.phitheta, I.GxEm.thetatheta)
      ##########################################
      ## (II) marginal G test
      ##########################################
      ## source("/Users/jytzeng/Research/HSregression/Software/hsreg.fun.r")
      ## source("/Users/jytzeng/Research/HSregression/Software/hsreg.conversion.r")
      source("hsreg.fun.r")
      source("hsreg.conversion.r")
      ## need to get sigma estimate
      pmat.data     = hsreg.pmat.fun(y=y, x=x.adj, trait.type=trait.type)
      P.mat       = pmat.data$pmat
      P.mat.class = pmat.data$pmat.class
      newy          = pmat.data$newy
      Vy            = pmat.data$v
      ## debug:
      ## tmp1 = (  diag(1,ngeno) - cbind(1, x.adj) %*% solve(t(cbind(1,x.adj) )%*% cbind(1,x.adj))%*%t(cbind(1,x.adj))  ) /P0$vy
      ## tmp2 = hsreg.pmat.fun(y=yy, x=xx, trait.type="gaussian")
      ## range(pmat.data$pmat - tmp1)      
      ##--- get newgeno (i.e., geno from "hsreg.prepare.geno.data") ----
#CSCS      geno.impute = haplo.to.impute(geno)
#CSCS      newgeno     = hsreg.prepare.geno.data(geno.impute, "impute", miss=NULL)
#CSCS      n.loci      = nrow(newgeno)
#CSCS      tvc.tmp     = tcrossprod(newy, newgeno)
#CSCS      T.G         = 1/n.loci * tcrossprod(tvc.tmp, tvc.tmp) 
      T.G = 1/2 * newy %*% SSS %*% newy;
      PS          = P.mat %*% SSS ## <--CSCS ps.calc.2(newgeno, P.mat)
      ##---- Gamma (2-moment) approximation
      trace.ps    = sum(diag(PS))
      E.T.G       = 1/2 * trace.ps
      Itautau     = 1/2 * get.sum.diag.square(PS)
      ##Itautau   = 1/2 * sum(diag(PS %*% PS)) 
      Itausigma   = 1/2 * trace.ps / Vy
      Isigmasigma = 1/2 * sum(diag( P.mat )) / Vy
      V.T.G       = Itautau - Itausigma^2/Isigmasigma
      ## Get gamma.
      aa.G = E.T.G^2/V.T.G
      bb.G = E.T.G  /V.T.G
      pval.G = 1 - pgamma(T.G, aa.G, bb.G)
      ##########################################
      ## (III) Joint G and GxE test
      ##########################################
      PSigma       = P.mat %*% Sigma.mat
      trace.psigma = sum(diag(PSigma))
      T.GxEjoint  = 1/2 * t(y) %*% PSigma %*% P.mat %*% y
      ##-- get E(T.joint) and V(T.joint)
      E.T.GxEjoint = 1/2 * trace.psigma
      Iphiphi      = 1/2 * get.sum.diag(PSigma)
      Iphitau      = 1/2 * get.sum.diag(PSigma, PS)
      Iphisigma    = 1/2 * trace.psigma/Vy               ## = 1/2 * get.sum.diag(PSigma, P.mat)
      V.T.GxEjoint = Iphiphi - Iphisigma^2/Isigmasigma
      cov.Tg.Tgxe  =  (Iphitau - Iphisigma *Itausigma)/Isigmasigma
      ##--------
      wG=wGxE=1
      T.joint      = wGxE   *   T.GxEjoint + wG   *   T.G  # same as       1/2 * t(y) %*% P.mat %*% (Sigma.mat+SSS) %*% P.mat %*% y
      E.T.joint    = wGxE   * E.T.GxEjoint + wG   * E.T.G
      V.T.joint    = wGxE^2 * V.T.GxEjoint + wG^2 * V.T.G + 2 * wG*wGxE * cov.Tg.Tgxe
      ##-- Gamma (2-moment) approximation
      aa.joint = E.T.joint^2/V.T.joint
      bb.joint = E.T.joint  /V.T.joint
      pval.joint = 1 - pgamma(T.joint, aa.joint, bb.joint)
      #################################################
      ## (IV) Joint-weighted G and GxE test
      ##----------------------------------------------
      ##      T_{jointwt} = w_G * T_G + w_GxE * T_GxE
      ##      where w_i = E_i/V_i
      #################################################
      wG        = bb.G
      wGxE      = E.T.GxEjoint / V.T.GxEjoint
      T.jointwt   = wGxE   *   T.GxEjoint + wG   *   T.G  # same as       1/2 * t(y) %*% P.mat %*% (Sigma.mat+SSS) %*% P.mat %*% y
      E.T.jointwt = wGxE   * E.T.GxEjoint + wG   * E.T.G
      V.T.jointwt = wGxE^2 * V.T.GxEjoint + wG^2 * V.T.G + 2 * wG*wGxE * cov.Tg.Tgxe
      ##-- Gamma (2-moment) approximation
      aa.jointwt   = E.T.jointwt^2/V.T.jointwt
      bb.jointwt   = E.T.jointwt  /V.T.jointwt
      pval.jointwt = 1 - pgamma(T.jointwt, aa.jointwt, bb.jointwt)
      return(list(
                  "pval"=
                  c("tau"= tau.eme,    "sigma2"=sigma2.eme,
                    
                    "T.GxEm"  =T.GxEm,   "pval.GxEm"=pval.GxEm,
                    "E.T.GxEm"=E.T.GxEm, "V.T.GxEm" =V.T.GxEm,
                    "aa.GxEm" =aa.GxEm,  "bb.GxEm"  =bb.GxEm,
                    
                    "T.G"  =T.G,       "pval.G"    =pval.G,
                    "E.T.G"=E.T.G,     "V.T.G"     =V.T.G,
                    "aa.G" =aa.G,      "bb.G"      =bb.G,
                    
                    "T.joint"  =T.joint,    "pval.joint"=pval.joint,
                    "E.T.joint"=E.T.joint,  "V.T.joint" =V.T.joint,                    
                    "aa.joint" =aa.joint,   "bb.joint"  =bb.joint,

                    "T.jointwt"  =T.jointwt,    "pval.jointwt"=pval.jointwt,
                    "E.T.jointwt"=E.T.jointwt,  "V.T.jointwt" =V.T.jointwt,                    
                    "aa.jointwt" =aa.jointwt,   "bb.jointwt"  =bb.jointwt)
)
##       
##       "S.mat" =SSS,
##       "P.mat.GxEm"=P.mat.GxEm,
##       "Sigma.mat"=        Sigma.mat,         
##       "PSigma.GxEm"=      PSigma.GxEm,        
##       "PS.GxEm"=          PS.GxEm,            
##       "I.GxEm.phiphi"=    I.GxEm.phiphi,      
##       "I.GxEm.phitau"=    I.GxEm.phitau,      
##       "I.GxEm.phisigma"=  I.GxEm.phisigma,    
##       "I.GxEm.tautau"=    I.GxEm.tautau,      
##       "I.GxEm.tausigma"=  I.GxEm.tausigma,    
##       "I.GxEm.sigmasigma"=I.GxEm.sigmasigma,  
##       "I.GxEm.phitheta"=  I.GxEm.phitheta,    
##       "I.GxEm.thetatheta"=I.GxEm.thetatheta,
##       "P.mat"=      P.mat,
##       "PS"=         PS, 
##       "PSigma"=     PSigma,
##       "Iphiphi"=    Iphiphi,      
##       "Iphitau"=    Iphitau,      
##       "Iphisigma"=  Iphisigma,    
##       "Itautau"=    Itautau,      
##       "Itausigma"=  Itausigma,    
##       "Isigmasigma"=Isigmasigma)
             )
    }

