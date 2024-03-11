require(GJRM)

##############TEST CODE
source("common_functions.R")

# a. Simulation parameters
set.seed(1000);options(scipen=999);
#dist="NO";a=1; b=2; c=0.75; mu1=1; mu2=2; n=1000
#dist="GA";a=.25; b=1.75; c=NA; mu1=10; mu2=12; n=1000
dist="PO";a=NA; b=NA; c=2; mu1=1; mu2=2; n=1000
#a=1; b=1; c=0.75; mu1=1; mu2=2; n=1000

#dataset <- generateBivDist(a=.25, b=1.75, c=NA, mu1=10, mu2=12, n=1000,dist)
#dataset <- generateBivDist(a=1, b=2, c=0.75, mu1=1, mu2=2, n=1000,dist)
dataset <- generateBivDist(n,a,b,c,mu1,mu2,dist)

eq.mu.1 <- formula(random_variable~1)
eq.mu.2 <- formula(random_variable.1~1)
fl <- list(eq.mu.1, eq.mu.2)

gamma_c_mu1<-dataset[dataset$time==0,]
gamma_c_mu2<-dataset[dataset$time==1,]

if(dist=="NO"){margin_dist="N"}
if(dist=="GA"){margin_dist="GA"}
if(dist=="PO"){margin_dist="PO"}

model_copula_n<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
model_copula_n_l<-  gjrm_l(fl, margins = c(margin_dist,margin_dist) , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")

summary(model_copula_n)

#####################GJRM_L
gjrm_l <- function (formula, data = list(), weights = NULL, subset = NULL, 
                    copula = "N", copula2 = "N", margins, model, dof = 3, dof2 = 3, 
                    ordinal = FALSE, surv = FALSE, cens1 = NULL, cens2 = NULL, 
                    cens3 = NULL, dep.cens = FALSE, upperBt1 = NULL, upperBt2 = NULL, 
                    gamlssfit = FALSE, fp = FALSE, infl.fac = 1, rinit = 1, rmax = 100, 
                    iterlimsp = 50, tolsp = 0.0000001, gc.l = FALSE, parscale, 
                    extra.regI = "t", k1.tvc = 0, k2.tvc = 0, knots = NULL, penCor = "unpen", 
                    sp.penCor = 3, Chol = FALSE, gamma = 1, w.alasso = NULL, 
                    drop.unused.levels = TRUE, min.dn = 0.0000000000000000000000000000000000000001, 
                    min.pr = 0.0000000000000001, max.pr = 0.999999) 
{
  BivD <- copula
  BivD2 <- copula2
  Model <- model
  if (missing(margins)) 
    stop("You must choose the margins' values.")
  if (missing(Model)) 
    stop("You must choose a model type.")
  if (margins[1] == "PH" && surv == TRUE) 
    margins[1] <- "cloglog"
  if (margins[1] == "PO" && surv == TRUE) 
    margins[1] <- "logit"
  if (margins[2] == "PH" && surv == TRUE) 
    margins[2] <- "cloglog"
  if (margins[2] == "PO" && surv == TRUE) 
    margins[2] <- "logit"
  bl <- c("probit", "logit", "cloglog")
  v.rB1 <- upperBt1
  v.rB2 <- upperBt2
  if (Model == "ROY") {
    L <- eval(substitute(SemiParROY(formula, data, weights, 
                                    subset, BivD1 = BivD, BivD2, margins, dof1 = dof, 
                                    dof2, gamlssfit, fp, infl.fac, rinit, rmax, iterlimsp, 
                                    tolsp, gc.l, parscale, extra.regI, knots = knots, 
                                    drop.unused.levels = drop.unused.levels, min.dn = min.dn, 
                                    min.pr = min.pr, max.pr = max.pr), list(weights = weights)))
  }
  else {
    #################################MAIN FUNCTION OF INTEREST###########################################################
    if (surv == FALSE && ordinal == FALSE) {
      if ((margins[1] %in% bl && margins[2] %in% bl && 
           is.na(margins[3])) || (margins[1] %in% bl && 
                                  !(margins[2] %in% bl) && Model == "B" && is.na(margins[3]))) {
        L <- eval(substitute(SemiParBIV_l(formula, data, 
                                        weights, subset, Model, BivD, margins, dof, 
                                        gamlssfit, fp, hess = TRUE, infl.fac, rinit, 
                                        rmax, iterlimsp, tolsp, gc.l, parscale, extra.regI, 
                                        intf = TRUE, theta.fx = NULL, knots = knots, 
                                        drop.unused.levels = drop.unused.levels, min.dn = min.dn, 
                                        min.pr = min.pr, max.pr = max.pr), list(weights = weights)))
      }
    }
    ####################################################################################
    if (surv == FALSE && ordinal == TRUE) {
      if ((margins[1] %in% bl && margins[2] %in% bl && 
           is.na(margins[3])) || (margins[1] %in% bl && 
                                  !(margins[2] %in% bl) && is.na(margins[3]))) {
        L <- eval(substitute(CopulaCLM(formula, data, 
                                       weights, subset, Model, BivD, margins, dof, 
                                       gamlssfit, fp, hess = TRUE, infl.fac, rinit, 
                                       rmax, iterlimsp, tolsp, gc.l, parscale, extra.regI, 
                                       intf = TRUE, theta.fx = NULL, knots = knots, 
                                       drop.unused.levels = drop.unused.levels, min.dn = min.dn, 
                                       min.pr = min.pr, max.pr = max.pr), list(weights = weights)))
      }
    }
    if (margins[1] %in% bl && !(margins[2] %in% bl) && surv == 
        FALSE && is.na(margins[3]) && Model == "BSS" && ordinal == 
        FALSE) {
      L <- eval(substitute(copulaSampleSel(formula, data, 
                                           weights, subset, BivD, margins, dof, fp, infl.fac, 
                                           rinit, rmax, iterlimsp, tolsp, gc.l, parscale, 
                                           extra.regI, knots, drop.unused.levels = drop.unused.levels, 
                                           min.dn = min.dn, min.pr = min.pr, max.pr = max.pr), 
                           list(weights = weights)))
    }
    if (!is.na(margins[3])) {
      if (margins[1] %in% bl && margins[2] %in% bl && margins[3] %in% 
          bl && surv == FALSE && ordinal == FALSE) {
        L <- eval(substitute(SemiParTRIV(formula, data, 
                                         weights, subset, Model, margins, penCor, sp.penCor, 
                                         approx = FALSE, Chol, infl.fac, gamma, w.alasso, 
                                         rinit, rmax, iterlimsp, tolsp, gc.l, parscale, 
                                         extra.regI, knots, drop.unused.levels = drop.unused.levels, 
                                         min.dn = min.dn, min.pr = min.pr, max.pr = max.pr), 
                             list(weights = weights)))
      }
    }
    if ((!(margins[1] %in% bl) || surv == TRUE) && ordinal == 
        FALSE) {
      robust <- FALSE
      t.c = 3
      sp <- qu.mag <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss1 <- gamlss2 <- gam1 <- gam2 <- y1m <- y2m <- indexTeq1 <- indexTeq2 <- NULL
      i.rho <- log.sig2.2 <- log.nu.2 <- log.nu.1 <- log.sig2.1 <- dof.st <- NULL
      end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- l.sp9 <- 0
      sp1 <- sp2 <- NULL
      sp3 <- gp3 <- gam3 <- X3 <- sp4 <- gp4 <- gam4 <- X4 <- sp5 <- gp5 <- gam5 <- X5 <- gam9 <- NULL
      sp6 <- gp6 <- gam6 <- X6 <- sp7 <- gp7 <- gam7 <- X7 <- sp8 <- gp8 <- gam8 <- X8 <- sp9 <- NULL
      c11 <- c10 <- c01 <- c00 <- NA
      cens1Mix <- cens2Mix <- NULL
      Sl.sf <- NULL
      sp.method <- "perf"
      Xd1 <- Xd2 <- mono.sm.pos1 <- mono.sm.pos2 <- mono.sm.pos <- NULL
      surv.flex <- FALSE
      Deq1 <- pos.pbeq1 <- Deq2 <- pos.pbeq2 <- list()
      BivD2 <- c("C0C90", "C0C270", "C180C90", "C180C270", 
                 "J0J90", "J0J270", "J180J90", "J180J270", "G0G90", 
                 "G0G270", "G180G90", "G180G270", "GAL0GAL90", 
                 "GAL0GAL270", "GAL180GAL90", "GAL180GAL270")
      opc <- c("N", "C0", "C90", "C180", "C270", "J0", 
               "J90", "J180", "J270", "G0", "G90", "G180", "G270", 
               "F", "AMH", "FGM", "T", "PL", "HO", "GAL0", "GAL90", 
               "GAL180", "GAL270")
      scc <- c("C0", "C180", "GAL0", "GAL180", "J0", "J180", 
               "G0", "G180", BivD2)
      sccn <- c("C90", "C270", "GAL90", "GAL270", "J90", 
                "J270", "G90", "G270")
      m2 <- c("N", "GU", "rGU", "LO", "LN", "WEI", "iG", 
              "GA", "BE", "FISK", "GP", "GPII", "GPo")
      m3 <- c("DAGUM", "SM", "TW")
      m1d <- c("PO", "ZTP", "DGP0")
      m2d <- c("NBI", "NBII", "PIG", "DGP", "DGPII")
      m3d <- c("DEL", "SICHEL")
      ct <- data.frame(c(opc), c(1:14, 55, 56, 57, 60, 
                                 61, 62:65))
      cta <- data.frame(c(opc), c(1, 3, 23, 13, 33, 6, 
                                  26, 16, 36, 4, 24, 14, 34, 5, 55, 56, 2, 60, 
                                  61, 62:65))
      if (BivD %in% BivD2) {
        if (BivD %in% BivD2[1:4]) 
          BivDt <- "C0"
        if (BivD %in% BivD2[5:12]) 
          BivDt <- "J0"
        if (BivD %in% BivD2[13:16]) 
          BivDt <- "C0"
        nC <- ct[which(ct[, 1] == BivDt), 2]
        nCa <- cta[which(cta[, 1] == BivDt), 2]
      }
      if (!(BivD %in% BivD2)) {
        nC <- ct[which(ct[, 1] == BivD), 2]
        nCa <- cta[which(cta[, 1] == BivD), 2]
      }
      if (!is.list(formula)) 
        stop("You must specify a list of equations.")
      l.flist <- length(formula)
      GJRM:::form.check(formula, l.flist)
      cl <- match.call()
      mf <- match.call(expand.dots = FALSE)
      pred.varR <- GJRM:::pred.var(formula, l.flist)
      v1 <- pred.varR$v1
      v2 <- pred.varR$v2
      pred.n <- pred.varR$pred.n
      if (!is.null(v.rB1)) 
        pred.n <- c(pred.n, v.rB1)
      if (!is.null(v.rB2)) 
        pred.n <- c(pred.n, v.rB2)
      fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + "))
      environment(fake.formula) <- environment(formula[[1]])
      mf$formula <- fake.formula
      mf$upperBt1 <- mf$upperBt2 <- mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$dep.cens <- mf$ordinal <- mf$Model <- mf$model <- mf$knots <- mf$k1.tvc <- mf$k2.tvc <- mf$surv <- mf$BivD <- mf$copula <- mf$copula2 <- mf$margins <- mf$fp <- mf$dof <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL
      mf$drop.unused.levels <- drop.unused.levels
      if (surv == TRUE) 
        mf$na.action <- na.pass
      mf[[1]] <- as.name("model.frame")
      data <- eval(mf, parent.frame())
      if (surv == TRUE) {
        if (!("(cens1)" %in% names(data)) && margins[1] %in% 
            bl) 
          stop("You must provide both censoring indicators.")
        if (!("(cens2)" %in% names(data)) && margins[2] %in% 
            bl) 
          stop("You must provide both censoring indicators.")
      }
      if (gc.l == TRUE) 
        gc()
      if (!("(weights)" %in% names(data))) {
        weights <- rep(1, dim(data)[1])
        data$weights <- weights
        names(data)[length(names(data))] <- "(weights)"
      }
      else weights <- data[, "(weights)"]
      if (!("(cens1)" %in% names(data))) {
        cens1 <- rep(0, dim(data)[1])
        data$cens1 <- cens1
        names(data)[length(names(data))] <- "(cens1)"
      }
      else cens1 <- data[, "(cens1)"]
      if (!("(cens2)" %in% names(data))) {
        cens2 <- rep(0, dim(data)[1])
        data$cens2 <- cens2
        names(data)[length(names(data))] <- "(cens2)"
      }
      else cens2 <- data[, "(cens2)"]
      if (!("(cens3)" %in% names(data))) {
        cens3 <- rep(0, dim(data)[1])
        data$cens3 <- cens3
        names(data)[length(names(data))] <- "(cens3)"
      }
      else cens3 <- data[, "(cens3)"]
      if (surv == TRUE) {
        if (is.factor(cens1) && !is.factor(cens2)) 
          stop("Both censoring indicators have to be factor variables for mixed censoring case.")
        if (!is.factor(cens1) && is.factor(cens2)) 
          stop("Both censoring indicators have to be factor variables for mixed censoring case.")
      }
      if (surv == TRUE && is.factor(cens1) && !is.null(v.rB1)) 
        data[!(cens1 == "I"), v.rB1] <- data[!(cens1 == 
                                                 "I"), v1[1]]
      if (surv == TRUE && is.factor(cens2) && !is.null(v.rB2)) 
        data[!(cens2 == "I"), v.rB2] <- data[!(cens2 == 
                                                 "I"), v2[1]]
      if (surv == TRUE) {
        if (any(is.na(data[, v1[1]]) | is.na(data[, v2[1]]))) 
          stop("Time to event with NA's. Please check your time covariates.")
        actual.NAs = as.numeric(which(apply(apply(data, 
                                                  1, is.na), 2, any)))
        data <- na.omit(data)
        if (length(actual.NAs) > 0) {
          cens1 <- cens1[-actual.NAs]
          cens2 <- cens2[-actual.NAs]
        }
      }
      n <- dim(data)[1]
      if (surv == TRUE && is.factor(cens1) && is.factor(cens2)) {
        cens1Mix <- cens1
        cens2Mix <- cens2
        cens1 <- cens2 <- rep(1, n)
      }
      M <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, 
                m3d = m3d, BivD = BivD, bl = bl, robust = robust, 
                opc = opc, extra.regI = extra.regI, margins = margins, 
                BivD2 = BivD2, dof = dof, surv = surv, c1 = cens1, 
                c2 = cens2, c3 = cens3, dep.cens = dep.cens)
      M$K1 <- NULL
      M$type.cens1 <- M$type.cens2 <- "R"
      pream.wm(formula, margins, M, l.flist)
      formula.eq1 <- formula[[1]]
      formula.eq2 <- formula[[2]]
      form.eq12R <- form.eq12(formula.eq1, data, v1, margins[1], 
                              m1d, m2d)
      formula.eq1 <- form.eq12R$formula.eq1
      formula.eq1r <- form.eq12R$formula.eq1r
      y1 <- form.eq12R$y1
      y1.test <- form.eq12R$y1.test
      y1m <- form.eq12R$y1m
      if (surv == FALSE) 
        gam1 <- eval(substitute(gam(formula.eq1, gamma = infl.fac, 
                                    weights = weights, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), list(weights = weights)))
      if (surv == TRUE && margins[1] %in% c(m2, m3) && 
          margins[2] %in% bl) 
        gam1 <- eval(substitute(gam(formula.eq1, gamma = infl.fac, 
                                    weights = weights, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), list(weights = weights)))
      else {
        if (surv == TRUE && !(margins[1] %in% bl)) 
          gam1 <- eval(substitute(gam(formula.eq1, gamma = infl.fac, 
                                      weights = weights * cens1, data = data, knots = knots, 
                                      drop.unused.levels = drop.unused.levels), 
                                  list(weights = weights, cens1 = cens1)))
      }
      if (surv == TRUE && margins[1] %in% bl) {
        surv.flex <- TRUE
        f.eq1 <- form.eq12R$f.eq1
        data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
        tempb <- eval(substitute(gam(f.eq1, family = cox.ph(), 
                                     data = data, weights = cens1, drop.unused.levels = drop.unused.levels), 
                                 list(cens1 = cens1)))
        data$Sh <- as.vector(mm(predict(tempb, type = "response"), 
                                min.pr = min.pr, max.pr = max.pr))
        cens11 <- ifelse(cens1 == 0, 0.0000001, cens1)
        gam1 <- eval(substitute(scam(formula.eq1, gamma = infl.fac, 
                                     weights = weights * cens11, data = data), list(weights = weights, 
                                                                                    cens11 = cens11)))
        lsgam1 <- length(gam1$smooth)
        if (lsgam1 == 0) 
          stop("You must use at least a monotonic smooth function of time in the first equation.")
        clsm <- ggr <- NA
        for (i in 1:lsgam1) {
          clsm[i] <- class(gam1$smooth[[i]])[1]
        }
        if (sum(as.numeric(clsm[1] %in% c("mpi.smooth"))) == 
            0) 
          stop("You must have a monotonic smooth of time and it has to be the first to be included.")
        l.sp1 <- length(gam1$sp)
        if (l.sp1 != 0) 
          sp1 <- gam1$sp
        sp1[1] <- 1
        gam.call <- gam1$call
        gam.call$sp <- sp1
        gam1 <- eval(gam.call)
        j <- 1
        for (i in 1:lsgam1) {
          if (max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$term))) != 
              0 && clsm[i] == "mpi.smooth") 
            mono.sm.pos1 <- c(mono.sm.pos1, c(gam1$smooth[[i]]$first.para:gam1$smooth[[i]]$last.para))
        }
        X1 <- predict(gam1, type = "lpmatrix")
        if (!is.null(indexTeq1) && k1.tvc != 0) {
          if (range(X1[, indexTeq1])[1] < 0) 
            stop("Check design matrix for smooth(s) of tvc term(s) in eq. 1.")
        }
        Xd1 <- Xdpred(gam1, data, v1[1])
        gam1$y <- data[, v1[1]]
        st.v1 <- c(gam1$coefficients)
        if (!is.null(indexTeq1)) {
          st.v1[mono.sm.pos1] <- exp(st.v1[mono.sm.pos1])
          while (range(Xd1 %*% st.v1)[1] < 0) st.v1[indexTeq1] <- 0.999 * 
              st.v1[indexTeq1]
          gam1$coefficients <- gam1$coefficients.t <- st.v1
          gam1$coefficients.t[mono.sm.pos1] <- exp(gam1$coefficients.t[mono.sm.pos1])
        }
      }
      gam1$formula <- formula.eq1r
      lsgam1 <- length(gam1$smooth)
      y1 <- y1.test
      if (margins[1] %in% c("LN")) 
        y1 <- log(y1)
      attr(data, "terms") <- NULL
      if (!(surv == TRUE && margins[1] %in% bl)) {
        names(gam1$model)[1] <- as.character(formula.eq1r[2])
        X1 <- predict(gam1, type = "lpmatrix")
        l.sp1 <- length(gam1$sp)
        sp1 <- gam1$sp
      }
      gp1 <- gam1$nsdf
      X1.d2 <- dim(X1)[2]
      form.eq12R <- form.eq12(formula.eq2, data, v2, margins[2], 
                              m1d, m2d)
      formula.eq2 <- form.eq12R$formula.eq1
      formula.eq2r <- form.eq12R$formula.eq1r
      y2 <- form.eq12R$y1
      y2.test <- form.eq12R$y1.test
      y2m <- form.eq12R$y1m
      if (surv == FALSE) 
        gam2 <- eval(substitute(gam(formula.eq2, gamma = infl.fac, 
                                    weights = weights, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), list(weights = weights)))
      if (surv == TRUE && !(margins[2] %in% bl)) 
        gam2 <- eval(substitute(gam(formula.eq2, gamma = infl.fac, 
                                    weights = weights * cens2, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), list(weights = weights, 
                                                                                   cens2 = cens2)))
      if (surv == TRUE && margins[2] %in% bl) {
        surv.flex <- TRUE
        f.eq2 <- form.eq12R$f.eq1
        data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
        tempb <- eval(substitute(gam(f.eq2, family = cox.ph(), 
                                     data = data, weights = cens2, drop.unused.levels = drop.unused.levels), 
                                 list(cens2 = cens2)))
        data$Sh <- as.vector(mm(predict(tempb, type = "response"), 
                                min.pr = min.pr, max.pr = max.pr))
        cens22 <- ifelse(cens2 == 0, 0.0000001, cens2)
        gam2 <- eval(substitute(scam(formula.eq2, gamma = infl.fac, 
                                     weights = weights * cens22, data = data), list(weights = weights, 
                                                                                    cens22 = cens22)))
        lsgam2 <- length(gam2$smooth)
        if (lsgam2 == 0) 
          stop("You must use at least a monotonic smooth function of time in the second equation.")
        clsm <- ggr <- NA
        for (i in 1:lsgam2) {
          clsm[i] <- class(gam2$smooth[[i]])[1]
        }
        if (sum(as.numeric(clsm[1] %in% c("mpi.smooth"))) == 
            0) 
          stop("You must have a monotonic smooth of time and it has to be the first to be included.")
        l.sp2 <- length(gam2$sp)
        if (l.sp2 != 0) 
          sp2 <- gam2$sp
        sp2[1] <- 1
        gam.call <- gam2$call
        gam.call$sp <- sp2
        gam2 <- eval(gam.call)
        j <- 1
        for (i in 1:lsgam2) {
          if (max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$term))) != 
              0 && clsm[i] == "mpi.smooth") 
            mono.sm.pos2 <- c(mono.sm.pos2, c(gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para))
        }
        X2 <- predict(gam2, type = "lpmatrix")
        if (!is.null(indexTeq2) && k2.tvc != 0) {
          if (range(X2[, indexTeq2])[1] < 0) 
            stop("Check design matrix for smooth(s) of tvc term(s) in eq. 2.")
        }
        Xd2 <- Xdpred(gam2, data, v2[1])
        gam2$y <- data[, v2[1]]
        st.v2 <- c(gam2$coefficients)
        if (!is.null(indexTeq2)) {
          st.v2[mono.sm.pos2] <- exp(st.v2[mono.sm.pos2])
          while (range(Xd2 %*% st.v2)[1] < 0) st.v2[indexTeq2] <- 0.999 * 
              st.v2[indexTeq2]
          gam2$coefficients <- gam2$coefficients.t <- st.v2
          gam2$coefficients.t[mono.sm.pos2] <- exp(gam2$coefficients.t[mono.sm.pos2])
        }
      }
      gam2$formula <- formula.eq2r
      lsgam2 <- length(gam2$smooth)
      y2 <- y2.test
      if (margins[2] %in% c("LN")) 
        y2 <- log(y2)
      attr(data, "terms") <- NULL
      if (!(surv == TRUE && margins[2] %in% bl)) {
        names(gam2$model)[1] <- as.character(formula.eq2r[2])
        X2 <- predict(gam2, type = "lpmatrix")
        l.sp2 <- length(gam2$sp)
        sp2 <- gam2$sp
      }
      gp2 <- gam2$nsdf
      X2.d2 <- dim(X2)[2]
      res1 <- residuals(gam1)
      res2 <- residuals(gam2)
      ass.s <- cor(res1, res2, method = "kendall")
      ass.s <- sign(ass.s) * ifelse(abs(ass.s) > 0.9, 0.9, 
                                    abs(ass.s))
      i.rho <- ass.dp(ass.s, BivD, scc, sccn, nCa)
      dof.st <- log(dof - 2)
      names(dof.st) <- "dof.star"
      if (!(margins[1] %in% c(m1d, bl))) {
        start.snR <- startsn(margins[1], y1)
        log.sig2.1 <- start.snR$log.sig2.1
        names(log.sig2.1) <- "sigma1.star"
        if (margins[1] %in% c(m3)) {
          log.nu.1 <- start.snR$log.nu.1
          names(log.nu.1) <- "nu.1.star"
        }
      }
      if (!(margins[2] %in% c(m1d, bl))) {
        start.snR <- startsn(margins[2], y2)
        log.sig2.2 <- start.snR$log.sig2.1
        names(log.sig2.2) <- "sigma2.star"
        if (margins[2] %in% c(m3)) {
          log.nu.2 <- start.snR$log.nu.1
          names(log.nu.2) <- "nu.2.star"
        }
      }
      vo <- list(gam1 = gam1, gam2 = gam2, i.rho = i.rho, 
                 log.sig2.2 = log.sig2.2, log.nu.2 = log.nu.2, 
                 log.nu.1 = log.nu.1, log.sig2.1 = log.sig2.1, 
                 dof.st = dof.st, n = n, drop.unused.levels = drop.unused.levels)
      start.v <- overall.sv(margins, M, vo)
      if (l.flist > 2) {
        overall.svGR <- overall.svG(formula, data, ngc = 2, 
                                    margins, M, vo, gam1, gam2, knots = knots)
        start.v = overall.svGR$start.v
        X3 = overall.svGR$X3
        X4 = overall.svGR$X4
        X5 = overall.svGR$X5
        X6 = overall.svGR$X6
        X7 = overall.svGR$X7
        X8 = overall.svGR$X8
        X3.d2 = overall.svGR$X3.d2
        X4.d2 = overall.svGR$X4.d2
        X5.d2 = overall.svGR$X5.d2
        X6.d2 = overall.svGR$X6.d2
        X7.d2 = overall.svGR$X7.d2
        X8.d2 = overall.svGR$X8.d2
        gp3 = overall.svGR$gp3
        gp4 = overall.svGR$gp4
        gp5 = overall.svGR$gp5
        gp6 = overall.svGR$gp6
        gp7 = overall.svGR$gp7
        gp8 = overall.svGR$gp8
        gam3 = overall.svGR$gam3
        gam4 = overall.svGR$gam4
        gam5 = overall.svGR$gam5
        gam6 = overall.svGR$gam6
        gam7 = overall.svGR$gam7
        gam8 = overall.svGR$gam8
        l.sp3 = overall.svGR$l.sp3
        l.sp4 = overall.svGR$l.sp4
        l.sp5 = overall.svGR$l.sp5
        l.sp6 = overall.svGR$l.sp6
        l.sp7 = overall.svGR$l.sp7
        l.sp8 = overall.svGR$l.sp8
        sp3 = overall.svGR$sp3
        sp4 = overall.svGR$sp4
        sp5 = overall.svGR$sp5
        sp6 = overall.svGR$sp6
        sp7 = overall.svGR$sp7
        sp8 = overall.svGR$sp8
      }
      GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, 
                  gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7, 
                  gam8 = gam8, gam9 = gam9)
      if ((l.sp1 != 0 || l.sp2 != 0 || l.sp3 != 0 || l.sp4 != 
           0 || l.sp5 != 0 || l.sp6 != 0 || l.sp7 != 0 || 
           l.sp8 != 0) && fp == FALSE) {
        L.GAM <- list(l.gam1 = length(gam1$coefficients), 
                      l.gam2 = length(gam2$coefficients), l.gam3 = length(gam3$coefficients), 
                      l.gam4 = length(gam4$coefficients), l.gam5 = length(gam5$coefficients), 
                      l.gam6 = length(gam6$coefficients), l.gam7 = length(gam7$coefficients), 
                      l.gam8 = length(gam8$coefficients), l.gam9 = 0)
        L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                     l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
                     l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = l.sp9)
        sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, 
                sp9)
        qu.mag <- S.m(GAM, L.SP, L.GAM)
      }
      if (missing(parscale)) 
        parscale <- 1
      respvec <- respvec2 <- respvec3 <- list(y1 = y1, 
                                              y2 = y2, y1.y2 = NULL, y1.cy2 = NULL, cy1.y2 = NULL, 
                                              cy1.cy2 = NULL, cy1 = NULL, cy = NULL, univ = 0)
      my.env <- new.env()
      my.env$signind <- 1
      lsgam3 <- length(gam3$smooth)
      lsgam4 <- length(gam4$smooth)
      lsgam5 <- length(gam5$smooth)
      lsgam6 <- length(gam6$smooth)
      lsgam7 <- length(gam7$smooth)
      lsgam8 <- length(gam8$smooth)
      lsgam9 <- length(gam9$smooth)
      indUR <- indUL <- indUI <- indUU <- indRR <- indRL <- indRI <- indRU <- indLR <- indLL <- indLI <- indLU <- indIR <- indIL <- indII <- indIU <- rep(0, 
                                                                                                                                                          n)
      if (surv == TRUE && dep.cens == FALSE) {
        if ((surv == TRUE && margins[1] %in% bl && margins[2] %in% 
             bl && !is.factor(cens1) && !is.factor(cens2)) || 
            (surv == TRUE && margins[1] %in% m2 && margins[2] %in% 
             m2)) {
          c11 <- cens1 * cens2
          c10 <- cens1 * (1 - cens2)
          c01 <- (1 - cens1) * cens2
          c00 <- (1 - cens1) * (1 - cens2)
        }
        if (surv == TRUE && margins[1] %in% c(m2, m3) && 
            margins[2] %in% bl) {
          c11 <- cens2
          c10 <- 1 - cens2
          c01 <- NULL
          c00 <- NULL
        }
        if (!is.null(cens1Mix) && !is.null(cens2Mix)) {
          if (surv == TRUE && margins[1] %in% bl && margins[2] %in% 
              bl && is.factor(cens1Mix) && is.factor(cens2Mix)) {
            gamlssfit <- TRUE
            indUR <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "R")
            indUL <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "L")
            indUI <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "I")
            indUU <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "U")
            indRR <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "R")
            indRL <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "L")
            indRI <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "I")
            indRU <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "U")
            indLR <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "R")
            indLL <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "L")
            indLI <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "I")
            indLU <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "U")
            indIR <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "R")
            indIL <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "L")
            indII <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "I")
            indIU <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "U")
          }
        }
      }
      if (surv == TRUE && dep.cens == TRUE) {
        c11 <- NULL
        c10 <- cens1
        c01 <- cens2
        c00 <- cens3
      }
      my.env$k1 <- k1.tvc
      my.env$k2 <- k2.tvc
      VC <- list(lsgam1 = lsgam1, indexTeq1 = indexTeq1, 
                 indexTeq2 = indexTeq2, lsgam2 = lsgam2, Deq1 = Deq1, 
                 pos.pbeq1 = pos.pbeq1, Deq2 = Deq2, pos.pbeq2 = pos.pbeq2, 
                 lsgam3 = lsgam3, robust = FALSE, sp.fixed = NULL, 
                 lsgam4 = lsgam4, Sl.sf = Sl.sf, sp.method = sp.method, 
                 lsgam5 = lsgam5, K1 = NULL, lsgam6 = lsgam6, 
                 lsgam7 = lsgam7, lsgam8 = lsgam8, lsgam9 = lsgam9, 
                 X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, 
                 X6 = X6, X7 = X7, X8 = X8, X1.d2 = X1.d2, X2.d2 = X2.d2, 
                 X3.d2 = X3.d2, X4.d2 = X4.d2, X5.d2 = X5.d2, 
                 X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2, 
                 gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, 
                 gp6 = gp6, gp7 = gp7, gp8 = gp8, l.sp1 = l.sp1, 
                 l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, 
                 l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, 
                 l.sp8 = l.sp8, l.sp9 = 0, my.env = my.env, infl.fac = infl.fac, 
                 weights = weights, fp = fp, gamlssfit = gamlssfit, 
                 hess = NULL, Model = "CC", univ.gamls = FALSE, 
                 model = model, end = end, BivD = BivD, nCa = nCa, 
                 copula = copula, copula2 = copula2, nC = nC, 
                 gc.l = gc.l, n = n, extra.regI = extra.regI, 
                 parscale = parscale, margins = margins, Cont = "YES", 
                 ccss = "no", m2 = m2, m3 = m3, m1d = m1d, m2d = m2d, 
                 m3d = m3d, bl = bl, triv = FALSE, y1m = y1m, 
                 y2m = y2m, tc = t.c, i.rho = i.rho, dof = dof, 
                 dof.st = dof.st, BivD2 = BivD2, cta = cta, ct = ct, 
                 zerov = -10, c11 = c11, c10 = c10, c01 = c01, 
                 c00 = c00, indUR = indUR, indUL = indUL, indUI = indUI, 
                 indUU = indUU, indRR = indRR, indRL = indRL, 
                 indRI = indRI, indRU = indRU, indLR = indLR, 
                 indLL = indLL, indLI = indLI, indLU = indLU, 
                 indIR = indIR, indIL = indIL, indII = indII, 
                 indIU = indIU, surv = surv, Xd1 = Xd1, Xd2 = Xd2, 
                 mono.sm.pos1 = mono.sm.pos1, mono.sm.pos2 = mono.sm.pos2, 
                 surv.flex = surv.flex, mono.sm.pos = mono.sm.pos, 
                 gp2.inf = NULL, informative = "no", zero.tol = 0.01, 
                 min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
      if (gc.l == TRUE) 
        gc()
      if (gamlssfit == TRUE) {
        type.cens1 <- type.cens2 <- "R"
        surv1 <- surv2 <- surv
        form.gamlR <- form.gaml(formula, l.flist, M)
        if (surv == TRUE && margins[1] %in% c(m2, m3) && 
            margins[2] %in% bl) 
          surv1 <- FALSE
        if (surv == TRUE && margins[1] %in% bl && margins[2] %in% 
            bl && is.factor(cens1Mix) && is.factor(cens2Mix)) {
          cens1 <- cens1Mix
          cens2 <- cens2Mix
          type.cens1 <- type.cens2 <- "mixed"
          M$type.cens1 = type.cens1
          M$type.cens2 = type.cens2
        }
        gamlss1 <- eval(substitute(gamlss(form.gamlR$formula.gamlss1, 
                                          data = data, weights = weights, subset = subset, 
                                          margin = margins[1], surv = surv1, cens = cens1, 
                                          type.cens = type.cens1, upperB = upperBt1, 
                                          infl.fac = infl.fac, rinit = rinit, rmax = rmax, 
                                          iterlimsp = iterlimsp, tolsp = tolsp, gc.l = gc.l, 
                                          parscale = 1, extra.regI = extra.regI, k.tvc = k1.tvc, 
                                          drop.unused.levels = drop.unused.levels), list(weights = weights, 
                                                                                         cens1 = cens1)))
        gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, 
                                          data = data, weights = weights, subset = subset, 
                                          margin = margins[2], surv = surv2, cens = cens2, 
                                          type.cens = type.cens2, upperB = upperBt2, 
                                          infl.fac = infl.fac, rinit = rinit, rmax = rmax, 
                                          iterlimsp = iterlimsp, tolsp = tolsp, gc.l = gc.l, 
                                          parscale = 1, extra.regI = extra.regI, k.tvc = k2.tvc, 
                                          drop.unused.levels = drop.unused.levels), list(weights = weights, 
                                                                                         cens2 = cens2)))
        SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4, 
                   sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8)
        gamls.upsvR <- gamls.upsv(gamlss1, gamlss2, margins, 
                                  M, l.flist, nstv = names(start.v), VC, GAM, 
                                  SP)
        sp <- gamls.upsvR$sp
        start.v <- gamls.upsvR$start.v
        VC$X1 <- gamlss1$VC$X1
        VC$Xd1 <- gamlss1$VC$Xd1
        VC$X1.2 <- gamlss1$VC$X2
        VC$X2 <- gamlss2$VC$X1
        VC$Xd2 <- gamlss2$VC$Xd1
        VC$X2.2 <- gamlss2$VC$X2
        rangeSurv1 <- gamlss1$rangeSurv
        rangeSurv2 <- gamlss2$rangeSurv
      }
      func.opt <- func.OPT_l(margins, M)
      SemiParFit <- SemiParBIV.fit(func.opt = func.opt, 
                                   start.v = start.v, rinit = rinit, rmax = rmax, 
                                   iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp, 
                                   respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag)
      SemiParFit.p <- copulaReg.fit.post(SemiParFit = SemiParFit, 
                                         VC = VC, GAM)
      y1.m <- y1
      if (margins[1] == "LN") 
        y1.m <- exp(y1)
      y2.m <- y2
      if (margins[2] == "LN") 
        y2.m <- exp(y2)
      SemiParFit <- SemiParFit.p$SemiParFit
      if (gc.l == TRUE) 
        gc()
      GJRM:::cov.c(SemiParFit)
      gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data
      L <- list(fit = SemiParFit$fit, dataset = NULL, n = n, 
                gamlss1 = gamlss1, gamlss2 = gamlss2, formula = formula, 
                robust = FALSE, edf11 = SemiParFit.p$edf11, surv = surv, 
                gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
                gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8, 
                coefficients = SemiParFit$fit$argument, coef.t = SemiParFit.p$coef.t, 
                iterlimsp = iterlimsp, weights = weights, cens1 = cens1, 
                cens2 = cens2, cens3 = cens3, sp = SemiParFit.p$sp, 
                iter.sp = SemiParFit$iter.sp, l.sp1 = l.sp1, 
                l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, 
                l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, 
                l.sp8 = l.sp8, bl = bl, l.sp9 = l.sp9, gam9 = gam9, 
                fp = fp, iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner, 
                theta = SemiParFit.p$theta, theta.a = SemiParFit.p$theta.a, 
                sigma21 = SemiParFit.p$sigma21, sigma22 = SemiParFit.p$sigma22, 
                sigma21.a = SemiParFit.p$sigma21.a, sigma22.a = SemiParFit.p$sigma22.a, 
                sigma1 = SemiParFit.p$sigma21, sigma2 = SemiParFit.p$sigma22, 
                sigma1.a = SemiParFit.p$sigma21.a, sigma2.a = SemiParFit.p$sigma22.a, 
                nu1 = SemiParFit.p$nu1, nu2 = SemiParFit.p$nu2, 
                nu1.a = SemiParFit.p$nu1.a, nu2.a = SemiParFit.p$nu2.a, 
                dof.a = SemiParFit.p$dof.a, dof = SemiParFit.p$dof, 
                X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, 
                X6 = X6, X7 = X7, X8 = X8, X1.d2 = X1.d2, X2.d2 = X2.d2, 
                X3.d2 = X3.d2, X4.d2 = X4.d2, X5.d2 = X5.d2, 
                X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2, 
                He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, 
                Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, F = SemiParFit.p$F, 
                F1 = SemiParFit.p$F1, t.edf = SemiParFit.p$t.edf, 
                edf = SemiParFit.p$edf, edf1 = SemiParFit.p$edf1, 
                edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3, 
                edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, 
                edf6 = SemiParFit.p$edf6, edf7 = SemiParFit.p$edf7, 
                edf8 = SemiParFit.p$edf8, edf1.1 = SemiParFit.p$edf1.1, 
                edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3, 
                edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, 
                edf1.6 = SemiParFit.p$edf1.6, edf1.7 = SemiParFit.p$edf1.7, 
                edf1.8 = SemiParFit.p$edf1.8, R = SemiParFit.p$R, 
                bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
                wor.c = SemiParFit$wor.c, eta1 = SemiParFit$fit$eta1, 
                eta2 = SemiParFit$fit$eta2, etad = SemiParFit$fit$etad, 
                etas1 = SemiParFit$fit$etas1, etas2 = SemiParFit$fit$etas2, 
                y1 = y1.m, y2 = y2.m, BivD = BivD, margins = margins, 
                copula = copula, copula2 = copula2, logLik = SemiParFit.p$logLik, 
                nC = nC, respvec = respvec, hess = TRUE, qu.mag = qu.mag, 
                gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, 
                gp6 = gp6, gp7 = gp7, gp8 = gp8, VC = VC, magpp = SemiParFit$magpp, 
                gamlssfit = gamlssfit, Cont = "YES", tau = SemiParFit.p$tau, 
                tau.a = SemiParFit.p$tau.a, l.flist = l.flist, 
                v1 = v1, v2 = v2, triv = FALSE, univar.gamlss = FALSE, 
                BivD2 = BivD2, call = cl, surv = surv, surv.flex = surv.flex, 
                Vb.t = SemiParFit.p$Vb.t, coef.t = SemiParFit.p$coef.t, 
                Model = "CC", model = model)
      if (BivD %in% BivD2) {
        L$teta1 <- SemiParFit$fit$teta1
        L$teta.ind1 <- SemiParFit$fit$teta.ind1
        L$teta2 <- SemiParFit$fit$teta2
        L$teta.ind2 <- SemiParFit$fit$teta.ind2
        L$Cop1 <- SemiParFit$fit$Cop1
        L$Cop2 <- SemiParFit$fit$Cop2
      }
      class(L) <- c("gjrm", "SemiParBIV")
    }
  }
  L
}


SemiParBIV_l <- function (formula, data = list(), weights = NULL, subset = NULL, 
                          Model = "B", BivD = "N", margins = c("probit", "probit"), 
                          dof = 3, gamlssfit = FALSE, fp = FALSE, hess = TRUE, infl.fac = 1, 
                          theta.fx = NULL, rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 0.0000001, 
                          gc.l = FALSE, parscale, extra.regI = "t", intf = FALSE, knots = NULL, 
                          drop.unused.levels = TRUE, min.dn = 0.0000000000000000000000000000000000000001, 
                          min.pr = 0.0000000000000001, max.pr = 0.999999) 
{
  i.rho <- sp <- qu.mag <- n.sel <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- inde <- y2m <- NULL
  y00 <- y10 <- y0p <- y1p <- gam2TW <- NULL
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- l.sp9 <- i.rho <- 0
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- gam8 <- gam9 <- gamlss2 <- dof.st <- NULL
  gamlss2 <- NULL
  sp1 <- sp2 <- c.gam2 <- X2s <- X3s <- NULL
  sp3 <- gp3 <- gam3 <- X3 <- NULL
  sp4 <- gp4 <- gam4 <- X4 <- NULL
  sp5 <- gp5 <- gam5 <- X5 <- NULL
  sp6 <- gp6 <- gam6 <- X6 <- NULL
  sp7 <- gp7 <- gam7 <- X7 <- NULL
  sp8 <- gp8 <- gam8 <- X8 <- NULL
  log.sig2 <- log.nu <- NULL
  Sl.sf <- NULL
  sp.method <- "perf"
  BivD2 <- c("C0C90", "C0C270", "C180C90", "C180C270", "J0J90", 
             "J0J270", "J180J90", "J180J270", "G0G90", "G0G270", "G180G90", 
             "G180G270", "GAL0GAL90", "GAL0GAL270", "GAL180GAL90", 
             "GAL180GAL270")
  opc <- c("N", "C0", "C90", "C180", "C270", "J0", "J90", "J180", 
           "J270", "G0", "G90", "G180", "G270", "F", "AMH", "FGM", 
           "T", "PL", "HO", "GAL0", "GAL90", "GAL180", "GAL270")
  scc <- c("C0", "C180", "GAL0", "GAL180", "J0", "J180", "G0", 
           "G180", BivD2)
  sccn <- c("C90", "C270", "GAL90", "GAL270", "J90", "J270", 
            "G90", "G270")
  mb <- c("B", "BSS", "BPO", "BPO0")
  m2 <- c("N", "GU", "rGU", "LO", "LN", "WEI", "iG", "GA", 
          "BE", "FISK", "GP", "GPII", "GPo")
  m3 <- c("DAGUM", "SM", "TW")
  m1d <- c("PO", "ZTP", "DGP0")
  m2d <- c("NBI", "NBII", "PIG", "DGP", "DGPII")
  bl <- c("probit", "logit", "cloglog")
  M <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, BivD = BivD, 
            opc = opc, extra.regI = extra.regI, margins = margins, 
            bl = bl, intf = intf, theta.fx = theta.fx, Model = Model, 
            mb = mb, BivD2 = BivD2, dof = dof)
  surv.flex <- FALSE
  M$K1 <- NULL
  ct <- data.frame(c(opc), c(1:14, 55, 56, 57, 60, 61, 62:65))
  cta <- data.frame(c(opc), c(1, 3, 23, 13, 33, 6, 26, 16, 
                              36, 4, 24, 14, 34, 5, 55, 56, 2, 60, 61, 62:65))
  if (BivD %in% BivD2) {
    if (BivD %in% BivD2[1:4]) 
      BivDt <- "C0"
    if (BivD %in% BivD2[5:12]) 
      BivDt <- "J0"
    if (BivD %in% BivD2[13:16]) 
      BivDt <- "C0"
    nC <- ct[which(ct[, 1] == BivDt), 2]
    nCa <- cta[which(cta[, 1] == BivDt), 2]
  }
  if (!(BivD %in% BivD2)) {
    nC <- ct[which(ct[, 1] == BivD), 2]
    nCa <- cta[which(cta[, 1] == BivD), 2]
  }
  if (!is.list(formula)) 
    stop("You must specify a list of equations.")
  M$l.flist <- l.flist <- length(formula)
  pream.wm(formula, margins, M, l.flist, type = "biv")
  GJRM:::form.check(formula, l.flist)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  pred.varR <- GJRM:::pred.var(formula, l.flist)
  v1 <- pred.varR$v1
  v2 <- pred.varR$v2
  pred.n <- pred.varR$pred.n
  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + "))
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula
  mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$ordinal <- mf$knots <- mf$dof <- mf$intf <- mf$theta.fx <- mf$Model <- mf$BivD <- mf$margins <- mf$fp <- mf$hess <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL
  mf$drop.unused.levels <- drop.unused.levels
  if (Model == "BSS") 
    mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  if (gc.l == TRUE) 
    gc()
  if (Model == "BSS") {
    data[is.na(data[, v1[1]]), v1[1]] <- 0
    indS <- data[, v1[1]]
    indS[is.na(indS)] <- 0
    indS <- as.logical(indS)
    data[indS == FALSE, v2[1]] <- 0
    data <- na.omit(data)
  }
  if (!("(weights)" %in% names(data))) {
    weights <- rep(1, dim(data)[1])
    data$weights <- weights
    names(data)[length(names(data))] <- "(weights)"
  }
  else weights <- data[, "(weights)"]
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]]
  if (Model == "B") {
    if (v1[1] %in% v2[-1]) 
      end <- 1
    if (v2[1] %in% v1[-1]) 
      end <- 2
  }
  gam1 <- eval(substitute(gam(formula.eq1, binomial(link = margins[1]), 
                              gamma = infl.fac, weights = weights, data = data, knots = knots, 
                              drop.unused.levels = drop.unused.levels), list(weights = weights)))
  X1 <- model.matrix(gam1)
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp)
  y1 <- gam1$y
  n <- length(y1)
  if (l.sp1 != 0) 
    sp1 <- gam1$sp
  gp1 <- gam1$nsdf
  inde <- rep(TRUE, n)
  if ((Model == "B" || Model == "BPO" || Model == "BPO0") && 
      margins[2] %in% bl) {
    gam2 <- eval(substitute(gam(formula.eq2, binomial(link = margins[2]), 
                                gamma = infl.fac, weights = weights, data = data, 
                                knots = knots, drop.unused.levels = drop.unused.levels), 
                            list(weights = weights)))
    X2 <- model.matrix(gam2)
    X2.d2 <- dim(X2)[2]
    l.sp2 <- length(gam2$sp)
    y2 <- gam2$y
    if (l.sp2 != 0) 
      sp2 <- gam2$sp
    if (Model == "B") {
      y1.y2 <- y1 * y2
      y1.cy2 <- y1 * (1 - y2)
      cy1.y2 <- (1 - y1) * y2
      cy1.cy2 <- (1 - y1) * (1 - y2)
    }
    if (Model == "BPO" || Model == "BPO0") 
      cy <- 1 - y1
  }
  if (Model == "B" && !(margins[2] %in% bl)) {
    form.eq12R <- form.eq12(formula.eq2, data, v2, margins[2], 
                            m1d, m2d)
    formula.eq2 <- form.eq12R$formula.eq1
    formula.eq2r <- form.eq12R$formula.eq1r
    y2 <- form.eq12R$y1
    y2.test <- form.eq12R$y1.test
    y2m <- form.eq12R$y1m
    gam2 <- eval(substitute(gam(formula.eq2, gamma = infl.fac, 
                                weights = weights, data = data, knots = knots, drop.unused.levels = drop.unused.levels), 
                            list(weights = weights)))
    gam2$formula <- formula.eq2r
    names(gam2$model)[1] <- as.character(formula.eq2r[2])
    y2 <- y2.test
    if (margins[2] %in% c("LN")) 
      y2 <- log(y2)
    X2 <- model.matrix(gam2)
    X2.d2 <- dim(X2)[2]
    l.sp2 <- length(gam2$sp)
    if (l.sp2 != 0) 
      sp2 <- gam2$sp
    if (margins[2] != "TW") 
      cy <- 1 - y1
    if (margins[2] == "TW") {
      formulamgcv <- formula[-c(1, 5)]
      formulamgcv[[2]] <- formula[-c(1, 5)][[3]]
      formulamgcv[[3]] <- formula[-c(1, 5)][[2]]
      gam2TW <- eval(substitute(gam(formulamgcv, gamma = infl.fac, 
                                    weights = weights, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels, family = twlss()), 
                                list(weights = weights)))
      y2b <- y2 > 0
      y00 <- (1 - y1) * (1 - y2b)
      y10 <- y1 * (1 - y2b)
      y0p <- (1 - y1) * y2b
      y1p <- y1 * y2b
    }
  }
  if (Model == "BSS") {
    inde <- as.logical(y1)
    gam2 <- eval(substitute(gam(formula.eq2, binomial(link = margins[2]), 
                                gamma = infl.fac, weights = weights, data = data, 
                                subset = inde, knots = knots, drop.unused.levels = drop.unused.levels), 
                            list(weights = weights, inde = inde)))
    X2s <- try(predict.gam(gam2, newdata = data[, -dim(data)[2]], 
                           type = "lpmatrix"), silent = TRUE)
    if (any(class(X2s) == "try-error")) 
      stop("Check that the factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.")
    X2.d2 <- length(gam2$coefficients)
    X2 <- model.matrix(gam2)
    y2 <- gam2$y
    n.sel <- sum(as.numeric(inde))
    l.sp2 <- length(gam2$sp)
    if (l.sp2 != 0) 
      sp2 <- gam2$sp
    cy1 <- (1 - y1)
    y1.y2 <- y1[inde] * y2
    y1.cy2 <- y1[inde] * (1 - y2)
    form.eq2imr <- update.formula(formula.eq2, ~. + imrGUANN)
    p.g1 <- predict.gam(gam1)
    imrGUANN <- data$imrGUANN <- dnorm(p.g1)/pnorm(p.g1)
    gam2.1 <- eval(substitute(gam(form.eq2imr, gamma = infl.fac, 
                                  binomial(link = margins[2]), weights = weights, data = data, 
                                  subset = inde, knots = knots, drop.unused.levels = drop.unused.levels), 
                              list(weights = weights, inde = inde)))
    pimr <- which(names(gam2.1$coefficients) == "imrGUANN")
    c.gam2 <- gam2.1$coefficients[-pimr]
    sia <- sqrt(mean(residuals(gam2.1, type = "deviance")^2) + 
                  mean(imrGUANN[inde] * (imrGUANN[inde] + p.g1[inde])) * 
                  gam2.1$coefficients["imrGUANN"]^2)[[1]]
    co <- (gam2.1$coefficients["imrGUANN"]/sia)[[1]]
    ass.s <- co
    ass.s <- sign(ass.s) * ifelse(abs(ass.s) > 0.2, 0.2, 
                                  abs(ass.s))
    if (l.sp2 != 0) 
      sp2 <- gam2.1$sp
  }
  gp2 <- gam2$nsdf
  if (is.null(theta.fx)) {
    if (!(Model %in% c("BPO0"))) {
      if (Model == "B") {
        res1 <- residuals(gam1)
        res2 <- residuals(gam2)
        ass.s <- cor(res1, res2, method = "kendall")
        ass.s <- sign(ass.s) * ifelse(abs(ass.s) > 0.9, 
                                      0.9, abs(ass.s))
      }
      if (Model == "BPO") 
        ass.s <- 0.01
      i.rho <- ass.dp(ass.s, BivD, scc, sccn, nCa)
    }
    if (Model == "BPO0") 
      i.rho <- 0
  }
  names(i.rho) <- "theta.star"
  if (margins[2] != "TW") {
    if (margins[1] %in% bl && margins[2] %in% c(bl, m1d) && 
        Model %in% c("B", "BPO") && is.null(theta.fx)) 
      start.v <- c(gam1$coefficients, gam2$coefficients, 
                   i.rho)
    if (Model == "BSS") 
      start.v <- c(gam1$coefficients, c.gam2, i.rho)
    if (Model == "BPO0") 
      start.v <- c(gam1$coefficients, gam2$coefficients)
    if (!is.null(theta.fx) && Model == "B" && margins[2] %in% 
        bl) 
      start.v <- c(gam1$coefficients, gam2$coefficients)
    if (margins[1] %in% bl && margins[2] %in% c(m2, m3, m2d)) {
      start.snR <- startsn(margins[2], y2)
      log.sig2 <- start.snR$log.sig2.1
      names(log.sig2) <- "sigma.star"
      if (margins[2] %in% c(m3)) {
        log.nu <- start.snR$log.nu.1
        names(log.nu) <- "nu.star"
      }
      if (margins[2] %in% c(m2, m2d)) 
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     log.sig2, i.rho)
      if (margins[2] %in% m3) 
        start.v <- c(gam1$coefficients, gam2$coefficients, 
                     log.sig2, log.nu, i.rho)
    }
  }
  if (l.flist > 2) {
    if (margins[2] == "TW") 
      log.nu <- log.sig2 <- 0.1
    vo <- list(gam1 = gam1, gam2 = gam2, i.rho = i.rho, log.sig2 = log.sig2, 
               log.nu = log.nu, n = n, drop.unused.levels = drop.unused.levels)
    overall.svGR <- overall.svG(formula, data, ngc = 2, margins, 
                                M, vo, gam1, gam2, type = "biv", inde = inde, c.gam2 = c.gam2, 
                                knots = knots)
    start.v = overall.svGR$start.v
    X3 = overall.svGR$X3
    X4 = overall.svGR$X4
    X5 = overall.svGR$X5
    X6 = overall.svGR$X6
    X7 = overall.svGR$X7
    X8 = overall.svGR$X8
    X3.d2 = overall.svGR$X3.d2
    X4.d2 = overall.svGR$X4.d2
    X5.d2 = overall.svGR$X5.d2
    X6.d2 = overall.svGR$X6.d2
    X7.d2 = overall.svGR$X7.d2
    X8.d2 = overall.svGR$X8.d2
    gp3 = overall.svGR$gp3
    gp4 = overall.svGR$gp4
    gp5 = overall.svGR$gp5
    gp6 = overall.svGR$gp6
    gp7 = overall.svGR$gp7
    gp8 = overall.svGR$gp8
    gam3 = overall.svGR$gam3
    gam4 = overall.svGR$gam4
    gam5 = overall.svGR$gam5
    gam6 = overall.svGR$gam6
    gam7 = overall.svGR$gam7
    gam8 = overall.svGR$gam8
    l.sp3 = overall.svGR$l.sp3
    l.sp4 = overall.svGR$l.sp4
    l.sp5 = overall.svGR$l.sp5
    l.sp6 = overall.svGR$l.sp6
    l.sp7 = overall.svGR$l.sp7
    l.sp8 = overall.svGR$l.sp8
    sp3 = overall.svGR$sp3
    sp4 = overall.svGR$sp4
    sp5 = overall.svGR$sp5
    sp6 = overall.svGR$sp6
    sp7 = overall.svGR$sp7
    sp8 = overall.svGR$sp8
    X3s = overall.svGR$X3s
    X4s = overall.svGR$X4s
    if (margins[2] == "TW") {
      gamlssfit <- FALSE
      nams <- names(start.v)
      start.v1TW <- start.v1 <- gam2TW$coefficients
      X2.d2mgcv <- X2.d2
      X3.d2mgcv <- X4.d2
      X4.d2mgcv <- X3.d2
      start.v1TW[1:X2.d2] <- start.v1[1:X2.d2mgcv]
      start.v1TW[(X2.d2 + 1):(X2.d2 + X3.d2)] <- start.v1[(X2.d2mgcv + 
                                                             X3.d2mgcv + 1):(X2.d2mgcv + X3.d2mgcv + X4.d2mgcv)]
      start.v1TW[(X2.d2 + X3.d2 + 1):(X2.d2 + X3.d2 + X4.d2)] <- start.v1[(X2.d2mgcv + 
                                                                             1):(X2.d2mgcv + X3.d2mgcv)]
      start.v1 <- start.v1TW
      start.v <- c(gam1$coefficients, start.v1, gam5$coefficients)
      names(start.v) <- nams
      rm(start.v1TW)
      spmgcv <- gam2TW$sp
      l.sp2mgcv <- l.sp2
      l.sp3mgcv <- l.sp4
      l.sp4mgcv <- l.sp3
      if (l.sp2 != 0) {
        sp2 <- spmgcv[1:l.sp2mgcv]
        names(sp2) <- names(gam2$sp)
      }
      if (l.sp3 != 0) {
        sp3 <- spmgcv[(l.sp2mgcv + l.sp3mgcv + 1):(l.sp2mgcv + 
                                                     l.sp3mgcv + l.sp4mgcv)]
        names(sp3) <- names(gam3$sp)
      }
      if (l.sp4 != 0) {
        sp4 <- spmgcv[(l.sp2mgcv + 1):(l.sp2mgcv + l.sp3mgcv)]
        names(sp4) <- names(gam4$sp)
      }
    }
  }
  GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
              gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8, gam9 = gam9)
  if ((l.sp1 != 0 || l.sp2 != 0 || l.sp3 != 0 || l.sp4 != 0 || 
       l.sp5 != 0 || l.sp6 != 0 || l.sp7 != 0 || l.sp8 != 0) && 
      fp == FALSE) {
    L.GAM <- list(l.gam1 = length(gam1$coefficients), l.gam2 = length(gam2$coefficients), 
                  l.gam3 = length(gam3$coefficients), l.gam4 = length(gam4$coefficients), 
                  l.gam5 = length(gam5$coefficients), l.gam6 = length(gam6$coefficients), 
                  l.gam7 = length(gam7$coefficients), l.gam8 = length(gam8$coefficients), 
                  l.gam9 = 0)
    L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                 l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, 
                 l.sp8 = l.sp8, l.sp9 = l.sp9)
    sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8)
    qu.mag <- S.m(GAM, L.SP, L.GAM)
  }
  if (missing(parscale)) 
    parscale <- 1
  respvec <- list(y1 = y1, y2 = y2, y1.y2 = y1.y2, y1.cy2 = y1.cy2, 
                  cy1.y2 = cy1.y2, cy1.cy2 = cy1.cy2, cy1 = cy1, cy = cy, 
                  univ = 0, y00 = y00, y10 = y10, y0p = y0p, y1p = y1p)
  my.env <- new.env()
  my.env$signind <- 1
  lsgam1 <- length(gam1$smooth)
  lsgam2 <- length(gam2$smooth)
  lsgam3 <- length(gam3$smooth)
  lsgam4 <- length(gam4$smooth)
  lsgam5 <- length(gam5$smooth)
  lsgam6 <- length(gam6$smooth)
  lsgam7 <- length(gam7$smooth)
  lsgam8 <- length(gam8$smooth)
  lsgam9 <- length(gam9$smooth)
  VC <- list(lsgam1 = lsgam1, robust = FALSE, sp.fixed = NULL, 
             K1 = NULL, lsgam2 = lsgam2, Sl.sf = Sl.sf, sp.method = sp.method, 
             lsgam3 = lsgam3, lsgam4 = lsgam4, lsgam5 = lsgam5, lsgam6 = lsgam6, 
             lsgam7 = lsgam7, lsgam8 = lsgam8, lsgam9 = lsgam9, X1 = X1, 
             inde = inde, my.env = my.env, X2 = X2, X3 = X3, X4 = X4, 
             X5 = X5, X6 = X6, X7 = X7, X8 = X8, X1.d2 = X1.d2, X2.d2 = X2.d2, 
             X3.d2 = X3.d2, X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, 
             X7.d2 = X7.d2, X8.d2 = X8.d2, gp1 = gp1, gp2 = gp2, gp3 = gp3, 
             gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, 
             l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, 
             l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8, 
             l.sp9 = 0, infl.fac = infl.fac, weights = weights, fp = fp, 
             univ.gamls = FALSE, hess = hess, nCa = nCa, Model = Model, 
             gamlssfit = gamlssfit, end = end, BivD = BivD, dof.st = log(dof - 
                                                                           2), dof = dof, nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI, 
             parscale = parscale, margins = margins, Cont = "NO", 
             ccss = "no", m2 = m2, m3 = m3, m2d = m2d, m1d = m1d, 
             m3d = NULL, bl = bl, X2s = X2s, X3s = X3s, triv = FALSE, 
             y2m = y2m, theta.fx = theta.fx, i.rho = i.rho, BivD2 = BivD2, 
             cta = cta, ct = ct, zerov = -10, surv.flex = surv.flex, 
             gp2.inf = NULL, informative = "no", zero.tol = 0.01, 
             min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, y00 = y00, 
             y10 = y10, y0p = y0p, y1p = y1p)
  if (gc.l == TRUE) 
    gc()
  if (gamlssfit == TRUE && !(margins[2] %in% bl)) {
    form.gamlR <- form.gaml(formula, l.flist, M, type = "biv")
    gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, 
                                      data = data, weights = weights, subset = subset, 
                                      margin = margins[2], infl.fac = infl.fac, rinit = rinit, 
                                      rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp, 
                                      gc.l = gc.l, parscale = 1, extra.regI = extra.regI, 
                                      drop.unused.levels = drop.unused.levels), list(weights = weights)))
    SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4, 
               sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8)
    gamls.upsvR <- gamls.upsv(gamlss1 = NULL, gamlss2, margins, 
                              M, l.flist, nstv = NULL, VC, GAM, SP, type = "biv")
    sp <- gamls.upsvR$sp
    start.v <- gamls.upsvR$start.v
  }
  if (Model %in% c("BSS", "BPO", "BPO0") || (Model == "B" && 
                                             margins[2] %in% bl)) 
    gamlssfit <- VC$gamlssfit <- TRUE
  func.opt <- func.OPT_l(margins, M, type = "biv")
  SemiParFit <- SemiParBIV.fit_l(func.opt = func.opt, start.v = start.v, 
                               rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, 
                               tolsp = tolsp, respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag)
  SemiParFit.p <- SemiParBIV.fit.post(SemiParFit = SemiParFit, 
                                      Model = Model, VC = VC, GAM)
  SemiParFit <- SemiParFit.p$SemiParFit
  y2.m <- y2
  if (margins[2] == "LN") 
    y2.m <- exp(y2)
  if (gc.l == TRUE) 
    gc()
  cov.c(SemiParFit)
  gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data
  if (!(Model == "B" && !(margins[2] %in% bl) && end == 2)) {
    dataset <- NULL
    rm(data)
  }
  else {
    attr(data, "terms") <- NULL
    dataset <- data
    rm(data)
  }
  formula.aux <- NULL
  if (Model == "BSS") {
    formula.aux <- formula
    formula.aux[[1]] <- reformulate(attr(terms(formula.aux[[1]]), 
                                         "term.labels"), response = "ry")
    environment(formula.aux[[1]]) <- environment(formula[[1]])
    formula.aux[[2]] <- reformulate(attr(terms(formula.aux[[2]]), 
                                         "term.labels"), response = "y")
    environment(formula.aux[[2]]) <- environment(formula[[2]])
  }
  L <- list(fit = SemiParFit$fit, dataset = dataset, formula = formula, 
            SemiParFit = SemiParFit, mice.formula = formula.aux, 
            gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, 
            gam6 = gam6, robust = FALSE, gam7 = gam7, gam8 = gam8, 
            gam9 = gam9, gam2TW = gam2TW, coefficients = SemiParFit$fit$argument, 
            coef.t = NULL, iterlimsp = iterlimsp, weights = weights, 
            sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, l.sp1 = l.sp1, 
            l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, l.sp5 = l.sp5, 
            l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = l.sp9, 
            bl = bl, fp = fp, iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner, 
            theta = SemiParFit.p$theta, theta.a = SemiParFit.p$theta.a, 
            OR = SemiParFit.p$OR, GM = SemiParFit.p$GM, n = n, n.sel = n.sel, 
            X1 = X1, X2 = X2, X3 = X3, X1.d2 = X1.d2, X2.d2 = X2.d2, 
            X3.d2 = X3.d2, X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8, 
            X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, X7.d2 = X7.d2, 
            X8.d2 = X8.d2, He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, 
            Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, F = SemiParFit.p$F, 
            F1 = SemiParFit.p$F1, t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
            edf11 = SemiParFit.p$edf11, edf1 = SemiParFit.p$edf1, 
            edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3, edf4 = SemiParFit.p$edf4, 
            edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6, edf7 = SemiParFit.p$edf7, 
            edf8 = SemiParFit.p$edf8, edf1.1 = SemiParFit.p$edf1.1, 
            edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3, 
            edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, 
            edf1.6 = SemiParFit.p$edf1.6, edf1.7 = SemiParFit.p$edf1.7, 
            edf1.8 = SemiParFit.p$edf1.8, R = SemiParFit.p$R, bs.mgfit = SemiParFit$bs.mgfit, 
            conv.sp = SemiParFit$conv.sp, wor.c = SemiParFit$wor.c, 
            p11 = SemiParFit$fit$p11, p10 = SemiParFit$fit$p10, p01 = SemiParFit$fit$p01, 
            p00 = SemiParFit$fit$p00, p1 = SemiParFit$fit$p1, p2 = SemiParFit$fit$p2, 
            eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, 
            etad = SemiParFit$fit$etad, etas = SemiParFit$fit$etas, 
            etan = SemiParFit$fit$etan, y1 = y1, y2 = y2.m, BivD = BivD, 
            margins = margins, logLik = SemiParFit.p$logLik, nC = nC, 
            hess = hess, respvec = respvec, inde = inde, qu.mag = qu.mag, 
            sigma2 = SemiParFit.p$sigma2, sigma2.a = SemiParFit.p$sigma2.a, 
            sigma = SemiParFit.p$sigma2, sigma.a = SemiParFit.p$sigma2.a, 
            nu = SemiParFit.p$nu, nu.a = SemiParFit.p$nu.a, Vb.t = SemiParFit.p$Vb.t, 
            gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, 
            gp6 = gp6, gp7 = gp7, gp8 = gp8, X2s = X2s, X3s = X3s, 
            p1n = SemiParFit.p$p1n, p2n = SemiParFit.p$p2n, VC = VC, 
            Model = Model, magpp = SemiParFit$magpp, gamlssfit = gamlssfit, 
            Cont = "NO", tau = SemiParFit.p$tau, tau.a = SemiParFit.p$tau.a, 
            l.flist = l.flist, v1 = v1, v2 = v2, triv = FALSE, univar.gamlss = FALSE, 
            gamlss = gamlss2, BivD2 = BivD2, dof = dof, dof.a = dof, 
            call = cl, surv = FALSE, surv.flex = surv.flex)
  if (BivD %in% BivD2) {
    L$teta1 <- SemiParFit$fit$teta1
    L$teta.ind1 <- SemiParFit$fit$teta.ind1
    L$teta2 <- SemiParFit$fit$teta2
    L$teta.ind2 <- SemiParFit$fit$teta.ind2
    L$Cop1 <- SemiParFit$fit$Cop1
    L$Cop2 <- SemiParFit$fit$Cop2
  }
  class(L) <- c("SemiParBIV", "gjrm")
  L
}



SemiParBIV_l.fit <- function (func.opt, start.v, rinit, rmax, iterlim, iterlimsp, 
          tolsp, respvec, VC, sp = NULL, qu.mag = NULL) 
{
  l.sp1 <- VC$l.sp1
  l.sp2 <- VC$l.sp2
  l.sp3 <- VC$l.sp3
  l.sp4 <- VC$l.sp4
  l.sp5 <- VC$l.sp5
  l.sp6 <- VC$l.sp6
  l.sp7 <- VC$l.sp7
  l.sp8 <- VC$l.sp8
  l.sp9 <- VC$l.sp9
  score.hist <- rp <- D <- L <- Sl.sfTemp <- St <- NULL
  l.splist <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                   l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, 
                   l.sp8 = l.sp8, l.sp9 = l.sp9)
  if (!is.null(VC$sp.fixed)) 
    sp <- VC$sp.fixed
  if ((l.sp1 == 0 && l.sp2 == 0 && l.sp3 == 0 && l.sp4 == 0 && 
       l.sp5 == 0 && l.sp6 == 0 && l.sp7 == 0 && l.sp8 == 0 && 
       l.sp9 == 0) || VC$fp == TRUE) 
    ps <- ps1 <- list(S.h = 0, S.h1 = 0, S.h2 = 0, qu.mag = NULL)
  else ps <- ps1 <- pen(qu.mag, sp, VC, univ = respvec$univ, 
                        l.splist)
  if (VC$triv == TRUE) {
    if (VC$penCor == "ridge") 
      qu.mag <- ps$qu.mag
    if (VC$penCor %in% c("lasso", "alasso")) 
      VC$sp <- sp
  }
  parsc <- rep(VC$parscale, length(start.v))
  sc <- TRUE
  fit <- fit1 <- try(trust(func.opt, start.v, rinit = rinit, 
                           rmax = rmax, parscale = parsc, respvec = respvec, VC = VC, 
                           ps = ps, blather = TRUE, iterlim = iterlim), silent = sc)
  if (inherits(fit, "try-error") || is.null(fit$l)) {
    fit <- fit1 <- try(trust(func.opt, start.v, rinit = rinit, 
                             rmax = rmax, parscale = parsc, respvec = respvec, 
                             VC = VC, ps = ps, blather = TRUE, iterlim = iterlim/4), 
                       silent = sc)
    if (inherits(fit, "try-error") || is.null(fit$l)) {
      fit <- fit1 <- try(trust(func.opt, start.v, rinit = rinit, 
                               rmax = rmax, parscale = parsc, respvec = respvec, 
                               VC = VC, ps = ps, blather = TRUE, iterlim = iterlim/10), 
                         silent = sc)
      if ((inherits(fit, "try-error") || is.null(fit$l)) && 
          VC$gamlssfit == FALSE) 
        stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
      if ((inherits(fit, "try-error") || is.null(fit$l)) && 
          VC$gamlssfit == TRUE) 
        stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")
    }
  }
  iter.if <- fit$iterations
  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- magpp <- NULL
  conv.sp <- TRUE
  if ((VC$fp == FALSE && is.null(VC$sp.fixed) && (l.sp1 != 
                                                  0 || l.sp2 != 0 || l.sp3 != 0 || l.sp4 != 0 || l.sp5 != 
                                                  0 || l.sp6 != 0 || l.sp7 != 0 || l.sp8 != 0 || l.sp9 != 
                                                  0))) {
    if (VC$sp.method == "perf") {
      stoprule.SP <- 1
      conv.sp <- TRUE
      iter.inner <- iter.sp <- 0
      while (stoprule.SP > tolsp) {
        fito <- fit$l
        o.ests <- c(fit$argument)
        spo <- sp
        wor.c <- working.comp(fit)
        if (VC$triv == TRUE && VC$penCor %in% c("lasso", 
                                                "alasso")) 
          qu.mag <- fit$qu.mag
        bs.mgfit <- try(magic(y = wor.c$Z, X = wor.c$X, 
                              sp = sp, S = qu.mag$Ss, off = qu.mag$off, rank = qu.mag$rank, 
                              gcv = FALSE, gamma = VC$infl.fac), silent = sc)
        if (inherits(bs.mgfit, "try-error")) {
          conv.sp <- FALSE
          break
        }
        if (any(is.na(bs.mgfit$sp)) == TRUE) {
          conv.sp <- FALSE
          break
        }
        sp <- bs.mgfit$sp
        if (!is.null(VC$sp.fixed)) 
          sp <- VC$sp.fixed
        iter.sp <- iter.sp + 1
        names(sp) <- names(spo)
        if (VC$triv == TRUE && VC$penCor %in% c("lasso", 
                                                "alasso")) 
          VC$sp <- sp
        ps <- pen(qu.mag, sp, VC, univ = respvec$univ, 
                  l.splist)
        fit <- try(trust(func.opt, o.ests, rinit = rinit, 
                         rmax = rmax, parscale = parsc, respvec = respvec, 
                         VC = VC, ps = ps, blather = TRUE, iterlim = iterlim), 
                   silent = sc)
        if (inherits(fit, "try-error") || is.null(fit$l)) {
          conv.sp <- FALSE
          ps <- ps1
          fit <- try(trust(func.opt, c(fit1$argument), 
                           rinit = rinit, rmax = rmax, parscale = parsc, 
                           respvec = respvec, VC = VC, ps = ps, blather = TRUE, 
                           iterlim = iterlim), silent = sc)
          if ((inherits(fit, "try-error") || is.null(fit$l)) && 
              VC$gamlssfit == FALSE) 
            stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
          if ((inherits(fit, "try-error") || is.null(fit$l)) && 
              VC$gamlssfit == TRUE) 
            stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")
        }
        iter.inner <- iter.inner + fit$iterations
        if (iter.sp >= iterlimsp) {
          conv.sp <- FALSE
          break
        }
        stoprule.SP <- abs(fit$l - fito)/(0.1 + abs(fit$l))
      }
      if (VC$gc.l == TRUE) 
        gc()
      magpp <- magic.post.proc(wor.c$X, bs.mgfit)
    }
    if (VC$sp.method == "efs") {
      LDfun <- function(Hp, eigen.fix) {
        rank <- dim(Hp)[1]
        D <- diag(Hp)
        if (sum(!is.finite(D)) > 0) 
          stop("non finite values in Hessian")
        if (min(D) < 0) {
          Dthresh <- max(D) * sqrt(.Machine$double.eps)
          if (-min(D) < Dthresh) {
            indefinite <- FALSE
            D[D < Dthresh] <- Dthresh
          }
          else indefinite <- TRUE
        }
        else indefinite <- FALSE
        if (indefinite) {
          if (eigen.fix) {
            eh <- eigen(Hp, symmetric = TRUE)
            ev <- abs(eh$values)
            Hp <- eh$vectors %*% (ev * t(eh$vectors))
          }
          else {
            Ib <- diag(rank) * abs(min(D))
            Ip <- diag(rank) * abs(max(D) * .Machine$double.eps^0.5)
            Hp <- Hp + Ip + Ib
          }
          D <- rep(1, ncol(Hp))
          indefinite <- TRUE
        }
        else {
          D <- D^-0.5
          Hp <- D * t(D * Hp)
          Ip <- diag(rank) * .Machine$double.eps^0.5
        }
        L <- suppressWarnings(chol(Hp, pivot = TRUE))
        mult <- 1
        while (attr(L, "rank") < rank) {
          if (eigen.fix) {
            eh <- eigen(Hp, symmetric = TRUE)
            ev <- eh$values
            thresh <- max(min(ev[ev > 0]), max(ev) * 
                            1e-06) * mult
            mult <- mult * 10
            ev[ev < thresh] <- thresh
            Hp <- eh$vectors %*% (ev * t(eh$vectors))
            L <- suppressWarnings(chol(Hp, pivot = TRUE))
          }
          else {
            L <- suppressWarnings(chol(Hp + Ip, pivot = TRUE))
            Ip <- Ip * 100
          }
          indefinite <- TRUE
        }
        list(L = L, D = D)
      }
      stoprule.SP <- 1
      conv.sp <- TRUE
      iter.inner <- iter.sp <- 0
      controlEFS <- list(efs.lspmax = 15, eps = 1e-07, 
                         tol = 1e-06, tiny = .Machine$double.eps^0.5, 
                         efs.tol = 0.1)
      score.hist <- rep(0, 200)
      mult <- 1
      lsp <- log(sp)
      gamma <- 1
      Mp <- -1
      eigen.fix <- FALSE
      Sl.termMult <- getFromNamespace("Sl.termMult", "mgcv")
      ldetS <- getFromNamespace("ldetS", "mgcv")
      Mp <- ncol(totalPenaltySpace(qu.mag$Ss, NULL, qu.mag$off, 
                                   length(fit$argument))$Z)
      Sl.sfTemp <- VC$Sl.sf
      for (i in 1:length(Sl.sfTemp)) Sl.sfTemp[[i]]$D <- solve(Sl.sfTemp[[i]]$D)
      for (iter in 1:200) {
        o.ests <- c(fit$argument)
        rp <- ldetS(VC$Sl.sf, rho = lsp, fixed = rep(FALSE, 
                                                     length(lsp)), np = length(fit$argument), root = TRUE)
        o.estsStar <- Sl.initial.repara(Sl.sfTemp, o.ests, 
                                        inverse = TRUE)
        Vb <- Sl.initial.repara(Sl.sfTemp, PDef(fit$hessian)$res.inv, 
                                inverse = TRUE)
        LD <- LDfun(PDef(Vb)$res.inv, eigen.fix)
        L <- LD$L
        D <- LD$D
        ipiv <- piv <- attr(L, "pivot")
        p <- length(piv)
        ipiv[piv] <- 1:p
        Vb <- crossprod(forwardsolve(t(L), diag(D, nrow = p)[piv, 
                                                             , drop = FALSE])[ipiv, , drop = FALSE])
        Vb <- Sl.repara(rp$rp, Vb, inverse = TRUE)
        SVb <- Sl.termMult(VC$Sl.sf, Vb)
        trVS <- rep(0, length(SVb))
        for (i in 1:length(SVb)) {
          ind <- attr(SVb[[i]], "ind")
          trVS[i] <- sum(diag(SVb[[i]][, ind]))
        }
        start <- Sl.repara(rp$rp, o.estsStar)
        Sb <- Sl.termMult(VC$Sl.sf, start, full = TRUE)
        bSb <- rep(0, length(Sb))
        for (i in 1:length(Sb)) bSb[i] <- sum(start * 
                                                Sb[[i]])
        S1 <- rp$ldet1
        a <- pmax(controlEFS$tiny, S1 * exp(-lsp) - trVS)
        r <- a/pmax(controlEFS$tiny, bSb)
        r[a == 0 & bSb == 0] <- 1
        r[!is.finite(r)] <- 1e+06
        lsp1 <- pmin(lsp + log(r) * mult, controlEFS$efs.lspmax)
        max.step <- max(abs(lsp1 - lsp))
        ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
        old.reml <- -as.numeric((-fit$l - drop(t(fit$argument) %*% 
                                                 fit$S.h %*% fit$argument)/2)/gamma + rp$ldetS/2 - 
                                  ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)
        sp1 <- exp(lsp1)
        names(sp1) <- names(lsp1)
        if (!is.null(VC$sp.fixed)) 
          sp1 <- VC$sp.fixed
        if (VC$triv == TRUE && VC$penCor %in% c("lasso", 
                                                "alasso")) 
          VC$sp <- sp1
        ps <- pen(qu.mag, sp1, VC, univ = respvec$univ, 
                  l.splist)
        fit <- try(trust(func.opt, o.ests, rinit = rinit, 
                         rmax = rmax, parscale = parsc, respvec = respvec, 
                         VC = VC, ps = ps, blather = TRUE, iterlim = iterlim), 
                   silent = sc)
        if (inherits(fit, "try-error") || is.null(fit$l)) {
          conv.sp <- FALSE
          ps <- ps1
          fit <- try(trust(func.opt, c(fit1$argument), 
                           rinit = rinit, rmax = rmax, parscale = parsc, 
                           respvec = respvec, VC = VC, ps = ps, blather = TRUE, 
                           iterlim = iterlim), silent = sc)
          if ((inherits(fit, "try-error") || is.null(fit$l)) && 
              VC$gamlssfit == FALSE) 
            stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
          if ((inherits(fit, "try-error") || is.null(fit$l)) && 
              VC$gamlssfit == TRUE) 
            stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")
        }
        iter.inner <- iter.inner + fit$iterations
        rp <- ldetS(VC$Sl.sf, rho = lsp1, fixed = rep(FALSE, 
                                                      length(lsp1)), np = length(fit$argument), root = TRUE)
        Vb <- Sl.initial.repara(Sl.sfTemp, PDef(fit$hessian)$res.inv, 
                                inverse = TRUE)
        LD <- LDfun(PDef(Vb)$res.inv, eigen.fix)
        L <- LD$L
        D <- LD$D
        ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
        fit$REML <- -as.numeric((-fit$l - drop(t(fit$argument) %*% 
                                                 fit$S.h %*% fit$argument)/2)/gamma + rp$ldetS/2 - 
                                  ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)
        if (fit$REML <= old.reml) {
          if (max.step < 0.05) {
            lsp2 <- pmin(lsp + log(r) * mult * 2, 12)
            sp2 <- exp(lsp2)
            names(sp2) <- names(lsp2)
            if (!is.null(VC$sp.fixed)) 
              sp2 <- VC$sp.fixed
            if (VC$triv == TRUE && VC$penCor %in% c("lasso", 
                                                    "alasso")) 
              VC$sp <- sp2
            ps <- pen(qu.mag, sp2, VC, univ = respvec$univ, 
                      l.splist)
            fit2 <- try(trust(func.opt, o.ests, rinit = rinit, 
                              rmax = rmax, parscale = parsc, respvec = respvec, 
                              VC = VC, ps = ps, blather = TRUE, iterlim = iterlim), 
                        silent = sc)
            if (inherits(fit2, "try-error") || is.null(fit2$l)) {
              conv.sp <- FALSE
              ps <- ps1
              fit2 <- try(trust(func.opt, c(fit1$argument), 
                                rinit = rinit, rmax = rmax, parscale = parsc, 
                                respvec = respvec, VC = VC, ps = ps, 
                                blather = TRUE, iterlim = iterlim), silent = sc)
              if ((inherits(fit2, "try-error") || is.null(fit2$l)) && 
                  VC$gamlssfit == FALSE) 
                stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
              if ((inherits(fit2, "try-error") || is.null(fit2$l)) && 
                  VC$gamlssfit == TRUE) 
                stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")
            }
            iter.inner <- iter.inner + fit2$iterations
            rp <- ldetS(VC$Sl.sf, rho = lsp2, fixed = rep(FALSE, 
                                                          length(lsp2)), np = length(fit$argument), 
                        root = TRUE)
            Vb <- Sl.initial.repara(Sl.sfTemp, PDef(fit2$hessian)$res.inv, 
                                    inverse = TRUE)
            LD <- LDfun(PDef(Vb)$res.inv, eigen.fix)
            L <- LD$L
            D <- LD$D
            ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
            fit2$REML <- -as.numeric((-fit2$l - drop(t(fit2$argument) %*% 
                                                       fit2$S.h %*% fit2$argument)/2)/gamma + 
                                       rp$ldetS/2 - ldetHp/2 + Mp * (log(2 * pi)/2) - 
                                       log(gamma)/2)
            if (fit2$REML < fit$REML) {
              fit <- fit2
              lsp <- lsp2
              mult <- mult * 2
            }
            else {
              lsp <- lsp1
            }
          }
          else lsp <- lsp1
        }
        else {
          while (fit$REML > old.reml && mult > 1) {
            mult <- mult/2
            lsp1 <- pmin(lsp + log(r) * mult, controlEFS$efs.lspmax)
            sp1 <- exp(lsp1)
            names(sp1) <- names(lsp1)
            if (!is.null(VC$sp.fixed)) 
              sp1 <- VC$sp.fixed
            if (VC$triv == TRUE && VC$penCor %in% c("lasso", 
                                                    "alasso")) 
              VC$sp <- sp1
            ps <- pen(qu.mag, sp1, VC, univ = respvec$univ, 
                      l.splist)
            fit <- try(trust(func.opt, o.ests, rinit = rinit, 
                             rmax = rmax, parscale = parsc, respvec = respvec, 
                             VC = VC, ps = ps, blather = TRUE, iterlim = iterlim), 
                       silent = sc)
            if (inherits(fit, "try-error") || is.null(fit$l)) {
              conv.sp <- FALSE
              ps <- ps1
              fit <- try(trust(func.opt, c(fit1$argument), 
                               rinit = rinit, rmax = rmax, parscale = parsc, 
                               respvec = respvec, VC = VC, ps = ps, 
                               blather = TRUE, iterlim = iterlim), silent = sc)
              if ((inherits(fit, "try-error") || is.null(fit$l)) && 
                  VC$gamlssfit == FALSE) 
                stop("It is not possible to fit the model.\nTry re-fitting the model and setting gamlssfit = TRUE if allowed.\nAlso, read the WARNINGS section of ?gjrm.")
              if ((inherits(fit, "try-error") || is.null(fit$l)) && 
                  VC$gamlssfit == TRUE) 
                stop("It is not possible to fit the model.\nRead the WARNINGS section of ?gjrm.")
            }
            iter.inner <- iter.inner + fit$iterations
            rp <- ldetS(VC$Sl.sf, rho = lsp1, fixed = rep(FALSE, 
                                                          length(lsp1)), np = length(fit$argument), 
                        root = TRUE)
            Vb <- Sl.initial.repara(Sl.sfTemp, PDef(fit$hessian)$res.inv, 
                                    inverse = TRUE)
            LD <- LDfun(PDef(Vb)$res.inv, eigen.fix)
            L <- LD$L
            D <- LD$D
            ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
            fit$REML <- -as.numeric((-fit$l - drop(t(fit$argument) %*% 
                                                     fit$S.h %*% fit$argument)/2)/gamma + rp$ldetS/2 - 
                                      ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)
          }
          lsp <- lsp1
          if (mult < 1) 
            mult <- 1
        }
        score.hist[iter] <- fit$REML
        if (iter > 3 && max.step < 0.05 && max(abs(diff(score.hist[(iter - 
                                                                    3):iter]))) < controlEFS$efs.tol) 
          break
        if (iter == 1) 
          old.ll <- fit$l
        else {
          if (abs(old.ll - fit$l) < 100 * controlEFS$eps * 
              abs(fit$l)) 
            break
          old.ll <- fit$l
        }
      }
      sp <- exp(lsp)
      iter.sp <- iter
      if (iter > 200) 
        conv.sp <- FALSE
      else conv.sp <- TRUE
      if (VC$gc.l == TRUE) 
        gc()
      St <- crossprod(rp$E)
    }
  }
  else {
    wor.c <- working.comp(fit)
    bs.mgfit <- magic(wor.c$Z, wor.c$X, numeric(0), list(), 
                      numeric(0))
    magpp <- magic.post.proc(wor.c$X, bs.mgfit)
  }
  rm(fit1, ps1)
  list(fit = fit, score.hist = score.hist, iter.if = iter.if, 
       sp.method = VC$sp.method, conv.sp = conv.sp, iter.sp = iter.sp, 
       iter.inner = iter.inner, bs.mgfit = bs.mgfit, wor.c = wor.c, 
       sp = sp, magpp = magpp, rp = rp$rp, Sl = VC$Sl.sf, D = D, 
       L = L, Sl.sfTemp = Sl.sfTemp, St = St)
}

func.OPT_l <- function (margins, M, type = "copR") 
{
  if (type == "ROY") {
    if (margins[2] %in% M$m1d && margins[3] %in% M$m1d) 
      func.opt <- bprobgHsDiscr1ROY
    if (margins[2] %in% M$m2d && margins[3] %in% M$m2d) 
      func.opt <- bprobgHsDiscr2ROY
    if (margins[2] %in% M$m2 && margins[3] %in% M$m2) 
      func.opt <- bprobgHsCont2ROY
    if (margins[2] %in% M$m3 && margins[3] %in% M$m3) 
      func.opt <- bprobgHsCont3ROY
    if (margins[2] %in% M$bl && margins[3] %in% M$bl) 
      func.opt <- bprobgHsBinROY
  }
  if (type == "biv") {
    if (M$Model == "B" && margins[2] %in% M$bl && is.null(M$theta.fx) && 
        is.null(M$K2)) 
      func.opt <- bprobgHs
    if (M$Model == "B" && margins[2] %in% M$bl && !is.null(M$theta.fx) && 
        is.null(M$K2)) 
      func.opt <- bprobgHsFixTheta
    if (M$Model == "B" && margins[2] %in% M$m1d) {
      func.opt <- bprobgHsDiscr1
      func.optUniv <- bprobgHsContUniv
    }
    if (M$Model == "B" && margins[2] %in% M$m2d) {
      func.opt <- bprobgHsDiscr2
      func.optUniv <- bprobgHsContUniv
    }
    if (M$Model == "BPO") 
      func.opt <- bprobgHsPO
    if (M$Model == "BPO0") 
      func.opt <- bprobgHsPO0
    if (M$Model == "BSS") 
      func.opt <- bprobgHsSS
    if (M$Model == "B" && margins[2] %in% M$m2 && is.null(M$K1)) {
      func.opt <- bprobgHsCont
      func.optUniv <- bprobgHsContUniv
    }
    if (M$Model == "B" && margins[2] %in% M$m3 && is.null(M$K1)) {
      func.opt <- bprobgHsCont3
      func.optUniv <- bprobgHsContUniv3
    }
    if (M$Model == "B" && margins[2] == "TW" && is.null(M$K1)) {
      func.opt <- bprobgHsCont3binTW
      func.optUniv <- bprobgHsContUniv3
    }
    if (M$Model == "B" && margins[2] %in% M$m2 && !is.null(M$K1)) {
      func.opt <- bCopulaCLMgHsCont
    }
    if (M$Model == "B" && margins[2] %in% M$bl && !is.null(M$K2)) {
      func.opt <- bCopulaCLMgHsOrd
    }
  }
  if (type == "copR") {
    if (margins[1] %in% M$m1d && margins[2] %in% M$m1d) 
      func.opt <- bdiscrdiscr11
    if (margins[1] %in% M$m1d && margins[2] %in% M$m2d) 
      func.opt <- bdiscrdiscr12
    if (margins[1] %in% M$m2d && margins[2] %in% M$m2d) 
      func.opt <- bdiscrdiscr
    if (margins[1] %in% M$m2 && margins[2] %in% M$m2 && M$BivD != 
        "T") {
      if (M$robust == FALSE) 
        func.opt <- bcont
      else func.opt <- bcontROB
    }
    if (margins[1] %in% M$m2 && margins[2] %in% M$m2 && M$BivD == 
        "T") 
      func.opt <- bconttwoParC
    if (margins[1] %in% M$m3 && margins[2] %in% M$m3 && M$BivD != 
        "T") 
      func.opt <- bcont3
    if (margins[1] %in% M$m3 && margins[2] %in% M$m3 && M$BivD == 
        "T") 
      func.opt <- bcont3twoParC
    if (margins[1] %in% M$m2 && margins[2] %in% M$m3 && M$BivD != 
        "T") 
      func.opt <- bcont23
    if (margins[1] %in% M$m2 && margins[2] %in% M$m3 && M$BivD == 
        "T") 
      func.opt <- bcont23twoParC
    if (margins[1] %in% M$m3 && margins[2] %in% M$m2 && M$BivD != 
        "T") 
      func.opt <- bcont32
    if (margins[1] %in% M$m3 && margins[2] %in% M$m2 && M$BivD == 
        "T") 
      func.opt <- bcont32twoParC
    if (margins[1] %in% M$m2 && margins[2] %in% M$m2 && M$surv == 
        TRUE) 
      func.opt <- bcontSurv
    if (margins[1] %in% M$bl && margins[2] %in% M$bl && M$surv == 
        TRUE && M$dep.cens == FALSE) 
      func.opt <- bcontSurvG
    if (margins[1] %in% M$bl && margins[2] %in% M$bl && M$surv == 
        TRUE && M$dep.cens == FALSE && M$type.cens1 == "mixed" && 
        M$type.cens2 == "mixed") 
      func.opt <- bcontSurvG_extended
    if (margins[1] %in% M$bl && margins[2] %in% M$bl && M$surv == 
        TRUE && M$dep.cens == TRUE && is.null(M$c3)) 
      func.opt <- bcontSurvGDep
    if (margins[1] %in% M$bl && margins[2] %in% M$bl && M$surv == 
        TRUE && M$dep.cens == TRUE && !is.null(M$c3)) 
      func.opt <- bcontSurvGDepA
  }
  if (type == "copSS") {
    if (margins[2] %in% M$m2) 
      func.opt <- bprobgHsContSS
    if (margins[2] %in% M$m3) 
      func.opt <- bprobgHsCont3SS
    if (margins[2] %in% M$m1d) 
      func.opt <- bprobgHsDiscr1SS
    if (margins[2] %in% M$m2d) 
      func.opt <- bprobgHsDiscr2SS
  }
  func.opt
}



function (params, respvec, VC, ps, AT = FALSE) 
{
  p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA
  eta1 <- VC$X1 %*% params[1:VC$X1.d2]
  eta2 <- VC$X2 %*% params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
  etad <- etas1 <- etas2 <- l.ln <- NULL
  if (is.null(VC$X3)) {
    sigma21.st <- etas1 <- params[(VC$X1.d2 + VC$X2.d2 + 
                                     1)]
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 
                                     2)]
    teta.st <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 3)]
  }
  if (!is.null(VC$X3)) {
    sigma21.st <- etas1 <- VC$X3 %*% params[(VC$X1.d2 + VC$X2.d2 + 
                                               1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    sigma22.st <- etas2 <- VC$X4 %*% params[(VC$X1.d2 + VC$X2.d2 + 
                                               VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    teta.st <- etad <- VC$X5 %*% params[(VC$X1.d2 + VC$X2.d2 + 
                                           VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 
                                                                       VC$X4.d2 + VC$X5.d2)]
  }
  sstr1 <- esp.tr(sigma21.st, VC$margins[1])
  sstr2 <- esp.tr(sigma22.st, VC$margins[2])
  sigma21.st <- sstr1$vrb.st
  sigma22.st <- sstr2$vrb.st
  sigma21 <- sstr1$vrb
  sigma22 <- sstr2$vrb
  eta1 <- eta.tr(eta1, VC$margins[1])
  eta2 <- eta.tr(eta2, VC$margins[2])
  resT <- teta.tr(VC, teta.st)
  teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
  teta1 <- teta2 <- teta <- resT$teta
  Cop1 <- Cop2 <- VC$BivD
  teta.ind1 <- as.logical(c(1, 0, round(runif(VC$n - 2))))
  teta.ind2 <- teta.ind1 == FALSE
  if (!(VC$BivD %in% VC$BivD2) && length(teta.st) > 1) {
    teta.st1 <- teta.st[teta.ind1]
    teta.st2 <- teta.st[teta.ind2]
    teta1 <- teta[teta.ind1]
    teta2 <- teta[teta.ind2]
  }
  if (VC$BivD %in% VC$BivD2) {
    if (VC$BivD %in% VC$BivD2[c(1:4, 13:16)]) 
      teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov), 
                          TRUE, FALSE)
    if (VC$BivD %in% VC$BivD2[5:12]) 
      teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov) + 
                            1, TRUE, FALSE)
    teta.ind2 <- teta.ind1 == FALSE
    VC$my.env$signind <- ifelse(teta.ind1 == TRUE, 1, -1)
    teta1 <- teta[teta.ind1]
    teta2 <- -teta[teta.ind2]
    teta.st1 <- teta.st[teta.ind1]
    teta.st2 <- teta.st[teta.ind2]
    if (length(teta) == 1) 
      teta.ind2 <- teta.ind1 <- rep(TRUE, VC$n)
    Cop1Cop2R <- Cop1Cop2(VC$BivD)
    Cop1 <- Cop1Cop2R$Cop1
    Cop2 <- Cop1Cop2R$Cop2
  }
  dHs1 <- distrHs(respvec$y1, eta1, sigma21, sigma21.st, nu = 1, 
                  nu.st = 1, margin2 = VC$margins[1], naive = FALSE, min.dn = VC$min.dn, 
                  min.pr = VC$min.pr, max.pr = VC$max.pr)
  dHs2 <- distrHs(respvec$y2, eta2, sigma22, sigma22.st, nu = 1, 
                  nu.st = 1, margin2 = VC$margins[2], naive = FALSE, min.dn = VC$min.dn, 
                  min.pr = VC$min.pr, max.pr = VC$max.pr)
  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2
  p1 <- dHs1$p2
  p2 <- dHs2$p2
  if (length(teta1) != 0) 
    dH1 <- copgHsAT(p1[teta.ind1], p2[teta.ind1], teta1, 
                    Cop1, Ln = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                    max.pr = VC$max.pr)
  if (length(teta2) != 0) 
    dH2 <- copgHsAT(p1[teta.ind2], p2[teta.ind2], teta2, 
                    Cop2, Ln = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                    max.pr = VC$max.pr)
  c.copula2.be1be2 <- NA
  if (length(teta1) != 0) 
    c.copula2.be1be2[teta.ind1] <- dH1$c.copula2.be1be2
  if (length(teta2) != 0) 
    c.copula2.be1be2[teta.ind2] <- dH2$c.copula2.be1be2
  l.par <- VC$weights * (log(pdf1) + log(pdf2) + log(c.copula2.be1be2))
  derpdf1.dereta1 <- dHs1$derpdf2.dereta2
  derpdf1.dersigma21.st <- dHs1$derpdf2.dersigma2.st
  derpdf2.dereta2 <- dHs2$derpdf2.dereta2
  derpdf2.dersigma22.st <- dHs2$derpdf2.dersigma2.st
  derp1.dereta1 <- dHs1$derp2.dereta2
  derp1.dersigma21.st <- dHs1$derp2.dersigma.st
  derp2.dereta2 <- dHs2$derp2.dereta2
  derp2.dersigma22.st <- dHs2$derp2.dersigma.st
  if (length(teta1) != 0) 
    BITS1 <- copgHsCont(p1[teta.ind1], p2[teta.ind1], teta1, 
                        teta.st1, Cop1, Cont = TRUE)
  if (length(teta2) != 0) 
    BITS2 <- copgHsCont(p1[teta.ind2], p2[teta.ind2], teta2, 
                        teta.st2, Cop2, Cont = TRUE)
  der2h.derp1p1 <- NA
  if (length(teta1) != 0) 
    der2h.derp1p1[teta.ind1] <- BITS1$der2h.derp1p1
  if (length(teta2) != 0) 
    der2h.derp1p1[teta.ind2] <- BITS2$der2h.derp1p1
  derc.dereta1 <- der2h.derp1p1 * derp1.dereta1
  derc.dersigma21.st <- der2h.derp1p1 * derp1.dersigma21.st
  der2h.derp1p2 <- NA
  if (length(teta1) != 0) 
    der2h.derp1p2[teta.ind1] <- BITS1$der2h.derp1p2
  if (length(teta2) != 0) 
    der2h.derp1p2[teta.ind2] <- BITS2$der2h.derp1p2
  derc.dereta2 <- der2h.derp1p2 * derp2.dereta2
  derc.dersigma22.st <- der2h.derp1p2 * derp2.dersigma22.st
  der2h.derp1teta <- NA
  derteta.derteta.st <- NA
  if (length(teta1) != 0) 
    der2h.derp1teta[teta.ind1] <- BITS1$der2h.derp1teta
  if (length(teta2) != 0) 
    der2h.derp1teta[teta.ind2] <- BITS2$der2h.derp1teta
  if (length(teta1) != 0) 
    derteta.derteta.st[teta.ind1] <- BITS1$derteta.derteta.st
  if (length(teta2) != 0) 
    derteta.derteta.st[teta.ind2] <- BITS2$derteta.derteta.st
  der2h.derp1teta.st <- der2h.derp1teta * derteta.derteta.st
  dl.dbe1 <- VC$weights * (derpdf1.dereta1/pdf1 + derc.dereta1/c.copula2.be1be2)
  dl.dbe2 <- VC$weights * (derpdf2.dereta2/pdf2 + derc.dereta2/c.copula2.be1be2)
  dl.dsigma21.st <- VC$weights * (derpdf1.dersigma21.st/pdf1 + 
                                    derc.dersigma21.st/c.copula2.be1be2)
  dl.dsigma22.st <- VC$weights * (derpdf2.dersigma22.st/pdf2 + 
                                    derc.dersigma22.st/c.copula2.be1be2)
  dl.dteta.st <- VC$weights * (der2h.derp1teta.st/c.copula2.be1be2)
  der2c.derrho.derrho <- NA
  der2c.derp1.derp1 <- NA
  der2c.derp2.derp2 <- NA
  der2c.derp1.derp2 <- NA
  der2c.derp1.derrho <- NA
  der2c.derp2.derrho <- NA
  der2teta.derteta.stteta.st <- NA
  if (length(teta1) != 0) {
    der2c.derrho.derrho[teta.ind1] <- BITS1$der2c.derrho.derrho
    der2c.derp1.derp1[teta.ind1] <- BITS1$der2c.derp1.derp1
    der2c.derp2.derp2[teta.ind1] <- BITS1$der2c.derp2.derp2
    der2c.derp1.derp2[teta.ind1] <- BITS1$der2c.derp1.derp2
    der2c.derp1.derrho[teta.ind1] <- BITS1$der2c.derp1.derrho
    der2c.derp2.derrho[teta.ind1] <- BITS1$der2c.derp2.derrho
  }
  if (length(teta2) != 0) {
    der2c.derrho.derrho[teta.ind2] <- BITS2$der2c.derrho.derrho
    der2c.derp1.derp1[teta.ind2] <- BITS2$der2c.derp1.derp1
    der2c.derp2.derp2[teta.ind2] <- BITS2$der2c.derp2.derp2
    der2c.derp1.derp2[teta.ind2] <- BITS2$der2c.derp1.derp2
    der2c.derp1.derrho[teta.ind2] <- BITS2$der2c.derp1.derrho
    der2c.derp2.derrho[teta.ind2] <- BITS2$der2c.derp2.derrho
  }
  if (length(teta1) != 0) 
    der2teta.derteta.stteta.st[teta.ind1] <- BITS1$der2teta.derteta.stteta.st
  if (length(teta2) != 0) 
    der2teta.derteta.stteta.st[teta.ind2] <- BITS2$der2teta.derteta.stteta.st
  der2pdf1.dereta1 <- dHs1$der2pdf2.dereta2
  der2pdf2.dereta2 <- dHs2$der2pdf2.dereta2
  der2pdf1.dersigma21.st2 <- dHs1$der2pdf2.dersigma2.st2
  der2pdf2.dersigma22.st2 <- dHs2$der2pdf2.dersigma2.st2
  der2p1.dereta1eta1 <- dHs1$der2p2.dereta2eta2
  der2p2.dereta2eta2 <- dHs2$der2p2.dereta2eta2
  der2p1.dersigma21.st2 <- dHs1$der2p2.dersigma2.st2
  der2p2.dersigma22.st2 <- dHs2$der2p2.dersigma2.st2
  der2pdf1.dereta1dersigma21.st <- dHs1$der2pdf2.dereta2dersigma2.st
  der2pdf2.dereta2dersigma22.st <- dHs2$der2pdf2.dereta2dersigma2.st
  der2p1.dereta1dersigma21.st <- dHs1$der2p2.dereta2dersigma2.st
  der2p2.dereta2dersigma22.st <- dHs2$der2p2.dereta2dersigma2.st
  d2l.be1.be1 <- -VC$weights * ((der2pdf1.dereta1 * pdf1 - 
                                   derpdf1.dereta1^2)/pdf1^2 + ((der2c.derp1.derp1 * derp1.dereta1^2 + 
                                                                   der2h.derp1p1 * der2p1.dereta1eta1) * c.copula2.be1be2 - 
                                                                  derc.dereta1^2)/c.copula2.be1be2^2)
  d2l.be2.be2 <- -VC$weights * ((der2pdf2.dereta2 * pdf2 - 
                                   derpdf2.dereta2^2)/pdf2^2 + ((der2c.derp2.derp2 * derp2.dereta2^2 + 
                                                                   der2h.derp1p2 * der2p2.dereta2eta2) * c.copula2.be1be2 - 
                                                                  derc.dereta2^2)/c.copula2.be1be2^2)
  d2l.rho.rho <- -VC$weights * (((der2c.derrho.derrho * derteta.derteta.st^2 + 
                                    der2h.derp1teta * der2teta.derteta.stteta.st) * c.copula2.be1be2 - 
                                   der2h.derp1teta.st^2)/c.copula2.be1be2^2)
  d2l.sigma21.sigma21 <- -VC$weights * ((der2pdf1.dersigma21.st2 * 
                                           pdf1 - derpdf1.dersigma21.st^2)/pdf1^2 + ((der2c.derp1.derp1 * 
                                                                                        derp1.dersigma21.st^2 + der2h.derp1p1 * der2p1.dersigma21.st2) * 
                                                                                       c.copula2.be1be2 - derc.dersigma21.st^2)/c.copula2.be1be2^2)
  d2l.sigma22.sigma22 <- -VC$weights * ((der2pdf2.dersigma22.st2 * 
                                           pdf2 - derpdf2.dersigma22.st^2)/pdf2^2 + ((der2c.derp2.derp2 * 
                                                                                        derp2.dersigma22.st^2 + der2h.derp1p2 * der2p2.dersigma22.st2) * 
                                                                                       c.copula2.be1be2 - derc.dersigma22.st^2)/c.copula2.be1be2^2)
  d2l.be1.be2 <- -VC$weights * ((der2c.derp1.derp2 * derp1.dereta1 * 
                                   derp2.dereta2 * c.copula2.be1be2 - derc.dereta1 * derc.dereta2)/c.copula2.be1be2^2)
  d2l.be1.rho <- -VC$weights * ((der2c.derp1.derrho * derp1.dereta1 * 
                                   derteta.derteta.st * c.copula2.be1be2 - derc.dereta1 * 
                                   der2h.derp1teta * derteta.derteta.st)/c.copula2.be1be2^2)
  d2l.be2.rho <- -VC$weights * ((der2c.derp2.derrho * derp2.dereta2 * 
                                   derteta.derteta.st * c.copula2.be1be2 - derc.dereta2 * 
                                   der2h.derp1teta * derteta.derteta.st)/c.copula2.be1be2^2)
  d2l.be1.sigma21 <- -VC$weights * ((der2pdf1.dereta1dersigma21.st * 
                                       pdf1 - derpdf1.dereta1 * derpdf1.dersigma21.st)/pdf1^2 + 
                                      ((der2c.derp1.derp1 * derp1.dereta1 * derp1.dersigma21.st + 
                                          der2h.derp1p1 * der2p1.dereta1dersigma21.st) * c.copula2.be1be2 - 
                                         derc.dereta1 * derc.dersigma21.st)/c.copula2.be1be2^2)
  d2l.be2.sigma22 <- -VC$weights * ((der2pdf2.dereta2dersigma22.st * 
                                       pdf2 - derpdf2.dereta2 * derpdf2.dersigma22.st)/pdf2^2 + 
                                      ((der2c.derp2.derp2 * derp2.dereta2 * derp2.dersigma22.st + 
                                          der2h.derp1p2 * der2p2.dereta2dersigma22.st) * c.copula2.be1be2 - 
                                         derc.dereta2 * derc.dersigma22.st)/c.copula2.be1be2^2)
  d2l.be2.sigma21 <- -VC$weights * ((der2c.derp1.derp2 * derp2.dereta2 * 
                                       derp1.dersigma21.st * c.copula2.be1be2 - derc.dereta2 * 
                                       derc.dersigma21.st)/c.copula2.be1be2^2)
  d2l.be1.sigma22 <- -VC$weights * ((der2c.derp1.derp2 * derp1.dereta1 * 
                                       derp2.dersigma22.st * c.copula2.be1be2 - derc.dereta1 * 
                                       derc.dersigma22.st)/c.copula2.be1be2^2)
  d2l.rho.sigma21 <- -VC$weights * ((der2c.derp1.derrho * derp1.dersigma21.st * 
                                       derteta.derteta.st * c.copula2.be1be2 - derc.dersigma21.st * 
                                       der2h.derp1teta * derteta.derteta.st)/c.copula2.be1be2^2)
  d2l.rho.sigma22 <- -VC$weights * ((der2c.derp2.derrho * derp2.dersigma22.st * 
                                       derteta.derteta.st * c.copula2.be1be2 - derc.dersigma22.st * 
                                       der2h.derp1teta * derteta.derteta.st)/c.copula2.be1be2^2)
  d2l.sigma21.sigma22 <- -VC$weights * ((der2c.derp1.derp2 * 
                                           derp1.dersigma21.st * derp2.dersigma22.st * c.copula2.be1be2 - 
                                           derc.dersigma21.st * derc.dersigma22.st)/c.copula2.be1be2^2)
  if (is.null(VC$X3)) {
    G <- -c(colSums(c(dl.dbe1) * VC$X1), colSums(c(dl.dbe2) * 
                                                   VC$X2), sum(dl.dsigma21.st), sum(dl.dsigma22.st), 
            sum(dl.dteta.st))
    be1.be1 <- crossprod(VC$X1 * c(d2l.be1.be1), VC$X1)
    be2.be2 <- crossprod(VC$X2 * c(d2l.be2.be2), VC$X2)
    be1.be2 <- crossprod(VC$X1 * c(d2l.be1.be2), VC$X2)
    be1.rho <- t(t(rowSums(t(VC$X1 * c(d2l.be1.rho)))))
    be1.sigma21 <- t(t(rowSums(t(VC$X1 * c(d2l.be1.sigma21)))))
    be1.sigma22 <- t(t(rowSums(t(VC$X1 * c(d2l.be1.sigma22)))))
    be2.rho <- t(t(rowSums(t(VC$X2 * c(d2l.be2.rho)))))
    be2.sigma21 <- t(t(rowSums(t(VC$X2 * c(d2l.be2.sigma21)))))
    be2.sigma22 <- t(t(rowSums(t(VC$X2 * c(d2l.be2.sigma22)))))
    H <- rbind(cbind(be1.be1, be1.be2, be1.sigma21, be1.sigma22, 
                     be1.rho), cbind(t(be1.be2), be2.be2, be2.sigma21, 
                                     be2.sigma22, be2.rho), cbind(t(be1.sigma21), t(be2.sigma21), 
                                                                  sum(d2l.sigma21.sigma21), sum(d2l.sigma21.sigma22), 
                                                                  sum(d2l.rho.sigma21)), cbind(t(be1.sigma22), t(be2.sigma22), 
                                                                                               sum(d2l.sigma21.sigma22), sum(d2l.sigma22.sigma22), 
                                                                                               sum(d2l.rho.sigma22)), cbind(t(be1.rho), t(be2.rho), 
                                                                                                                            sum(d2l.rho.sigma21), sum(d2l.rho.sigma22), sum(d2l.rho.rho)))
  }
  if (!is.null(VC$X3)) {
    G <- -c(colSums(c(dl.dbe1) * VC$X1), colSums(c(dl.dbe2) * 
                                                   VC$X2), colSums(c(dl.dsigma21.st) * VC$X3), colSums(c(dl.dsigma22.st) * 
                                                                                                         VC$X4), colSums(c(dl.dteta.st) * VC$X5))
    be1.be1 <- crossprod(VC$X1 * c(d2l.be1.be1), VC$X1)
    be2.be2 <- crossprod(VC$X2 * c(d2l.be2.be2), VC$X2)
    be1.be2 <- crossprod(VC$X1 * c(d2l.be1.be2), VC$X2)
    be1.rho <- crossprod(VC$X1 * c(d2l.be1.rho), VC$X5)
    be2.rho <- crossprod(VC$X2 * c(d2l.be2.rho), VC$X5)
    be1.sigma21 <- crossprod(VC$X1 * c(d2l.be1.sigma21), 
                             VC$X3)
    be1.sigma22 <- crossprod(VC$X1 * c(d2l.be1.sigma22), 
                             VC$X4)
    be2.sigma21 <- crossprod(VC$X2 * c(d2l.be2.sigma21), 
                             VC$X3)
    be2.sigma22 <- crossprod(VC$X2 * c(d2l.be2.sigma22), 
                             VC$X4)
    sigma21.sigma21 <- crossprod(VC$X3 * c(d2l.sigma21.sigma21), 
                                 VC$X3)
    sigma21.sigma22 <- crossprod(VC$X3 * c(d2l.sigma21.sigma22), 
                                 VC$X4)
    rho.sigma21 <- crossprod(VC$X3 * c(d2l.rho.sigma21), 
                             VC$X5)
    sigma22.sigma22 <- crossprod(VC$X4 * c(d2l.sigma22.sigma22), 
                                 VC$X4)
    rho.sigma22 <- crossprod(VC$X4 * c(d2l.rho.sigma22), 
                             VC$X5)
    rho.rho <- crossprod(VC$X5 * c(d2l.rho.rho), VC$X5)
    H <- rbind(cbind(be1.be1, be1.be2, be1.sigma21, be1.sigma22, 
                     be1.rho), cbind(t(be1.be2), be2.be2, be2.sigma21, 
                                     be2.sigma22, be2.rho), cbind(t(be1.sigma21), t(be2.sigma21), 
                                                                  sigma21.sigma21, sigma21.sigma22, rho.sigma21), cbind(t(be1.sigma22), 
                                                                                                                        t(be2.sigma22), t(sigma21.sigma22), sigma22.sigma22, 
                                                                                                                        rho.sigma22), cbind(t(be1.rho), t(be2.rho), t(rho.sigma21), 
                                                                                                                                            t(rho.sigma22), rho.rho))
  }
  res <- -sum(l.par)
  if (VC$extra.regI == "pC") 
    H <- regH(H, type = 1)
  S.h <- ps$S.h
  if (length(S.h) != 1) {
    S.h1 <- 0.5 * crossprod(params, S.h) %*% params
    S.h2 <- S.h %*% params
  }
  else S.h <- S.h1 <- S.h2 <- 0
  S.res <- res
  res <- S.res + S.h1
  G <- G + S.h2
  H <- H + S.h
  if (VC$extra.regI == "sED") 
    H <- regH(H, type = 2)
  if (VC$margins[1] == "LN" || VC$margins[2] == "LN") {
    if (VC$margins[1] == "LN") 
      dHs1 <- distrHsAT(exp(respvec$y1), eta1, sigma21, 
                        1, margin2 = VC$margins[1], min.dn = VC$min.dn, 
                        min.pr = VC$min.pr, max.pr = VC$max.pr)
    if (VC$margins[2] == "LN") 
      dHs2 <- distrHsAT(exp(respvec$y2), eta2, sigma22, 
                        1, margin2 = VC$margins[2], min.dn = VC$min.dn, 
                        min.pr = VC$min.pr, max.pr = VC$max.pr)
    pdf1 <- dHs1$pdf2
    pdf2 <- dHs2$pdf2
    p1 <- dHs1$p2
    p2 <- dHs2$p2
    if (length(teta1) != 0) 
      dH1 <- copgHsAT(p1[teta.ind1], p2[teta.ind1], teta1, 
                      Cop1, Ln = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                      max.pr = VC$max.pr)
    if (length(teta2) != 0) 
      dH2 <- copgHsAT(p1[teta.ind2], p2[teta.ind2], teta2, 
                      Cop2, Ln = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                      max.pr = VC$max.pr)
    c.copula2.be1be2 <- NA
    if (length(teta1) != 0) 
      c.copula2.be1be2[teta.ind1] <- dH1$c.copula2.be1be2
    if (length(teta2) != 0) 
      c.copula2.be1be2[teta.ind2] <- dH2$c.copula2.be1be2
    l.ln <- -sum(VC$weights * (log(pdf1) + log(pdf2) + log(c.copula2.be1be2)))
  }
  list(value = res, gradient = G, hessian = H, S.h = S.h, S.h1 = S.h1, 
       S.h2 = S.h2, l = S.res, l.ln = l.ln, l.par = l.par, ps = ps, 
       eta1 = eta1, eta2 = eta2, etad = etad, etas1 = etas1, 
       etas2 = etas2, BivD = VC$BivD, p1 = p1, p2 = p2, pdf1 = pdf1, 
       pdf2 = pdf2, c.copula.be2 = c.copula.be2, c.copula.be1 = c.copula.be1, 
       c.copula2.be1be2 = c.copula2.be1be2, dl.dbe1 = dl.dbe1, 
       dl.dbe2 = dl.dbe2, dl.dsigma21.st = dl.dsigma21.st, dl.dsigma22.st = dl.dsigma22.st, 
       dl.dteta.st = dl.dteta.st, teta.ind2 = teta.ind2, teta.ind1 = teta.ind1, 
       Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2)
}

