makeFormula<-function(x.var, y.var){
    formula<-as.formula(paste(y.var, paste(x.var, collapse="+"), sep="~"))
    return(formula)
}
mdl.cv<-function(dataSet, x.var, y.var, model="lm", step=FALSE, criteria=NULL, split=12, ncomp=NULL){
    segment<-split(1:nrow(dataSet), ceiling(1:nrow(dataSet)/split))
    formula=makeFormula(x.var, y.var)
    mdl<-list()
    predVec<-rep(NA, nrow(dataSet))
    errVec<-rep(NA, nrow(dataSet))
    
    for(i in seq_along(segment)){
        dataset<-dataSet[-segment[[i]],]
        testset<-dataSet[segment[[i]],]
        if(step & model=="lm"){
            if(!criteria %in% c("AIC", "BIC", "Cp", "R2adj", "forward", "backward")){
                stop("Please! enter the correct criteria")
            }else{
                require(leaps)
                if(criteria=="Cp"){
                    ## Model selected by Mallows Cp Criteria
                    cp.leaps<-leaps(x=dataset[,x.var],
                                    y=dataset[,y.var],
                                    method="Cp", nbest = 1, names = x.var)
                    # Model fitting
                    cp.which<-names(which(cp.leaps$which[which.min(cp.leaps$Cp),]))
                    formula<-makeFormula(cp.which, y.var)
                    mdl[[i]]<-lm(formula, data=dataset)
                }else if(criteria=="R2adj"){
                    ## Model selected by R2adj Criteria
                    r2adj.leaps<-leaps(x=dataset[,x.var],
                                       y=dataset[,y.var],
                                       method="adjr2", nbest = 1, names=x.var)
                    # Model fitting
                    r2.which<-names(which(r2adj.leaps$which[which.max(r2adj.leaps$adjr2),]))
                    formula<-makeFormula(r2.which, y.var)
                    mdl[[i]]<-lm(formula, data=dataset)
                }else if(criteria=="AIC" | criteria=="BIC"){
                    lmBstSetSmry <- summary(regsubsets(dataset[,x.var],
                                                       dataset[,y.var], 
                                                       nbest = 1, nvmax = length(x.var)))
                    nvars<-apply(lmBstSetSmry$which, 1, sum)
                    bic.vec<-lmBstSetSmry$bic
                    aic.vec<-bic.vec-nvars*log(sum(train))+nvars
                    
                    ## Fitting selected linear model
                    aic.which<-names(which(lmBstSetSmry$which[which.min(aic.vec),]))[-1]
                    bic.which<-names(which(lmBstSetSmry$which[which.min(bic.vec),]))[-1]
                    if(criteria=="AIC"){
                        formula<-makeFormula(aic.which, y.var)
                        mdl[[i]]<-lm(formula, data=dataset)
                    }else if(criteria=="BIC"){
                        formula<-makeFormula(bic.which, y.var)
                        mdl[[i]]<-lm(formula, data=dataset)
                    }
                }else if(criteria=="forward"){
                    require(mixlm)
                    fm.log<-capture.output({
                        mdl[[i]]<- forward(do.call(lm, list(formula, dataset)), alpha = 0.05, full = FALSE)
                    })
                }else if(criteria=="backward"){
                    require(mixlm)
                    fm.log<-capture.output({
                        mdl[[i]]<- backward(do.call(lm, list(formula, dataset)), alpha = 0.05, full = FALSE)
                    })
                }
            }
        }else if(step & model!='lm'){
            stop("Stepwise can only be performed using Linear Model, Please input 'lm' in the model.")
        }else if(model=='lm'){
            mdl[[i]]<-lm(formula, dataset)
        }else if(model=='ridge'){
            require(ridge)
            mdl[[i]]<- linearRidge(formula, dataset)
        }else if(model=="pls" | model=="pcr"){
            require(pls)
            if(model=="pls"){
                mdl[[i]]<-plsr(formula, data=dataset, scale=TRUE)
            }else{
                mdl[[i]]<-pcr(formula, data=dataset, scale=TRUE)
            }
            if(!is.null(ncomp) & is.numeric(ncomp)){
                predVec[segment[[i]]]<-predict(mdl[[i]], newdata=testset[,x.var], ncomp=ncomp)[,,]
                errVec[segment[[i]]]<-testset[,y.var]-predVec[segment[[i]]]
                next
            }
            else{
                stop("`ncomp' is needed for PLS and PCR prediction and it should be numeric.")
            }
        }else{
            stop("Model can take 'lm' or 'ridge' value.")
        }
        predVec[segment[[i]]]<-predict(mdl[[i]], newdata=testset[,x.var])
        errVec[segment[[i]]]<-testset[,y.var]-predVec[segment[[i]]]
    }
    rmse.cv<-sqrt(1/nrow(dataSet)*sum(errVec^2))
    r2pred<-1-sum(errVec^2)/sum((predVec-mean(dataSet[,y.var]))^2)
    invisible(list(Model=mdl, Predicted=predVec, Error=errVec, rmsep=rmse.cv, r2pred=r2pred))
}