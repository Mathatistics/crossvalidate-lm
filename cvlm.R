makeFormula<-function(x.var, y.var){
    formula<-paste(y.var, paste(x.var, collapse="+"), sep="~")
    return(formula)
}
mdl.cv<-function(dataSet, x.var, y.var, step=FALSE, criteria=NULL, split=12){
    segment<-split(1:nrow(dataSet), ceiling(1:nrow(dataSet)/split))
    formula=makeFormula(x.var, y.var)
    mdl<-list()
    predVec<-rep(NA, nrow(dataSet))
    errVec<-rep(NA, nrow(dataSet))
    
    for(i in seq_along(segment)){
        dataset<-dataSet[-segment[[i]],]
        if(step){
            if(!criteria %in% c("AIC", "BIC", "Cp", "R2adj", "forward", "backward")){
                stop("Please! enter the correct criteria")
            }else{
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
                    fm.log<-capture.output({
                        mdl[[i]]<- forward(lm(formula, data=dataset), alpha = 0.05, full = FALSE)
                    })
                }else if(criteria=="backward"){
                    fm.log<-capture.output({
                        mdl[[i]]<- backward(lm(formula, data=dataset), alpha = 0.05, full = FALSE)
                    })
                }
            }
        }else{
            mdl[i]<-lm(formula, dataset)
        }
        predVec[segment[[i]]]<-predict(mdl[[i]], newdata=dataSet[segment[[i]], x.var])
        errVec[segment[[i]]]<-dataSet[segment[[i]],y.var]-predict(mdl[[i]], newdata=dataSet[segment[[i]], x.var])
    }
    rmse.cv<-sqrt(1/nrow(dataset)*sum(errVec^2))
    invisible(list(Model=mdl, Predicted=predVec, Error=errVec, rmsep=rmse.cv))
}
