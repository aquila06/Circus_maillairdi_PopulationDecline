
compute_distri<-function(model, dfpred, time_var="year", factor_name=NA, nreplicates=1000,...){
			
			
			#Check that the model was run with parameters estimated through REML
			if(!model$method %in% c("REML", "ML")) print("Predictions computed through this function are using the unconditionnal=T option of the vcov function: this requires models to be fitted with method=REML or ML")
			#Some information on the variables to properly index the resulting matrix
			
			#Need to check if the factor_name isn't included in the model parameters -> this will mess up the results matrix (this could happen when the objective is to predict over a grid for example (though X-Y), as one needs to keep track of the cell id
			
			
			Xp<-predict.gam(model, newdata=dfpred, type="lpmatrix") #
			betas_mod<-coef(model)
			vcov_mod<-vcov(model, unconditionnal=T) #unconditionnal = T allows parameter estimates to be treated as uncertain at their estimates, but this requires the model to be estimated through REML or ML--> add some line of code to check that beforehand and throw a warning if needed) 
			
			#sample MVN values of beta
			betas_rep<-rmvn(n=nreplicates, betas_mod, vcov_mod)
			
			#Product 
			lp<-t(exp(Xp %*% t(betas_rep))) #
			
			
			
			#If there is a structuring factor, store the results in a list
			if(!is.na(factor_name)){
			#A list to store the results
			list_distri<-list()
			
			#Convert grouping factor if it is not
			if(!is.factor(dfpred[,factor_name])){
					dfpred[,factor_name]<-factor(dfpred[,factor_name])
					}
			
			#return the ids of each 
			index_fact<-lapply(levels(dfpred[,factor_name]), function(x) which(dfpred[,factor_name]==x))
			
			for (i in 1:length(index_fact)){
						list_distri[[i]]<-lp[,index_fact[[i]]]
						#colnames(list_distri[[i]])<-dfpred[,time_var][index_fact[[i]]]
						colnames(list_distri[[i]])<-dfpred[,time_var][index_fact[[i]]]
						}
			# names for each element of the list
			names(list_distri)<-levels(dfpred[,factor_name])
						
			#this list is returned as the basis to compute mean trend and confidene intervals
			return(list_distri)
			}
			else{ #if no structuring factor, then just name lp (with the year index) and return lp
			colnames(lp)<-dfpred[,time_var]
			return(lp)
			}
	
	}

#A function	to compute trends from a table with a grouping factor
compute_trends<-function(list_distri, factor_name=NA, ystart, yend, alpha = 0.05,...){
		
		ystart<-as.character(ystart)
		yend<-as.character(yend)
		perChange<-function(x) 100 * (x[,yend] - x[,ystart])/x[,ystart]
		absChange<-function(x) x[,yend] - x[,ystart]
		
		if(!is.data.frame(list_distri)  & is.list(list_distri) & is.na(factor_name)){
			
			#merge distributions of values for the ystart & yend in 2 vectors
			l_sta_end<-lapply(list_distri, function(x) {
								x[,c(ystart,yend)] 
								})
			
			df_sta_end<-do.call("rbind", l_sta_end)			
					trT<- perChange(df_sta_end)
					CIT<-quantile(trT, probs = c(alpha/2, 1 - alpha/2), na.rm=T)
					mean_trT<-mean(trT)		
					resT<-data.frame(mean_trend=mean_trT, CI_low=CIT[1], CI_high=CIT[2])
					
					}else if(!is.na(factor_name)){ #if there is a grouping factor
					resT<-lapply(list_distri, function(y) { #apply on the list
					trT<- perChange(y)
					CIT<-quantile(trT, probs = c(alpha/2, 1 - alpha/2), na.rm=T)
					mean_trT<-mean(trT)		
					resT<-data.frame(mean_trend=mean_trT, CI_low=CIT[1], CI_high=CIT[2])
					
					})
		resT<-do.call("rbind", resT)
		resT$id<-names(list_distri)
		resT<-resT[,c(4,1:3)]

			}else{
			trT<- perChange(list_distri)
			CIT<-quantile(trT, probs = c(alpha/2, 1 - alpha/2))
			mean_trT<-mean(trT)		
			resT<-cbind(mean_trend=mean_trT, CI_low=CIT[1], CI_high=CIT[2])	
			
			
	}	
		return(resT)
}



compute_abs_diff<-function(list_distri, factor_name=NA, ystart, yend, alpha = 0.05,...){
		
		ystart<-as.character(ystart)
		yend<-as.character(yend)
		absChange<-function(x) x[,yend] - x[,ystart]
		
		if(!is.na(factor_name)){ #if there is a grouping factor
		resT<-lapply(list_distri, function(y) {
			trT<- absChange(y)
			CIT<-quantile(trT, probs = c(alpha/2, 1 - alpha/2))
			mean_trT<-mean(trT)		
			resT<-data.frame(mean_diff=mean_trT, CI_low=CIT[1], CI_high=CIT[2])
			
			})
		resT<-do.call("rbind", resT)
		resT$id<-names(list_distri)
		resT<-resT[,c(4,1:3)]

			}else{
			trT<- absChange(list_distri)
			CIT<-quantile(trT, probs = c(alpha/2, 1 - alpha/2))
			mean_trT<-mean(trT)		
			resT<-cbind(mean_diff=mean_trT, CI_low=CIT[1], CI_high=CIT[2])	
			
			
	}	
		return(resT)
}



#abundance to ggplot: transform the lists from the computed distri function
abund2ggplot<-function(list_distrib, alpha=0.05,...){
	
		
		temp<-lapply(list_distrib, function(x) {
		 
		 
		 lT<-apply(x, 2, function(y){
				
			mean<-mean(y)
			CIs<-quantile(y, probs = c(alpha/2, 1 - alpha/2))
		    return(rbind(mean, CIs[1], CIs[2]))
			}
		)
		}
		)
		tempdf<-lapply(temp, function(z) {
		 dfT<-data.frame(t(z))
		return(dfT)
		}
		)
		
		tempr<-do.call("rbind", tempdf)
		tempr$year<-as.numeric(dimnames(list_distrib[[1]])[[2]])
		tempr$factor_id<- rep(names(list_distrib), each=dim(list_distrib[[1]])[2])
		names(tempr)[1:3]<-c("mean", "CI_low", "CI_high")
		return(tempr)
		
		
		
}