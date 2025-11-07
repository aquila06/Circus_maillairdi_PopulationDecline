require(mgcv)
require(ggplot2)
require(ncf)
require(dplyr)
require(tidyverse)
require(flextable)
require(ggplot2)
require(ggpubr)
require(gratia)
require(sf)
require(viridis)
require(ggridges)

#Set the working directory
#setwd("/Path/to/the/folder")
setwd("P:\\MS\\MaillardH_Comptages\\Submission_bioRxiv\\")
dirout<-if(dir.exists()"./CodeOutputs/"


# A home-made function to convert the ncf::spline_correlog() object for ggplot2 plotting
source("./Scripts/ncf_sp_correl2ggplot.r")



#ggplot spatial theme for thee figures
theme_sp<-theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          #panel.grid.major=element_blank(),
		  plot.background=element_rect(fill="white"),
		  panel.background = element_rect(fill="transparent", colour="black"),
		  text=element_text(size=20) )


#Labels for the facet_wrap 
new_labels1<-c("1"="Census 1", "2"="Census 2", "3"="Census 3")
new_labels2<-c("1"="1998-2000", "2"="2009-2010", "3"="2017-2019")
new_labels3<-c("1"="1998", "2"="2009", "3"="2017")



# A path to store the elements produced by the code

#Load the data
load("./Data/Data_MS_CircusMaillairdi.RData")


#Run the models


#Keep the Tweedie distribution (as it includes the Poisson, even if runtime is a bit longer)

#Step: searching for a sufficient number of knots for the spatial term

m_full01 <-gam(ptot~surveyF + s(yearF, bs="re") + s(id_newF, bs="re") + s(log_length, k=5) + s(jd_s, bs="cc", k=5) +  s(X,Y, k=50), family=tw, data=df_pts_full, method="REML")
gam.check(m_full01)
# smoother for X,Y doesn't have enough dof... increase it s(X,Y) : edf ~ 19

m_full02 <-gam(ptot~surveyF + s(yearF, bs="re") + s(id_newF, bs="re") + s(log_length, k=5) + s(jd_s, bs="cc", k=5) +  s(X,Y, k=80), family=tw, data=df_pts_full, method="REML")# 
gam.check(m_full02)
# s(X,Y) : edf ~ 22


#Model GI: A spatial smoother common to all survey and a survey factor specific intercept (can be modelled as a RE given it only has 3 levels)
m_full03 <-gam(ptot~surveyF + s(yearF, bs="re") + s(id_newF, bs="re") + s(log_length, k=5) + s(jd_s, bs="cc", k=5) +  s(X,Y, k=100), family=tw, data=df_pts_full, method="REML")# 
gam.check(m_full03)
appraise(m_full03)

# A common smoother for the 3 surveys plus a separate (independant, no constrain on it) smoother for each level.
m_full05<-gam(ptot~surveyF +  s(yearF, bs="re") + s(id_newF, bs="re") + s(log_length, k=5) + s(jd_s, bs="cc", k=5) +  s(X,Y, k=100) + s(X,Y, by=surveyF, m=1, k=100), family=tw, data=df_pts_full, method="REML")
summary(m_full05)
appraise(m_full05, point_col="steelblue", point_alpha=0.4)

# Model bs="fs" --> a spatial term s(X,Y) as well as group-level smooth terms for each survey similar to RE on slopes in GLMM. No need to add the factor survey in the formula, each group will has its own intercept as part of the penalty construction (see Pedersen et al. 2019)
m_full06<-gam(ptot~  s(yearF, bs="re") + s(id_newF, bs="re") + s(log_length, k=5) + s(jd_s, bs="cc", k=5) +  s(X,Y, k=100) + s(X,Y, surveyF, bs="fs", k=100), family=tw, data=df_pts_full, method="REML")
summary(m_full06) 
appraise(m_full06, point_col="steelblue", point_alpha=0.4)

#Best Model 
m_full04<-gam(ptot~  s(yearF, bs="re") + s(id_newF, bs="re") + s(log_length, k=5) + s(jd_s, bs="cc", k=5) +  s(X,Y, k=100) + s(X,Y,surveyF, bs="sz", k=100), family=tw, data=df_pts_full, method="REML")# 
summary(m_full04)
m04_fit<-appraise(m_full04)

ggsave(paste0(dirout, "/Best_model_fit.jpeg"), plot=m04_fit,  units="cm",  width=18, height=18, dpi="retina")	

#Save as a docx table
sum_bestm_ft<-flextable::as_flextable(m_full04) 
sum_bestm_ft %>%
	colformat_double(., j=3:6, digits=2)%>%
		flextable::save_as_docx(path=paste0(dirout, "Summary_GAM_papangue.docx" ))
		%>%
			flextable::save_as_image((., paste0(dirout, "TableSummaryGAM.png"), expand = 10, res = 300)

#Assess if residuals are problematic, spatially speaking

# A spatial plot of the residuals
resid_df<-m_full04$model
resid_df$resids<-resid(m_full04)


resid_xy_gam<-ggplot(resid_df) + geom_point(aes(x=X, y=Y, size=abs(resids), fill=resids>0, colour=resids>0), alpha=0.9, shape=21) + 
scale_colour_manual(name = "Sign of the residual", values = setNames(c('blue','red'),c(T, F)), labels=c('<0', '>0')) + 
scale_fill_manual(name = "Sign of the residual", values = setNames(c('blue','red'),c(T, F)), labels=c('<0', '>0')) +  geom_jitter(aes(x=X, y=Y, size=abs(resids), fill=resids>0, colour=resids>0), width=200, height=200) + 
scale_size_continuous(name="abs(residual)") + facet_wrap(~surveyF, labeller=labeller(surveyF=new_labels2)) + coord_fixed() +
geom_sf(data=lim_reu, fill="transparent") + 
theme_sp + theme(legend.position="bottom", legend.key=element_rect(colour = NA, fill = "transparent"))

ggsave(paste0(dirout, "/Spatial_residuals.jpeg"), plot=resid_xy_gam,  units="cm",  width=28, height=12, dpi="retina")	
	
#Check the residuals of the different models (using the above resids_df)

sp_correl_surv1<-spline.correlog(resid_df$X[resid_df$surveyF=="1"], resid_df$Y[resid_df$surveyF=="1"], resid_df$resids[resid_df$surveyF=="1"], xmax=20000)
sp_correl_surv2<-spline.correlog(resid_df$X[resid_df$surveyF=="2"], resid_df$Y[resid_df$surveyF=="2"], resid_df$resids[resid_df$surveyF=="2"], xmax=20000)
sp_correl_surv3<-spline.correlog(resid_df$X[resid_df$surveyF=="3"], resid_df$Y[resid_df$surveyF=="3"], resid_df$resids[resid_df$surveyF=="3"], xmax=20000)
sp_correl_survall<-spline.correlog(resid_df$X, resid_df$Y, resid_df$resids, xmax=20000)

#Convert for plotting with ggplot2
df_surv1<-extract_CI_sp_correl(sp_correl_surv1)
df_surv1$survey<-"Survey 1"
df_surv2<-extract_CI_sp_correl(sp_correl_surv2)
df_surv2$survey<-"Survey 2"
df_surv3<-extract_CI_sp_correl(sp_correl_surv3)
df_surv3$survey<-"Survey 3"
df_surv_all<-extract_CI_sp_correl(sp_correl_survall)
df_surv_all$survey<-"All surveys"

df_correl_ggplot<-rbind(df_surv1, df_surv2, df_surv3, df_surv_all)

df_correl_ggplot$ic_high[which(df_correl_ggplot$ic_high>1)]<-1
					
												
					
#Create a plot the correlograms of the residuals
residual_correlograms<-ggplot(df_correl_ggplot, aes(x=x, y=y)) + 	geom_hline(yintercept=0, linewidth=0.25, alpha=0.5) +
			geom_line(color="red") + geom_ribbon(aes(ymin=ic_low, ymax=ic_high), alpha=0.3, fill="red") + 
			ylab("Correlation") + xlab("Distance (in meters)") + ylim(c(-1,1)) + facet_wrap(~survey, ncol=2)
ggsave(paste0(dirout, "/Correlograms.jpeg"), plot=residual_correlograms,  units="cm",  width=17, height=10, dpi="retina")	
			
#Assess overdisperion
sum(residuals(m_full04, type = "pearson")^2) / df.residual(m_full04)


#Save a RData for later (or to generate the figures) (~200Mo !!!)		
save.image(paste0(dirout, "ModelsRun_CircusMaillardi.RData"))
