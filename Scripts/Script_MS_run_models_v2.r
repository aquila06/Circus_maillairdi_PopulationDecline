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
require(terra)
require(viridis)
require(ggridges)

#Set the working directory
#setwd("/Path/to/the/folder")
setwd("P:\\MS\\MaillardH_Comptages\\Submission_bioRxiv\\v2\\")
dirout<-"./CodeOutputs/"
if(!dir.exists(dirout)) dir.create(dirout)


# A home-made function to convert the ncf::spline_correlog() object for ggplot2 plotting purposes
source("./Scripts/ncf_sp_correl2ggplot.r")

#Read the limits of the islands + areas (leeward and windward) as a csv
ew_lim<-read.table("./Data/coords_areas.csv", sep=";", h=T) %>%
			st_as_sf(coords = c("X", "Y"), crs = 2975) %>%
				group_by(REGION) %>%
					summarise(geometry = st_combine(geometry)) %>%
						st_cast("POLYGON")	

# Can be done directly from the provided shapefile
# ew_lim<-st_read("./Data/Windward_leeward_reunion.shp")


#The study area
lim_reu<- ew_lim %>%
			st_union()

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

#For R users, load the RData
#load("./Data/Data_MS_CircusMaillairdi.RData")

# Else read the csv file and convert columns to appropriate format
df_pts_full<-read.table("./Data/df_pts_full.csv", sep=";", h=T)

#convert to factors
df_pts_full$surveyF<-as.factor(df_pts_full$surveyF)
df_pts_full$id_newF<-as.factor(df_pts_full$id_newF)
df_pts_full$yearF<-as.factor(df_pts_full$yearF)

############################
##### Run the models #######
############################

#Keep the Tweedie distribution (as it includes the Poisson, even if runtime is a bit longer)

#Step: searching for a sufficient number of knots for the spatial term

set.seed(1234)

if(!file.exists("./Data/ModelsRun_CircusMaillardi.RData")){




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

#Best Model (see text in the MS)
m_full04<-gam(ptot~  s(yearF, bs="re") + s(id_newF, bs="re") + s(log_length, k=5) + s(jd_s, bs="cc", k=5) +  s(X,Y, k=100) + s(X,Y,surveyF, bs="sz", k=100), family=tw, data=df_pts_full, method="REML")# 
summary(m_full04)
m04_fit<-appraise(m_full04, point_col="steelblue", point_alpha=0.4)


save.image("./Data/ModelsRun_CircusMaillardi_fall_surv.RData")
} else {


	load("./Data/ModelsRun_CircusMaillardi.RData")

}


# a diag plot of the model
plot_appraise<-appraise(m_full04)
ggsave(paste0(dirout, "Best_mod_diag_plot.jpeg" ), plot=plot_appraise,  units="cm",  width=17, height=17, dpi="retina")



set.seed(1234)
#Save the summary as a docx table
sum_bestm_ft<-flextable::as_flextable(m_full04) 
sum_bestm_ft %>%
	colformat_double(., j=3:6, digits=2)%>%
		flextable::save_as_docx(path=paste0(dirout, "Table2_Summary_GAM_papangue.docx" ))

sum_bestm_ft %>%
	colformat_double(., j=3:6, digits=2)%>%			
			flextable::save_as_image(.,paste0(dirout, "Table2_TableSummaryGAM.png"), expand = 10, res = 300)

#Assess if residuals are problematic, spatially speaking

# A spatial plot of the residuals
resid_df<-m_full04$model
resid_df$resids<-resid(m_full04)


# A figure to get a grasp of residual autocorrelation
resid_xy_gam<-ggplot(resid_df) + geom_point(aes(x=X, y=Y, size=abs(resids), fill=resids>0, colour=resids>0), alpha=0.9, shape=21) + 
scale_colour_manual(name = "Sign of the residual", values = setNames(c('blue','red'),c(T, F)), labels=c('<0', '>0')) + 
scale_fill_manual(name = "Sign of the residual", values = setNames(c('blue','red'),c(T, F)), labels=c('<0', '>0')) +  geom_jitter(aes(x=X, y=Y, size=abs(resids), fill=resids>0, colour=resids>0), width=200, height=200) + 
scale_size_continuous(name="abs(residual)") + facet_wrap(~surveyF, labeller=labeller(surveyF=new_labels2)) + coord_fixed() +
geom_sf(data=lim_reu, fill="transparent") + 
theme_sp + theme(legend.position="bottom", legend.key=element_rect(colour = NA, fill = "transparent"))

ggsave(paste0(dirout, "/Spatial_residuals.jpeg"), plot=resid_xy_gam,  units="cm",  width=28, height=12, dpi="retina")	
	
#Compute the residuals for each survey

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
					
												
					
# A plot with the correlograms of the residuals for the different surveys
residual_correlograms<-ggplot(df_correl_ggplot, aes(x=x, y=y)) + 	geom_hline(yintercept=0, linewidth=0.25, alpha=0.5) +
			geom_line(color="red") + geom_ribbon(aes(ymin=ic_low, ymax=ic_high), alpha=0.3, fill="red") + 
			ylab("Correlation") + xlab("Distance (in meters)") + ylim(c(-1,1)) + facet_wrap(~survey, ncol=2)
ggsave(paste0(dirout, "/Correlograms.jpeg"), plot=residual_correlograms,  units="cm",  width=17, height=10, dpi="retina")	
			
#Assess overdisperion
sum(residuals(m_full04, type = "pearson")^2) / df.residual(m_full04)


#Eventually save a RData for later (or to generate the figures) 	
#save.image(paste0(dirout, "ModelsRun_CircusMaillardi.RData"))

#Read the data

# a figure 2 rows of 3 panels showing the sampling scheme in space for the complete data set and the points surveyed during the breeding season.
df_all <- df_pts_full %>%
  mutate(Data = "All data")

# Bloc 2 : uniquement avril-juin
df_amj <- df_pts_full %>%
  filter(month %in% 4:6) %>%
  mutate(Data ="Breeding period (Apr-June)")

df_restofyear <- df_pts_full %>%
  filter(!month %in% 4:6) %>%
  mutate(Data = "Outside the main breeding period")


# Empilement des deux blocs
df_combined <- bind_rows(df_all, df_amj, df_restofyear) %>%
  mutate(Data = factor(Data, levels = c("All data", "Breeding period (Apr-June)", "Outside the main breeding period")))

pl_sp_sampl_supp_mat<-ggplot(df_combined, aes(x = X, y = Y)) +
  geom_point(alpha = 0.6, size = 1.5) +
	facet_grid(Data ~ survey_name) +
		coord_equal() +
		geom_contour(data=mnt_reu_df, aes(x=X, y=Y, z=value), binwidth=100, colour="black", alpha=0.2, linewidth=0.2) + 
			theme_void()


ggsave(paste0(dirout, "Spatial_sampling_supp_mat.jpeg"), plot=pl_sp_sampl_supp_mat,  units="cm",  width=24, height=24, dpi="retina")


		

#MNT for plotting purpose with ggplot2 (though the
mnt_reu<-rast("./Data/reu_mnt_tiles_res100m.tif")
mnt_reu_df<-as.data.frame(mnt_reu, xy=T)
colnames(mnt_reu_df)<-c("X", "Y", "value")

#A grid of 1 km x 1 km over which to compute the predictions
reu_grd_ref<-rast(extent=vect(lim_reu), resolution=c(1000,1000))
reu_grid<-rasterize(vect(lim_reu),  reu_grd_ref)

#The base df for predictions
reu_df_pred<-as.data.frame(reu_grid, xy=T)
colnames(reu_df_pred)[c(1,2)]<-c("X", "Y")
reu_df_pred$id<-as.character(1:nrow(reu_df_pred))

years<-c(1998:2000, 2009:2010, 2017:2019)

surveys<-levels(df_pts_full$surveyF)

#Rbind so to get a prediction for each survey over the whole island

all_surveys_pred<-do.call("rbind", replicate(length(surveys), reu_df_pred, simplify=F))
all_surveys_pred$log_length<-mean(df_pts_full$log_length)
all_surveys_pred$jd_s<-0
all_surveys_pred$id_newF<-"550"#levels(df_pts_full$id_newF)[1]
all_surveys_pred$yearF<-"2023"#levels(m_full04$model$yearF)[2]
all_surveys_pred$surveyN<-rep(surveys, each=nrow(reu_df_pred))
all_surveys_pred$surveyF<-factor(all_surveys_pred$surveyN)
all_surveys_pred$id<-rep(1:nrow(reu_df_pred), length(surveys))

#Regionalise the trend bas to compare trends in the east and in the west
l_pix_areas<-ew_lim %>%		
					st_intersects(st_as_sf(reu_df_pred, coords=c("X","Y"), crs=2975), .)	
df_pix_areas<-data.frame(id=as.numeric(attributes(l_pix_areas)["region.id"]$region.id), area=as.numeric(unlist(l_pix_areas)))			
					
	
# The home made functions to predict in space and time (from a spatio-temporal model)
source("./Scripts/compute_spatio_temp_trends_gam_vcov.r")


#Predicions for plotting purpose
pred_abund<-predict(m_full04, newdata=all_surveys_pred, type="response", exclude=c("s(yearF)", "s(id_newF)", "s(jd_s)", "s(log_length)"), se.fit=T)
all_surveys_pred$pred_abund04<-as.numeric(pred_abund$fit)
all_surveys_pred$pred_abund04_se<-as.numeric(pred_abund$se.fit)
pred_sz<-predict(m_full04, newdata=all_surveys_pred, type="iterms", exclude=c("s(yearF)", "s(id_newF)", "s(jd_s)", "s(log_length)","s(X,Y)"), se.fit=T)
all_surveys_pred$sz_term<-as.numeric(pred_sz$fit)
all_surveys_pred$sz_term_se<-as.numeric(pred_sz$se.fit)
pred_sxy<-predict(m_full04, newdata=all_surveys_pred, type="iterms", exclude=c("s(yearF)", "s(id_newF)", "s(jd_s)", "s(log_length)","s(X,Y,surveyF)"), se.fit=T)
all_surveys_pred$sxy_alone<-as.numeric(pred_sxy$fit)
all_surveys_pred$sxy_alone_se<-as.numeric(pred_sxy$se.fit)

#The spatial term alone
plot_sxy_alone<-    ggplot(all_surveys_pred) + geom_raster(aes(x=X, y=Y, fill=sxy_alone)) + coord_fixed() + 
								scale_fill_viridis(name="Partial effect",  guide=guide_colorbar(barwidth=unit(5, "cm")) ) + 
									geom_sf(data=lim_reu, fill="transparent", size=0.5) +
									theme_void()+ 
									 theme(legend.position="bottom", plot.title=element_text(size=10, h=0.5,  face="bold", vjust=2)) + ggtitle("s(X,Y)")

plot_sxy_alone_se<-ggplot(all_surveys_pred) + geom_raster(aes(x=X, y=Y, fill=sxy_alone_se)) + coord_fixed() + 
								scale_fill_viridis(name="SE",  guide=guide_colorbar(barwidth=unit(5, "cm")) ) + 
									geom_sf(data=lim_reu, fill="transparent", size=0.5) +
									theme_void()+ 
									 theme(legend.position="bottom", plot.title=element_text(size=10, h=0.5,  face="bold", vjust=2)) + 
									 ggtitle("Standard Error of the s(X,Y) smoother")

#A plot of the different terms of the model
plot_sz_terms<-ggplot(all_surveys_pred) + geom_raster(aes(x=X, y=Y, fill=sz_term)) + facet_wrap(~surveyF, labeller=labeller(surveyF=new_labels2)) + coord_fixed() + 
							scale_fill_viridis(name="Partial effect", guide=guide_colorbar(barwidth=unit(5, "cm")) )  +
							geom_sf(data=lim_reu, fill="transparent", size=0.5) +
									 theme_void()+ 
									 
									 theme(legend.position="bottom", plot.title=element_text(size=10, h=0.5,  face="bold", vjust=2)) + ggtitle("s(X,Y, surveyF)")


plot_sz_terms_se<-ggplot(all_surveys_pred) + geom_raster(aes(x=X, y=Y, fill=sz_term_se)) + facet_wrap(~surveyF, labeller=labeller(surveyF=new_labels2)) + coord_fixed() + 
							scale_fill_viridis(name="SE", guide=guide_colorbar(barwidth=unit(5, "cm")) ) +  
							geom_sf(data=lim_reu, fill="transparent", size=0.5) +
									theme_void()+
									 theme(legend.position="bottom", plot.title=element_text(size=10, h=0.5,  face="bold", vjust=2)) + 
									 ggtitle("Standard Error of the s(X,Y, surveyF) smoother")


#Predicted abundance
plot_pred_abund<-ggplot(all_surveys_pred) + geom_raster(aes(x=X, y=Y, fill=pred_abund04)) + 
									facet_wrap(~surveyF, labeller=labeller(surveyF=new_labels2)) + coord_fixed() + 	
									scale_fill_viridis(name="relative abundance index",  guide=guide_colorbar(barwidth=unit(5, "cm"))) +
									geom_sf(data=lim_reu, fill="transparent", size=0.5) +
									 theme_void()+ 
									 theme(legend.position="bottom", plot.title=element_text(size=10, h=0.5,  face="bold", vjust=2)) + ggtitle("Predicted abundance")

	

#The plot of the SE of the overall response
plot_se_abund<-ggplot(all_surveys_pred) + geom_raster(aes(x=X, y=Y, fill=pred_abund04_se))  + coord_fixed() +  
									scale_fill_viridis(name="SE",   guide=guide_colorbar(barwidth=unit(5, "cm"))) +
									geom_sf(data=lim_reu, fill="transparent", size=0.5) +
									geom_point(data=m_full04$model, aes(x=X, y=Y), size=0.5, colour="red") +
									facet_wrap(~surveyF, labeller=labeller(surveyF=new_labels2)) +
									theme_void()+ 
									 theme(legend.position="bottom", 
										plot.title=element_text(size=10, h=0.5,  face="bold", vjust=2)) +
								
									ggtitle("Standard error of predicted abundance")


#A figure to assess the terms of the model

p_jd<-gratia::draw(m_full04, select="s(jd_s)", ci_col="dodgerblue1",  smooth_col= "dodgerblue1", resid_col="dodgerblue1") + xlab("Julian date (z-score)") + ggtitle("s(Date)")
p_length <-gratia::draw(m_full04, select="s(log_length)", ci_col="dodgerblue1",  smooth_col= "dodgerblue1", resid_col="dodgerblue1") + xlab("Duration of observation (log-transformed)") + ggtitle("s(Observation length)")


#Save the figures

#Non sptail terms
plot_mod_ass_nonsp<-ggarrange(p_jd, p_length, ncol=2,  labels=c("a", "b"))
ggsave(paste0(dirout, "/Non_sp_terms.jpeg"), plot=plot_mod_ass_nonsp,  units="cm",  width=17, height=10, dpi="retina")

# a figure with spatiel terms of the model and perdicted abundance index
plot_all_spterms<-ggarrange(plot_sxy_alone, plot_sz_terms, plot_pred_abund,  ncol=1, labels=c("a", "b", "c"))
ggsave(paste0(dirout,  "/All_spatial_terms_means.jpeg"), plot=plot_all_spterms,  units="cm",  width=16, height=18, dpi="retina")


# a figure for the supp mat 

plot_pred_abund_se<-ggarrange(plot_pred_abund,  plot_se_abund, ncol=1, labels=c("a", "b"))
ggsave(paste0(dirout,  "/Spatial_pred_abund_se.jpeg"), plot=plot_pred_abund_se,  units="cm",  width=16, height=10, dpi="retina")





#the figures with predicted / simultaed abundance and computed trends 

#Number of replicates
nreplicates<-1000
#Distribution of abundance value
l_pred_pix_surv<-compute_distri(model=m_full04, time_var="surveyF", dfpred=all_surveys_pred, factor_name="id", nreplicates=nreplicates)




#Combine distribution per pixels --> not accurate so average abundance per survey per simulation (rep_id)
df_pix_repl<-do.call("rbind", l_pred_pix_surv) %>%
								as.data.frame(.) %>%
									mutate(id=rep(1:nrow(reu_df_pred), each=nreplicates), rep_id=rep(1:nreplicates, nrow(reu_df_pred))) %>%
								pivot_longer(cols=1:3, names_to="surveyF", values_to="Abundance") %>%
								group_by(rep_id, surveyF) %>%
								 summarise(mean_abund=mean(Abundance)) %>%
								 mutate(area="Island")

stat_mean_abund<- df_pix_repl %>%
							group_by(surveyF) %>%
								summarise(mean_abund_surv=mean(mean_abund), IC_low=as.numeric(quantile(mean_abund, probs=0.025)),
											IC_high=as.numeric(quantile(mean_abund, probs=0.975)))

plot_isl_abund<-ggplot(df_pix_repl) + geom_line(aes(x=surveyF, y=mean_abund, group=rep_id), alpha=0.1, color="dodgerblue1") +
						geom_point(data=stat_mean_abund, aes(x=surveyF, y=mean_abund_surv), colour="dodgerblue3", size=4) + 
									geom_errorbar(data=stat_mean_abund, aes(x=surveyF, ymin=IC_low, ymax=IC_high), colour="dodgerblue3", width=0.1, linewidth=1) +
									ylab("Mean abundance index (island)") + xlab("Survey") +
									scale_x_discrete(labels=c("1998-2000", "2009-2010", "2017-2019"))
									#theme(panel.background = element_blank())

#back to a wide data.frame to test the trends with the home made function
mat_mean_abund<-df_pix_repl %>%
					pivot_wider(names_from=surveyF, values_from=mean_abund) %>%
					data.frame(.)

tr_1_3<-compute_trends(mat_mean_abund, ystart="X1", yend="X3")
tr_1_2<-compute_trends(mat_mean_abund, ystart="X1", yend="X2")
tr_2_3<-compute_trends(mat_mean_abund, ystart="X2", yend="X3")

#A table for that
tab_res_trends<-data.frame(cbind(Periods=c("1998-2009","1998-2017", "2009-2017" ),area="Island"), 
								rbind(tr_1_2, tr_1_3, tr_2_3))
rownames(tab_res_trends)<-NULL

#Represent each trends in separate graph

survey1_survey2<-100 * (mat_mean_abund[,"X2"]-mat_mean_abund[,"X1"])/mat_mean_abund[,"X1"]
survey1_survey3<-100 * (mat_mean_abund[,"X3"]-mat_mean_abund[,"X1"])/mat_mean_abund[,"X1"]
survey2_survey3<-100 * (mat_mean_abund[,"X3"]-mat_mean_abund[,"X2"])/mat_mean_abund[,"X2"]
plot_dist_trends<-data.frame(area="Island", Periods=rep(c("1998 - 2009","1998 - 2017", "2009 - 2017" ), each=nreplicates), Trends=c(survey1_survey2, survey1_survey3, survey2_survey3))

dist_trends_change<-ggplot(plot_dist_trends, aes(x=Trends, y=Periods)) + 
				geom_density_ridges(fill="dodgerblue1", colour="dodgerblue1", alpha=0.7, bandwidth = 5) + 
				xlab("% of change in abundance") +  ylab("") + 
				geom_vline(xintercept = 0, linetype="dotted", color = "black", linewidth=1.5)

#A plot combining the 2 graphs

abund_trends<-ggarrange(plot_isl_abund, dist_trends_change, ncol=2, labels=c("a", "b"))



#need to merge the region ID (east vs west for each pixel) in the dplyr below
df_pix_reg_repl<-do.call("rbind", l_pred_pix_surv) %>%
								as.data.frame(.) %>%
									mutate(id=rep(1:nrow(reu_df_pred), each=nreplicates), rep_id=rep(1:nreplicates, nrow(reu_df_pred))) %>%
											left_join(., df_pix_areas, by=join_by(id==id)) %>%
													pivot_longer(cols=1:3, names_to="surveyF", values_to="Abundance") %>%
															group_by(area, rep_id, surveyF) %>%
																summarise(mean_abund=mean(Abundance)) %>%
																	mutate(area=as.character(area)) %>%
																	as.data.frame(.)
#Change names
df_pix_reg_repl$area<-recode(df_pix_reg_repl$area, "1"="West", "2"="East")
																		
#regional trends
stat_mean_reg_abund<- df_pix_reg_repl %>%
							group_by(surveyF, area) %>%
								summarise(mean_abund_surv=mean(mean_abund), IC_low=as.numeric(quantile(mean_abund, probs=0.025)),
											IC_high=as.numeric(quantile(mean_abund, probs=0.975)))

area.labs<-c("West", "East")
names(area.labs)<-c(1,2)

plot_reg_abund<-ggplot() + geom_line(data=df_pix_reg_repl, aes(x=surveyF, y=mean_abund, group=rep_id), colour="dodgerblue1", alpha=0.05) + 
						geom_point(data=stat_mean_reg_abund, aes(x=surveyF, y=mean_abund_surv, group=factor(area)), colour="dodgerblue1", size=4) + 
									geom_errorbar(data=stat_mean_reg_abund, aes(x=surveyF, ymin=IC_low, ymax=IC_high, group=factor(area)), colour="dodgerblue1",  width=0.1, linewidth=1) +
									ylab("Mean abundance index (island)") + xlab("Survey") +
									scale_x_discrete(labels=c("1998-2000", "2009-2010", "2017-2019")) + 
									facet_wrap(~area, ncol=2, labeller=labeller(area=area.labs)) + 
									theme(plot.margin = unit(c(5.5, 5.5, 5.5, 35), units="pt"))

#back to a wide data.frame to test the trensd with the home made function
mat_mean_reg_abund<-df_pix_reg_repl %>%
					group_by(area) %>%
					pivot_wider(names_from=surveyF, values_from=mean_abund) %>%
					data.frame(.)

tr_reg1_1_3<-compute_trends(mat_mean_reg_abund[which(mat_mean_reg_abund$area=="West"),], ystart="X1", yend="X3")
tr_reg1_1_2<-compute_trends(mat_mean_reg_abund[which(mat_mean_reg_abund$area=="West"),], ystart="X1", yend="X2")
tr_reg1_2_3<-compute_trends(mat_mean_reg_abund[which(mat_mean_reg_abund$area=="West"),], ystart="X2", yend="X3")

tr_reg2_1_3<-compute_trends(mat_mean_reg_abund[which(mat_mean_reg_abund$area=="East"),], ystart="X1", yend="X3")
tr_reg2_1_2<-compute_trends(mat_mean_reg_abund[which(mat_mean_reg_abund$area=="East"),], ystart="X1", yend="X2")
tr_reg2_2_3<-compute_trends(mat_mean_reg_abund[which(mat_mean_reg_abund$area=="East"),], ystart="X2", yend="X3")



#A table for that (udes later to add trends on the final graph)
tab_reg_trends<-data.frame(cbind(Periods=rep(c("1998-2009","1998-2017", "2009-2017"),  2), area=rep(c("West", "East"), each=3)),
								rbind(tr_reg1_1_2, tr_reg1_1_3,tr_reg1_2_3, tr_reg2_1_2, tr_reg2_1_3, tr_reg2_2_3))
rownames(tab_res_trends)<-NULL

#Represent each trends in separate graph

#trends by region
df_reg1<-mat_mean_reg_abund[which(mat_mean_reg_abund$area=="West"),]
df_reg2<-mat_mean_reg_abund[which(mat_mean_reg_abund$area=="East"),]
reg1_survey1_survey2<-100 * (df_reg1[,"X2"]-df_reg1[,"X1"])/df_reg1[,"X1"]
reg1_survey1_survey3<-100 * (df_reg1[,"X3"]-df_reg1[,"X1"])/df_reg1[,"X1"]
reg1_survey2_survey3<-100 * (df_reg1[,"X3"]-df_reg1[,"X2"])/df_reg1[,"X2"]

reg2_survey1_survey2<-100 * (df_reg2[,"X2"]-df_reg2[,"X1"])/df_reg2[,"X1"]
reg2_survey1_survey3<-100 * (df_reg2[,"X3"]-df_reg2[,"X1"])/df_reg2[,"X1"]
reg2_survey2_survey3<-100 * (df_reg2[,"X3"]-df_reg2[,"X2"])/df_reg2[,"X2"]


df_dist_reg1_trends<-data.frame(area = "West",  Periods=rep(c("1998 - 2009","1998 - 2017", "2009 - 2017"), each=nreplicates), Trends=c(reg1_survey1_survey2, reg1_survey1_survey3, reg1_survey2_survey3)) 
df_dist_reg2_trends<-data.frame(area = "East",  Periods=rep(c("1998 - 2009","1998 - 2017", "2009 - 2017" ), each=nreplicates), Trends=c(reg2_survey1_survey2, reg2_survey1_survey3, reg2_survey2_survey3)) 
df_dist_reg_trends<-rbind(df_dist_reg1_trends, df_dist_reg2_trends)

#Compute lambdas
lamb_1_3_isl<-1+(tr_1_3[1]/100)^1/21
lamb_1_3_east<-1+(tr_reg2_1_3[1]/100)^1/21
lamb_1_3_west<-1+(tr_reg1_1_3[1]/100)^1/21
lamb_1_2_isl<-1+(tr_1_2[1]/100)^1/12
lamb_1_2_east<-1+(tr_reg2_1_2[1]/100)^1/12
lamb_1_2_west<-1+(tr_reg1_1_2[1]/100)^1/12
lamb_2_3_isl<-1+(tr_2_3[1]/100)^1/9
lamb_2_3_east<-1+(tr_reg2_2_3[1]/100)^1/9
lamb_2_3_west<-1+(tr_reg1_2_3[1]/100)^1/9

df_lambda<-data.frame(Area= rep(c("Island", "East", "West"), 3), 
					Period=rep(c("1998-2009","1998-2017", "2009-2017"), each=3),
					Lambda=c(lamb_1_2_isl, lamb_1_2_east, lamb_1_2_west, lamb_1_3_isl, lamb_1_3_east, lamb_1_3_west, lamb_2_3_isl, lamb_2_3_east, lamb_2_3_west), row.names=NULL)
df_lambda$Yearly_pc<-(df_lambda$Lambda-1)*100
colnames(df_lambda)<-c("Area", "Period", "Yearly lambda", "Yearly %")

#export as a latex table
df_lambda %>%
flextable::as_flextable(.) %>%
	flextable::save_as_image(., paste0(dirout, "Lambdas.png"), expand = 10, res = 300)



# Now a figure for mean and simulated abundances

#Combine the df of abundance used to plot changes in abundance

df_all<-rbind(df_pix_repl, df_pix_reg_repl)
df_all$area<-factor(df_all$area, levels=c("Island", "West", "East"))

stat_all<- df_all %>%
			group_by(surveyF, area) %>%
								summarise(mean_abund_surv=mean(mean_abund), IC_low=as.numeric(quantile(mean_abund, probs=0.025)),
											IC_high=as.numeric(quantile(mean_abund, probs=0.975)))



#the a plot 
plot_all_abund<-ggplot() + geom_line(data=df_all, aes(x=surveyF, y=mean_abund, group=rep_id), colour="black", alpha=0.01) + 
						geom_point(data=stat_all, aes(x=surveyF, y=mean_abund_surv, group=factor(area)), colour="dodgerblue1", size=4) + 
									geom_errorbar(data=stat_all, aes(x=surveyF, ymin=IC_low, ymax=IC_high, group=factor(area)), colour="dodgerblue1",  width=0.1, linewidth=1) +
									ylab("Abundance index") + xlab("Survey") +
									scale_x_discrete(labels=c("1998-2000", "2009-2010", "2017-2019")) + 
									facet_wrap(~area, ncol=3, labeller=labeller(area=area.labs)) + 
									theme(plot.margin = unit(c(5.5, 5.5, 5.5, 35), units="pt"),
											axis.title.y = element_text(vjust=8))




# Combine estimated trends over the island or by regions
df_all_trends<-rbind( plot_dist_trends, df_dist_reg_trends)
df_all_trends$area<-factor(df_all_trends$area, levels=c("Island", "West", "East"))

# A data. frame with the values for the mean IC95% 
tab_res_all<-rbind(tab_res_trends, tab_reg_trends) %>%
				mutate_at(c("mean_trend","CI_low", "CI_high"), as.numeric) 

tab_res_all$tr_char<-paste0(round(tab_res_all$mean_trend,1),"% \n[", round(tab_res_all$CI_low,1),"%;",round(tab_res_all$CI_high,1),"%]")
tab_res_all$area<-factor(tab_res_all$area, levels=c("Island", "West", "East"))

#A plot of the trends by region(Island, West, East)
plot_all_trends<-ggplot() + 
						geom_density_ridges(data=df_all_trends, aes(x=Trends, y=Periods), fill="dodgerblue1", colour="dodgerblue1", alpha=0.7, bandwidth=5) + 
									xlab("% of change in abundance") +  ylab("") + 
								geom_vline(xintercept = 0, linetype="dotted", color = "black", linewidth=1) +
									facet_wrap(~area,  ncol=3, labeller=labeller(area=area.labs))	+
										theme(axis.text.y = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
											geom_text(data = tab_res_all, aes(x=50, y=as.numeric(as.factor(Periods))+0.5), label=tab_res_all$tr_char, size=3)


#ggplot(data=df_all_trends, aes(x=Trends)) + geom_histogram() + facet_wrap(~area)


abund_trends_all<-ggarrange(plot_all_abund, plot_all_trends, nrow=2, labels=c("a", "b"))	

ggsave(paste0(dirout, "Abundance_trends_region.jpeg"), plot=abund_trends_all,  units="cm",  width=24, height=17, dpi="retina")



df_plot_pres<-df_pts_full %>% 
		group_by(surveyF, id_newF) %>%
				summarise(sum_obs=sum(ptot),X=X[1], Y=Y[1], max_obs=max(ptot)) %>%
					mutate(pres_sp=ifelse(sum_obs==0,0,1)) %>%
					mutate(mean_obsCat=cut(max_obs, c(-1,0,1,2,3,4), labels=c("0", "1", "2", "3", "4"))) %>%
						ungroup()

plot_sp_abraw<-ggplot(df_plot_pres) +  #geom_sf(data=lim_reu, colour="black", fill="white") + 
				geom_contour(data=mnt_reu_df, aes(x=X, y=Y, z=value), binwidth=100, colour="black", alpha=0.2, linewidth=0.2) + 
				geom_sf(data=east_west, colour="black", fill="transparent",linewidth=0.5, alpha=0) + 
				#scale_colour_manual(name="Region", values=c("red3","royalblue1"), labels=c("East", "West")) +
				geom_point(aes(x=X, y=Y, size=factor(mean_obsCat), colour=factor(mean_obsCat)))+
				scale_size_manual(name="Maximum n° of breeding pair / survey", breaks=0:4, values=c(0.25, 0.75, 1.25,  1.75, 2.25))+ 
				scale_colour_manual(name="Maximum n° of breeding pair / survey", breaks=0:4, values=c("blue", rep("red", 4))) +
				facet_grid(~surveyF, labeller=labeller(surveyF=new_labels2)) +  theme_sp +
				theme(legend.position="bottom", 
				legend.text=element_text(size=8),
				legend.title=element_text(size=9),
						legend.title.align = 0.5,
						legend.background=element_rect(color="transparent")
					 #legend.box.just = "center",
					 #legend.box = "vertical"
					 )#, label.vjust = 0.5

ggsave(paste0(dirout, "Map_max_abund_per_survey.jpeg"), plot=plot_sp_abraw , units="cm",  width=24, height=10, dpi="retina" )



########################################################################################################
####################################### Comparing trends ###############################################
##################################### for data collected at  ###########################################
####################################### different periods ##############################################

# Else read the csv file and convert columns to appropriate format
df_pts_fullR<-read.table("./Data/df_pts_full.csv", sep=";", h=T)

df_pts_full<-df_pts_fullR
df_pts_full$surveyF<-as.factor(df_pts_full$surveyF)
df_pts_full$id_newF<-as.factor(df_pts_full$id_newF)
df_pts_full$yearF<-as.factor(df_pts_full$yearF)


# A data set  for the points common to all surveys

id_visited_all <- df_pts_fullR %>%
  distinct(id_newF, surveyF) %>%
  count(id_newF, name = "n_survey") %>%
  filter(n_survey == n_distinct(df_pts_fullR$surveyF)) %>%
  pull(id_newF)


#convert to factors
df_pts_comm<-subset(df_pts_fullR, df_pts_fullR$id_newF %in% id_visited_all)
df_pts_comm$surveyF<-as.factor(df_pts_comm$surveyF)
df_pts_comm$id_newF<-as.factor(df_pts_comm$id_newF)
df_pts_comm$yearF<-as.factor(df_pts_comm$yearF)


# A data set for all the points surveyed in Apr-May-June
df_pts_fall<-subset(df_pts_fullR, df_pts_fullR$jd > 100 & df_pts_fullR$jd < 200)
df_pts_fall$surveyF<-as.factor(df_pts_fall$surveyF)
df_pts_fall$id_newF<-as.factor(df_pts_fall$id_newF)
df_pts_fall$yearF<-as.factor(df_pts_fall$yearF)

#A data set just for the points common to all surveys, surveyed in Apr-May-June
df_pts_comm_fall<-subset(df_pts_fullR, df_pts_fullR$id_newF %in% id_visited_all & (df_pts_fullR$jd > 100 & df_pts_fullR$jd < 200))
df_pts_comm_fall$surveyF<-as.factor(df_pts_comm_fall$surveyF)
df_pts_comm_fall$id_newF<-as.factor(df_pts_comm_fall$id_newF)
df_pts_comm_fall$yearF<-as.factor(df_pts_comm_fall$yearF)

# We test the same model structure
datasets<-c("df_pts_full", "df_pts_comm", "df_pts_fall", "df_pts_comm_fall")

df_pred<-data.frame(surveyF=levels(df_pts_full$surveyF))
df_pred$log_length<-mean(df_pts_full$log_length)
df_pred$jd_s<-0
df_pred$id_newF<-"550"#levels(df_pts_full$id_newF)[1]
df_pred$yearF<-"2023"


nreplicates<-1000

df_res<-NULL


set.seed(1234)

for (i in 1:length(datasets)){
		
		mtemp<-gam(ptot~  s(yearF, bs="re") + s(id_newF, bs="re") + s(log_length, k=5) + s(jd_s, bs="cc", k=5) +  surveyF, family=tw, data=get(datasets[i]), method="REML")# 
		distri_temp<-compute_distri(model=mtemp, time_var="surveyF", dfpred=df_pred, nreplicates=nreplicates)
		tr_1_3<-compute_trends(distri_temp, ystart="1", yend="3")
		tr_1_2<-compute_trends(distri_temp, ystart="1", yend="2")
		tr_2_3<-compute_trends(distri_temp, ystart="2", yend="3")
		
		#Combine trends in 1 object
		trAllT<-rbind(tr_1_2, tr_1_3, tr_2_3)
		
		df_resT<-data.frame(data_set=datasets[i], diff_surv=c("1_3", "1_2", "2_3"), mean_trend=NA, IC_low=NA, IC_high=NA)

		df_resT[1:3, c("mean_trend", "IC_low", "IC_high")]<-trAllT
		df_res<-rbind(df_res, df_resT)
}	

df_res$Periods<-rep(c("1998 - 2009","1998 - 2017", "2009 - 2017" ), length(datasets))
df_res$Data_sets<-rep(c("All points", "Points common to the 3 surveys", "Points surveyed in Apr-June", "Points common to the 3 surveys (Apr-June)"), each=3)

plot_comp_tr<-ggplot(df_res) + geom_point(aes(x=Data_sets, y=mean_trend, colour=Data_sets), size=4) + geom_errorbar(aes(x=Data_sets, ymin=IC_low, ymax=IC_high, colour=Data_sets), width = 0.2, linewidth=1) + facet_wrap(~Periods) + 
			geom_hline(yintercept=0, linetype="dashed", colour="black") + ylab("Population trend (in %)")+
			scale_colour_discrete(name="Data set") + 
			 theme(axis.title.x=element_blank(),
					axis.text.x=element_blank(),
						axis.ticks.x=element_blank(),
						legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
						 guides(color=guide_legend(nrow=2, byrow=TRUE)) 

ggsave(paste0(dirout, "/Comp_trends_data_sets.jpeg"), plot=plot_comp_tr,  units="cm",  width=17, height=10, dpi=300)	
