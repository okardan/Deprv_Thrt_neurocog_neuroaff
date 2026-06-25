# Scripts for "Cognitive and affective neurodevelopment in youth exposed to deprivation and threat"
# Contact Omid Kardan omidk@med.umich.edu
# This script is for compiling the variables used in the neuroaffective analyses

library(dplyr)
library(ggplot2)
library(tidyverse)
library(ppcor)
library(fastDummies)
library(sjPlot)
library(corrplot)
library(lme4)
library(lmerTest)
library(plyr)
library(mediation)
library(effectsize)

setwd('study directory')

####################            Run the Kardan_dprv_thrt_variables_compilation_cog R script first                  ##########################
####################            ###########################
###################                                 ##########################
df_filled_demog_withcog_Y0 <- read.csv('~/df_filled_demog_withcog_Y0_more.csv')
df_filled_demog_withcog_Y0$subid <- df_filled_demog_withcog_Y0$subid.x
df_filled_demog_withcog_Y0$picvocab <- df_filled_demog_withcog_Y0$nihtbx_picvocab_uncorrected
df_filled_demog_withcog_Y0$read <- df_filled_demog_withcog_Y0$nihtbx_reading_uncorrected
df_filled_demog_withcog_Y0$picture <- df_filled_demog_withcog_Y0$nihtbx_picture_uncorrected
df_filled_demog_withcog_Y0$flanker <- df_filled_demog_withcog_Y0$nihtbx_flanker_uncorrected
df_filled_demog_withcog_Y0$pattern <- df_filled_demog_withcog_Y0$nihtbx_pattern_uncorrected

df_filled_demog_withcog_Y2 <- read.csv('~/df_filled_demog_withcog_Y2_more.csv')
df_filled_demog_withcog_Y2$subid <- df_filled_demog_withcog_Y2$subid.x
df_filled_demog_withcog_Y2$rel_family_id <- df_filled_demog_withcog_Y0$rel_family_id
df_filled_demog_withcog_Y2$Achieva <- df_filled_demog_withcog_Y0$Achieva
df_filled_demog_withcog_Y2$Discovery <- df_filled_demog_withcog_Y0$Discovery
df_filled_demog_withcog_Y2$Ingenia <- df_filled_demog_withcog_Y0$Ingenia
df_filled_demog_withcog_Y2$Orchestra <- df_filled_demog_withcog_Y0$Orchestra
df_filled_demog_withcog_Y2$Prisma <- df_filled_demog_withcog_Y0$Prisma
df_filled_demog_withcog_Y2$Pfit <- df_filled_demog_withcog_Y0$Pfit
df_filled_demog_withcog_Y2$Premier <- df_filled_demog_withcog_Y0$Premier
df_filled_demog_withcog_Y2$UHP <- df_filled_demog_withcog_Y0$UHP
df_filled_demog_withcog_Y2$picvocab <- df_filled_demog_withcog_Y2$nihtbx_picvocab_uncorrected
df_filled_demog_withcog_Y2$read <- df_filled_demog_withcog_Y2$nihtbx_reading_uncorrected
df_filled_demog_withcog_Y2$picture <- df_filled_demog_withcog_Y2$nihtbx_picture_uncorrected
df_filled_demog_withcog_Y2$flanker <- df_filled_demog_withcog_Y2$nihtbx_flanker_uncorrected
df_filled_demog_withcog_Y2$pattern <- df_filled_demog_withcog_Y2$nihtbx_pattern_uncorrected

# frame displacement and quality-based exclusions for the MID and EN-back task fMRI
wbdiff_mot <- read.csv('~/mri_y_qc_motion.csv')
wbdiff_a<- read.csv('~/mri_y_qc_incl.csv')
wbdiff_b <- merge(wbdiff_mot[wbdiff_mot$eventname == 'baseline_year_1_arm_1',c('subkey','tfmri_mid_all_meanmotion','tfmri_nback_all_meanmotion')],
                  wbdiff_a[wbdiff_a$eventname == 'baseline_year_1_arm_1',c('subkey','imgincl_mid_include','imgincl_nback_include')],
                  by.x = 'subkey', by.y = 'subkey', all.y = TRUE)
wbdiff_all_b <- merge(wbdiff_b,
                      df_filled_demog_withcog_Y0[,c('subid.x','interview_age','Income','HighestEd',
                                                             'Male_bin','White','Black','Hispanic','Asian','Other',
                                                             'Achieva','Discovery','Ingenia','Orchestra','Prisma','Pfit','Premier','UHP',
                                                             'site_id_l','rel_family_id','anthroheightcalc','pds_ss',
                                                             'picvocab','read','picture','flanker','pattern',
                                                             'nback0_acc','nback2_acc','PF10_INT_lavaan','PF10_EXT_lavaan','PF10_lavaan','reshist_addr1_popdensity',
                                                             'nih_abcd_cpm','tfmri_meanFD','tfmri_numruns',
                                                             'FamilyConflict', 'famhx_ss_momdad_dg_p' , 'famhx_ss_momdad_alc_p',
                                                             'reshist_addr1_adi_wsum','NeighCrime','NeighSafety',
                                                             'ace_y0','ace_y1y2','no_to_yes_ace','yes_to_more_ace','cumul_ace')],
                      by.x = 'subkey', by.y = 'subid.x', all.x = TRUE)  # 

wbdiff_all_b <- wbdiff_all_b[wbdiff_all_b$imgincl_mid_include==1 & wbdiff_all_b$imgincl_nback_include==1,] 

wbdiff_y2 <- merge(wbdiff_mot[wbdiff_mot$eventname == '2_year_follow_up_y_arm_1',c('subkey','tfmri_mid_all_meanmotion','tfmri_nback_all_meanmotion')],
                   wbdiff_a[wbdiff_a$eventname == '2_year_follow_up_y_arm_1',c('subkey','imgincl_mid_include','imgincl_nback_include')],
                   by.x = 'subkey', by.y = 'subkey', all.y = TRUE)
wbdiff_all_y2 <- merge(wbdiff_y2,df_filled_demog_withcog_Y2[,c('subid.x','interview_age','Income','HighestEd',
                                                               'Male_bin','White','Black','Hispanic','Asian','Other',
                                                               'Achieva','Discovery','Ingenia','Orchestra','Prisma','Pfit','Premier','UHP',
                                                               'site_id_l','rel_family_id','anthroheightcalc','pds_ss',
                                                               'picvocab','read','picture','flanker','pattern',
                                                               'nback0_acc','nback2_acc','PF10_INT_lavaan','PF10_EXT_lavaan','PF10_lavaan','reshist_addr1_popdensity',
                                                               'nih_abcd_cpm','tfmri_meanFD','tfmri_numruns',
                                                               'FamilyConflict', 'famhx_ss_momdad_dg_p' , 'famhx_ss_momdad_alc_p',
                                                               'reshist_addr1_adi_wsum','NeighCrime','NeighSafety',
                                                               'ace_y0','ace_y1y2','no_to_yes_ace','yes_to_more_ace','cumul_ace')],
                      by.x = 'subkey', by.y = 'subid.x', all.x = TRUE)  # 

wbdiff_all_y2 <- wbdiff_all_y2[wbdiff_all_y2$imgincl_mid_include==1 & wbdiff_all_y2$imgincl_nback_include==1,] 


dat10 <- wbdiff_all_b %>% mutate(age = interview_age/12, age_q = (scale(interview_age))^2)
dat20 <- wbdiff_all_y2 %>% mutate( age = interview_age/12, age_q = (scale(interview_age))^2)

# NAc and Caudate reward anticipation
rwd_subc <- read.csv('~/mri_y_tfmr_mid_alrvn_aseg.csv')
rwd_subc_0 <- rwd_subc[rwd_subc$eventname == 'baseline_year_1_arm_1',]
rwd_subc_2 <- rwd_subc[rwd_subc$eventname == '2_year_follow_up_y_arm_1',]

dat100 <- merge(dat10,rwd_subc_0[,c('subkey','tfmri_ma_alrvn_b_scs_aalh','tfmri_ma_alrvn_b_scs_aarh','tfmri_ma_alrvn_b_scs_cdlh','tfmri_ma_alrvn_b_scs_cdrh')],
                by.x = 'subkey', by.y = 'subkey', all.x = TRUE)  # 
dat200 <- merge(dat20,rwd_subc_2[,c('subkey','tfmri_ma_alrvn_b_scs_aalh','tfmri_ma_alrvn_b_scs_aarh','tfmri_ma_alrvn_b_scs_cdlh','tfmri_ma_alrvn_b_scs_cdrh')],
                by.x = 'subkey', by.y = 'subkey', all.x = TRUE)  # 

# Amygdala and Insula negative emotional face

aff_subc <- read.csv('~/mri_y_tfmr_nback_ngfvntf_aseg.csv')
aff_cort <- read.csv('~/mri_y_tfmr_nback_ngfvntf_dsk.csv')
aff_subc_0 <- aff_subc[aff_subc$eventname == 'baseline_year_1_arm_1',] 
aff_subc_2 <- aff_subc[aff_subc$eventname == '2_year_follow_up_y_arm_1',]
aff_cort_0 <- aff_cort[aff_cort$eventname == 'baseline_year_1_arm_1',] 
aff_cort_2 <- aff_cort[aff_cort$eventname == '2_year_follow_up_y_arm_1',]
# names in Release 6 are different: tfmri_nback_all_224 and 238 are left and right amygdala (mr_y__nback__ngfvntf__aseg__ag__rh_beta)
# names in Release 6 are different: tfmri_nback_all_780 and 814 are left and right insula (mr_y__nback__ngfvntf__dsk__ins__rh_beta)
dat101 <- merge(dat100,aff_subc_0[,c('subkey','tfmri_nback_all_224','tfmri_nback_all_238')],
                by.x = 'subkey', by.y = 'subkey', all.x = TRUE)  # 
dat202 <- merge(dat200,aff_subc_2[,c('subkey','tfmri_nback_all_224','tfmri_nback_all_238')],
                by.x = 'subkey', by.y = 'subkey', all.x = TRUE)  # 
dat1000 <- merge(dat101,aff_cort_0[,c('subkey','tfmri_nback_all_780','tfmri_nback_all_814')],
                 by.x = 'subkey', by.y = 'subkey', all.x = TRUE)  # 
dat2000 <- merge(dat202,aff_cort_2[,c('subkey','tfmri_nback_all_780','tfmri_nback_all_814')],
                 by.x = 'subkey', by.y = 'subkey', all.x = TRUE)  # 

dat1000 <- dat1000 %>% mutate(ng_nt_face_amyg = (tfmri_nback_all_224 + tfmri_nback_all_238)/2, ng_nt_face_insu = (tfmri_nback_all_780 + tfmri_nback_all_814)/2, mid_NAc = (tfmri_ma_alrvn_b_scs_aalh + tfmri_ma_alrvn_b_scs_aarh)/2, mid_Caud = (tfmri_ma_alrvn_b_scs_cdlh + tfmri_ma_alrvn_b_scs_cdrh)/2)
dat2000 <- dat2000 %>% mutate(ng_nt_face_amyg = (tfmri_nback_all_224 + tfmri_nback_all_238)/2, ng_nt_face_insu = (tfmri_nback_all_780 + tfmri_nback_all_814)/2, mid_NAc = (tfmri_ma_alrvn_b_scs_aalh + tfmri_ma_alrvn_b_scs_aarh)/2, mid_Caud = (tfmri_ma_alrvn_b_scs_cdlh + tfmri_ma_alrvn_b_scs_cdrh)/2)

##########################################
write.csv(dat1000,'~/df_filled_beh_Y0_more_neuroaff.csv')
write.csv(dat2000,'~/df_filled_beh_Y2_more_neuroaff.csv')