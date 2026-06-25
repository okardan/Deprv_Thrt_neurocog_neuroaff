% # Scripts for "Cognitive and affective neurodevelopment in youth exposed to deprivation and threat"
% # Contact Omid Kardan omidk@med.umich.edu
%%% This script is for running PLS based on neuroaffective measures and
%%% also makes the Figure 1 (at the end of the script)
% Requirements: 
% 1) Uses csv file generated from the R script Kardan_dprv_thrt_variables_compilation_cog.R
% % 2) Uses the PLS scripts from https://www.rotman-baycrest.on.ca/index.php?section=345
% Download plscmd and plsgui and place in Pls_folder
% 3) Uses export_fig from https://github.com/altmany/export_fig
%%  ######################################### Deprivation only #####################################################################################
%%  %% Create matched control group for the depriv_only participants
clear all
T_0 = read_abcd_beh_Aff('~/df_filled_beh_Y0_more_neuroaff.csv');
T_2 =read_abcd_beh_Aff('~/df_filled_beh_Y2_more_neuroaff.csv');
s_table = readtable("~/df_filled_demog_withcog_Y0_more.csv");
aff_long = innerjoin(T_0,T_2,"Keys","subkey");
with_2b_aff = find(ismember(string(s_table.subid_x), string(aff_long.subkey)));


usable_subs = aff_long(~isnan(aff_long.mid_NAc_T_2 + aff_long.ng_nt_face_insu_T_2 + ...
    aff_long.mid_NAc_T_0 + aff_long.ng_nt_face_insu_T_0),:); % save as usable_subs_aff.mat

% this section runs iteratively to make the matched group which is then
% saved as num_matches_subs_deprv_only_aff.mat
depv_high1 = find(usable_subs.cumul_ace_T_2 == 1      ); % 1 is depriv only

%load(['~/num_matches_subs_deprv_only_aff.mat']);%uncomment this after the first
% time running the script which generates the list


depv_high = depv_high1;
% depv_high =depv_high1(num_matches_subs_depv_only_aff >1);  %uncomment this after the first
% time running the script which generates the list

depv_high_n = length(depv_high);
depv_low = find(usable_subs.cumul_ace_T_0 == 0 ); % & ~isnan(usable_subs.logM_T_0+usable_subs.logM_T_2)

usable_subs.raceth = usable_subs.White_T_0;
usable_subs.raceth(usable_subs.Black_T_0 == 1) = 2;
usable_subs.raceth(usable_subs.Hispanic_T_0 == 1) = 3;
usable_subs.raceth(usable_subs.Asian_T_0 == 1) = 4;
usable_subs.raceth(usable_subs.Other_T_0 == 1) = 5;
usable_subs.interview_age(isnan(usable_subs.interview_age_T_0)) = nanmean(usable_subs.interview_age_T_0);
usable_subs.agey = round(usable_subs.interview_age_T_0/6);
usable_subs.HighestEd_T_0(isnan(usable_subs.HighestEd_T_0)) = nanmean(usable_subs.HighestEd_T_0);
usable_subs.edu = round(usable_subs.HighestEd_T_0/5);

income_matches = cell(depv_high_n,1);
age_matches = cell(depv_high_n,1);
sex_matches = cell(depv_high_n,1);
raceth_matches = cell(depv_high_n,1);
edu_matches = cell(depv_high_n,1);

for k =1:depv_high_n
    
    if ~isnan(usable_subs.Income_T_0(depv_high(k)) )
        a = find(round(usable_subs.Income_T_0) == round(usable_subs.Income_T_0(depv_high(k))));
    else a = depv_low;
    end
        income_matches{k} = intersect(depv_low,a);
        
    if ~isnan(usable_subs.agey(depv_high(k)))
        a = find(usable_subs.agey == usable_subs.agey(depv_high(k)));
    else a = depv_low;
    end
        age_matches{k} = intersect(depv_low,a);
    
    if ~isnan(usable_subs.Male_bin_T_0(depv_high(k)))
        a = find(usable_subs.Male_bin_T_0 == usable_subs.Male_bin_T_0(depv_high(k)));
    else a = depv_low;
    end
         sex_matches{k} = intersect(depv_low,a);
         
    if ~isnan(usable_subs.raceth(depv_high(k)))
        a = find(usable_subs.raceth == usable_subs.raceth(depv_high(k)));
    else a = depv_low;
    end
         raceth_matches{k} = intersect(depv_low,a);
    
    if ~isnan(usable_subs.edu(depv_high(k)))
        a = find(usable_subs.edu == usable_subs.edu(depv_high(k)));
    else a = depv_low;
    end
         edu_matches{k} = intersect(depv_low,a);
            
end
cont_mj = cell(depv_high_n,1); cont_mj_full = cell(depv_high_n,1); cont_mj_bare = cell(depv_high_n,1);
for k = 1:depv_high_n
    gg = intersect( intersect(sex_matches{k},income_matches{k}), intersect(raceth_matches{k},edu_matches{k}));
    cont_mj_full{k} = setdiff(intersect( gg, age_matches{k}  ),  depv_high);
    cont_mj{k} = setdiff( gg,  depv_high);
    gg2 = intersect(income_matches{k}, raceth_matches{k});
    cont_mj_bare{k} = setdiff( intersect(gg2,edu_matches{k}),  depv_high);

end

num_matches_subs_depv_only_aff =[];
for l=1:764
num_matches_subs_depv_only_aff = [num_matches_subs_depv_only_aff; length(cont_mj_bare{l})];
end
% save list as num_matches_subs_deprv_only_aff.mat

HCss = []; llls =[];
for ww = 1:10000
    HCs =NaN(depv_high_n,1);
    for k=1:depv_high_n
        if ~isempty(cont_mj_full{k})
            if length(cont_mj_full{k}) >1
                for jjk = 1: length(cont_mj_full{k})
                    dd = randperm(length(cont_mj_full{k}),1);
                    if ~ismember(cont_mj_full{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj_full{k}(dd);
            
            
        elseif ~isempty(cont_mj{k})
            if length(cont_mj{k}) >1
                for jjk = 1: length(cont_mj{k})
                    dd = randperm(length(cont_mj{k}),1);
                    if ~ismember(cont_mj{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj{k}(dd);
            
            
        elseif ~isempty(cont_mj_bare{k})
            if length(cont_mj_bare{k}) >1
                for jjk = 1: length(cont_mj_bare{k})
                    dd = randperm(length(cont_mj_bare{k}),1);
                    if ~ismember(cont_mj_bare{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj_bare{k}(dd);
        end
    end
    llls = [llls; length(unique(HCs))];
        
        HCss = [HCss  HCs];

end

 matches = HCss(:,1);
for j=1:10000
    [v, w] = unique( matches, 'stable' );
    duplicate_indices = setdiff( 1:numel(matches), w );
    if length(duplicate_indices)==0
        break
    end
    matches(duplicate_indices) = HCss(duplicate_indices,j);
end
        
matches_603 = unique(matches); 
% check if cohen's d is samll for group differences in each category

T_ = usable_subs;
cohd = (nanmean(T_.age_T_0(depv_high))-nanmean(T_.age_T_0(setdiff(1:4344,depv_high))) )/ sqrt(nanvar(T_.age_T_0(depv_high))+nanvar(T_.age_T_0(setdiff(1:4344,depv_high))) )
cohd = (nanmean(T_.age_T_0(depv_high))-nanmean(T_.age_T_0(depv_low)) )/ sqrt(nanvar(T_.age_T_0(depv_high))+nanvar(T_.age_T_0(depv_low)) )
cohd = (nanmean(T_.age_T_0(depv_high))-nanmean(T_.age_T_0(matches_603)) )/ sqrt(nanvar(T_.age_T_0(depv_high))+nanvar(T_.age_T_0(matches_603)) )

cohd = (nanmean(T_.Income_T_0(depv_high))-nanmean(T_.Income_T_0(setdiff(1:4344,depv_high))) )/ sqrt(nanvar(T_.Income_T_0(depv_high))+nanvar(T_.Income_T_0(setdiff(1:4344,depv_high))) )
cohd = (nanmean(T_.Income_T_0(depv_high))-nanmean(T_.Income_T_0(depv_low)) )/ sqrt(nanvar(T_.Income_T_0(depv_high))+nanvar(T_.Income_T_0(depv_low)) )
cohd = (nanmean(T_.Income_T_0(depv_high))-nanmean(T_.Income_T_0(matches_603)) )/ sqrt(nanvar(T_.Income_T_0(depv_high))+nanvar(T_.Income_T_0(matches_603)) )


cohd = (nanmean(T_.HighestEd_T_0(depv_high))-nanmean(T_.HighestEd_T_0(setdiff(1:4344,depv_high))) )/ sqrt(nanvar(T_.HighestEd_T_0(depv_high))+nanvar(T_.HighestEd_T_0(setdiff(1:4344,depv_high))) )
cohd = (nanmean(T_.HighestEd_T_0(depv_high))-nanmean(T_.HighestEd_T_0(depv_low)) )/ sqrt(nanvar(T_.HighestEd_T_0(depv_high))+nanvar(T_.HighestEd_T_0(depv_low)) )
cohd = (nanmean(T_.HighestEd_T_0(depv_high))-nanmean(T_.HighestEd_T_0(matches_603)) )/ sqrt(nanvar(T_.HighestEd_T_0(depv_high))+nanvar(T_.HighestEd_T_0(matches_603)) )

cohd = (nanmean(T_.Black_T_0(depv_high))-nanmean(T_.Black_T_0(setdiff(1:4344,depv_high))) )/ sqrt(nanvar(T_.Black_T_0(depv_high))+nanvar(T_.Black_T_0(setdiff(1:4344,depv_high))) )
cohd = (nanmean(T_.Black_T_0(depv_high))-nanmean(T_.Black_T_0(depv_low)) )/ sqrt(nanvar(T_.Black_T_0(depv_high))+nanvar(T_.Black_T_0(depv_low)) )
cohd = (nanmean(T_.Black_T_0(depv_high))-nanmean(T_.Black_T_0(matches_603)) )/ sqrt(nanvar(T_.Black_T_0(depv_high))+nanvar(T_.Black_T_0(matches_603)) )

cohd = (nanmean(T_.White_T_0(depv_high))-nanmean(T_.White_T_0(setdiff(1:4344,depv_high))) )/ sqrt(nanvar(T_.White_T_0(depv_high))+nanvar(T_.White_T_0(setdiff(1:4344,depv_high))) )
cohd = (nanmean(T_.White_T_0(depv_high))-nanmean(T_.White_T_0(depv_low)) )/ sqrt(nanvar(T_.White_T_0(depv_high))+nanvar(T_.White_T_0(depv_low)) )
cohd = (nanmean(T_.White_T_0(depv_high))-nanmean(T_.White_T_0(matches_603)) )/ sqrt(nanvar(T_.White_T_0(depv_high))+nanvar(T_.White_T_0(matches_603)) )

cohd = (nanmean(T_.Hispanic_T_0(depv_high))-nanmean(T_.Hispanic_T_0(setdiff(1:4344,depv_high))) )/ sqrt(nanvar(T_.Hispanic_T_0(depv_high))+nanvar(T_.Hispanic_T_0(setdiff(1:4344,depv_high))) )
cohd = (nanmean(T_.Hispanic_T_0(depv_high))-nanmean(T_.Hispanic_T_0(depv_low)) )/ sqrt(nanvar(T_.Hispanic_T_0(depv_high))+nanvar(T_.Hispanic_T_0(depv_low)) )
cohd = (nanmean(T_.Hispanic_T_0(depv_high))-nanmean(T_.Hispanic_T_0(matches_603)) )/ sqrt(nanvar(T_.Hispanic_T_0(depv_high))+nanvar(T_.Hispanic_T_0(matches_603)) )


matches_depv_603_aff = unique(matches); % save these as your matched controls
depv_721_aff = depv_high; % save these as your exposure group


%% %%%%% PLS
%% neuroaff latent variable

clear all
T_0 = read_abcd_beh_Aff('~/df_filled_beh_Y0_more_neuroaff.csv');
T_2 =read_abcd_beh_Aff('~/df_filled_beh_Y2_more_neuroaff.csv');


load('~/usable_subs_aff.mat')
T_ = usable_subs;

load('matches_depv_603_aff.mat')
load('depv_721_aff.mat')

addpath(genpath('~/Pls_folder/Pls'));
rng('default');
groups{1} = [matches_depv_603_aff];
groups{2} = [depv_721_aff];
nTimes = 2;

datamat_lst = cell(length(groups),1); 
% ####*****************####

for g = 1:numel(groups)
    for c = 1:nTimes
        for s = 1:numel(groups{g})
            sub_id = char(T_.subkey(groups{g}(s)));
            
            vec2 =[];
            if c==1
                locb = find(T_0.subkey == sub_id);
                if length(locb)==0
                    sub_id
                end
                vec2 = [T_0.mid_NAc(locb)   T_0.mid_Caud(locb)   T_0.ng_nt_face_amyg(locb) T_0.ng_nt_face_insu(locb)]; % T_0.logM(locb) T_0.logK(locb)
                                                    
            end
            
            if c==2
                loct = find(T_2.subkey == sub_id);
                if length(loct)==0
                    sub_id
                end
               vec2 = [T_2.mid_NAc(loct)   T_2.mid_Caud(loct)   T_2.ng_nt_face_amyg(loct) T_2.ng_nt_face_insu(loct)]; % T_2.logM(loct) T_2.logK(loct)
                  

            end
            
            datamat_lst{g} = [datamat_lst{g}; vec2];
            s;
            
        end
    end
end

num_subj = [length(groups{1})  length(groups{2})];
num_cond = nTimes;
option.method = 1; 
option.num_boot = 2000;
option.num_perm = 2000;
option.meancentering_type=[2]; 
result = pls_analysis(datamat_lst, num_subj, num_cond, option);
result.datamat_lst = datamat_lst;
% save as PLS_result_603low_721hdepv_y0_y2_aff.mat
%% plotting the results (figure 2 panel B and table 1 bottom half)
load('PLS_result_603low_721hdepv_y0_y2_aff.mat');  
ps = result.perm_result.sprob
load('~/usable_subs_aff.mat')
T_ = usable_subs;

load('matches_depv_603_aff.mat')
load('depv_721_aff.mat')
matches = matches_depv_603_aff;
depvhigh = depv_721_aff;

% *** table 1 bottom half
depv_low = find(usable_subs.cumul_ace_T_2 == 0 );
T_.AsiOth = T_.Asian_T_0 + T_.Other_T_0;

[nanmean(T_.HighestEd_T_0(depvhigh)) sqrt(nanvar(T_.HighestEd_T_0(depvhigh)))]
[nanmean(T_.HighestEd_T_0(depv_low)) sqrt(nanvar(T_.HighestEd_T_0(depv_low)))]
[nanmean(T_.HighestEd_T_0(matches)) sqrt(nanvar(T_.HighestEd_T_0(matches)))]
cohd = (nanmean(T_.HighestEd_T_0(depvhigh))-nanmean(T_.HighestEd_T_0(depv_low)) )/ sqrt(nanvar(T_.HighestEd_T_0(depvhigh))+nanvar(T_.HighestEd_T_0(depv_low)) )
cohd = (nanmean(T_.HighestEd_T_0(depvhigh))-nanmean(T_.HighestEd_T_0(matches)) )/ sqrt(nanvar(T_.HighestEd_T_0(depvhigh))+nanvar(T_.HighestEd_T_0(matches)) )

      % levels= 1:10, labels = c("5000", "8500", "14000", "20500", "30000","42500", 7:"62500", 8:"87500", "150000", "200000") )
nanmean(T_.Income_T_0(depvhigh)); 
42500 + (6.7895-6)*( 62500 - 42500)
nanmean(T_.Income_T_0(depv_low)) ;
87500 + (8.308-8)*(150000 - 87500)
nanmean(T_.Income_T_0(matches)); 
62500 + (7.1663-7)*(87500 - 62500)
cohd = (nanmean(T_.Income_T_0(depvhigh))-nanmean(T_.Income_T_0(depv_low)) )/ sqrt(nanvar(T_.Income_T_0(depvhigh))+nanvar(T_.Income_T_0(depv_low)) )
cohd = (nanmean(T_.Income_T_0(depvhigh))-nanmean(T_.Income_T_0(matches)) )/ sqrt(nanvar(T_.Income_T_0(depvhigh))+nanvar(T_.Income_T_0(matches)) )

[length(find(T_.White_T_0(depvhigh))) length(find(T_.White_T_0(depvhigh)))/length(depvhigh)]
[length(find(T_.White_T_0(depv_low))) length(find(T_.White_T_0(depv_low)))/length(depv_low)]
[length(find(T_.White_T_0(matches))) length(find(T_.White_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.White_T_0(depvhigh))-nanmean(T_.White_T_0(depv_low)) )/ sqrt(nanvar(T_.White_T_0(depvhigh))+nanvar(T_.White_T_0(depv_low)) )
cohd = (nanmean(T_.White_T_0(depvhigh))-nanmean(T_.White_T_0(matches)) )/ sqrt(nanvar(T_.White_T_0(depvhigh))+nanvar(T_.White_T_0(matches)) )


[length(find(T_.Hispanic_T_0(depvhigh))) length(find(T_.Hispanic_T_0(depvhigh)))/length(depvhigh)]
[length(find(T_.Hispanic_T_0(depv_low))) length(find(T_.Hispanic_T_0(depv_low)))/length(depv_low)]
[length(find(T_.Hispanic_T_0(matches))) length(find(T_.Hispanic_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.Hispanic_T_0(depvhigh))-nanmean(T_.Hispanic_T_0(depv_low)) )/ sqrt(nanvar(T_.Hispanic_T_0(depvhigh))+nanvar(T_.Hispanic_T_0(depv_low)) )
cohd = (nanmean(T_.Hispanic_T_0(depvhigh))-nanmean(T_.Hispanic_T_0(matches)) )/ sqrt(nanvar(T_.Hispanic_T_0(depvhigh))+nanvar(T_.Hispanic_T_0(matches)) )


[length(find(T_.Black_T_0(depvhigh))) length(find(T_.Black_T_0(depvhigh)))/length(depvhigh)]
[length(find(T_.Black_T_0(depv_low))) length(find(T_.Black_T_0(depv_low)))/length(depv_low)]
[length(find(T_.Black_T_0(matches))) length(find(T_.Black_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.Black_T_0(depvhigh))-nanmean(T_.Black_T_0(depv_low)) )/ sqrt(nanvar(T_.Black_T_0(depvhigh))+nanvar(T_.Black_T_0(depv_low)) )
cohd = (nanmean(T_.Black_T_0(depvhigh))-nanmean(T_.Black_T_0(matches)) )/ sqrt(nanvar(T_.Black_T_0(depvhigh))+nanvar(T_.Black_T_0(matches)) )


[length(find(T_.AsiOth(depvhigh))) length(find(T_.AsiOth(depvhigh)))/length(depvhigh)]
[length(find(T_.AsiOth(depv_low))) length(find(T_.AsiOth(depv_low)))/length(depv_low)]
[length(find(T_.AsiOth(matches))) length(find(T_.AsiOth(matches)))/length(matches)]
cohd = (nanmean(T_.AsiOth(depvhigh))-nanmean(T_.AsiOth(depv_low)) )/ sqrt(nanvar(T_.AsiOth(depvhigh))+nanvar(T_.AsiOth(depv_low)) )
cohd = (nanmean(T_.AsiOth(depvhigh))-nanmean(T_.AsiOth(matches)) )/ sqrt(nanvar(T_.AsiOth(depvhigh))+nanvar(T_.AsiOth(matches)) )

[length(find(~T_.Male_bin_T_0(depvhigh))) length(find(~T_.Male_bin_T_0(depvhigh)))/length(depvhigh)]
[length(find(~T_.Male_bin_T_0(depv_low))) length(find(~T_.Male_bin_T_0(depv_low)))/length(depv_low)]
[length(find(~T_.Male_bin_T_0(matches))) length(find(~T_.Male_bin_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.Male_bin_T_0(depvhigh))-nanmean(T_.Male_bin_T_0(depv_low)) )/ sqrt(nanvar(T_.Male_bin_T_0(depvhigh))+nanvar(T_.Male_bin_T_0(depv_low)) )
cohd = (nanmean(T_.Male_bin_T_0(depvhigh))-nanmean(T_.Male_bin_T_0(matches)) )/ sqrt(nanvar(T_.Male_bin_T_0(depvhigh))+nanvar(T_.Male_bin_T_0(matches)) )


[nanmean(T_.age_T_0(depvhigh)) sqrt(nanvar(T_.age_T_0(depvhigh)))]
[nanmean(T_.age_T_0(depv_low)) sqrt(nanvar(T_.age_T_0(depv_low)))]
[nanmean(T_.age_T_0(matches)) sqrt(nanvar(T_.age_T_0(matches)))]
cohd = (nanmean(T_.age_T_0(depvhigh))-nanmean(T_.age_T_0(depv_low)) )/ sqrt(nanvar(T_.age_T_0(depvhigh))+nanvar(T_.age_T_0(depv_low)) )
cohd = (nanmean(T_.age_T_0(depvhigh))-nanmean(T_.age_T_0(matches)) )/ sqrt(nanvar(T_.age_T_0(depvhigh))+nanvar(T_.age_T_0(matches)) )

[nanmean(T_.age_T_2(depvhigh)) sqrt(nanvar(T_.age_T_2(depvhigh)))]
[nanmean(T_.age_T_2(depv_low)) sqrt(nanvar(T_.age_T_2(depv_low)))]
[nanmean(T_.age_T_2(matches)) sqrt(nanvar(T_.age_T_2(matches)))]
cohd = (nanmean(T_.age_T_2(depvhigh))-nanmean(T_.age_T_2(depv_low)) )/ sqrt(nanvar(T_.age_T_2(depvhigh))+nanvar(T_.age_T_2(depv_low)) )
cohd = (nanmean(T_.age_T_2(depvhigh))-nanmean(T_.age_T_2(matches)) )/ sqrt(nanvar(T_.age_T_2(depvhigh))+nanvar(T_.age_T_2(matches)) )

% ***



Sites = zeros(3244,22);
for K=1:22
    Sites((T_.site_id_l_T_0 == K),K) = 1;
end

lv=1;   %~#

lg1 = length(matches); lg2 = length(depvhigh);
adjy0cont = zeros(2*(lg1+lg2),1);
adjy0cont(1:lg1) = 1;

adjy2cont = zeros(2*(lg1+lg2),1);
adjy2cont(lg1+1:2*lg1) = 1;

adjy0pm = zeros(2*(lg1+lg2),1);
adjy0pm(2*lg1+1:2*lg1+lg2) = 1;

adjy2pm = zeros(2*(lg1+lg2),1);
adjy2pm(2*lg1+lg2+1:2*(lg1+lg2)) = 1;


covariates_y0cont = [  T_.tfmri_mid_all_meanmotion_T_0(matches)   T_.tfmri_nback_all_meanmotion_T_0(matches)   Sites(matches,:)  T_.White_T_0(matches) T_.Black_T_0(matches) T_.Hispanic_T_0(matches) T_.HighestEd_T_0(matches) T_.Income_T_0(matches) T_.Male_bin_T_0(matches)  T_.interview_age_T_0(matches)  T_.Ingenia_T_0(matches) T_.Achieva_T_0(matches) T_.Discovery_T_0(matches)  T_.Pfit_T_0(matches) T_.reshist_addr1_adi_wsum_T_0(matches)  T_.pds_ss_T_0(matches)]; % T_.PF10_lavaan(matches)
covariates_y2cont = [ T_.tfmri_mid_all_meanmotion_T_2(matches)   T_.tfmri_nback_all_meanmotion_T_2(matches)   Sites(matches,:)  T_.White_T_0(matches) T_.Black_T_0(matches) T_.Hispanic_T_0(matches) T_.HighestEd_T_0(matches) T_.Income_T_0(matches) T_.Male_bin_T_0(matches)  T_.interview_age_T_2(matches)   T_.Ingenia_T_0(matches) T_.Achieva_T_0(matches) T_.Discovery_T_0(matches)  T_.Pfit_T_0(matches)  T_.reshist_addr1_adi_wsum_T_0(matches) T_.pds_ss_T_2(matches)];
covariates_y0pm = [   T_.tfmri_mid_all_meanmotion_T_0(depvhigh)   T_.tfmri_nback_all_meanmotion_T_0(depvhigh)   Sites(depvhigh,:)  T_.White_T_0(depvhigh) T_.Black_T_0(depvhigh) T_.Hispanic_T_0(depvhigh) T_.HighestEd_T_0(depvhigh) T_.Income_T_0(depvhigh) T_.Male_bin_T_0(depvhigh)  T_.interview_age_T_0(depvhigh)  T_.Ingenia_T_0(depvhigh) T_.Achieva_T_0(depvhigh) T_.Discovery_T_0(depvhigh)  T_.Pfit_T_0(depvhigh)  T_.reshist_addr1_adi_wsum_T_0(depvhigh)  T_.pds_ss_T_0(depvhigh)];
covariates_y2pm = [ T_.tfmri_mid_all_meanmotion_T_2(depvhigh)   T_.tfmri_nback_all_meanmotion_T_2(depvhigh)   Sites(depvhigh,:)  T_.White_T_0(depvhigh) T_.Black_T_0(depvhigh) T_.Hispanic_T_0(depvhigh) T_.HighestEd_T_0(depvhigh) T_.Income_T_0(depvhigh) T_.Male_bin_T_0(depvhigh)  T_.interview_age_T_2(depvhigh)   T_.Ingenia_T_0(depvhigh) T_.Achieva_T_0(depvhigh) T_.Discovery_T_0(depvhigh)  T_.Pfit_T_0(depvhigh)  T_.reshist_addr1_adi_wsum_T_0(depvhigh) T_.pds_ss_T_2(depvhigh)];


covariates = [covariates_y0cont; covariates_y2cont; covariates_y0pm; covariates_y2pm];
% 24,25,26, are race/ethnic
% exclude_race_as_cov = [1:23,27:37];
exclude_race_as_cov = [1:size(covariates,2)];
[a1 ] = partialcorr([result.usc(:,lv),adjy0cont],covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a2 ] = partialcorr([result.usc(:,lv),adjy2cont], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a3 ] = partialcorr([result.usc(:,lv),adjy0pm], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a4 ] = partialcorr([result.usc(:,lv),adjy2pm], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');

a1s =[]; a2s =[]; a3s =[]; a4s =[]; b1s=[];b2s=[];b3s=[]; b4s=[];
dats = [result.datamat_lst{1, 1}; result.datamat_lst{2, 1}];

nboot = 2000;
for jj=1:nboot
   rand_inds1 = randi([1,lg1],lg1,1); rand_inds11 = randi([1,lg1],lg1,1); 
   rand_inds2 = randi([1,lg2],lg2,1); rand_inds22 = randi([1,lg2],lg2,1);
   inds = [rand_inds1; lg1+ rand_inds11; 2*lg1+ rand_inds2; 2*lg1+lg2+ rand_inds22];
   a11 = partialcorr([result.usc(inds,lv),adjy0cont(inds)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
   a1s = [a1s; a11(1,2)];
   a22 = partialcorr([result.usc(inds,lv),adjy2cont(inds)], covariates(inds,:), 'Rows', 'pairwise');
a2s = [a2s; a22(1,2)];
a33 = partialcorr([result.usc(inds,lv),adjy0pm(inds)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
a3s = [a3s; a33(1,2)];
a44 = partialcorr([result.usc(inds,lv),adjy2pm(inds)] , covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
a4s = [a4s; a44(1,2)];


   b11 = partialcorr([result.usc(inds,lv), dats(inds,1)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
   b1s = [b1s; b11(1,2)];
   b22 = partialcorr([result.usc(inds,lv),dats(inds,2)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
b2s = [b2s; b22(1,2)];
b33 = partialcorr([result.usc(inds,lv),dats(inds,3)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
b3s = [b3s; b33(1,2)];

b44 = partialcorr([result.usc(inds,lv),dats(inds,4)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
b4s = [b4s; b44(1,2)];

end


 gt_contrast2 = (+1).*[a1s, a2s, a3s, a4s];
 gt_contrast2b = (+1).*[b1s, b2s, b3s b4s];

 p1 = 1 - length(find(  (gt_contrast2(:,2)-gt_contrast2(:,1)) > (mean(gt_contrast2(:,4)) - mean(gt_contrast2(:,3)) ) ) )/nboot
p2 = 1 - (length(find(  (gt_contrast2(:,4)>gt_contrast2(:,3))  ) )/nboot);
p3 = 1 - (length(find(  (gt_contrast2(:,2)>gt_contrast2(:,1))  ) )/nboot);
 p4 = 1 - length(find(  (gt_contrast2(:,2)+gt_contrast2(:,1)) > (mean(gt_contrast2(:,4) + (gt_contrast2(:,3)))  )))/nboot
ef = mean(gt_contrast2(:,2))-mean(gt_contrast2(:,1))

 addpath(genpath('~\export_fig-master'));
pmcol = [1 .7 .9];
f1 = figure;
hold on
xss = [1,1.8,3.2,4];
sqlist = [1,2,3,4]; % order of grouping columns
for sq=1:4
if (sq==1 | sq==2); col = [.45 .75 0.75];  % color for columns 1 and 2 (HC)
else
    col = pmcol;
end
k = sqlist(sq);
   mm = mean(gt_contrast2(:,k));
    negs = mm - prctile(gt_contrast2(:,k),5);
     pos = prctile(gt_contrast2(:,k),95) - mm;
     bar(xss(sq),mean(gt_contrast2(:,k)),'FaceColor',col,'EdgeColor',col);
   errorbar(xss(sq),mean(gt_contrast2(:,k)),negs,pos,'.k');
end
set(gca,'XTick',xss,'XTickLabels',{'Baseline (9-10 y)','Follow-up (11-12 y)','Baseline (9-10 y)','Follow-up (11-12 y)'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS group-by-time loading (partial corr)'); ylim([-.2,.22]);
f1.Position = [100 100 500 900];
hold off
export_fig -m5 -transparent depv_pls_Aff_lv1_gtbars_partial_nop.jpg

f2 = figure;
hold on
locs = [.8,1.4,2.6,3.2];
for k=1:4
  col = [.74 .65 .73];
 mmb = mean(gt_contrast2b(:,k));
   negsb = mmb - prctile(gt_contrast2b(:,k),5);
   posb = prctile(gt_contrast2b(:,k),95) - mmb;
   bar(locs(k),mean(gt_contrast2b(:,k)),'FaceColor',col,'EdgeColor',col,'BarWidth',.45);
   errorbar(locs(k),mean(gt_contrast2b(:,k)),negsb,posb,'.k');
end
set(gca,'XTick',locs,'XTickLabels',{'NAc (rew vs nt)','Caud (rew vs nt)','Amyg (neg vs nt face)','Insula (neg vs nt face)'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS neuroaffective loading (partial corr)'); ylim([-1,1]); %axis square
f2.Position = [100 100 500 900];
hold off  
export_fig -m5 -transparent depv_pls_Aff_lv1_ncbars_partial_nop.jpg


%%  ######################################### Threat only #####################################################################################
%%  %% Create matched control group for the threat_only participants
clear all
T_0 = read_abcd_beh_Aff('~/df_filled_beh_Y0_more_neuroaff.csv');
T_2 =read_abcd_beh_Aff('~/df_filled_beh_Y2_more_neuroaff.csv');
s_table = readtable("~/df_filled_demog_withcog_Y0_more.csv");
aff_long = innerjoin(T_0,T_2,"Keys","subkey");
with_2b_aff = find(ismember(string(s_table.subid_x), string(aff_long.subkey)));


usable_subs = aff_long(~isnan(aff_long.mid_NAc_T_2 + aff_long.ng_nt_face_insu_T_2 + ...
    aff_long.mid_NAc_T_0 + aff_long.ng_nt_face_insu_T_0),:); % save as usable_subs_aff.mat

% this section runs iteratively to make the matched group which is then
% saved as num_matches_subs_thrt_only_aff.mat
thrt_high1 = find(usable_subs.cumul_ace_T_0 == 2  )    ; % 2 is threat only 

%load(['~/num_matches_subs_thrt_only_aff.mat']); %uncomment this after the first
% time running the script which generates the list

thrt_high = thrt_high1;
% thrt_high =thrt_high1(num_matches_subs_thrt_aff >1); %uncomment this after the first
% time running the script which generates the list

thrt_high_n = length(thrt_high);
thrt_low = find(usable_subs.cumul_ace_T_0 == 0 ); % 

usable_subs.raceth = usable_subs.White_T_0;
usable_subs.raceth(usable_subs.Black_T_0 == 1) = 2;
usable_subs.raceth(usable_subs.Hispanic_T_0 == 1) = 3;
usable_subs.raceth(usable_subs.Asian_T_0 == 1) = 4;
usable_subs.raceth(usable_subs.Other_T_0 == 1) = 5;
usable_subs.interview_age(isnan(usable_subs.interview_age_T_0)) = nanmean(usable_subs.interview_age_T_0);
usable_subs.agey = round(usable_subs.interview_age_T_0/6);
usable_subs.HighestEd_T_0(isnan(usable_subs.HighestEd_T_0)) = nanmean(usable_subs.HighestEd_T_0);
usable_subs.edu = round(usable_subs.HighestEd_T_0/5);

income_matches = cell(thrt_high_n,1);
age_matches = cell(thrt_high_n,1);
sex_matches = cell(thrt_high_n,1);
raceth_matches = cell(thrt_high_n,1);
edu_matches = cell(thrt_high_n,1);

for k =1:thrt_high_n
    
    if ~isnan(usable_subs.Income_T_0(thrt_high(k)) )
        a = find(round(usable_subs.Income_T_0) == round(usable_subs.Income_T_0(thrt_high(k))));
    else a = thrt_low;
    end
        income_matches{k} = intersect(thrt_low,a);
        
    if ~isnan(usable_subs.agey(thrt_high(k)))
        a = find(usable_subs.agey == usable_subs.agey(thrt_high(k)));
    else a = thrt_low;
    end
        age_matches{k} = intersect(thrt_low,a);
    
    if ~isnan(usable_subs.Male_bin_T_0(thrt_high(k)))
        a = find(usable_subs.Male_bin_T_0 == usable_subs.Male_bin_T_0(thrt_high(k)));
    else a = thrt_low;
    end
         sex_matches{k} = intersect(thrt_low,a);
         
    if ~isnan(usable_subs.raceth(thrt_high(k)))
        a = find(usable_subs.raceth == usable_subs.raceth(thrt_high(k)));
    else a = thrt_low;
    end
         raceth_matches{k} = intersect(thrt_low,a);
    
    if ~isnan(usable_subs.edu(thrt_high(k)))
        a = find(usable_subs.edu == usable_subs.edu(thrt_high(k)));
    else a = thrt_low;
    end
         edu_matches{k} = intersect(thrt_low,a);
            
end
cont_mj = cell(thrt_high_n,1); cont_mj_full = cell(thrt_high_n,1); cont_mj_bare = cell(thrt_high_n,1);
for k = 1:thrt_high_n
    gg = intersect( intersect(sex_matches{k},income_matches{k}), intersect(raceth_matches{k},edu_matches{k}));
    cont_mj_full{k} = setdiff(intersect( gg, age_matches{k}  ),  thrt_high);
    cont_mj{k} = setdiff( gg,  thrt_high);
    gg2 = intersect(income_matches{k}, raceth_matches{k});
    cont_mj_bare{k} = setdiff( intersect(gg2,edu_matches{k}),  thrt_high);

end

num_matches_subs_thrt_aff =[];
for l=1:239
num_matches_subs_thrt_aff = [num_matches_subs_thrt_aff; length(cont_mj_bare{l})];
end
% save list as num_matches_subs_thrt_only_aff.mat

HCss = []; llls =[];
for ww = 1:10000
    HCs =NaN(thrt_high_n,1);
    for k=1:thrt_high_n
        if ~isempty(cont_mj_full{k})
            if length(cont_mj_full{k}) >1
                for jjk = 1: length(cont_mj_full{k})
                    dd = randperm(length(cont_mj_full{k}),1);
                    if ~ismember(cont_mj_full{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj_full{k}(dd);
            
            
        elseif ~isempty(cont_mj{k})
            if length(cont_mj{k}) >1
                for jjk = 1: length(cont_mj{k})
                    dd = randperm(length(cont_mj{k}),1);
                    if ~ismember(cont_mj{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj{k}(dd);
            
            
        elseif ~isempty(cont_mj_bare{k})
            if length(cont_mj_bare{k}) >1
                for jjk = 1: length(cont_mj_bare{k})
                    dd = randperm(length(cont_mj_bare{k}),1);
                    if ~ismember(cont_mj_bare{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj_bare{k}(dd);
        end
    end
    llls = [llls; length(unique(HCs))];
        
        HCss = [HCss  HCs];

end

 matches = HCss(:,1);
for j=1:10000
    [v, w] = unique( matches, 'stable' );
    duplicate_indices = setdiff( 1:numel(matches), w );
    if length(duplicate_indices)==0
        break
    end
    matches(duplicate_indices) = HCss(duplicate_indices,j);
end
        
matches_228 = unique(matches); 
% check if cohen's d is samll for group differences in each category
T_ = usable_subs;
cohd = (nanmean(T_.age_T_0(thrt_high))-nanmean(T_.age_T_0(setdiff(1:4344,thrt_high))) )/ sqrt(nanvar(T_.age_T_0(thrt_high))+nanvar(T_.age_T_0(setdiff(1:4344,thrt_high))) )
cohd = (nanmean(T_.age_T_0(thrt_high))-nanmean(T_.age_T_0(thrt_low)) )/ sqrt(nanvar(T_.age_T_0(thrt_high))+nanvar(T_.age_T_0(thrt_low)) )
cohd = (nanmean(T_.age_T_0(thrt_high))-nanmean(T_.age_T_0(matches_228)) )/ sqrt(nanvar(T_.age_T_0(thrt_high))+nanvar(T_.age_T_0(matches_228)) )

cohd = (nanmean(T_.Income_T_0(thrt_high))-nanmean(T_.Income_T_0(setdiff(1:4344,thrt_high))) )/ sqrt(nanvar(T_.Income_T_0(thrt_high))+nanvar(T_.Income_T_0(setdiff(1:4344,thrt_high))) )
cohd = (nanmean(T_.Income_T_0(thrt_high))-nanmean(T_.Income_T_0(thrt_low)) )/ sqrt(nanvar(T_.Income_T_0(thrt_high))+nanvar(T_.Income_T_0(thrt_low)) )
cohd = (nanmean(T_.Income_T_0(thrt_high))-nanmean(T_.Income_T_0(matches_228)) )/ sqrt(nanvar(T_.Income_T_0(thrt_high))+nanvar(T_.Income_T_0(matches_228)) )


cohd = (nanmean(T_.HighestEd_T_0(thrt_high))-nanmean(T_.HighestEd_T_0(setdiff(1:4344,thrt_high))) )/ sqrt(nanvar(T_.HighestEd_T_0(thrt_high))+nanvar(T_.HighestEd_T_0(setdiff(1:4344,thrt_high))) )
cohd = (nanmean(T_.HighestEd_T_0(thrt_high))-nanmean(T_.HighestEd_T_0(thrt_low)) )/ sqrt(nanvar(T_.HighestEd_T_0(thrt_high))+nanvar(T_.HighestEd_T_0(thrt_low)) )
cohd = (nanmean(T_.HighestEd_T_0(thrt_high))-nanmean(T_.HighestEd_T_0(matches_228)) )/ sqrt(nanvar(T_.HighestEd_T_0(thrt_high))+nanvar(T_.HighestEd_T_0(matches_228)) )

cohd = (nanmean(T_.Black_T_0(thrt_high))-nanmean(T_.Black_T_0(setdiff(1:4344,thrt_high))) )/ sqrt(nanvar(T_.Black_T_0(thrt_high))+nanvar(T_.Black_T_0(setdiff(1:4344,thrt_high))) )
cohd = (nanmean(T_.Black_T_0(thrt_high))-nanmean(T_.Black_T_0(thrt_low)) )/ sqrt(nanvar(T_.Black_T_0(thrt_high))+nanvar(T_.Black_T_0(thrt_low)) )
cohd = (nanmean(T_.Black_T_0(thrt_high))-nanmean(T_.Black_T_0(matches_228)) )/ sqrt(nanvar(T_.Black_T_0(thrt_high))+nanvar(T_.Black_T_0(matches_228)) )

cohd = (nanmean(T_.White_T_0(thrt_high))-nanmean(T_.White_T_0(setdiff(1:4344,thrt_high))) )/ sqrt(nanvar(T_.White_T_0(thrt_high))+nanvar(T_.White_T_0(setdiff(1:4344,thrt_high))) )
cohd = (nanmean(T_.White_T_0(thrt_high))-nanmean(T_.White_T_0(thrt_low)) )/ sqrt(nanvar(T_.White_T_0(thrt_high))+nanvar(T_.White_T_0(thrt_low)) )
cohd = (nanmean(T_.White_T_0(thrt_high))-nanmean(T_.White_T_0(matches_228)) )/ sqrt(nanvar(T_.White_T_0(thrt_high))+nanvar(T_.White_T_0(matches_228)) )

cohd = (nanmean(T_.Hispanic_T_0(thrt_high))-nanmean(T_.Hispanic_T_0(setdiff(1:4344,thrt_high))) )/ sqrt(nanvar(T_.Hispanic_T_0(thrt_high))+nanvar(T_.Hispanic_T_0(setdiff(1:4344,thrt_high))) )
cohd = (nanmean(T_.Hispanic_T_0(thrt_high))-nanmean(T_.Hispanic_T_0(thrt_low)) )/ sqrt(nanvar(T_.Hispanic_T_0(thrt_high))+nanvar(T_.Hispanic_T_0(thrt_low)) )
cohd = (nanmean(T_.Hispanic_T_0(thrt_high))-nanmean(T_.Hispanic_T_0(matches_228)) )/ sqrt(nanvar(T_.Hispanic_T_0(thrt_high))+nanvar(T_.Hispanic_T_0(matches_228)) )


matches_thrt_228_aff = unique(matches); % save these as your matched controls
thrt_232_aff = thrt_high; % save these as your exposure group


%% %%%%% PLS
%% neuroaff latent variable

clear all
T_0 = read_abcd_beh_Aff('~/df_filled_beh_Y0_more_neuroaff.csv');
T_2 =read_abcd_beh_Aff('~/df_filled_beh_Y2_more_neuroaff.csv');


load('~/usable_subs_aff.mat')
T_ = usable_subs;

load('matches_thrt_228_aff.mat')
load('thrt_232_aff.mat')

addpath(genpath('Z:\PM2.5\New Manuscript\Pls'));
rng('default');
groups{1} = [matches_thrt_228_aff];
groups{2} = [thrt_232_aff];
nTimes = 2;

datamat_lst = cell(length(groups),1); 
% ####*****************####

for g = 1:numel(groups)
    for c = 1:nTimes
        for s = 1:numel(groups{g})
            sub_id = char(T_.subkey(groups{g}(s)));
            
            vec2 =[];
            if c==1
                locb = find(T_0.subkey == sub_id);
                if length(locb)==0
                    sub_id
                end
                vec2 = [T_0.mid_NAc(locb)   T_0.mid_Caud(locb)  T_0.ng_nt_face_amyg(locb) T_0.ng_nt_face_insu(locb)];
                                                    
            end
            
            if c==2
                loct = find(T_2.subkey == sub_id);
                if length(loct)==0
                    sub_id
                end
               vec2 = [T_2.mid_NAc(loct)   T_2.mid_Caud(loct)  T_2.ng_nt_face_amyg(loct) T_2.ng_nt_face_insu(loct)];
                  

            end
            
            datamat_lst{g} = [datamat_lst{g}; vec2];
            s;
            
        end
    end
end

num_subj = [length(groups{1})  length(groups{2})];
num_cond = nTimes;
option.method = 1; 
option.num_boot = 1000;
option.num_perm = 1000;
option.meancentering_type=[2]; 
result = pls_analysis(datamat_lst, num_subj, num_cond, option);
result.datamat_lst = datamat_lst;
% save as PLS_result_228low_232hthrt_y0_y2_aff.mat
%% plotting the results (Table 2 bottom half and Figure 3C)
clear all

load('PLS_result_228low_232hthrt_y0_y2_aff.mat');  
ps = result.perm_result.sprob
load('~/usable_subs_aff.mat')
T_ = usable_subs;

load('matches_thrt_228_aff.mat')
load('thrt_232_aff.mat')
matches = matches_thrt_228_aff;
thrthigh = thrt_232_aff;

% *** table 2 bottom half
thrt_low = find(usable_subs.cumul_ace_T_2 == 0 );
T_.AsiOth = T_.Asian_T_0 + T_.Other_T_0;

[nanmean(T_.HighestEd_T_0(thrthigh)) sqrt(nanvar(T_.HighestEd_T_0(thrthigh)))]
[nanmean(T_.HighestEd_T_0(thrt_low)) sqrt(nanvar(T_.HighestEd_T_0(thrt_low)))]
[nanmean(T_.HighestEd_T_0(matches)) sqrt(nanvar(T_.HighestEd_T_0(matches)))]
cohd = (nanmean(T_.HighestEd_T_0(thrthigh))-nanmean(T_.HighestEd_T_0(thrt_low)) )/ sqrt(nanvar(T_.HighestEd_T_0(thrthigh))+nanvar(T_.HighestEd_T_0(thrt_low)) )
cohd = (nanmean(T_.HighestEd_T_0(thrthigh))-nanmean(T_.HighestEd_T_0(matches)) )/ sqrt(nanvar(T_.HighestEd_T_0(thrthigh))+nanvar(T_.HighestEd_T_0(matches)) )

      % levels= 1:10, labels = c("5000", "8500", "14000", "20500", "30000","42500", 7:"62500", 8:"87500", "150000", "200000") )
nanmean(T_.Income_T_0(thrthigh)); 
62500 + (7.8333-7)*(87500 - 62500)
nanmean(T_.Income_T_0(thrt_low)) ;
87500 + (8.308-8)*(150000 - 87500)
nanmean(T_.Income_T_0(matches)); 
62500 + (7.9039-7)*(87500 - 62500)
cohd = (nanmean(T_.Income_T_0(thrthigh))-nanmean(T_.Income_T_0(thrt_low)) )/ sqrt(nanvar(T_.Income_T_0(thrthigh))+nanvar(T_.Income_T_0(thrt_low)) )
cohd = (nanmean(T_.Income_T_0(thrthigh))-nanmean(T_.Income_T_0(matches)) )/ sqrt(nanvar(T_.Income_T_0(thrthigh))+nanvar(T_.Income_T_0(matches)) )

[length(find(T_.White_T_0(thrthigh))) length(find(T_.White_T_0(thrthigh)))/length(thrthigh)]
[length(find(T_.White_T_0(thrt_low))) length(find(T_.White_T_0(thrt_low)))/length(thrt_low)]
[length(find(T_.White_T_0(matches))) length(find(T_.White_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.White_T_0(thrthigh))-nanmean(T_.White_T_0(thrt_low)) )/ sqrt(nanvar(T_.White_T_0(thrthigh))+nanvar(T_.White_T_0(thrt_low)) )
cohd = (nanmean(T_.White_T_0(thrthigh))-nanmean(T_.White_T_0(matches)) )/ sqrt(nanvar(T_.White_T_0(thrthigh))+nanvar(T_.White_T_0(matches)) )


[length(find(T_.Hispanic_T_0(thrthigh))) length(find(T_.Hispanic_T_0(thrthigh)))/length(thrthigh)]
[length(find(T_.Hispanic_T_0(thrt_low))) length(find(T_.Hispanic_T_0(thrt_low)))/length(thrt_low)]
[length(find(T_.Hispanic_T_0(matches))) length(find(T_.Hispanic_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.Hispanic_T_0(thrthigh))-nanmean(T_.Hispanic_T_0(thrt_low)) )/ sqrt(nanvar(T_.Hispanic_T_0(thrthigh))+nanvar(T_.Hispanic_T_0(thrt_low)) )
cohd = (nanmean(T_.Hispanic_T_0(thrthigh))-nanmean(T_.Hispanic_T_0(matches)) )/ sqrt(nanvar(T_.Hispanic_T_0(thrthigh))+nanvar(T_.Hispanic_T_0(matches)) )


[length(find(T_.Black_T_0(thrthigh))) length(find(T_.Black_T_0(thrthigh)))/length(thrthigh)]
[length(find(T_.Black_T_0(thrt_low))) length(find(T_.Black_T_0(thrt_low)))/length(thrt_low)]
[length(find(T_.Black_T_0(matches))) length(find(T_.Black_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.Black_T_0(thrthigh))-nanmean(T_.Black_T_0(thrt_low)) )/ sqrt(nanvar(T_.Black_T_0(thrthigh))+nanvar(T_.Black_T_0(thrt_low)) )
cohd = (nanmean(T_.Black_T_0(thrthigh))-nanmean(T_.Black_T_0(matches)) )/ sqrt(nanvar(T_.Black_T_0(thrthigh))+nanvar(T_.Black_T_0(matches)) )


[length(find(T_.AsiOth(thrthigh))) length(find(T_.AsiOth(thrthigh)))/length(thrthigh)]
[length(find(T_.AsiOth(thrt_low))) length(find(T_.AsiOth(thrt_low)))/length(thrt_low)]
[length(find(T_.AsiOth(matches))) length(find(T_.AsiOth(matches)))/length(matches)]
cohd = (nanmean(T_.AsiOth(thrthigh))-nanmean(T_.AsiOth(thrt_low)) )/ sqrt(nanvar(T_.AsiOth(thrthigh))+nanvar(T_.AsiOth(thrt_low)) )
cohd = (nanmean(T_.AsiOth(thrthigh))-nanmean(T_.AsiOth(matches)) )/ sqrt(nanvar(T_.AsiOth(thrthigh))+nanvar(T_.AsiOth(matches)) )

[length(find(~T_.Male_bin_T_0(thrthigh))) length(find(~T_.Male_bin_T_0(thrthigh)))/length(thrthigh)]
[length(find(~T_.Male_bin_T_0(thrt_low))) length(find(~T_.Male_bin_T_0(thrt_low)))/length(thrt_low)]
[length(find(~T_.Male_bin_T_0(matches))) length(find(~T_.Male_bin_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.Male_bin_T_0(thrthigh))-nanmean(T_.Male_bin_T_0(thrt_low)) )/ sqrt(nanvar(T_.Male_bin_T_0(thrthigh))+nanvar(T_.Male_bin_T_0(thrt_low)) )
cohd = (nanmean(T_.Male_bin_T_0(thrthigh))-nanmean(T_.Male_bin_T_0(matches)) )/ sqrt(nanvar(T_.Male_bin_T_0(thrthigh))+nanvar(T_.Male_bin_T_0(matches)) )


[nanmean(T_.age_T_0(thrthigh)) sqrt(nanvar(T_.age_T_0(thrthigh)))]
[nanmean(T_.age_T_0(thrt_low)) sqrt(nanvar(T_.age_T_0(thrt_low)))]
[nanmean(T_.age_T_0(matches)) sqrt(nanvar(T_.age_T_0(matches)))]
cohd = (nanmean(T_.age_T_0(thrthigh))-nanmean(T_.age_T_0(thrt_low)) )/ sqrt(nanvar(T_.age_T_0(thrthigh))+nanvar(T_.age_T_0(thrt_low)) )
cohd = (nanmean(T_.age_T_0(thrthigh))-nanmean(T_.age_T_0(matches)) )/ sqrt(nanvar(T_.age_T_0(thrthigh))+nanvar(T_.age_T_0(matches)) )

[nanmean(T_.age_T_2(thrthigh)) sqrt(nanvar(T_.age_T_2(thrthigh)))]
[nanmean(T_.age_T_2(thrt_low)) sqrt(nanvar(T_.age_T_2(thrt_low)))]
[nanmean(T_.age_T_2(matches)) sqrt(nanvar(T_.age_T_2(matches)))]
cohd = (nanmean(T_.age_T_2(thrthigh))-nanmean(T_.age_T_2(thrt_low)) )/ sqrt(nanvar(T_.age_T_2(thrthigh))+nanvar(T_.age_T_2(thrt_low)) )
cohd = (nanmean(T_.age_T_2(thrthigh))-nanmean(T_.age_T_2(matches)) )/ sqrt(nanvar(T_.age_T_2(thrthigh))+nanvar(T_.age_T_2(matches)) )

% ***


Sites = zeros(3244,22);
for K=1:22
    Sites((T_.site_id_l_T_0 == K),K) = 1;
end

lv=1;   %~#

lg1 = length(matches); lg2 = length(thrthigh);
adjy0cont = zeros(2*(lg1+lg2),1);
adjy0cont(1:lg1) = 1;

adjy2cont = zeros(2*(lg1+lg2),1);
adjy2cont(lg1+1:2*lg1) = 1;

adjy0pm = zeros(2*(lg1+lg2),1);
adjy0pm(2*lg1+1:2*lg1+lg2) = 1;

adjy2pm = zeros(2*(lg1+lg2),1);
adjy2pm(2*lg1+lg2+1:2*(lg1+lg2)) = 1;


covariates_y0cont = [  T_.tfmri_mid_all_meanmotion_T_0(matches)   T_.tfmri_nback_all_meanmotion_T_0(matches)   Sites(matches,:)  T_.White_T_0(matches) T_.Black_T_0(matches) T_.Hispanic_T_0(matches) T_.HighestEd_T_0(matches) T_.Income_T_0(matches) T_.Male_bin_T_0(matches)  T_.interview_age_T_0(matches)  T_.Ingenia_T_0(matches) T_.Achieva_T_0(matches) T_.Discovery_T_0(matches)  T_.Pfit_T_0(matches) T_.reshist_addr1_adi_wsum_T_0(matches)  T_.pds_ss_T_0(matches)]; % T_.PF10_lavaan(matches)
covariates_y2cont = [ T_.tfmri_mid_all_meanmotion_T_2(matches)   T_.tfmri_nback_all_meanmotion_T_2(matches)   Sites(matches,:)  T_.White_T_0(matches) T_.Black_T_0(matches) T_.Hispanic_T_0(matches) T_.HighestEd_T_0(matches) T_.Income_T_0(matches) T_.Male_bin_T_0(matches)  T_.interview_age_T_2(matches)   T_.Ingenia_T_0(matches) T_.Achieva_T_0(matches) T_.Discovery_T_0(matches)  T_.Pfit_T_0(matches)  T_.reshist_addr1_adi_wsum_T_0(matches) T_.pds_ss_T_2(matches)];
covariates_y0pm = [   T_.tfmri_mid_all_meanmotion_T_0(thrthigh)   T_.tfmri_nback_all_meanmotion_T_0(thrthigh)   Sites(thrthigh,:)  T_.White_T_0(thrthigh) T_.Black_T_0(thrthigh) T_.Hispanic_T_0(thrthigh) T_.HighestEd_T_0(thrthigh) T_.Income_T_0(thrthigh) T_.Male_bin_T_0(thrthigh)  T_.interview_age_T_0(thrthigh)  T_.Ingenia_T_0(thrthigh) T_.Achieva_T_0(thrthigh) T_.Discovery_T_0(thrthigh)  T_.Pfit_T_0(thrthigh)  T_.reshist_addr1_adi_wsum_T_0(thrthigh)  T_.pds_ss_T_0(thrthigh)];
covariates_y2pm = [ T_.tfmri_mid_all_meanmotion_T_2(thrthigh)   T_.tfmri_nback_all_meanmotion_T_2(thrthigh)   Sites(thrthigh,:)  T_.White_T_0(thrthigh) T_.Black_T_0(thrthigh) T_.Hispanic_T_0(thrthigh) T_.HighestEd_T_0(thrthigh) T_.Income_T_0(thrthigh) T_.Male_bin_T_0(thrthigh)  T_.interview_age_T_2(thrthigh)   T_.Ingenia_T_0(thrthigh) T_.Achieva_T_0(thrthigh) T_.Discovery_T_0(thrthigh)  T_.Pfit_T_0(thrthigh)  T_.reshist_addr1_adi_wsum_T_0(thrthigh) T_.pds_ss_T_2(thrthigh)];



covariates = [covariates_y0cont; covariates_y2cont; covariates_y0pm; covariates_y2pm];
% 24,25,26, are race/ethnic
% exclude_race_as_cov = [1:23,27:37];
exclude_race_as_cov = [1:size(covariates,2)];
[a1 ] = partialcorr([result.usc(:,lv),adjy0cont],covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a2 ] = partialcorr([result.usc(:,lv),adjy2cont], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a3 ] = partialcorr([result.usc(:,lv),adjy0pm], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a4 ] = partialcorr([result.usc(:,lv),adjy2pm], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');

a1s =[]; a2s =[]; a3s =[]; a4s =[]; b1s=[];b2s=[];b3s=[]; b4s=[];
dats = [result.datamat_lst{1, 1}; result.datamat_lst{2, 1}];

nboot = 5000;
for jj=1:nboot
   rand_inds1 = randi([1,lg1],lg1,1); rand_inds11 = randi([1,lg1],lg1,1); 
   rand_inds2 = randi([1,lg2],lg2,1); rand_inds22 = randi([1,lg2],lg2,1);
   inds = [rand_inds1; lg1+ rand_inds11; 2*lg1+ rand_inds2; 2*lg1+lg2+ rand_inds22];
   a11 = partialcorr([result.usc(inds,lv),adjy0cont(inds)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
   a1s = [a1s; a11(1,2)];
   a22 = partialcorr([result.usc(inds,lv),adjy2cont(inds)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
a2s = [a2s; a22(1,2)];
a33 = partialcorr([result.usc(inds,lv),adjy0pm(inds)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
a3s = [a3s; a33(1,2)];
a44 = partialcorr([result.usc(inds,lv),adjy2pm(inds)] , covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
a4s = [a4s; a44(1,2)];


   b11 = partialcorr([result.usc(inds,lv), dats(inds,1)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
   b1s = [b1s; b11(1,2)];
   b22 = partialcorr([result.usc(inds,lv),dats(inds,2)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
b2s = [b2s; b22(1,2)];
b33 = partialcorr([result.usc(inds,lv),dats(inds,3)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
b3s = [b3s; b33(1,2)];

b44 = partialcorr([result.usc(inds,lv),dats(inds,4)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
b4s = [b4s; b44(1,2)];

end

 
 gt_contrast2 = (-1).*[a1s, a2s, a3s, a4s];
 gt_contrast2b = (-1).*[b1s, b2s, b3s b4s];

 p1 = 1 - length(find(  (gt_contrast2(:,2)-gt_contrast2(:,1)) > (mean(gt_contrast2(:,4)) - mean(gt_contrast2(:,3)) ) ) )/nboot
p2 = 1 - (length(find(  (gt_contrast2(:,4)>gt_contrast2(:,3))  ) )/nboot);
p3 = 1 - (length(find(  (gt_contrast2(:,2)>gt_contrast2(:,1))  ) )/nboot);
 p4 = 1 - length(find(  (gt_contrast2(:,2)+gt_contrast2(:,1)) > (mean(gt_contrast2(:,4) + (gt_contrast2(:,3)))  )))/nboot
ef = mean(gt_contrast2(:,2))-mean(gt_contrast2(:,1))

 addpath(genpath('Z:\MatlabGraphics\export_fig-master'));
pmcol = [.97 .57 .31];
f1 = figure;
hold on
xss = [1,1.8,3.2,4];
sqlist = [1,2,3,4]; % order of grouping columns
for sq=1:4
if (sq==1 | sq==2); col = [.45 .75 0.75];  % color for columns 1 and 2 (HC)
else
    col = pmcol;
end
k = sqlist(sq);
   mm = mean(gt_contrast2(:,k));
 negs = mm - prctile(gt_contrast2(:,k),5);
  pos = prctile(gt_contrast2(:,k),95) - mm;
  bar(xss(sq),mean(gt_contrast2(:,k)),'FaceColor',col,'EdgeColor',col);
   errorbar(xss(sq),mean(gt_contrast2(:,k)),negs,pos,'.k');
end
set(gca,'XTick',xss,'XTickLabels',{'Baseline (9-10 y)','Follow-up (11-12 y)','Baseline (9-10 y)','Follow-up (11-12 y)'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS group-by-time loading (partial corr)'); ylim([-.2,.22]);
f1.Position = [100 100 500 900];
hold off
export_fig -m5 -transparent thrt_pls_Aff_lv1_gtbars_partial_nop.jpg

f2 = figure;
hold on
locs = [.8,1.4,2.6,3.2];
for k=1:4
    col = [.74 .65 .73];
 mmb = mean(gt_contrast2b(:,k));
    negsb = mmb - prctile(gt_contrast2b(:,k),5);
     posb = prctile(gt_contrast2b(:,k),95) - mmb;
     bar(locs(k),mean(gt_contrast2b(:,k)),'FaceColor',col,'EdgeColor',col,'BarWidth',.45);
     errorbar(locs(k),mean(gt_contrast2b(:,k)),negsb,posb,'.k');
end
set(gca,'XTick',locs,'XTickLabels',{'NAc (rew vs nt)','Caud (rew vs nt)','Amyg (neg vs nt face)','Insula (neg vs nt face)'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS neuroaffective loading (partial corr)'); ylim([-1,1]); %axis square
f2.Position = [100 100 500 900];
hold off
export_fig -m5 -transparent thrt_pls_Aff_lv1_ncbars_partial_nop.jpg

%%  ######################################### Threat plus #####################################################################################
%%  %% Create matched control group for the threat_plus participants
clear all
T_0 = read_abcd_beh_Aff('~/df_filled_beh_Y0_more_neuroaff.csv');
T_2 =read_abcd_beh_Aff('~/df_filled_beh_Y2_more_neuroaff.csv');
s_table = readtable("~/df_filled_demog_withcog_Y0_more.csv");
aff_long = innerjoin(T_0,T_2,"Keys","subkey");
with_2b_aff = find(ismember(string(s_table.subid_x), string(aff_long.subkey)));


usable_subs = aff_long(~isnan(aff_long.mid_NAc_T_2 + aff_long.ng_nt_face_insu_T_2 + ...
    aff_long.mid_NAc_T_0 + aff_long.ng_nt_face_insu_T_0),:); % save as usable_subs_aff.mat


% this section runs iteratively to make the maatched group which is then
% saved as num_matches_subs_thplus_aff.mat
thplus_high1 = find(usable_subs.cumul_ace_T_0 > 1  )    ; 

%load(['~/num_matches_subs_thplus_aff.mat']);%uncomment this after the first
% time running the script which generates the list

thplus_high = thplus_high1;
% thplus_high =thplus_high1(num_matches_subs_thplus_aff >1); %uncomment this after the first
% time running the script which generates the list
 
thplus_high_n = length(thplus_high);
thplus_low = find(usable_subs.cumul_ace_T_0 == 0 ); 

usable_subs.raceth = usable_subs.White_T_0;
usable_subs.raceth(usable_subs.Black_T_0 == 1) = 2;
usable_subs.raceth(usable_subs.Hispanic_T_0 == 1) = 3;
usable_subs.raceth(usable_subs.Asian_T_0 == 1) = 4;
usable_subs.raceth(usable_subs.Other_T_0 == 1) = 5;
usable_subs.interview_age(isnan(usable_subs.interview_age_T_0)) = nanmean(usable_subs.interview_age_T_0);
usable_subs.agey = round(usable_subs.interview_age_T_0/6);
usable_subs.HighestEd_T_0(isnan(usable_subs.HighestEd_T_0)) = nanmean(usable_subs.HighestEd_T_0);
usable_subs.edu = round(usable_subs.HighestEd_T_0/5);

income_matches = cell(thplus_high_n,1);
age_matches = cell(thplus_high_n,1);
sex_matches = cell(thplus_high_n,1);
raceth_matches = cell(thplus_high_n,1);
edu_matches = cell(thplus_high_n,1);

for k =1:thplus_high_n
    
    if ~isnan(usable_subs.Income_T_0(thplus_high(k)) )
        a = find(round(usable_subs.Income_T_0) == round(usable_subs.Income_T_0(thplus_high(k))));
    else a = thplus_low;
    end
        income_matches{k} = intersect(thplus_low,a);
        
    if ~isnan(usable_subs.agey(thplus_high(k)))
        a = find(usable_subs.agey == usable_subs.agey(thplus_high(k)));
    else a = thplus_low;
    end
        age_matches{k} = intersect(thplus_low,a);
    
    if ~isnan(usable_subs.Male_bin_T_0(thplus_high(k)))
        a = find(usable_subs.Male_bin_T_0 == usable_subs.Male_bin_T_0(thplus_high(k)));
    else a = thplus_low;
    end
         sex_matches{k} = intersect(thplus_low,a);
         
    if ~isnan(usable_subs.raceth(thplus_high(k)))
        a = find(usable_subs.raceth == usable_subs.raceth(thplus_high(k)));
    else a = thplus_low;
    end
         raceth_matches{k} = intersect(thplus_low,a);
    
    if ~isnan(usable_subs.edu(thplus_high(k)))
        a = find(usable_subs.edu == usable_subs.edu(thplus_high(k)));
    else a = thplus_low;
    end
         edu_matches{k} = intersect(thplus_low,a);
            
end
cont_mj = cell(thplus_high_n,1); cont_mj_full = cell(thplus_high_n,1); cont_mj_bare = cell(thplus_high_n,1);
for k = 1:thplus_high_n
    gg = intersect( intersect(sex_matches{k},income_matches{k}), intersect(raceth_matches{k},edu_matches{k}));
    cont_mj_full{k} = setdiff(intersect( gg, age_matches{k}  ),  thplus_high);
    cont_mj{k} = setdiff( gg,  thplus_high);
    gg2 = intersect(income_matches{k}, raceth_matches{k});
    cont_mj_bare{k} = setdiff( intersect(gg2,edu_matches{k}),  thplus_high);

end

num_matches_subs_thplus_aff =[];
for l=1:494
num_matches_subs_thplus_aff = [num_matches_subs_thplus_aff; length(cont_mj_bare{l})];
end
% save as num_matches_subs_thplus_aff.mat

HCss = []; llls =[];
for ww = 1:10000
    HCs =NaN(thplus_high_n,1);
    for k=1:thplus_high_n
        if ~isempty(cont_mj_full{k})
            if length(cont_mj_full{k}) >1
                for jjk = 1: length(cont_mj_full{k})
                    dd = randperm(length(cont_mj_full{k}),1);
                    if ~ismember(cont_mj_full{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj_full{k}(dd);
            
            
        elseif ~isempty(cont_mj{k})
            if length(cont_mj{k}) >1
                for jjk = 1: length(cont_mj{k})
                    dd = randperm(length(cont_mj{k}),1);
                    if ~ismember(cont_mj{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj{k}(dd);
            
            
        elseif ~isempty(cont_mj_bare{k})
            if length(cont_mj_bare{k}) >1
                for jjk = 1: length(cont_mj_bare{k})
                    dd = randperm(length(cont_mj_bare{k}),1);
                    if ~ismember(cont_mj_bare{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj_bare{k}(dd);
        end
    end
    llls = [llls; length(unique(HCs))];
        
        HCss = [HCss  HCs];

end

 matches = HCss(:,1);
for j=1:10000
    [v, w] = unique( matches, 'stable' );
    duplicate_indices = setdiff( 1:numel(matches), w );
    if length(duplicate_indices)==0
        break
    end
    matches(duplicate_indices) = HCss(duplicate_indices,j);
end
% check if cohen's d is samll for group differences in each category         
matches_423 = unique(matches); 

T_ = usable_subs;
cohd = (nanmean(T_.age_T_0(thplus_high))-nanmean(T_.age_T_0(setdiff(1:4344,thplus_high))) )/ sqrt(nanvar(T_.age_T_0(thplus_high))+nanvar(T_.age_T_0(setdiff(1:4344,thplus_high))) )
cohd = (nanmean(T_.age_T_0(thplus_high))-nanmean(T_.age_T_0(thplus_low)) )/ sqrt(nanvar(T_.age_T_0(thplus_high))+nanvar(T_.age_T_0(thplus_low)) )
cohd = (nanmean(T_.age_T_0(thplus_high))-nanmean(T_.age_T_0(matches_423)) )/ sqrt(nanvar(T_.age_T_0(thplus_high))+nanvar(T_.age_T_0(matches_423)) )

cohd = (nanmean(T_.Income_T_0(thplus_high))-nanmean(T_.Income_T_0(setdiff(1:4344,thplus_high))) )/ sqrt(nanvar(T_.Income_T_0(thplus_high))+nanvar(T_.Income_T_0(setdiff(1:4344,thplus_high))) )
cohd = (nanmean(T_.Income_T_0(thplus_high))-nanmean(T_.Income_T_0(thplus_low)) )/ sqrt(nanvar(T_.Income_T_0(thplus_high))+nanvar(T_.Income_T_0(thplus_low)) )
cohd = (nanmean(T_.Income_T_0(thplus_high))-nanmean(T_.Income_T_0(matches_423)) )/ sqrt(nanvar(T_.Income_T_0(thplus_high))+nanvar(T_.Income_T_0(matches_423)) )


cohd = (nanmean(T_.HighestEd_T_0(thplus_high))-nanmean(T_.HighestEd_T_0(setdiff(1:4344,thplus_high))) )/ sqrt(nanvar(T_.HighestEd_T_0(thplus_high))+nanvar(T_.HighestEd_T_0(setdiff(1:4344,thplus_high))) )
cohd = (nanmean(T_.HighestEd_T_0(thplus_high))-nanmean(T_.HighestEd_T_0(thplus_low)) )/ sqrt(nanvar(T_.HighestEd_T_0(thplus_high))+nanvar(T_.HighestEd_T_0(thplus_low)) )
cohd = (nanmean(T_.HighestEd_T_0(thplus_high))-nanmean(T_.HighestEd_T_0(matches_423)) )/ sqrt(nanvar(T_.HighestEd_T_0(thplus_high))+nanvar(T_.HighestEd_T_0(matches_423)) )

cohd = (nanmean(T_.Black_T_0(thplus_high))-nanmean(T_.Black_T_0(setdiff(1:4344,thplus_high))) )/ sqrt(nanvar(T_.Black_T_0(thplus_high))+nanvar(T_.Black_T_0(setdiff(1:4344,thplus_high))) )
cohd = (nanmean(T_.Black_T_0(thplus_high))-nanmean(T_.Black_T_0(thplus_low)) )/ sqrt(nanvar(T_.Black_T_0(thplus_high))+nanvar(T_.Black_T_0(thplus_low)) )
cohd = (nanmean(T_.Black_T_0(thplus_high))-nanmean(T_.Black_T_0(matches_423)) )/ sqrt(nanvar(T_.Black_T_0(thplus_high))+nanvar(T_.Black_T_0(matches_423)) )

cohd = (nanmean(T_.White_T_0(thplus_high))-nanmean(T_.White_T_0(setdiff(1:4344,thplus_high))) )/ sqrt(nanvar(T_.White_T_0(thplus_high))+nanvar(T_.White_T_0(setdiff(1:4344,thplus_high))) )
cohd = (nanmean(T_.White_T_0(thplus_high))-nanmean(T_.White_T_0(thplus_low)) )/ sqrt(nanvar(T_.White_T_0(thplus_high))+nanvar(T_.White_T_0(thplus_low)) )
cohd = (nanmean(T_.White_T_0(thplus_high))-nanmean(T_.White_T_0(matches_423)) )/ sqrt(nanvar(T_.White_T_0(thplus_high))+nanvar(T_.White_T_0(matches_423)) )

cohd = (nanmean(T_.Hispanic_T_0(thplus_high))-nanmean(T_.Hispanic_T_0(setdiff(1:4344,thplus_high))) )/ sqrt(nanvar(T_.Hispanic_T_0(thplus_high))+nanvar(T_.Hispanic_T_0(setdiff(1:4344,thplus_high))) )
cohd = (nanmean(T_.Hispanic_T_0(thplus_high))-nanmean(T_.Hispanic_T_0(thplus_low)) )/ sqrt(nanvar(T_.Hispanic_T_0(thplus_high))+nanvar(T_.Hispanic_T_0(thplus_low)) )
cohd = (nanmean(T_.Hispanic_T_0(thplus_high))-nanmean(T_.Hispanic_T_0(matches_423)) )/ sqrt(nanvar(T_.Hispanic_T_0(thplus_high))+nanvar(T_.Hispanic_T_0(matches_423)) )


matches_thplus_423_aff = unique(matches); % save these as your matched controls
thplus_464_aff = thplus_high; % save these as your exposure 


%% %%%%% PLS
%% neuroaff latent variable

clear all
T_0 = read_abcd_beh_Aff('~/df_filled_beh_Y0_more_neuroaff.csv');
T_2 =read_abcd_beh_Aff('~/df_filled_beh_Y2_more_neuroaff.csv');


load('~/usable_subs_aff.mat')
T_ = usable_subs;

load('matches_thplus_423_aff.mat')
load('thplus_464_aff.mat')

addpath(genpath('~\PLS_folder\Pls'));
rng('default');
groups{1} = [matches_thplus_423_aff];
groups{2} = [thplus_464_aff];
nTimes = 2;

datamat_lst = cell(length(groups),1); 
% ####*****************####

for g = 1:numel(groups)
    for c = 1:nTimes
        for s = 1:numel(groups{g})
            sub_id = char(T_.subkey(groups{g}(s)));
            
            vec2 =[];
            if c==1
                locb = find(T_0.subkey == sub_id);
                if length(locb)==0
                    sub_id
                end
                vec2 = [T_0.mid_NAc(locb)   T_0.mid_Caud(locb)  T_0.ng_nt_face_amyg(locb) T_0.ng_nt_face_insu(locb)];
                                                    
            end
            
            if c==2
                loct = find(T_2.subkey == sub_id);
                if length(loct)==0
                    sub_id
                end
               vec2 = [T_2.mid_NAc(loct)   T_2.mid_Caud(loct)  T_2.ng_nt_face_amyg(loct) T_2.ng_nt_face_insu(loct)];
                  

            end
            
            datamat_lst{g} = [datamat_lst{g}; vec2];
            s;
            
        end
    end
end

num_subj = [length(groups{1})  length(groups{2})];
num_cond = nTimes;
option.method = 1; 
option.num_boot = 1000;
option.num_perm = 1000;
option.meancentering_type=[2];
result = pls_analysis(datamat_lst, num_subj, num_cond, option);
result.datamat_lst = datamat_lst;
% save as PLS_result_423low_464hthplus_y0_y2_aff.mat

%% plotting the results (Table 3 bottom half and Figure 3D)
clear all

load('PLS_result_423low_464hthplus_y0_y2_aff.mat');  
ps = result.perm_result.sprob
load('~/usable_subs_aff.mat')
T_ = usable_subs;

load('matches_thplus_423_aff.mat')
load('thplus_464_aff.mat')
matches = matches_thplus_423_aff;
thplus = thplus_464_aff;

% *** table 3 bottom half
thrt_low = find(usable_subs.cumul_ace_T_2 == 0 );
T_.AsiOth = T_.Asian_T_0 + T_.Other_T_0;

[nanmean(T_.HighestEd_T_0(thplus)) sqrt(nanvar(T_.HighestEd_T_0(thplus)))]
[nanmean(T_.HighestEd_T_0(thrt_low)) sqrt(nanvar(T_.HighestEd_T_0(thrt_low)))]
[nanmean(T_.HighestEd_T_0(matches)) sqrt(nanvar(T_.HighestEd_T_0(matches)))]
cohd = (nanmean(T_.HighestEd_T_0(thplus))-nanmean(T_.HighestEd_T_0(thrt_low)) )/ sqrt(nanvar(T_.HighestEd_T_0(thplus))+nanvar(T_.HighestEd_T_0(thrt_low)) )
cohd = (nanmean(T_.HighestEd_T_0(thplus))-nanmean(T_.HighestEd_T_0(matches)) )/ sqrt(nanvar(T_.HighestEd_T_0(thplus))+nanvar(T_.HighestEd_T_0(matches)) )

      % levels= 1:10, labels = c("5000", "8500", "14000", "20500", "30000","42500", 7:"62500", 8:"87500", "150000", "200000") )
nanmean(T_.Income_T_0(thplus)); 
42500 + (6.9199-6)*(62500 - 42500)
nanmean(T_.Income_T_0(thrt_low)) ;
87500 + (8.308-8)*(150000 - 87500)
nanmean(T_.Income_T_0(matches)); 
62500 + ( 7.1351-7)*(87500 - 62500)
cohd = (nanmean(T_.Income_T_0(thplus))-nanmean(T_.Income_T_0(thrt_low)) )/ sqrt(nanvar(T_.Income_T_0(thplus))+nanvar(T_.Income_T_0(thrt_low)) )
cohd = (nanmean(T_.Income_T_0(thplus))-nanmean(T_.Income_T_0(matches)) )/ sqrt(nanvar(T_.Income_T_0(thplus))+nanvar(T_.Income_T_0(matches)) )

[length(find(T_.White_T_0(thplus))) length(find(T_.White_T_0(thplus)))/length(thplus)]
[length(find(T_.White_T_0(thrt_low))) length(find(T_.White_T_0(thrt_low)))/length(thrt_low)]
[length(find(T_.White_T_0(matches))) length(find(T_.White_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.White_T_0(thplus))-nanmean(T_.White_T_0(thrt_low)) )/ sqrt(nanvar(T_.White_T_0(thplus))+nanvar(T_.White_T_0(thrt_low)) )
cohd = (nanmean(T_.White_T_0(thplus))-nanmean(T_.White_T_0(matches)) )/ sqrt(nanvar(T_.White_T_0(thplus))+nanvar(T_.White_T_0(matches)) )


[length(find(T_.Hispanic_T_0(thplus))) length(find(T_.Hispanic_T_0(thplus)))/length(thplus)]
[length(find(T_.Hispanic_T_0(thrt_low))) length(find(T_.Hispanic_T_0(thrt_low)))/length(thrt_low)]
[length(find(T_.Hispanic_T_0(matches))) length(find(T_.Hispanic_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.Hispanic_T_0(thplus))-nanmean(T_.Hispanic_T_0(thrt_low)) )/ sqrt(nanvar(T_.Hispanic_T_0(thplus))+nanvar(T_.Hispanic_T_0(thrt_low)) )
cohd = (nanmean(T_.Hispanic_T_0(thplus))-nanmean(T_.Hispanic_T_0(matches)) )/ sqrt(nanvar(T_.Hispanic_T_0(thplus))+nanvar(T_.Hispanic_T_0(matches)) )


[length(find(T_.Black_T_0(thplus))) length(find(T_.Black_T_0(thplus)))/length(thplus)]
[length(find(T_.Black_T_0(thrt_low))) length(find(T_.Black_T_0(thrt_low)))/length(thrt_low)]
[length(find(T_.Black_T_0(matches))) length(find(T_.Black_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.Black_T_0(thplus))-nanmean(T_.Black_T_0(thrt_low)) )/ sqrt(nanvar(T_.Black_T_0(thplus))+nanvar(T_.Black_T_0(thrt_low)) )
cohd = (nanmean(T_.Black_T_0(thplus))-nanmean(T_.Black_T_0(matches)) )/ sqrt(nanvar(T_.Black_T_0(thplus))+nanvar(T_.Black_T_0(matches)) )


[length(find(T_.AsiOth(thplus))) length(find(T_.AsiOth(thplus)))/length(thplus)]
[length(find(T_.AsiOth(thrt_low))) length(find(T_.AsiOth(thrt_low)))/length(thrt_low)]
[length(find(T_.AsiOth(matches))) length(find(T_.AsiOth(matches)))/length(matches)]
cohd = (nanmean(T_.AsiOth(thplus))-nanmean(T_.AsiOth(thrt_low)) )/ sqrt(nanvar(T_.AsiOth(thplus))+nanvar(T_.AsiOth(thrt_low)) )
cohd = (nanmean(T_.AsiOth(thplus))-nanmean(T_.AsiOth(matches)) )/ sqrt(nanvar(T_.AsiOth(thplus))+nanvar(T_.AsiOth(matches)) )

[length(find(~T_.Male_bin_T_0(thplus))) length(find(~T_.Male_bin_T_0(thplus)))/length(thplus)]
[length(find(~T_.Male_bin_T_0(thrt_low))) length(find(~T_.Male_bin_T_0(thrt_low)))/length(thrt_low)]
[length(find(~T_.Male_bin_T_0(matches))) length(find(~T_.Male_bin_T_0(matches)))/length(matches)]
cohd = (nanmean(T_.Male_bin_T_0(thplus))-nanmean(T_.Male_bin_T_0(thrt_low)) )/ sqrt(nanvar(T_.Male_bin_T_0(thplus))+nanvar(T_.Male_bin_T_0(thrt_low)) )
cohd = (nanmean(T_.Male_bin_T_0(thplus))-nanmean(T_.Male_bin_T_0(matches)) )/ sqrt(nanvar(T_.Male_bin_T_0(thplus))+nanvar(T_.Male_bin_T_0(matches)) )


[nanmean(T_.age_T_0(thplus)) sqrt(nanvar(T_.age_T_0(thplus)))]
[nanmean(T_.age_T_0(thrt_low)) sqrt(nanvar(T_.age_T_0(thrt_low)))]
[nanmean(T_.age_T_0(matches)) sqrt(nanvar(T_.age_T_0(matches)))]
cohd = (nanmean(T_.age_T_0(thplus))-nanmean(T_.age_T_0(thrt_low)) )/ sqrt(nanvar(T_.age_T_0(thplus))+nanvar(T_.age_T_0(thrt_low)) )
cohd = (nanmean(T_.age_T_0(thplus))-nanmean(T_.age_T_0(matches)) )/ sqrt(nanvar(T_.age_T_0(thplus))+nanvar(T_.age_T_0(matches)) )

[nanmean(T_.age_T_2(thplus)) sqrt(nanvar(T_.age_T_2(thplus)))]
[nanmean(T_.age_T_2(thrt_low)) sqrt(nanvar(T_.age_T_2(thrt_low)))]
[nanmean(T_.age_T_2(matches)) sqrt(nanvar(T_.age_T_2(matches)))]
cohd = (nanmean(T_.age_T_2(thplus))-nanmean(T_.age_T_2(thrt_low)) )/ sqrt(nanvar(T_.age_T_2(thplus))+nanvar(T_.age_T_2(thrt_low)) )
cohd = (nanmean(T_.age_T_2(thplus))-nanmean(T_.age_T_2(matches)) )/ sqrt(nanvar(T_.age_T_2(thplus))+nanvar(T_.age_T_2(matches)) )

% ***


Sites = zeros(4344,22);
for K=1:22
    Sites((T_.site_id_l_T_0 == K),K) = 1;
end

lv=1;   %~#

lg1 = length(matches); lg2 = length(thplus);
adjy0cont = zeros(2*(lg1+lg2),1);
adjy0cont(1:lg1) = 1;

adjy2cont = zeros(2*(lg1+lg2),1);
adjy2cont(lg1+1:2*lg1) = 1;

adjy0pm = zeros(2*(lg1+lg2),1);
adjy0pm(2*lg1+1:2*lg1+lg2) = 1;

adjy2pm = zeros(2*(lg1+lg2),1);
adjy2pm(2*lg1+lg2+1:2*(lg1+lg2)) = 1;

T_.depv = zeros(4344,1);
T_.depv(find(T_.cumul_ace_T_0 == 1))=1;
covariates_y0cont = [T_.depv(matches)  T_.tfmri_mid_all_meanmotion_T_0(matches)   T_.tfmri_nback_all_meanmotion_T_0(matches)   Sites(matches,:)  T_.White_T_0(matches) T_.Black_T_0(matches) T_.Hispanic_T_0(matches) T_.HighestEd_T_0(matches) T_.Income_T_0(matches) T_.Male_bin_T_0(matches)  T_.interview_age_T_0(matches)  T_.Ingenia_T_0(matches) T_.Achieva_T_0(matches) T_.Discovery_T_0(matches)  T_.Pfit_T_0(matches) T_.reshist_addr1_adi_wsum_T_0(matches)  T_.pds_ss_T_0(matches)]; % T_.PF10_lavaan(matches)
covariates_y2cont = [T_.depv(matches) T_.tfmri_mid_all_meanmotion_T_2(matches)   T_.tfmri_nback_all_meanmotion_T_2(matches)   Sites(matches,:)  T_.White_T_0(matches) T_.Black_T_0(matches) T_.Hispanic_T_0(matches) T_.HighestEd_T_0(matches) T_.Income_T_0(matches) T_.Male_bin_T_0(matches)  T_.interview_age_T_2(matches)   T_.Ingenia_T_0(matches) T_.Achieva_T_0(matches) T_.Discovery_T_0(matches)  T_.Pfit_T_0(matches)  T_.reshist_addr1_adi_wsum_T_0(matches) T_.pds_ss_T_2(matches)];
covariates_y0pm = [T_.depv(thplus)   T_.tfmri_mid_all_meanmotion_T_0(thplus)   T_.tfmri_nback_all_meanmotion_T_0(thplus)   Sites(thplus,:)  T_.White_T_0(thplus) T_.Black_T_0(thplus) T_.Hispanic_T_0(thplus) T_.HighestEd_T_0(thplus) T_.Income_T_0(thplus) T_.Male_bin_T_0(thplus)  T_.interview_age_T_0(thplus)  T_.Ingenia_T_0(thplus) T_.Achieva_T_0(thplus) T_.Discovery_T_0(thplus)  T_.Pfit_T_0(thplus)  T_.reshist_addr1_adi_wsum_T_0(thplus)  T_.pds_ss_T_0(thplus)];
covariates_y2pm = [T_.depv(thplus)  T_.tfmri_mid_all_meanmotion_T_2(thplus)   T_.tfmri_nback_all_meanmotion_T_2(thplus)   Sites(thplus,:)  T_.White_T_0(thplus) T_.Black_T_0(thplus) T_.Hispanic_T_0(thplus) T_.HighestEd_T_0(thplus) T_.Income_T_0(thplus) T_.Male_bin_T_0(thplus)  T_.interview_age_T_2(thplus)   T_.Ingenia_T_0(thplus) T_.Achieva_T_0(thplus) T_.Discovery_T_0(thplus)  T_.Pfit_T_0(thplus)  T_.reshist_addr1_adi_wsum_T_0(thplus) T_.pds_ss_T_2(thplus)];



covariates = [covariates_y0cont; covariates_y2cont; covariates_y0pm; covariates_y2pm];
% 24,25,26, are race/ethnic
% exclude_race_as_cov = [1:23,27:37];
exclude_race_as_cov = [1:size(covariates,2)];
[a1 ] = partialcorr([result.usc(:,lv),adjy0cont],covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a2 ] = partialcorr([result.usc(:,lv),adjy2cont], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a3 ] = partialcorr([result.usc(:,lv),adjy0pm], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a4 ] = partialcorr([result.usc(:,lv),adjy2pm], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');

a1s =[]; a2s =[]; a3s =[]; a4s =[]; b1s=[];b2s=[];b3s=[]; b4s=[];
dats = [result.datamat_lst{1, 1}; result.datamat_lst{2, 1}];

nboot = 5000;
for jj=1:nboot
   rand_inds1 = randi([1,lg1],lg1,1); rand_inds11 = randi([1,lg1],lg1,1); 
   rand_inds2 = randi([1,lg2],lg2,1); rand_inds22 = randi([1,lg2],lg2,1);
   inds = [rand_inds1; lg1+ rand_inds11; 2*lg1+ rand_inds2; 2*lg1+lg2+ rand_inds22];
   a11 = partialcorr([result.usc(inds,lv),adjy0cont(inds)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
   a1s = [a1s; a11(1,2)];
   a22 = partialcorr([result.usc(inds,lv),adjy2cont(inds)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
a2s = [a2s; a22(1,2)];
a33 = partialcorr([result.usc(inds,lv),adjy0pm(inds)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
a3s = [a3s; a33(1,2)];
a44 = partialcorr([result.usc(inds,lv),adjy2pm(inds)] , covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
a4s = [a4s; a44(1,2)];


   b11 = partialcorr([result.usc(inds,lv), dats(inds,1)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
   b1s = [b1s; b11(1,2)];
   b22 = partialcorr([result.usc(inds,lv),dats(inds,2)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
b2s = [b2s; b22(1,2)];
b33 = partialcorr([result.usc(inds,lv),dats(inds,3)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
b3s = [b3s; b33(1,2)];

b44 = partialcorr([result.usc(inds,lv),dats(inds,4)], covariates(inds,exclude_race_as_cov), 'Rows', 'pairwise');
b4s = [b4s; b44(1,2)];

end

   
 gt_contrast2 = (-1).*[a1s, a2s, a3s, a4s];
 gt_contrast2b = (-1).*[b1s, b2s, b3s b4s];

 p1 = 1 - length(find(  (gt_contrast2(:,2)-gt_contrast2(:,1)) > (mean(gt_contrast2(:,4)) - mean(gt_contrast2(:,3)) ) ) )/nboot
p2 = 1 - (length(find(  (gt_contrast2(:,4)>gt_contrast2(:,3))  ) )/nboot);
p3 = 1 - (length(find(  (gt_contrast2(:,2)>gt_contrast2(:,1))  ) )/nboot);
 p4 = 1 - length(find(  (gt_contrast2(:,2)+gt_contrast2(:,1)) > (mean(gt_contrast2(:,4) + (gt_contrast2(:,3)))  )))/nboot
ef = mean(gt_contrast2(:,2))-mean(gt_contrast2(:,1))

 addpath(genpath('~\export_fig-master'));
pmcol = [.97 .57 .31];
f1 = figure;
hold on
xss = [1,1.8,3.2,4];
sqlist = [1,2,3,4]; % order of grouping columns
for sq=1:4
if (sq==1 | sq==2); col = [.45 .75 0.75];  % color for columns 1 and 2 (HC)
else
    col = pmcol;
end
k = sqlist(sq);
   mm = mean(gt_contrast2(:,k));
  negs = mm - prctile(gt_contrast2(:,k),5);
   pos = prctile(gt_contrast2(:,k),95) - mm;
   bar(xss(sq),mean(gt_contrast2(:,k)),'FaceColor',col,'EdgeColor',col);
   errorbar(xss(sq),mean(gt_contrast2(:,k)),negs,pos,'.k');
end
set(gca,'XTick',xss,'XTickLabels',{'Baseline (9-10 y)','Follow-up (11-12 y)','Baseline (9-10 y)','Follow-up (11-12 y)'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS group-by-time loading (partial corr)'); ylim([-.2,.22]);
f1.Position = [100 100 500 900];
hold off
export_fig -m5 -transparent thplus_pls_Aff_lv1_gtbars_partial_nop.jpg

f2 = figure;
hold on
locs = [.8,1.4,2.6,3.2];
for k=1:4

    col = [.74 .65 .73];

 mmb = mean(gt_contrast2b(:,k));
 
   negsb = mmb - prctile(gt_contrast2b(:,k),5);
  
   posb = prctile(gt_contrast2b(:,k),95) - mmb;
    bar(locs(k),mean(gt_contrast2b(:,k)),'FaceColor',col,'EdgeColor',col,'BarWidth',.45);
     errorbar(locs(k),mean(gt_contrast2b(:,k)),negsb,posb,'.k');
end
set(gca,'XTick',locs,'XTickLabels',{'NAc (rew vs nt)','Caud (rew vs nt)','Amyg (neg vs nt face)','Insula (neg vs nt face)'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS neuroaffective loading (partial corr)'); ylim([-1,1]); %axis square
f2.Position = [100 100 500 900];
hold off

export_fig -m5 -transparent thplus_pls_Aff_lv1_ncbars_partial_nop.jpg

%%
%############################## Figure 1 #######################
% requires:
% Jonas (2026). Violin Plots for plotting multiple distributions (distributionPlot.m) (https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m), MATLAB Central File Exchange. Retrieved June 25, 2026.
% and  Darik (2026). venn (https://www.mathworks.com/matlabcentral/fileexchange/22282-venn), MATLAB Central File Exchange. Retrieved June 5, 2026.
clear all

addpath(genpath('~\venn')) %

load('~\usable_subs_cog.mat')
T_c = usable_subs;

load('~s/usable_subs_aff.mat')
T_af = usable_subs;

addpath(genpath('~\distributionPlot'))
dist_stats_adj = [    T_c.mat_score, T_c.AFC_y2  ...
   nanmean([T_c.NBK T_c.NIH5./120],2), nanmean([T_c.NBK_y2 T_c.NIH_y2./120],2),...
   T_c.smri_thick_cdk_mean , T_c.thick_y2];
dist_aff = [    T_af.mid_NAc_T_0,T_af.mid_NAc_T_2  ...
   T_af.mid_Caud_T_0,T_af.mid_Caud_T_2 ,...
   T_af.ng_nt_face_amyg_T_0,T_af.ng_nt_face_amyg_T_2  ...
   T_af.ng_nt_face_insu_T_0,T_af.ng_nt_face_insu_T_2 ];


titles ={'Cort FC maturation','Cog Task Performance','Cortical GM Thickness',...
    'Reward anticipation (NAc)','Reward anticipation (Caud)', 'Fearful face (Amyg)','Fearful face (Insula)'};
xnames ={'Y0','Y2'};
figure;
subplot(4,4,[1])

[H1, S1] = venn([638+382-198,382], 382-198)
for i = 1:length(S1.ZoneCentroid)
text(S1.ZoneCentroid(i,1), S1.ZoneCentroid(i,2), num2str(S1.ZonePop(i)), 'FontSize',14,'HorizontalAlignment', 'center');
end
xlim([-35,45])

subplot(4,4,[5])
[H2, S2] = venn([721+464-232,464], 464-232)
for i = 1:length(S2.ZoneCentroid)
text(S2.ZoneCentroid(i,1), S2.ZoneCentroid(i,2), num2str(S2.ZonePop(i)), 'FontSize',14,'HorizontalAlignment', 'center');
end
xlim([-35,45]) % threat color =  [1 .4 .2] alpha=.8, deprv color =[1 .7 .9]

subplot(4,4,[2,6])
distributionPlot({dist_stats_adj(:,1)},...    
    'color',{[0.90,0.65,0.73]},'xyOri','normal','histOri','right','showMM',0,'xNames',{'Y0','Y2'})

distributionPlot({dist_stats_adj(:,2)},...
    'color',{[0.70,0.75,0.93]},'xyOri','normal','histOri','left','showMM',0,'xNames','Y2')

line([0.5 ,1],[nanmedian(dist_stats_adj(:,1)) ,nanmedian(dist_stats_adj(:,1))],'Color','red');
line([1 ,1.5],[nanmedian(dist_stats_adj(:,2)) ,nanmedian(dist_stats_adj(:,2))],'Color','blue');
ylabel('AFC score'); ylim([-.4,.65])
set(gca,'FontSize',16,'XTickLabelRotation',45,'XTick',[.5,1.5],'XTickLabel',{'Y0','Y2'})
title(titles{1});

subplot(4,4,[3,7])
distributionPlot({dist_stats_adj(:,3)},...    
    'color',{[0.90,0.65,0.73]},'xyOri','normal','histOri','right','showMM',0,'xNames','Y0')

distributionPlot({dist_stats_adj(:,4)},...
    'color',{[0.70,0.75,0.93]},'xyOri','normal','histOri','left','showMM',0,'xNames','Y2')
hold on

ylabel('Mean Accuracy');  ylim([.3,1])
set(gca,'FontSize',16,'XTickLabelRotation',45,'XTick',[.5,1.5],'XTickLabel',{'Y0','Y2'})
line([0.5 ,1],[nanmedian(dist_stats_adj(:,3)) ,nanmedian(dist_stats_adj(:,3))],'Color','red');
line([1,1.5],[nanmedian(dist_stats_adj(:,4)) ,nanmedian(dist_stats_adj(:,4))],'Color','blue');
title(titles{2});

subplot(4,4,[4,8])
distributionPlot({dist_stats_adj(:,5)},...    
    'color',{[0.90,0.65,0.73]},'xyOri','normal','histOri','right','showMM',0,'xNames','Y0')


distributionPlot({dist_stats_adj(:,6)},...
    'color',{[0.70,0.75,0.93]},'xyOri','normal','histOri','left','showMM',0,'xNames','Y2')

ylabel('Thickness (mm)');  ylim([2,3.1])
set(gca,'FontSize',16,'XTickLabelRotation',45,'XTick',[.5,1.5],'XTickLabel',{'Y0','Y2'})
hold on
line([0.5 ,1],[nanmedian(dist_stats_adj(:,5)) ,nanmedian(dist_stats_adj(:,5))],'Color','red');
line([1 ,1.5],[nanmedian(dist_stats_adj(:,6)) ,nanmedian(dist_stats_adj(:,6))],'Color','blue');
title(titles{3});


subplot(4,4,[9,13])
distributionPlot({dist_aff(:,1)},...    
    'color',{[0.90,0.65,0.73]},'xyOri','normal','histOri','right','showMM',0,'xNames','Y0')


distributionPlot({dist_aff(:,2)},...
    'color',{[0.70,0.75,0.93]},'xyOri','normal','histOri','left','showMM',0,'xNames','Y2')

ylabel('High reward vs. neutral contrast (/beta)');  ylim([-2.2,2.2])
set(gca,'FontSize',16,'XTickLabelRotation',45,'XTick',[.5,1.5],'XTickLabel',{'Baseline (Y0)','Follow-up (Y2)'})
hold on
line([0.5 ,1],[nanmedian(dist_aff(:,1)) ,nanmedian(dist_aff(:,1))],'Color','red');
line([1 ,1.5],[nanmedian(dist_aff(:,2)) ,nanmedian(dist_aff(:,2))],'Color','blue');
title(titles{4});

subplot(4,4,[10,14])
distributionPlot({dist_aff(:,3)},...    
    'color',{[0.90,0.65,0.73]},'xyOri','normal','histOri','right','showMM',0,'xNames','Y0')


distributionPlot({dist_aff(:,4)},...
    'color',{[0.70,0.75,0.93]},'xyOri','normal','histOri','left','showMM',0,'xNames','Y2')

ylabel('High reward vs. neutral contrast (/beta)');  ylim([-2.2,2.2])
set(gca,'FontSize',16,'XTickLabelRotation',45,'XTick',[.5,1.5],'XTickLabel',{'Baseline (9-10 y)','Follow-up (11-12 y)'})
hold on
line([0.5 ,1],[nanmedian(dist_aff(:,3)) ,nanmedian(dist_aff(:,3))],'Color','red');
line([1 ,1.5],[nanmedian(dist_aff(:,4)) ,nanmedian(dist_aff(:,4))],'Color','blue');
title(titles{5});

subplot(4,4,[11,15])
distributionPlot({dist_aff(:,5)},...    
    'color',{[0.90,0.65,0.73]},'xyOri','normal','histOri','right','showMM',0,'xNames','Y0')


distributionPlot({dist_aff(:,6)},...
    'color',{[0.70,0.75,0.93]},'xyOri','normal','histOri','left','showMM',0,'xNames','Y2')

ylabel('Fearful vs. neutral face contrast (/beta)');  ylim([-2.2,2.2])
set(gca,'FontSize',16,'XTickLabelRotation',45,'XTick',[.5,1.5],'XTickLabel',{'Baseline (9-10 y)','Follow-up (11-12 y)'})
hold on
line([0.5 ,1],[nanmedian(dist_aff(:,5)) ,nanmedian(dist_aff(:,5))],'Color','red');
line([1 ,1.5],[nanmedian(dist_aff(:,6)) ,nanmedian(dist_aff(:,6))],'Color','blue');
title(titles{6});


subplot(4,4,[12,16])
distributionPlot({dist_aff(:,7)},...    
    'color',{[0.90,0.65,0.73]},'xyOri','normal','histOri','right','showMM',0,'xNames','Y0')


distributionPlot({dist_aff(:,8)},...
    'color',{[0.70,0.75,0.93]},'xyOri','normal','histOri','left','showMM',0,'xNames','Y2')

ylabel('Fearful vs. neutral face contrast (/beta)');  ylim([-2.2,2.2])
set(gca,'FontSize',16,'XTickLabelRotation',45,'XTick',[.5,1.5],'XTickLabel',{'Baseline (9-10 y)','Follow-up (11-12 y)'})
hold on
line([0.5 ,1],[nanmedian(dist_aff(:,7)) ,nanmedian(dist_aff(:,7))],'Color','red');
line([1 ,1.5],[nanmedian(dist_aff(:,8)) ,nanmedian(dist_aff(:,8))],'Color','blue');
title(titles{7});

addpath(genpath('Z:\MatlabGraphics\export_fig-master'))
export_fig -m5 -transparent Figure1.jpg