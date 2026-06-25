% # Scripts for "Cognitive and affective neurodevelopment in youth exposed to deprivation and threat"
% # Contact Omid Kardan omidk@med.umich.edu
%%% This script is for running PLS based on neurocognitve measures
% Requirements: 
% 1) Uses csv file generated from the R script Kardan_dprv_thrt_variables_compilation_cog.R
% % 2) Uses the PLS scripts from https://www.rotman-baycrest.on.ca/index.php?section=345
% Download plscmd and plsgui and place in Pls_folder
% 3) Uses export_fig from https://github.com/altmany/export_fig
%%  ######################################### Deprivation only #####################################################################################
%%  %% Create matched control group for the depriv_only participants
clear all
T_0 = read_abcd_beh('~\df_filled_beh_Y0_more.csv');
T_2 =read_abcd_beh('~\df_filled_beh_Y2_more.csv');


load('~/with_2b_2024.mat') % list of partiicpants with sufficient rsfMRI data at two timepoints
s_table = readtable("~/df_filled_demog_withcog_Y0_more.csv");
sub_11868 = table(s_table.subid_x);
usable_subs0 = innerjoin(sub_11868(with_2b_2024,:),T_0,"LeftKeys","Var1","RightKeys","subids");
tfu = table(T_2.subids, T_2.NBK,  T_2.NIH5,  T_2.smri_thick_cdk_mean,'VariableNames',{'subids','NBK_y2' ,'NIH_y2' ,'thick_y2'});
usable_subs = innerjoin(usable_subs0,tfu,"LeftKeys","Var1","RightKeys","subids");
usable_subs.cog = nanmean([usable_subs.NBK usable_subs.NIH5./120],2) + nanmean([usable_subs.NBK_y2 usable_subs.NIH_y2./120],2);
usable_subs.neur = usable_subs.smri_thick_cdk_mean + usable_subs.thick_y2; % save as usable_subs_cog.mat

% this section runs iteratively to make the maatched group which is then
% saved as num_matches_subs_depv_only.mat
depv_high1 = find(usable_subs.cumul_ace == 1 & ~isnan(usable_subs.cog) & ~isnan(usable_subs.neur) ...
    ); % 1 is depriv only

% load(['~/num_matches_subs_depv_only.mat']); %uncomment this after the first
% time running the script which generates the list

depv_high = depv_high1;
% depv_high =depv_high1(num_matches_subs_depv_only >1); %uncomment this after the first
% time running the script which generates the list

depv_high_n = length(depv_high);
depv_low = find(usable_subs.cumul_ace == 0 & ~isnan(usable_subs.cog) & ~isnan(usable_subs.neur)); % 

usable_subs.raceth = usable_subs.White;
usable_subs.raceth(usable_subs.Black == 1) = 2;
usable_subs.raceth(usable_subs.Hispanic == 1) = 3;
usable_subs.raceth(usable_subs.Asian == 1) = 4;
usable_subs.raceth(usable_subs.Other == 1) = 5;
usable_subs.interview_age(isnan(usable_subs.interview_age)) = nanmean(usable_subs.interview_age);
usable_subs.agey = round(usable_subs.interview_age/6);
usable_subs.HighestEd(isnan(usable_subs.HighestEd)) = nanmean(usable_subs.HighestEd);
usable_subs.edu = round(usable_subs.HighestEd/5);

income_matches = cell(depv_high_n,1);
age_matches = cell(depv_high_n,1);
sex_matches = cell(depv_high_n,1);
raceth_matches = cell(depv_high_n,1);
edu_matches = cell(depv_high_n,1);

for k =1:depv_high_n
    
    if ~isnan(usable_subs.Income(depv_high(k)) )
        a = find(round(usable_subs.Income) == round(usable_subs.Income(depv_high(k))));
    else a = depv_low;
    end
        income_matches{k} = intersect(depv_low,a);
        
    if ~isnan(usable_subs.agey(depv_high(k)))
        a = find(usable_subs.agey == usable_subs.agey(depv_high(k)));
    else a = depv_low;
    end
        age_matches{k} = intersect(depv_low,a);
    
    if ~isnan(usable_subs.Male_bin(depv_high(k)))
        a = find(usable_subs.Male_bin == usable_subs.Male_bin(depv_high(k)));
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

num_matches_subs_depv_only =[];
for l=1:693
num_matches_subs_depv_only = [num_matches_subs_depv_only; length(cont_mj_bare{l})];
end
% save list as num_matches_subs_depv_only.mat.mat

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
        
matches_534 = unique(matches); 
% check if cohen's d is samll for group differences in each category
T_ = usable_subs;
cohd = (nanmean(T_.age(depv_high))-nanmean(T_.age(setdiff(1:3645,depv_high))) )/ sqrt(nanvar(T_.age(depv_high))+nanvar(T_.age(setdiff(1:3645,depv_high))) )
cohd = (nanmean(T_.age(depv_high))-nanmean(T_.age(depv_low)) )/ sqrt(nanvar(T_.age(depv_high))+nanvar(T_.age(depv_low)) )
cohd = (nanmean(T_.age(depv_high))-nanmean(T_.age(matches_534)) )/ sqrt(nanvar(T_.age(depv_high))+nanvar(T_.age(matches_534)) )

cohd = (nanmean(T_.Income(depv_high))-nanmean(T_.Income(setdiff(1:3645,depv_high))) )/ sqrt(nanvar(T_0.Income(depv_high))+nanvar(T_.Income(setdiff(1:3645,depv_high))) )
cohd = (nanmean(T_.Income(depv_high))-nanmean(T_.Income(depv_low)) )/ sqrt(nanvar(T_.Income(depv_high))+nanvar(T_.Income(depv_low)) )
cohd = (nanmean(T_.Income(depv_high))-nanmean(T_.Income(matches_534)) )/ sqrt(nanvar(T_.Income(depv_high))+nanvar(T_.Income(matches_534)) )


cohd = (nanmean(T_.HighestEd(depv_high))-nanmean(T_.HighestEd(setdiff(1:3645,depv_high))) )/ sqrt(nanvar(T_.HighestEd(depv_high))+nanvar(T_.HighestEd(setdiff(1:3645,depv_high))) )
cohd = (nanmean(T_.HighestEd(depv_high))-nanmean(T_.HighestEd(depv_low)) )/ sqrt(nanvar(T_.HighestEd(depv_high))+nanvar(T_.HighestEd(depv_low)) )
cohd = (nanmean(T_.HighestEd(depv_high))-nanmean(T_.HighestEd(matches_534)) )/ sqrt(nanvar(T_.HighestEd(depv_high))+nanvar(T_.HighestEd(matches_534)) )

cohd = (nanmean(T_.Black(depv_high))-nanmean(T_.Black(setdiff(1:3645,depv_high))) )/ sqrt(nanvar(T_.Black(depv_high))+nanvar(T_.Black(setdiff(1:3645,depv_high))) )
cohd = (nanmean(T_.Black(depv_high))-nanmean(T_.Black(depv_low)) )/ sqrt(nanvar(T_.Black(depv_high))+nanvar(T_.Black(depv_low)) )
cohd = (nanmean(T_.Black(depv_high))-nanmean(T_.Black(matches_534)) )/ sqrt(nanvar(T_.Black(depv_high))+nanvar(T_.Black(matches_534)) )

cohd = (nanmean(T_.White(depv_high))-nanmean(T_.White(setdiff(1:3645,depv_high))) )/ sqrt(nanvar(T_.White(depv_high))+nanvar(T_.White(setdiff(1:3645,depv_high))) )
cohd = (nanmean(T_.White(depv_high))-nanmean(T_.White(depv_low)) )/ sqrt(nanvar(T_.White(depv_high))+nanvar(T_.White(depv_low)) )
cohd = (nanmean(T_.White(depv_high))-nanmean(T_.White(matches_534)) )/ sqrt(nanvar(T_.White(depv_high))+nanvar(T_.White(matches_534)) )

cohd = (nanmean(T_.Hispanic(depv_high))-nanmean(T_.Hispanic(setdiff(1:3645,depv_high))) )/ sqrt(nanvar(T_.Hispanic(depv_high))+nanvar(T_.Hispanic(setdiff(1:3645,depv_high))) )
cohd = (nanmean(T_.Hispanic(depv_high))-nanmean(T_.Hispanic(depv_low)) )/ sqrt(nanvar(T_.Hispanic(depv_high))+nanvar(T_.Hispanic(depv_low)) )
cohd = (nanmean(T_.Hispanic(depv_high))-nanmean(T_.Hispanic(matches_534)) )/ sqrt(nanvar(T_.Hispanic(depv_high))+nanvar(T_.Hispanic(matches_534)) )

matches_depv_534 = unique(matches); % save these as your matched controls
depv_638 = depv_high; % save these as your exposure group


%% %%%%% PLS
%% neurocognitve latent variable

clear all
T_0 = read_abcd_beh('~\df_filled_beh_Y0_more.csv');
T_2 =read_abcd_beh('~\df_filled_beh_Y2_more.csv');


load('~/usable_subs_cog.mat')
T_ = usable_subs;

load('matches_depv_534.mat')
load('depv_638.mat')


addpath(genpath('~/Pls_folder/Pls'));
rng('default');
groups{1} = [matches_depv_534];
groups{2} = [depv_638];
nTimes = 2;

datamat_lst = cell(length(groups),1); 
% ####*****************####
T_0.cog = nanmean([T_0.NBK T_0.NIH5./120],2);
T_2.cog = nanmean([T_2.NBK T_2.NIH5./120],2);


for g = 1:numel(groups)
    for c = 1:nTimes
        for s = 1:numel(groups{g})
            sub_id = char(T_.Var1(groups{g}(s)));
            
            vec2 =[];
            if c==1
                locb = find(T_0.subids == sub_id);
                if length(locb)==0
                    sub_id
                end
                vec2 = [T_0.mat_score(locb)   T_0.cog(locb)  T_0.smri_thick_cdk_mean(locb)];
                                                    
            end
            
            if c==2
                loct = find(T_2.subids == sub_id);
                if length(loct)==0
                    sub_id
                end
                vec2 = [T_2.mat_score(loct)  T_2.cog(loct)  T_2.smri_thick_cdk_mean(loct)];
                  

            end
            
            datamat_lst{g} = [datamat_lst{g}; vec2];
            s;
            
        end
    end
end

num_subj = [length(groups{1})  length(groups{2})];
num_cond = nTimes;
option.method = 1; %
option.num_boot = 2000;
option.num_perm = 2000;
option.meancentering_type=[2]; % 
result = pls_analysis(datamat_lst, num_subj, num_cond, option);
result.datamat_lst = datamat_lst;
% save as PLS_result_534lowdepv_638hdepv_y0_y2.mat
%% plotting the results (figure 2 panel A and table 1 top half)
load('PLS_result_534lowdepv_638hdepv_y0_y2.mat');  
ps = result.perm_result.sprob
load('~/usable_subs_cog.mat')
T_q = usable_subs;

T_0 = read_abcd_beh('~\df_filled_beh_Y0_more.csv');
T_2 =read_abcd_beh('~\df_filled_beh_Y2_more.csv');
tfu = table(T_2.subids, T_2.interview_age,T_2.nruns,T_2.pds_ss,'VariableNames',{'subids','age_Y2','nruns_Y2','pds_ss_Y2'});
T_ = innerjoin(T_q,tfu,"LeftKeys","Var1","RightKeys","subids");

load('matches_depv_534.mat')
load('depv_638.mat')
matches = matches_depv_534;
thrthigh = depv_638;

% *** for Table 1 top half in the manuscript
thrt_low = find(usable_subs.cumul_ace == 0 & ~isnan(usable_subs.cog) & ~isnan(usable_subs.neur));
T_.AsiOth = T_.Asian + T_.Other;

[nanmean(T_.HighestEd(thrthigh)) sqrt(nanvar(T_.HighestEd(thrthigh)))]
[nanmean(T_.HighestEd(thrt_low)) sqrt(nanvar(T_.HighestEd(thrt_low)))]
[nanmean(T_.HighestEd(matches)) sqrt(nanvar(T_.HighestEd(matches)))]
cohd = (nanmean(T_.HighestEd(thrthigh))-nanmean(T_.HighestEd(thrt_low)) )/ sqrt(nanvar(T_.HighestEd(thrthigh))+nanvar(T_.HighestEd(thrt_low)) )
cohd = (nanmean(T_.HighestEd(thrthigh))-nanmean(T_.HighestEd(matches)) )/ sqrt(nanvar(T_.HighestEd(thrthigh))+nanvar(T_.HighestEd(matches)) )

      % levels= 1:10, labels = c("5000", "8500", "14000", "20500", "30000","42500", 7:"62500", 8:"87500", "150000", "200000") )
nanmean(T_.Income(thrthigh)); min(T_.Income(thrthigh)); max(T_.Income(thrthigh));
42500 + (6.8124-6)*(62500 - 42500)
nanmean(T_.Income(thrt_low)) ;
87500 + (8.2134-8)*(150000 - 87500)
nanmean(T_.Income(matches)); 
62500 + (7.17-7)*(87500 - 62500)
cohd = (nanmean(T_.Income(thrthigh))-nanmean(T_.Income(thrt_low)) )/ sqrt(nanvar(T_.Income(thrthigh))+nanvar(T_.Income(thrt_low)) )
cohd = (nanmean(T_.Income(thrthigh))-nanmean(T_.Income(matches)) )/ sqrt(nanvar(T_.Income(thrthigh))+nanvar(T_.Income(matches)) )

[length(find(T_.White(thrthigh))) length(find(T_.White(thrthigh)))/length(thrthigh)]
[length(find(T_.White(thrt_low))) length(find(T_.White(thrt_low)))/length(thrt_low)]
[length(find(T_.White(matches))) length(find(T_.White(matches)))/length(matches)]
cohd = (nanmean(T_.White(thrthigh))-nanmean(T_.White(thrt_low)) )/ sqrt(nanvar(T_.White(thrthigh))+nanvar(T_.White(thrt_low)) )
cohd = (nanmean(T_.White(thrthigh))-nanmean(T_.White(matches)) )/ sqrt(nanvar(T_.White(thrthigh))+nanvar(T_.White(matches)) )


[length(find(T_.Hispanic(thrthigh))) length(find(T_.Hispanic(thrthigh)))/length(thrthigh)]
[length(find(T_.Hispanic(thrt_low))) length(find(T_.Hispanic(thrt_low)))/length(thrt_low)]
[length(find(T_.Hispanic(matches))) length(find(T_.Hispanic(matches)))/length(matches)]
cohd = (nanmean(T_.Hispanic(thrthigh))-nanmean(T_.Hispanic(thrt_low)) )/ sqrt(nanvar(T_.Hispanic(thrthigh))+nanvar(T_.Hispanic(thrt_low)) )
cohd = (nanmean(T_.Hispanic(thrthigh))-nanmean(T_.Hispanic(matches)) )/ sqrt(nanvar(T_.Hispanic(thrthigh))+nanvar(T_.Hispanic(matches)) )


[length(find(T_.Black(thrthigh))) length(find(T_.Black(thrthigh)))/length(thrthigh)]
[length(find(T_.Black(thrt_low))) length(find(T_.Black(thrt_low)))/length(thrt_low)]
[length(find(T_.Black(matches))) length(find(T_.Black(matches)))/length(matches)]
cohd = (nanmean(T_.Black(thrthigh))-nanmean(T_.Black(thrt_low)) )/ sqrt(nanvar(T_.Black(thrthigh))+nanvar(T_.Black(thrt_low)) )
cohd = (nanmean(T_.Black(thrthigh))-nanmean(T_.Black(matches)) )/ sqrt(nanvar(T_.Black(thrthigh))+nanvar(T_.Black(matches)) )


[length(find(T_.AsiOth(thrthigh))) length(find(T_.AsiOth(thrthigh)))/length(thrthigh)]
[length(find(T_.AsiOth(thrt_low))) length(find(T_.AsiOth(thrt_low)))/length(thrt_low)]
[length(find(T_.AsiOth(matches))) length(find(T_.AsiOth(matches)))/length(matches)]
cohd = (nanmean(T_.AsiOth(thrthigh))-nanmean(T_.AsiOth(thrt_low)) )/ sqrt(nanvar(T_.AsiOth(thrthigh))+nanvar(T_.AsiOth(thrt_low)) )
cohd = (nanmean(T_.AsiOth(thrthigh))-nanmean(T_.AsiOth(matches)) )/ sqrt(nanvar(T_.AsiOth(thrthigh))+nanvar(T_.AsiOth(matches)) )

[length(find(~T_.Male_bin(thrthigh))) length(find(~T_.Male_bin(thrthigh)))/length(thrthigh)]
[length(find(~T_.Male_bin(thrt_low))) length(find(~T_.Male_bin(thrt_low)))/length(thrt_low)]
[length(find(~T_.Male_bin(matches))) length(find(~T_.Male_bin(matches)))/length(matches)]
cohd = (nanmean(T_.Male_bin(thrthigh))-nanmean(T_.Male_bin(thrt_low)) )/ sqrt(nanvar(T_.Male_bin(thrthigh))+nanvar(T_.Male_bin(thrt_low)) )
cohd = (nanmean(T_.Male_bin(thrthigh))-nanmean(T_.Male_bin(matches)) )/ sqrt(nanvar(T_.Male_bin(thrthigh))+nanvar(T_.Male_bin(matches)) )


[nanmean(T_.age(thrthigh)) sqrt(nanvar(T_.age(thrthigh)))]
[nanmean(T_.age(thrt_low)) sqrt(nanvar(T_.age(thrt_low)))]
[nanmean(T_.age(matches)) sqrt(nanvar(T_.age(matches)))]
cohd = (nanmean(T_.age(thrthigh))-nanmean(T_.age(thrt_low)) )/ sqrt(nanvar(T_.age(thrthigh))+nanvar(T_.age(thrt_low)) )
cohd = (nanmean(T_.age(thrthigh))-nanmean(T_.age(matches)) )/ sqrt(nanvar(T_.age(thrthigh))+nanvar(T_.age(matches)) )

[nanmean(T_.age_Y2(thrthigh)) sqrt(nanvar(T_.age_Y2(thrthigh)))]/12
[nanmean(T_.age_Y2(thrt_low)) sqrt(nanvar(T_.age_Y2(thrt_low)))]/12
[nanmean(T_.age_Y2(matches)) sqrt(nanvar(T_.age_Y2(matches)))]/12
cohd = (nanmean(T_.age_Y2(thrthigh))-nanmean(T_.age_Y2(thrt_low)) )/ sqrt(nanvar(T_.age_Y2(thrthigh))+nanvar(T_.age_Y2(thrt_low)) )
cohd = (nanmean(T_.age_Y2(thrthigh))-nanmean(T_.age_Y2(matches)) )/ sqrt(nanvar(T_.age_Y2(thrthigh))+nanvar(T_.age_Y2(matches)) )

% ***

Sites = zeros(3645,22);
for K=1:22
    Sites((T_.site_id_l == K),K) = 1;
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


covariates_y0cont = [T_.pFD_y0(matches)   Sites(matches,:)   T_.White(matches) T_.Black(matches) T_.Hispanic(matches) T_.HighestEd(matches) T_.Income(matches) T_.Male_bin(matches)  T_.interview_age(matches)  T_.Ingenia(matches) T_.Achieva(matches) T_.Discovery(matches)  T_.Pfit(matches) T_.nruns(matches) T_.reshist_addr1_adi_wsum(matches)  T_.pds_ss(matches) ]; % T_.PF10_lavaan(matches)
covariates_y2cont = [T_.pFD_y2(matches)   Sites(matches,:)  T_.White(matches) T_.Black(matches) T_.Hispanic(matches) T_.HighestEd(matches) T_.Income(matches) T_.Male_bin(matches)  T_.age_Y2(matches)   T_.Ingenia(matches) T_.Achieva(matches) T_.Discovery(matches)  T_.Pfit(matches) T_.nruns_Y2(matches) T_.reshist_addr1_adi_wsum(matches)  T_.pds_ss_Y2(matches) ]; % T_.PF10_lavaan(matches)
covariates_y0pm = [T_.pFD_y0(thrthigh)   Sites(thrthigh,:)  T_.White(thrthigh) T_.Black(thrthigh) T_.Hispanic(thrthigh) T_.HighestEd(thrthigh) T_.Income(thrthigh) T_.Male_bin(thrthigh)  T_.interview_age(thrthigh)   T_.Ingenia(thrthigh) T_.Achieva(thrthigh) T_.Discovery(thrthigh)  T_.Pfit(thrthigh) T_.nruns(thrthigh) T_.reshist_addr1_adi_wsum(thrthigh)  T_.pds_ss(thrthigh) ]; % T_.PF10_lavaan(thrthigh)
covariates_y2pm = [T_.pFD_y2(thrthigh)   Sites(thrthigh,:)  T_.White(thrthigh) T_.Black(thrthigh) T_.Hispanic(thrthigh) T_.HighestEd(thrthigh) T_.Income(thrthigh) T_.Male_bin(thrthigh)  T_.age_Y2(thrthigh)  T_.Ingenia(thrthigh) T_.Achieva(thrthigh) T_.Discovery(thrthigh)  T_.Pfit(thrthigh) T_.nruns_Y2(thrthigh) T_.reshist_addr1_adi_wsum(thrthigh)  T_.pds_ss_Y2(thrthigh) ]; %T_.PF10_lavaan(thrthigh)

covariates = [covariates_y0cont; covariates_y2cont; covariates_y0pm; covariates_y2pm];
% 24,25,26, are race/ethnic
% exclude_race_as_cov = [1:23,27:37];
exclude_race_as_cov = [1:size(covariates,2)];
[a1 ] = partialcorr([result.usc(:,lv),adjy0cont],covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a2 ] = partialcorr([result.usc(:,lv),adjy2cont], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a3 ] = partialcorr([result.usc(:,lv),adjy0pm], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');
[a4 ] = partialcorr([result.usc(:,lv),adjy2pm], covariates(:,exclude_race_as_cov), 'Rows', 'pairwise');

a1s =[]; a2s =[]; a3s =[]; a4s =[]; b1s=[];b2s=[];b3s=[]; b1sadj=[];b2sadj=[];b3sadj=[];
dats = [result.datamat_lst{1, 1}; result.datamat_lst{2, 1}];

nboot = 2000;
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


end

 gt_contrast2 = (1).*[a1s, a2s, a3s, a4s];
 gt_contrast2b = (1).*[b1s, b2s, b3s];

 p1 = 1 - length(find(  (gt_contrast2(:,2)-gt_contrast2(:,1)) > (mean(gt_contrast2(:,4)) - mean(gt_contrast2(:,3)) ) ) )/nboot
p2 = 1 - (length(find(  (gt_contrast2(:,4)>gt_contrast2(:,3))  ) )/nboot);
p3 = 1 - (length(find(  (gt_contrast2(:,2)>gt_contrast2(:,1))  ) )/nboot);
 p4 = 1 - length(find(  (gt_contrast2(:,2)+gt_contrast2(:,1)) > (mean(gt_contrast2(:,4) + (gt_contrast2(:,3)))  )))/nboot
 ef = mean(gt_contrast2(:,2))-mean(gt_contrast2(:,1))

ef = (mean(gt_contrast2(:,2))+mean(gt_contrast2(:,1)) )/2 - (mean(gt_contrast2(:,3))+mean(gt_contrast2(:,4)) )/2

 addpath(genpath('Z:\MatlabGraphics\export_fig-master'));
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
export_fig -m5 -transparent depv_pls_Cog_lv1_gtbars_nop.jpg

f2 = figure;
hold on
for k=1:3
    col = [.84 .75 .13];
 mmb = mean(gt_contrast2b(:,k));
     negsb = mmb - prctile(gt_contrast2b(:,k),5);
      posb = prctile(gt_contrast2b(:,k),95) - mmb;
    bar(k*.8,mean(gt_contrast2b(:,k)),'FaceColor',col,'EdgeColor',col,'BarWidth',.45);
      errorbar(k*.8,mean(gt_contrast2b(:,k)),negsb,posb,'.k');
end
set(gca,'XTick',[.8,1.6,2.4],'XTickLabels',{'Cort FC maturation','Cog task perf','Cort GM thickness'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS neurocognitve loading (partial corr)'); ylim([-1,1]); %axis square
f2.Position = [100 100 500 900];
hold off  
export_fig -m5 -transparent depv_pls_Cog_lv1_ncbars_nuis_nop.jpg

%%  ######################################### THREAT only #####################################################################################
%%  %% Create matched control group for the threat_only participants
clear all
T_0 = read_abcd_beh('~\df_filled_beh_Y0_more.csv');
T_2 =read_abcd_beh('~\df_filled_beh_Y2_more.csv');


load('~/with_2b_2024.mat')
s_table = readtable("~/df_filled_demog_withcog_Y0_more.csv");
sub_11868 = table(s_table.subid_x);
usable_subs0 = innerjoin(sub_11868(with_2b_2024,:),T_0,"LeftKeys","Var1","RightKeys","subids");
tfu = table(T_2.subids, T_2.NBK,  T_2.NIH5,  T_2.smri_thick_cdk_mean,T_2.mat_score, 'VariableNames',{'subids','NBK_y2' ,'NIH_y2' ,'thick_y2','AFC_y2'});
usable_subs = innerjoin(usable_subs0,tfu,"LeftKeys","Var1","RightKeys","subids");
usable_subs.cog = nanmean([usable_subs.NBK usable_subs.NIH5./120],2) + nanmean([usable_subs.NBK_y2 usable_subs.NIH_y2./120],2);
usable_subs.neur = usable_subs.smri_thick_cdk_mean + usable_subs.thick_y2;  % save as usable_subs_cog.mat

% this section runs iteratively to make the matched group which is then
% saved as num_matches_subs_thrt_only.mat
thrt_high1 = find(usable_subs.cumul_ace == 2 & ~isnan(usable_subs.cog) & ~isnan(usable_subs.neur) ...
    ); % 2 is threat only

%load(['~/num_matches_subs_thrt_only.mat']);%uncomment this after the first
% time running the script which generates the list

thrt_high = thrt_high1; % 
%thrt_high =thrt_high1(num_matches_subs_thrt_only >1); %uncomment this after the first
% time running the script which generates the list

thrt_high_n = length(thrt_high);
thrt_low = find(usable_subs.cumul_ace == 0 & ~isnan(usable_subs.cog) & ~isnan(usable_subs.neur)); % 

usable_subs.raceth = usable_subs.White;
usable_subs.raceth(usable_subs.Black == 1) = 2;
usable_subs.raceth(usable_subs.Hispanic == 1) = 3;
usable_subs.raceth(usable_subs.Asian == 1) = 4;
usable_subs.raceth(usable_subs.Other == 1) = 5;
usable_subs.interview_age(isnan(usable_subs.interview_age)) = nanmean(usable_subs.interview_age);
usable_subs.agey = round(usable_subs.interview_age/6);
usable_subs.HighestEd(isnan(usable_subs.HighestEd)) = nanmean(usable_subs.HighestEd);
usable_subs.edu = round(usable_subs.HighestEd/5);

income_matches = cell(thrt_high_n,1);
age_matches = cell(thrt_high_n,1);
sex_matches = cell(thrt_high_n,1);
raceth_matches = cell(thrt_high_n,1);
edu_matches = cell(thrt_high_n,1);

for k =1:thrt_high_n
    
    if ~isnan(usable_subs.Income(thrt_high(k)) )
        a = find(round(usable_subs.Income) == round(usable_subs.Income(thrt_high(k))));
    else a = thrt_low;
    end
        income_matches{k} = intersect(thrt_low,a);
        
    if ~isnan(usable_subs.agey(thrt_high(k)))
        a = find(usable_subs.agey == usable_subs.agey(thrt_high(k)));
    else a = thrt_low;
    end
        age_matches{k} = intersect(thrt_low,a);
    
    if ~isnan(usable_subs.Male_bin(thrt_high(k)))
        a = find(usable_subs.Male_bin == usable_subs.Male_bin(thrt_high(k)));
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

num_matches_subs_thrt_only =[];
for l=1:207
num_matches_subs_thrt_only = [num_matches_subs_thrt_only; length(cont_mj_bare{l})];
end
% save list as num_matches_subs_thrt_only.mat

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
        
matches_195 = unique(matches); 
% check if cohen's d is samll for group differences in each category
T_ = usable_subs;
cohd = (nanmean(T_.age(thrt_high))-nanmean(T_.age(setdiff(1:3645,thrt_high))) )/ sqrt(nanvar(T_.age(thrt_high))+nanvar(T_.age(setdiff(1:3645,thrt_high))) )
cohd = (nanmean(T_.age(thrt_high))-nanmean(T_.age(thrt_low)) )/ sqrt(nanvar(T_.age(thrt_high))+nanvar(T_.age(thrt_low)) )
cohd = (nanmean(T_.age(thrt_high))-nanmean(T_.age(matches_195)) )/ sqrt(nanvar(T_.age(thrt_high))+nanvar(T_.age(matches_195)) )

cohd = (nanmean(T_.Income(thrt_high))-nanmean(T_.Income(setdiff(1:3645,thrt_high))) )/ sqrt(nanvar(T_0.Income(thrt_high))+nanvar(T_.Income(setdiff(1:3645,thrt_high))) )
cohd = (nanmean(T_.Income(thrt_high))-nanmean(T_.Income(thrt_low)) )/ sqrt(nanvar(T_.Income(thrt_high))+nanvar(T_.Income(thrt_low)) )
cohd = (nanmean(T_.Income(thrt_high))-nanmean(T_.Income(matches_195)) )/ sqrt(nanvar(T_.Income(thrt_high))+nanvar(T_.Income(matches_195)) )


cohd = (nanmean(T_.HighestEd(thrt_high))-nanmean(T_.HighestEd(setdiff(1:3645,thrt_high))) )/ sqrt(nanvar(T_.HighestEd(thrt_high))+nanvar(T_.HighestEd(setdiff(1:3645,thrt_high))) )
cohd = (nanmean(T_.HighestEd(thrt_high))-nanmean(T_.HighestEd(thrt_low)) )/ sqrt(nanvar(T_.HighestEd(thrt_high))+nanvar(T_.HighestEd(thrt_low)) )
cohd = (nanmean(T_.HighestEd(thrt_high))-nanmean(T_.HighestEd(matches_195)) )/ sqrt(nanvar(T_.HighestEd(thrt_high))+nanvar(T_.HighestEd(matches_195)) )

cohd = (nanmean(T_.Black(thrt_high))-nanmean(T_.Black(setdiff(1:3645,thrt_high))) )/ sqrt(nanvar(T_.Black(thrt_high))+nanvar(T_.Black(setdiff(1:3645,thrt_high))) )
cohd = (nanmean(T_.Black(thrt_high))-nanmean(T_.Black(thrt_low)) )/ sqrt(nanvar(T_.Black(thrt_high))+nanvar(T_.Black(thrt_low)) )
cohd = (nanmean(T_.Black(thrt_high))-nanmean(T_.Black(matches_195)) )/ sqrt(nanvar(T_.Black(thrt_high))+nanvar(T_.Black(matches_195)) )

cohd = (nanmean(T_.White(thrt_high))-nanmean(T_.White(setdiff(1:3645,thrt_high))) )/ sqrt(nanvar(T_.White(thrt_high))+nanvar(T_.White(setdiff(1:3645,thrt_high))) )
cohd = (nanmean(T_.White(thrt_high))-nanmean(T_.White(thrt_low)) )/ sqrt(nanvar(T_.White(thrt_high))+nanvar(T_.White(thrt_low)) )
cohd = (nanmean(T_.White(thrt_high))-nanmean(T_.White(matches_195)) )/ sqrt(nanvar(T_.White(thrt_high))+nanvar(T_.White(matches_195)) )

cohd = (nanmean(T_.Hispanic(thrt_high))-nanmean(T_.Hispanic(setdiff(1:3645,thrt_high))) )/ sqrt(nanvar(T_.Hispanic(thrt_high))+nanvar(T_.Hispanic(setdiff(1:3645,thrt_high))) )
cohd = (nanmean(T_.Hispanic(thrt_high))-nanmean(T_.Hispanic(thrt_low)) )/ sqrt(nanvar(T_.Hispanic(thrt_high))+nanvar(T_.Hispanic(thrt_low)) )
cohd = (nanmean(T_.Hispanic(thrt_high))-nanmean(T_.Hispanic(matches_195)) )/ sqrt(nanvar(T_.Hispanic(thrt_high))+nanvar(T_.Hispanic(matches_195)) )


matches_thrt_195 = unique(matches); % save these as your matched controls
thrt_198 = thrt_high; % save these as your exposure group


%% %%%%% PLS
%% neurocognitve latent variable

clear all
T_0 = read_abcd_beh('~\df_filled_beh_Y0_more.csv');
T_2 =read_abcd_beh('~\df_filled_beh_Y2_more.csv');


load('Z:\PM2.5\New Manuscript/usable_subs.mat')
T_ = usable_subs;

load('matches_thrt_195.mat')
load('thrt_198.mat')

addpath(genpath('~\PLS_folder\Pls'));
rng('default');
groups{1} = [matches_thrt_195];
groups{2} = [thrt_198];
nTimes = 2;

datamat_lst = cell(length(groups),1); 
% ####*****************####
T_0.cog = nanmean([T_0.NBK T_0.NIH5./120],2);
T_2.cog = nanmean([T_2.NBK T_2.NIH5./120],2);


for g = 1:numel(groups)
    for c = 1:nTimes
        for s = 1:numel(groups{g})
            sub_id = char(T_.Var1(groups{g}(s)));
            
            vec2 =[];
            if c==1
                locb = find(T_0.subids == sub_id);
                if length(locb)==0
                    sub_id
                end
                vec2 = [T_0.mat_score(locb)   T_0.cog(locb)  T_0.smri_thick_cdk_mean(locb)];
                                                    
            end
            
            if c==2
                loct = find(T_2.subids == sub_id);
                if length(loct)==0
                    sub_id
                end
                vec2 = [T_2.mat_score(loct)  T_2.cog(loct)  T_2.smri_thick_cdk_mean(loct)];
                  

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
% save as PLS_result_195lowthrt_198hthrt_y0_y2.mat
%% plotting the results (Table 2 top half and Figure 3A)
clear all

load('PLS_result_195lowthrt_198hthrt_y0_y2.mat');  
ps = result.perm_result.sprob
load('~\usable_subs_cog.mat')
T_q = usable_subs;

T_0 = read_abcd_beh('~\df_filled_beh_Y0_more.csv');
T_2 =read_abcd_beh('~\df_filled_beh_Y2_more.csv');
tfu = table(T_2.subids, T_2.interview_age,T_2.nruns,T_2.pds_ss,'VariableNames',{'subids','age_Y2','nruns_Y2','pds_ss_Y2'});
T_ = innerjoin(T_q,tfu,"LeftKeys","Var1","RightKeys","subids");

load('matches_thrt_195.mat')
load('thrt_198.mat')
matches = matches_thrt_195;
thrthigh = thrt_198;

% *** table 2 top half
thrt_low = find(usable_subs.cumul_ace == 0 & ~isnan(usable_subs.cog) & ~isnan(usable_subs.neur));
T_.AsiOth = T_.Asian + T_.Other;

[nanmean(T_.HighestEd(thrthigh)) sqrt(nanvar(T_.HighestEd(thrthigh)))]
[nanmean(T_.HighestEd(thrt_low)) sqrt(nanvar(T_.HighestEd(thrt_low)))]
[nanmean(T_.HighestEd(matches)) sqrt(nanvar(T_.HighestEd(matches)))]
cohd = (nanmean(T_.HighestEd(thrthigh))-nanmean(T_.HighestEd(thrt_low)) )/ sqrt(nanvar(T_.HighestEd(thrthigh))+nanvar(T_.HighestEd(thrt_low)) )
cohd = (nanmean(T_.HighestEd(thrthigh))-nanmean(T_.HighestEd(matches)) )/ sqrt(nanvar(T_.HighestEd(thrthigh))+nanvar(T_.HighestEd(matches)) )

      % levels= 1:10, labels = c("5000", "8500", "14000", "20500", "30000","42500", 7:"62500", 8:"87500", "150000", "200000") )
nanmean(T_.Income(thrthigh)); 
62500 + (7.6899-7)*(87500 - 62500)
nanmean(T_.Income(thrt_low)) ;
87500 + (8.2134-8)*(150000 - 87500)
nanmean(T_.Income(matches)); 
62500 + (7.7611-7)*(87500 - 62500)
cohd = (nanmean(T_.Income(thrthigh))-nanmean(T_.Income(thrt_low)) )/ sqrt(nanvar(T_.Income(thrthigh))+nanvar(T_.Income(thrt_low)) )
cohd = (nanmean(T_.Income(thrthigh))-nanmean(T_.Income(matches)) )/ sqrt(nanvar(T_.Income(thrthigh))+nanvar(T_.Income(matches)) )

[length(find(T_.White(thrthigh))) length(find(T_.White(thrthigh)))/length(thrthigh)]
[length(find(T_.White(thrt_low))) length(find(T_.White(thrt_low)))/length(thrt_low)]
[length(find(T_.White(matches))) length(find(T_.White(matches)))/length(matches)]
cohd = (nanmean(T_.White(thrthigh))-nanmean(T_.White(thrt_low)) )/ sqrt(nanvar(T_.White(thrthigh))+nanvar(T_.White(thrt_low)) )
cohd = (nanmean(T_.White(thrthigh))-nanmean(T_.White(matches)) )/ sqrt(nanvar(T_.White(thrthigh))+nanvar(T_.White(matches)) )


[length(find(T_.Hispanic(thrthigh))) length(find(T_.Hispanic(thrthigh)))/length(thrthigh)]
[length(find(T_.Hispanic(thrt_low))) length(find(T_.Hispanic(thrt_low)))/length(thrt_low)]
[length(find(T_.Hispanic(matches))) length(find(T_.Hispanic(matches)))/length(matches)]
cohd = (nanmean(T_.Hispanic(thrthigh))-nanmean(T_.Hispanic(thrt_low)) )/ sqrt(nanvar(T_.Hispanic(thrthigh))+nanvar(T_.Hispanic(thrt_low)) )
cohd = (nanmean(T_.Hispanic(thrthigh))-nanmean(T_.Hispanic(matches)) )/ sqrt(nanvar(T_.Hispanic(thrthigh))+nanvar(T_.Hispanic(matches)) )


[length(find(T_.Black(thrthigh))) length(find(T_.Black(thrthigh)))/length(thrthigh)]
[length(find(T_.Black(thrt_low))) length(find(T_.Black(thrt_low)))/length(thrt_low)]
[length(find(T_.Black(matches))) length(find(T_.Black(matches)))/length(matches)]
cohd = (nanmean(T_.Black(thrthigh))-nanmean(T_.Black(thrt_low)) )/ sqrt(nanvar(T_.Black(thrthigh))+nanvar(T_.Black(thrt_low)) )
cohd = (nanmean(T_.Black(thrthigh))-nanmean(T_.Black(matches)) )/ sqrt(nanvar(T_.Black(thrthigh))+nanvar(T_.Black(matches)) )


[length(find(T_.AsiOth(thrthigh))) length(find(T_.AsiOth(thrthigh)))/length(thrthigh)]
[length(find(T_.AsiOth(thrt_low))) length(find(T_.AsiOth(thrt_low)))/length(thrt_low)]
[length(find(T_.AsiOth(matches))) length(find(T_.AsiOth(matches)))/length(matches)]
cohd = (nanmean(T_.AsiOth(thrthigh))-nanmean(T_.AsiOth(thrt_low)) )/ sqrt(nanvar(T_.AsiOth(thrthigh))+nanvar(T_.AsiOth(thrt_low)) )
cohd = (nanmean(T_.AsiOth(thrthigh))-nanmean(T_.AsiOth(matches)) )/ sqrt(nanvar(T_.AsiOth(thrthigh))+nanvar(T_.AsiOth(matches)) )

[length(find(T_.Male_bin(thrthigh))) length(find(T_.Male_bin(thrthigh)))/length(thrthigh)]
[length(find(T_.Male_bin(thrt_low))) length(find(T_.Male_bin(thrt_low)))/length(thrt_low)]
[length(find(T_.Male_bin(matches))) length(find(T_.Male_bin(matches)))/length(matches)]
cohd = (nanmean(T_.Male_bin(thrthigh))-nanmean(T_.Male_bin(thrt_low)) )/ sqrt(nanvar(T_.Male_bin(thrthigh))+nanvar(T_.Male_bin(thrt_low)) )
cohd = (nanmean(T_.Male_bin(thrthigh))-nanmean(T_.Male_bin(matches)) )/ sqrt(nanvar(T_.Male_bin(thrthigh))+nanvar(T_.Male_bin(matches)) )


[nanmean(T_.age(thrthigh)) sqrt(nanvar(T_.age(thrthigh)))]
[nanmean(T_.age(thrt_low)) sqrt(nanvar(T_.age(thrt_low)))]
[nanmean(T_.age(matches)) sqrt(nanvar(T_.age(matches)))]
cohd = (nanmean(T_.age(thrthigh))-nanmean(T_.age(thrt_low)) )/ sqrt(nanvar(T_.age(thrthigh))+nanvar(T_.age(thrt_low)) )
cohd = (nanmean(T_.age(thrthigh))-nanmean(T_.age(matches)) )/ sqrt(nanvar(T_.age(thrthigh))+nanvar(T_.age(matches)) )

[nanmean(T_.age_Y2(thrthigh)) sqrt(nanvar(T_.age_Y2(thrthigh)))]/12
[nanmean(T_.age_Y2(thrt_low)) sqrt(nanvar(T_.age_Y2(thrt_low)))]/12
[nanmean(T_.age_Y2(matches)) sqrt(nanvar(T_.age_Y2(matches)))]/12
cohd = (nanmean(T_.age_Y2(thrthigh))-nanmean(T_.age_Y2(thrt_low)) )/ sqrt(nanvar(T_.age_Y2(thrthigh))+nanvar(T_.age_Y2(thrt_low)) )
cohd = (nanmean(T_.age_Y2(thrthigh))-nanmean(T_.age_Y2(matches)) )/ sqrt(nanvar(T_.age_Y2(thrthigh))+nanvar(T_.age_Y2(matches)) )

% ***

Sites = zeros(3645,22);
for K=1:22
    Sites((T_.site_id_l == K),K) = 1;
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


covariates_y0cont = [T_.pFD_y0(matches)   Sites(matches,:)  T_.White(matches) T_.Black(matches) T_.Hispanic(matches) T_.HighestEd(matches) T_.Income(matches) T_.Male_bin(matches)    T_.Ingenia(matches) T_.Achieva(matches) T_.Discovery(matches)  T_.Pfit(matches) T_.nruns(matches) T_.reshist_addr1_adi_wsum(matches) T_.pds_ss(matches) T_.interview_age(matches)]; % T_.PF10_lavaan(matches)  
covariates_y2cont = [T_.pFD_y2(matches)   Sites(matches,:)  T_.White(matches) T_.Black(matches) T_.Hispanic(matches) T_.HighestEd(matches) T_.Income(matches) T_.Male_bin(matches)     T_.Ingenia(matches) T_.Achieva(matches) T_.Discovery(matches)  T_.Pfit(matches) T_.nruns_Y2(matches) T_.reshist_addr1_adi_wsum(matches) T_.pds_ss_Y2(matches) T_.age_Y2(matches)]; % T_.PF10_lavaan(matches)  
covariates_y0pm = [T_.pFD_y0(thrthigh)   Sites(thrthigh,:)  T_.White(thrthigh) T_.Black(thrthigh) T_.Hispanic(thrthigh) T_.HighestEd(thrthigh) T_.Income(thrthigh) T_.Male_bin(thrthigh)     T_.Ingenia(thrthigh) T_.Achieva(thrthigh) T_.Discovery(thrthigh)  T_.Pfit(thrthigh) T_.nruns(thrthigh) T_.reshist_addr1_adi_wsum(thrthigh)  T_.pds_ss(thrthigh)  T_.interview_age(thrthigh)]; % T_.PF10_lavaan(thrthigh)  
covariates_y2pm = [T_.pFD_y2(thrthigh)   Sites(thrthigh,:)  T_.White(thrthigh) T_.Black(thrthigh) T_.Hispanic(thrthigh) T_.HighestEd(thrthigh) T_.Income(thrthigh) T_.Male_bin(thrthigh)   T_.Ingenia(thrthigh) T_.Achieva(thrthigh) T_.Discovery(thrthigh)  T_.Pfit(thrthigh) T_.nruns_Y2(thrthigh) T_.reshist_addr1_adi_wsum(thrthigh)  T_.pds_ss_Y2(thrthigh)  T_.age_Y2(thrthigh)]; % T_.PF10_lavaan(thrthigh)   

covariates = [covariates_y0cont; covariates_y2cont; covariates_y0pm; covariates_y2pm];

[a1 ] = partialcorr([result.usc(:,lv),adjy0cont],covariates, 'Rows', 'pairwise');
[a2 ] = partialcorr([result.usc(:,lv),adjy2cont], covariates, 'Rows', 'pairwise');
[a3 ] = partialcorr([result.usc(:,lv),adjy0pm], covariates, 'Rows', 'pairwise');
[a4 ] = partialcorr([result.usc(:,lv),adjy2pm], covariates, 'Rows', 'pairwise');

a1s =[]; a2s =[]; a3s =[]; a4s =[]; b1s=[];b2s=[];b3s=[]; b1sadj=[];b2sadj=[];b3sadj=[];
dats = [result.datamat_lst{1, 1}; result.datamat_lst{2, 1}];

nboot = 5000;
for jj=1:nboot
   rand_inds1 = randi([1,lg1],lg1,1); rand_inds11 = randi([1,lg1],lg1,1); 
   rand_inds2 = randi([1,lg2],lg2,1); rand_inds22 = randi([1,lg2],lg2,1);
   inds = [rand_inds1; lg1+ rand_inds11; 2*lg1+ rand_inds2; 2*lg1+lg2+ rand_inds22];
   a11 = partialcorr([result.usc(inds,lv),adjy0cont(inds)], covariates(inds,:), 'Rows', 'pairwise');
   a1s = [a1s; a11(1,2)];
   a22 = partialcorr([result.usc(inds,lv),adjy2cont(inds)], covariates(inds,:), 'Rows', 'pairwise');
a2s = [a2s; a22(1,2)];
a33 = partialcorr([result.usc(inds,lv),adjy0pm(inds)], covariates(inds,:), 'Rows', 'pairwise');
a3s = [a3s; a33(1,2)];
a44 = partialcorr([result.usc(inds,lv),adjy2pm(inds)] , covariates(inds,:), 'Rows', 'pairwise');
a4s = [a4s; a44(1,2)];


   b11 = partialcorr([result.usc(inds,lv), dats(inds,1)], covariates(inds,:), 'Rows', 'pairwise');
   b1s = [b1s; b11(1,2)];
   b22 = partialcorr([result.usc(inds,lv),dats(inds,2)], covariates(inds,:), 'Rows', 'pairwise');
b2s = [b2s; b22(1,2)];
b33 = partialcorr([result.usc(inds,lv),dats(inds,3)], covariates(inds,:), 'Rows', 'pairwise');
b3s = [b3s; b33(1,2)];


end
   
 gt_contrast2 = (1).*[a1s, a2s, a3s, a4s];
 gt_contrast2b = (1).*[b1s, b2s, b3s];

 p1 = 1 - length(find(  (gt_contrast2(:,2)-gt_contrast2(:,1)) > (mean(gt_contrast2(:,4)) - mean(gt_contrast2(:,3)) ) ) )/nboot
p2 = 1 - (length(find(  (gt_contrast2(:,4)>gt_contrast2(:,3))  ) )/nboot);
p3 = 1 - (length(find(  (gt_contrast2(:,2)>gt_contrast2(:,1))  ) )/nboot);
 p4 = 1 - length(find(  (gt_contrast2(:,2)+gt_contrast2(:,1)) > (mean(gt_contrast2(:,4) + (gt_contrast2(:,3)))  )))/nboot
 ef = mean(gt_contrast2(:,2))-mean(gt_contrast2(:,1))

ef = (mean(gt_contrast2(:,2))+mean(gt_contrast2(:,1)) )/2 - (mean(gt_contrast2(:,3))+mean(gt_contrast2(:,4)) )/2


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
    'XTickLabelRotation',45,'FontSize',22,'YTick',[-.5,-.25,0,.25,.5]);
ylabel('PLS group-by-time loading (partial corr)'); ylim([-.5,.52]);
f1.Position = [100 100 500 900];
hold off
export_fig -m5 -transparent thrt_pls_Cog_lv1_gtbars_partial_nop.jpg

f2 = figure;
hold on
for k=1:3
    col = [.84 .75 .13];;
 mmb = mean(gt_contrast2b(:,k));
    negsb = mmb - prctile(gt_contrast2b(:,k),5);
    posb = prctile(gt_contrast2b(:,k),95) - mmb;
    bar(k*.8,mean(gt_contrast2b(:,k)),'FaceColor',col,'EdgeColor',col,'BarWidth',.45);
    errorbar(k*.8,mean(gt_contrast2b(:,k)),negsb,posb,'.k');
end
set(gca,'XTick',[.8,1.6,2.4],'XTickLabels',{'Cort FC maturation','Cog task perf','Cort GM thickness'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS neurocognitve loading (partial corr)'); ylim([-1,1]); %axis square
f2.Position = [100 100 500 900];
hold off  
export_fig -m5 -transparent thrt_pls_Cog_lv1_ncbars_partial_nop.jpg


%%

%%  ######################################### Threat+ #####################################################################################
%%  %% Create matched control group for the threat+ participants
clear all
T_0 = read_abcd_beh('~\df_filled_beh_Y0_more.csv');
T_2 =read_abcd_beh('~\df_filled_beh_Y2_more.csv');

load('~/with_2b_2024.mat')
s_table = readtable("~/df_filled_demog_withcog_Y0_more.csv");
sub_11868 = table(s_table.subid_x);
usable_subs0 = innerjoin(sub_11868(with_2b_2024,:),T_0,"LeftKeys","Var1","RightKeys","subids");
tfu = table(T_2.subids, T_2.NBK,  T_2.NIH5,  T_2.smri_thick_cdk_mean,'VariableNames',{'subids','NBK_y2' ,'NIH_y2' ,'thick_y2'});
usable_subs = innerjoin(usable_subs0,tfu,"LeftKeys","Var1","RightKeys","subids");
usable_subs.cog = nanmean([usable_subs.NBK usable_subs.NIH5./120],2) + nanmean([usable_subs.NBK_y2 usable_subs.NIH_y2./120],2);
usable_subs.neur = usable_subs.smri_thick_cdk_mean + usable_subs.thick_y2;% save as usable_subs_cog.mat

% this section runs iteratively to make the maatched group which is then
% saved as num_matches_subs_bothth.mat
both_high1 = find(usable_subs.cumul_ace >1 & ~isnan(usable_subs.cog) & ~isnan(usable_subs.neur) ...
    ); % 2 is threat exclusively and 3 is both

%load(['~/num_matches_subs_bothth.mat']);%uncomment this after the first
% time running the script which generates the list

both_high = both_high1;
% both_high =both_high1(num_matches_subs_both >1); %uncomment this after the first
% time running the script which generates the list

both_high_n = length(both_high);
depv_low = find(usable_subs.cumul_ace == 0 & ~isnan(usable_subs.cog) & ~isnan(usable_subs.neur)); % 

usable_subs.raceth = usable_subs.White;
usable_subs.raceth(usable_subs.Black == 1) = 2;
usable_subs.raceth(usable_subs.Hispanic == 1) = 3;
usable_subs.raceth(usable_subs.Asian == 1) = 4;
usable_subs.raceth(usable_subs.Other == 1) = 5;
usable_subs.interview_age(isnan(usable_subs.interview_age)) = nanmean(usable_subs.interview_age);
usable_subs.agey = round(usable_subs.interview_age/6);
usable_subs.HighestEd(isnan(usable_subs.HighestEd)) = nanmean(usable_subs.HighestEd);
usable_subs.edu = round(usable_subs.HighestEd/5);

income_matches = cell(both_high_n,1);
age_matches = cell(both_high_n,1);
sex_matches = cell(both_high_n,1);
raceth_matches = cell(both_high_n,1);
edu_matches = cell(both_high_n,1);

for k =1:both_high_n
    
    if ~isnan(usable_subs.Income(both_high(k)) )
        a = find(round(usable_subs.Income) == round(usable_subs.Income(both_high(k))));
    else a = depv_low;
    end
        income_matches{k} = intersect(depv_low,a);
        
    if ~isnan(usable_subs.agey(both_high(k)))
        a = find(usable_subs.agey == usable_subs.agey(both_high(k)));
    else a = depv_low;
    end
        age_matches{k} = intersect(depv_low,a);
    
    if ~isnan(usable_subs.Male_bin(both_high(k)))
        a = find(usable_subs.Male_bin == usable_subs.Male_bin(both_high(k)));
    else a = depv_low;
    end
         sex_matches{k} = intersect(depv_low,a);
         
    if ~isnan(usable_subs.raceth(both_high(k)))
        a = find(usable_subs.raceth == usable_subs.raceth(both_high(k)));
    else a = depv_low;
    end
         raceth_matches{k} = intersect(depv_low,a);
    
    if ~isnan(usable_subs.edu(both_high(k)))
        a = find(usable_subs.edu == usable_subs.edu(both_high(k)));
    else a = depv_low;
    end
         edu_matches{k} = intersect(depv_low,a);
            
end
cont_mj = cell(both_high_n,1); cont_mj_full = cell(both_high_n,1); cont_mj_bare = cell(both_high_n,1);
for k = 1:both_high_n
    gg = intersect( intersect(sex_matches{k},income_matches{k}), intersect(raceth_matches{k},edu_matches{k}));
    cont_mj_full{k} = setdiff(intersect( gg, age_matches{k}  ),  both_high);
    cont_mj{k} = setdiff( gg,  both_high);
    gg2 = intersect(income_matches{k}, raceth_matches{k});
    cont_mj_bare{k} = setdiff( intersect(gg2,edu_matches{k}),  both_high);

end

num_matches_subs_both =[];
for l=1:414
num_matches_subs_both = [num_matches_subs_both; length(cont_mj_bare{l})];
end
% save as num_matches_subs_bothth.mat

HCss = []; llls =[];
for ww = 1:10000
    HCs =NaN(both_high_n,1);
    for k=1:both_high_n
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
matches_346 = unique(matches); 

T_ = usable_subs;
cohd = (nanmean(T_.age(both_high))-nanmean(T_.age(setdiff(1:3645,both_high))) )/ sqrt(nanvar(T_.age(both_high))+nanvar(T_.age(setdiff(1:3645,both_high))) )
cohd = (nanmean(T_.age(both_high))-nanmean(T_.age(depv_low)) )/ sqrt(nanvar(T_.age(both_high))+nanvar(T_.age(depv_low)) )
cohd = (nanmean(T_.age(both_high))-nanmean(T_.age(matches_346)) )/ sqrt(nanvar(T_.age(both_high))+nanvar(T_.age(matches_346)) )

cohd = (nanmean(T_.Income(both_high))-nanmean(T_.Income(setdiff(1:3645,both_high))) )/ sqrt(nanvar(T_0.Income(both_high))+nanvar(T_.Income(setdiff(1:3645,both_high))) )
cohd = (nanmean(T_.Income(both_high))-nanmean(T_.Income(depv_low)) )/ sqrt(nanvar(T_.Income(both_high))+nanvar(T_.Income(depv_low)) )
cohd = (nanmean(T_.Income(both_high))-nanmean(T_.Income(matches_346)) )/ sqrt(nanvar(T_.Income(both_high))+nanvar(T_.Income(matches_346)) )


cohd = (nanmean(T_.HighestEd(both_high))-nanmean(T_.HighestEd(setdiff(1:3645,both_high))) )/ sqrt(nanvar(T_.HighestEd(both_high))+nanvar(T_.HighestEd(setdiff(1:3645,both_high))) )
cohd = (nanmean(T_.HighestEd(both_high))-nanmean(T_.HighestEd(depv_low)) )/ sqrt(nanvar(T_.HighestEd(both_high))+nanvar(T_.HighestEd(depv_low)) )
cohd = (nanmean(T_.HighestEd(both_high))-nanmean(T_.HighestEd(matches_346)) )/ sqrt(nanvar(T_.HighestEd(both_high))+nanvar(T_.HighestEd(matches_346)) )

cohd = (nanmean(T_.Black(both_high))-nanmean(T_.Black(setdiff(1:3645,both_high))) )/ sqrt(nanvar(T_.Black(both_high))+nanvar(T_.Black(setdiff(1:3645,both_high))) )
cohd = (nanmean(T_.Black(both_high))-nanmean(T_.Black(depv_low)) )/ sqrt(nanvar(T_.Black(both_high))+nanvar(T_.Black(depv_low)) )
cohd = (nanmean(T_.Black(both_high))-nanmean(T_.Black(matches_346)) )/ sqrt(nanvar(T_.Black(both_high))+nanvar(T_.Black(matches_346)) )

cohd = (nanmean(T_.White(both_high))-nanmean(T_.White(setdiff(1:3645,both_high))) )/ sqrt(nanvar(T_.White(both_high))+nanvar(T_.White(setdiff(1:3645,both_high))) )
cohd = (nanmean(T_.White(both_high))-nanmean(T_.White(depv_low)) )/ sqrt(nanvar(T_.White(both_high))+nanvar(T_.White(depv_low)) )
cohd = (nanmean(T_.White(both_high))-nanmean(T_.White(matches_346)) )/ sqrt(nanvar(T_.White(both_high))+nanvar(T_.White(matches_346)) )

cohd = (nanmean(T_.Hispanic(both_high))-nanmean(T_.Hispanic(setdiff(1:3645,both_high))) )/ sqrt(nanvar(T_.Hispanic(both_high))+nanvar(T_.Hispanic(setdiff(1:3645,both_high))) )
cohd = (nanmean(T_.Hispanic(both_high))-nanmean(T_.Hispanic(depv_low)) )/ sqrt(nanvar(T_.Hispanic(both_high))+nanvar(T_.Hispanic(depv_low)) )
cohd = (nanmean(T_.Hispanic(both_high))-nanmean(T_.Hispanic(matches_346)) )/ sqrt(nanvar(T_.Hispanic(both_high))+nanvar(T_.Hispanic(matches_346)) )


matches_both_346 = unique(matches); % save these as your matched controls
both_382 = both_high; % save these as your exposure group


%% %%%%% PLS
%% neurocognitve latent variable

clear all
T_0 = read_abcd_beh('~\df_filled_beh_Y0_more.csv');
T_2 =read_abcd_beh('~\df_filled_beh_Y2_more.csv');


load('~/usable_subs_cog.mat')
T_ = usable_subs;

load('matches_both_346.mat')
load('both_382.mat')

addpath(genpath('~\PLS_folder\Pls'));
rng('default');
groups{1} = [matches_both_346];
groups{2} = [both_382];
nTimes = 2;

datamat_lst = cell(length(groups),1); 
% ####*****************####
T_0.cog = nanmean([T_0.NBK T_0.NIH5./120],2);
T_2.cog = nanmean([T_2.NBK T_2.NIH5./120],2);


for g = 1:numel(groups)
    for c = 1:nTimes
        for s = 1:numel(groups{g})
            sub_id = char(T_.Var1(groups{g}(s)));
            
            vec2 =[];
            if c==1
                locb = find(T_0.subids == sub_id);
                if length(locb)==0
                    sub_id
                end
                vec2 = [T_0.mat_score(locb)   T_0.cog(locb)  T_0.smri_thick_cdk_mean(locb)];
                                                    
            end
            
            if c==2
                loct = find(T_2.subids == sub_id);
                if length(loct)==0
                    sub_id
                end
                vec2 = [T_2.mat_score(loct)  T_2.cog(loct)  T_2.smri_thick_cdk_mean(loct)];
                  

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
% save as PLS_result_346low_382hbothth_y0_y2.mat
%% plotting the results (Table 3 top half and Figure 2B)
load('PLS_result_346low_382hbothth_y0_y2.mat');  
ps = result.perm_result.sprob
load('~/usable_subs_cog.mat')
T_q = usable_subs;

T_0 = read_abcd_beh('~\df_filled_beh_Y0_more.csv');
T_2 =read_abcd_beh('~\df_filled_beh_Y2_more.csv');
tfu = table(T_2.subids, T_2.interview_age,T_2.nruns,T_2.pds_ss,'VariableNames',{'subids','age_Y2','nruns_Y2','pds_ss_Y2'});
T_ = innerjoin(T_q,tfu,"LeftKeys","Var1","RightKeys","subids");

load('matches_both_346.mat')
load('both_382.mat')
matches = matches_both_346;
thrthigh = both_382;

% *** table 3 top half
thrt_low = find(usable_subs.cumul_ace == 0 & ~isnan(usable_subs.cog) & ~isnan(usable_subs.neur));
T_.AsiOth = T_.Asian + T_.Other;

[nanmean(T_.HighestEd(thrthigh)) sqrt(nanvar(T_.HighestEd(thrthigh)))]
[nanmean(T_.HighestEd(thrt_low)) sqrt(nanvar(T_.HighestEd(thrt_low)))]
[nanmean(T_.HighestEd(matches)) sqrt(nanvar(T_.HighestEd(matches)))]
cohd = (nanmean(T_.HighestEd(thrthigh))-nanmean(T_.HighestEd(thrt_low)) )/ sqrt(nanvar(T_.HighestEd(thrthigh))+nanvar(T_.HighestEd(thrt_low)) )
cohd = (nanmean(T_.HighestEd(thrthigh))-nanmean(T_.HighestEd(matches)) )/ sqrt(nanvar(T_.HighestEd(thrthigh))+nanvar(T_.HighestEd(matches)) )

      % levels= 1:10, labels = c("5000", "8500", "14000", "20500", "30000","42500", 7:"62500", 8:"87500", "150000", "200000") )
nanmean(T_.Income(thrthigh)); 
42500 + (6.9725-6)*(62500 - 42500)
nanmean(T_.Income(thrt_low)) ;
87500 + (8.2134-8)*(150000 - 87500)
nanmean(T_.Income(matches)); 
62500 + (7.2052-7)*(87500 - 62500)
cohd = (nanmean(T_.Income(thrthigh))-nanmean(T_.Income(thrt_low)) )/ sqrt(nanvar(T_.Income(thrthigh))+nanvar(T_.Income(thrt_low)) )
cohd = (nanmean(T_.Income(thrthigh))-nanmean(T_.Income(matches)) )/ sqrt(nanvar(T_.Income(thrthigh))+nanvar(T_.Income(matches)) )

[length(find(T_.White(thrthigh))) length(find(T_.White(thrthigh)))/length(thrthigh)]
[length(find(T_.White(thrt_low))) length(find(T_.White(thrt_low)))/length(thrt_low)]
[length(find(T_.White(matches))) length(find(T_.White(matches)))/length(matches)]
cohd = (nanmean(T_.White(thrthigh))-nanmean(T_.White(thrt_low)) )/ sqrt(nanvar(T_.White(thrthigh))+nanvar(T_.White(thrt_low)) )
cohd = (nanmean(T_.White(thrthigh))-nanmean(T_.White(matches)) )/ sqrt(nanvar(T_.White(thrthigh))+nanvar(T_.White(matches)) )


[length(find(T_.Hispanic(thrthigh))) length(find(T_.Hispanic(thrthigh)))/length(thrthigh)]
[length(find(T_.Hispanic(thrt_low))) length(find(T_.Hispanic(thrt_low)))/length(thrt_low)]
[length(find(T_.Hispanic(matches))) length(find(T_.Hispanic(matches)))/length(matches)]
cohd = (nanmean(T_.Hispanic(thrthigh))-nanmean(T_.Hispanic(thrt_low)) )/ sqrt(nanvar(T_.Hispanic(thrthigh))+nanvar(T_.Hispanic(thrt_low)) )
cohd = (nanmean(T_.Hispanic(thrthigh))-nanmean(T_.Hispanic(matches)) )/ sqrt(nanvar(T_.Hispanic(thrthigh))+nanvar(T_.Hispanic(matches)) )


[length(find(T_.Black(thrthigh))) length(find(T_.Black(thrthigh)))/length(thrthigh)]
[length(find(T_.Black(thrt_low))) length(find(T_.Black(thrt_low)))/length(thrt_low)]
[length(find(T_.Black(matches))) length(find(T_.Black(matches)))/length(matches)]
cohd = (nanmean(T_.Black(thrthigh))-nanmean(T_.Black(thrt_low)) )/ sqrt(nanvar(T_.Black(thrthigh))+nanvar(T_.Black(thrt_low)) )
cohd = (nanmean(T_.Black(thrthigh))-nanmean(T_.Black(matches)) )/ sqrt(nanvar(T_.Black(thrthigh))+nanvar(T_.Black(matches)) )


[length(find(T_.AsiOth(thrthigh))) length(find(T_.AsiOth(thrthigh)))/length(thrthigh)]
[length(find(T_.AsiOth(thrt_low))) length(find(T_.AsiOth(thrt_low)))/length(thrt_low)]
[length(find(T_.AsiOth(matches))) length(find(T_.AsiOth(matches)))/length(matches)]
cohd = (nanmean(T_.AsiOth(thrthigh))-nanmean(T_.AsiOth(thrt_low)) )/ sqrt(nanvar(T_.AsiOth(thrthigh))+nanvar(T_.AsiOth(thrt_low)) )
cohd = (nanmean(T_.AsiOth(thrthigh))-nanmean(T_.AsiOth(matches)) )/ sqrt(nanvar(T_.AsiOth(thrthigh))+nanvar(T_.AsiOth(matches)) )

[length(find(~T_.Male_bin(thrthigh))) length(find(~T_.Male_bin(thrthigh)))/length(thrthigh)]
[length(find(~T_.Male_bin(thrt_low))) length(find(~T_.Male_bin(thrt_low)))/length(thrt_low)]
[length(find(~T_.Male_bin(matches))) length(find(~T_.Male_bin(matches)))/length(matches)]
cohd = (nanmean(T_.Male_bin(thrthigh))-nanmean(T_.Male_bin(thrt_low)) )/ sqrt(nanvar(T_.Male_bin(thrthigh))+nanvar(T_.Male_bin(thrt_low)) )
cohd = (nanmean(T_.Male_bin(thrthigh))-nanmean(T_.Male_bin(matches)) )/ sqrt(nanvar(T_.Male_bin(thrthigh))+nanvar(T_.Male_bin(matches)) )


[nanmean(T_.age(thrthigh)) sqrt(nanvar(T_.age(thrthigh)))]
[nanmean(T_.age(thrt_low)) sqrt(nanvar(T_.age(thrt_low)))]
[nanmean(T_.age(matches)) sqrt(nanvar(T_.age(matches)))]
cohd = (nanmean(T_.age(thrthigh))-nanmean(T_.age(thrt_low)) )/ sqrt(nanvar(T_.age(thrthigh))+nanvar(T_.age(thrt_low)) )
cohd = (nanmean(T_.age(thrthigh))-nanmean(T_.age(matches)) )/ sqrt(nanvar(T_.age(thrthigh))+nanvar(T_.age(matches)) )

[nanmean(T_.age_Y2(thrthigh)) sqrt(nanvar(T_.age_Y2(thrthigh)))]/12
[nanmean(T_.age_Y2(thrt_low)) sqrt(nanvar(T_.age_Y2(thrt_low)))]/12
[nanmean(T_.age_Y2(matches)) sqrt(nanvar(T_.age_Y2(matches)))]/12
cohd = (nanmean(T_.age_Y2(thrthigh))-nanmean(T_.age_Y2(thrt_low)) )/ sqrt(nanvar(T_.age_Y2(thrthigh))+nanvar(T_.age_Y2(thrt_low)) )
cohd = (nanmean(T_.age_Y2(thrthigh))-nanmean(T_.age_Y2(matches)) )/ sqrt(nanvar(T_.age_Y2(thrthigh))+nanvar(T_.age_Y2(matches)) )

% ***


Sites = zeros(3645,22);
for K=1:22
    Sites((T_.site_id_l == K),K) = 1;
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


T_.depv = zeros(3645,1);
T_.depv(find(T_.cumul_ace == 1))=1;
covariates_y0cont = [T_.depv(matches) T_.pFD_y0(matches)   Sites(matches,:)  T_.White(matches) T_.Black(matches) T_.Hispanic(matches) T_.HighestEd(matches) T_.Income(matches) T_.Male_bin(matches)  T_.interview_age(matches)  T_.Ingenia(matches) T_.Achieva(matches) T_.Discovery(matches)  T_.Pfit(matches) T_.nruns(matches) T_.reshist_addr1_adi_wsum(matches)  T_.pds_ss(matches)]; % T_.PF10_lavaan(matches)
covariates_y2cont = [T_.depv(matches) T_.pFD_y2(matches)   Sites(matches,:)  T_.White(matches) T_.Black(matches) T_.Hispanic(matches) T_.HighestEd(matches) T_.Income(matches) T_.Male_bin(matches)  T_.age_Y2(matches)   T_.Ingenia(matches) T_.Achieva(matches) T_.Discovery(matches)  T_.Pfit(matches) T_.nruns_Y2(matches) T_.reshist_addr1_adi_wsum(matches)  T_.pds_ss_Y2(matches)]; % T_.PF10_lavaan(matches)
covariates_y0pm = [T_.depv(thrthigh) T_.pFD_y0(thrthigh)   Sites(thrthigh,:)  T_.White(thrthigh) T_.Black(thrthigh) T_.Hispanic(thrthigh) T_.HighestEd(thrthigh) T_.Income(thrthigh) T_.Male_bin(thrthigh)  T_.interview_age(thrthigh)   T_.Ingenia(thrthigh) T_.Achieva(thrthigh) T_.Discovery(thrthigh)  T_.Pfit(thrthigh) T_.nruns(thrthigh) T_.reshist_addr1_adi_wsum(thrthigh)  T_.pds_ss(thrthigh)]; % T_.PF10_lavaan(thrthigh)
covariates_y2pm = [T_.depv(thrthigh) T_.pFD_y2(thrthigh)   Sites(thrthigh,:)  T_.White(thrthigh) T_.Black(thrthigh) T_.Hispanic(thrthigh) T_.HighestEd(thrthigh) T_.Income(thrthigh) T_.Male_bin(thrthigh)  T_.age_Y2(thrthigh)  T_.Ingenia(thrthigh) T_.Achieva(thrthigh) T_.Discovery(thrthigh)  T_.Pfit(thrthigh) T_.nruns_Y2(thrthigh) T_.reshist_addr1_adi_wsum(thrthigh)  T_.pds_ss_Y2(thrthigh)]; % T_.PF10_lavaan(thrthigh)


covariates = [covariates_y0cont; covariates_y2cont; covariates_y0pm; covariates_y2pm];

[a1 ] = partialcorr([result.usc(:,lv),adjy0cont],covariates, 'Rows', 'pairwise');
[a2 ] = partialcorr([result.usc(:,lv),adjy2cont], covariates, 'Rows', 'pairwise');
[a3 ] = partialcorr([result.usc(:,lv),adjy0pm], covariates, 'Rows', 'pairwise');
[a4 ] = partialcorr([result.usc(:,lv),adjy2pm], covariates, 'Rows', 'pairwise');

a1s =[]; a2s =[]; a3s =[]; a4s =[]; b1s=[];b2s=[];b3s=[]; b1sadj=[];b2sadj=[];b3sadj=[];
dats = [result.datamat_lst{1, 1}; result.datamat_lst{2, 1}];

nboot = 5000;
for jj=1:nboot
   rand_inds1 = randi([1,lg1],lg1,1); rand_inds11 = randi([1,lg1],lg1,1); 
   rand_inds2 = randi([1,lg2],lg2,1); rand_inds22 = randi([1,lg2],lg2,1);
   inds = [rand_inds1; lg1+ rand_inds11; 2*lg1+ rand_inds2; 2*lg1+lg2+ rand_inds22];
   a11 = partialcorr([result.usc(inds,lv),adjy0cont(inds)], covariates(inds,:), 'Rows', 'pairwise');
   a1s = [a1s; a11(1,2)];
   a22 = partialcorr([result.usc(inds,lv),adjy2cont(inds)], covariates(inds,:), 'Rows', 'pairwise');
a2s = [a2s; a22(1,2)];
a33 = partialcorr([result.usc(inds,lv),adjy0pm(inds)], covariates(inds,:), 'Rows', 'pairwise');
a3s = [a3s; a33(1,2)];
a44 = partialcorr([result.usc(inds,lv),adjy2pm(inds)] , covariates(inds,:), 'Rows', 'pairwise');
a4s = [a4s; a44(1,2)];


   b11 = partialcorr([result.usc(inds,lv), dats(inds,1)], covariates(inds,:), 'Rows', 'pairwise');
   b1s = [b1s; b11(1,2)];
   b22 = partialcorr([result.usc(inds,lv),dats(inds,2)], covariates(inds,:), 'Rows', 'pairwise');
b2s = [b2s; b22(1,2)];
b33 = partialcorr([result.usc(inds,lv),dats(inds,3)], covariates(inds,:), 'Rows', 'pairwise');
b3s = [b3s; b33(1,2)];

end

 gt_contrast2 = (1).*[a1s, a2s, a3s, a4s];
 gt_contrast2b = (1).*[b1s, b2s, b3s];

 p1 = 1 - length(find(  (gt_contrast2(:,2)-gt_contrast2(:,1)) > (mean(gt_contrast2(:,4)) - mean(gt_contrast2(:,3)) ) ) )/nboot
p2 = 1 - (length(find(  (gt_contrast2(:,4)>gt_contrast2(:,3))  ) )/nboot);
p3 = 1 - (length(find(  (gt_contrast2(:,2)>gt_contrast2(:,1))  ) )/nboot);
 p4 = 1 - length(find(  (gt_contrast2(:,2)+gt_contrast2(:,1)) > (mean(gt_contrast2(:,4) + (gt_contrast2(:,3)))  )))/nboot
ef = (mean(gt_contrast2(:,2)+gt_contrast2(:,1)) - mean(gt_contrast2(:,3)+gt_contrast2(:,4)))/2

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
export_fig -m5 -transparent bothth_pls_Cog_lv1_gtbars_partial_nop.jpg

f2 = figure;
hold on
for k=1:3
    col = [.84 .75 .13];

 mmb = mean(gt_contrast2b(:,k));
 
   negsb = mmb - prctile(gt_contrast2b(:,k),5);
     posb = prctile(gt_contrast2b(:,k),95) - mmb;
      bar(k*.8,mean(gt_contrast2b(:,k)),'FaceColor',col,'EdgeColor',col,'BarWidth',.45);
      errorbar(k*.8,mean(gt_contrast2b(:,k)),negsb,posb,'.k');
end
set(gca,'XTick',[.8,1.6,2.4],'XTickLabels',{'Cort FC maturation','Cog task perf','Cort GM thickness'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS neurocognitve loading (partial corr)'); ylim([-1,1]); %axis square
f2.Position = [100 100 500 900];
hold off  
export_fig -m5 -transparent bothth_pls_Cog_lv1_ncbars_partial_nop.jpg