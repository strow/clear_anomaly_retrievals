%% Override many settings and add covariance matrix, these are the settings I typically use, so don't touch
topts = [];

topts = struct;

topts.resetnorm2one = -1; %%% DEFAULT, keep my scaling
%topts.resetnorm2one = +1; %%% new, reset stuff to 1.00000000

%% topts.dataset = 3; %% AIRS 11 year
%% topts.dataset = 2; %% IASI2AIRS 11 year
topts.dataset = +1; %% AIRS 16 year DEFAULT, has     nu cal fixed
%% topts.dataset = -1; %% AIRS 16 year        , without nu cal 

%% topts.co2lays = 3; %% try 3 fat layers
topts.co2lays = 1; %% DEFAULT, co2 column DEFAULT

topts.numchan = 2645;
topts.chan_LW_SW =  0;  %% just LW/MW DEFAULT, keep strat T chans                       470 chans
%% topts.chan_LW_SW = -1;  %% all chans,LW/MW and SW, keep strat T chans <<<<<<<           769 chans
%% topts.chan_LW_SW = +1;  %% LW/MW but avoid deep 15 um                                   436 chans
%% topts.chan_LW_SW = +2;  %% LW/MW/SW but avoid deep 15 um not yet coded                  735 chans
%% topts.chan_LW_SW = -2;  %% SW only                                                      299 chans

%% topts.descORasc = -1;  %% ascending  new; note have not really re-done jacs
topts.descORasc = +1;   %% descending DEFAULT

%topts.addco2jacs = +1;   %% add co2 = 2.2 ppmv/yr to cal

%%%% <<<<<<< this is what I typically change >>>>>>>>>>>>>>   
%%%% <<<<<<< this is what I typically change >>>>>>>>>>>>>>   
%%%% <<<<<<< this is what I typically change >>>>>>>>>>>>>>   

topts.iFixTG_NoFit = -1;  %% fit all trace gases DEFAULT
topts.iFixTG_NoFit = [4]; %% fit all trace gases except CFC11 = gas 4

topts.iFixTz_NoFit = -1;  %% DEFAULT, fit for T
%% topts.iFixTz_NoFit = +1;  %% fix T and do NOT fit for it

topts.iFixO3_NoFit = -1;  %% DEFAULT, fit for O3
%% topts.iFixO3_NoFit = +1;  %% fix O3 and do NOT fit for it
%% topts.iFixO3_NoFit =  0;  %% fix O3 a-priori to zero and do NOT fit for it

topts.offsetrates = -1; %% DEFAULT do NOT offset!!!!!
%%topts.offsetrates = +1;  %% add in offset of 0.01 K/yr; note for anomaly, if topts.ocb_set = +1;topts.offsetrates = +1
                            %% then we LINEARLY add in 0,01 K/year, depending on timestep DEFAULT OFF

%%%%%%%%%% ANOM or RATES %%%%%%%%%%
topts.set_tracegas = +1;   %% set   a priori for trace gas as eg 2.2/1.0/4.5/0.0/0.0 for CO2/N2O/CH4/CFC11/CFC12; if anomaly, adjust according to timestep!!!!!  
                           %% DEFAULT for anom
topts.set_tracegas = -1;   %% leave a priori for trace gas as 0, ONLY WHEN DOING RATES; 
                           %% also set driver.i16daytimestep = -1;   in clust_run_retrieval_latbins_AIRS_iasitimespan_loop_anomaly.m
%% topts.set_tracegas = +2;   %% set   a priori for trace gas as eg 2.2/1.0/4.5/0.0/0.0 for CO2/N2O/CH4/CFC11/CFC12; if anomaly, adjust according to PREVIOUS timestep!!!!!  BOOTSTRAP
                              %% can only run with  sbatch --exclude=cnode[204,225,260,267] --array=1-40 sergio_matlab_jobB.sbatch 2    since depends on previous timestep
%%%%%%%%%% ANOM or RATES %%%%%%%%%%

topts.ocb_set = 0;        %% obs  DEFAULT >>>>>>>>
topts.ocb_set = +1;       %% cal

%% topts.iDoStrowFiniteJac = -1;       %% do not change the time varying anomaly tracegas jacs                            done for all anomaly timesteps
%% topts.iDoStrowFiniteJac = +1;       %% +1 stick to Sergio tracegas jacs = BT(1.001 X(t,latbin)) - BT(1.00 X(t,latbin)) interp in time
%% topts.iDoStrowFiniteJac = +2;       %% +2 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin))       interp in time
topts.iDoStrowFiniteJac = +3;       %% DEFAULT +3 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin))       done for all anomaly timesteps  DEFAULT >>>>>>>
%% topts.iDoStrowFiniteJac = +4;       %% +3 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin))       done for all anomaly timesteps  with AGE of AIR 
%% topts.iDoStrowFiniteJac = -1;       %% do not change the time varying anomaly tracegas jacs                            done for all anomaly timesteps works great with new kcarta jacs

%% topts.iXJac = 0; %% const geo kcarta jacs
%% topts.iXJac = 1; %% varying geo sarta jacs  
topts.iXJac = 2; %% varying geo kcarta jacs DEFAULT >>>>>>>

%% topts.iNlays_retrieve = 10;
%% topts.iNlays_retrieve = 40;
%% topts.iNlays_retrieve = 60;
topts.iNlays_retrieve = 20; %%% DEFAULT >>>>>>>>
%% topts.iNlays_retrieve = 96;

%% topts.iChSet = 1;  %% old chans, ORIG >>>>
topts.iChSet = 2;  %% new chans, yes CFC11, NEW DEFAULT >>>>>
topts.iChSet = 3;  %% new chans, no CFC11

%%%% <<<<<<< this is what I typically change >>>>>>>>>>>>>>   
%%%% <<<<<<< this is what I typically change >>>>>>>>>>>>>>
%%%% <<<<<<< this is what I typically change >>>>>>>>>>>>>>   

topts.obs_corr_matrix = -1; %% use diagnol obs uncertainty, gives nice 2.0 ppmv/yr CO2 for all latbins, after fixed the channel uncdrtainty DEFAULT
%% topts.obs_corr_matrix = +1; %% instead of diagnol obs cov, use full cov matrix

%% topts.tie_sst_lowestlayer = -1;  %% no  cross cov between SST and lowest layer
topts.tie_sst_lowestlayer = +1;  %% yes cross cov between SST and lowest layer DEFAULT

topts.invtype = 1; %% pinv, DEFAULT
%% topts.invtype = 0; %% inv
%% topts.invtype = 2; %% invillco
%% topts.invtype = 3; %% inverse, Texas A&M
%% topts.invtype = 4; %% Se ridge regression
%% topts.invtype = 5; %% inverse, minimum eigenvalue
