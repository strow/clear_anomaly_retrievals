function [driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,iNlays_retrieve,topts);

settings.dataset = 1;       % (1) AIRS 16 year dataset (2) AIRS 11 year (IASI) (3) IASI2AIRS 11 year
settings.co2lays = 1;       % assume column jac for CO2, or 3 lays (gnd-500,500-trop,trop-TOA

settings.ocb_set = 0;       % 0 = obs, 1 = cal, -1 = bias
settings.numchan = 2645;    % L1b = 2378, L1c = 2645
settings.chan_LW_SW = 0;    % 0 is 640 to 1640, 1 is 700 to 1640, -1 is 640 to 2740

settings.set_tracegas = -1;            %% do we leave apriori as 0 or set CO2/N2o/CH4/CFC to be 2.2, 1, 4.5, -1             
                                       %%   (-1 = no, +1 = yes)
                                       %%   note : if anomaly, adjust a priori for CO2,N2O,CH4,CFC
settings.offsetrates  = -1;            %% do we add a constant offset to the spectral rates (-1 = no, +1 = yes)
                                       %% if driver.rateset.ocb_set  == 'obs';
settings.addco2jacs   = -1;            %% do we add co2 jacs to the spectral rates (-1 = no, +1 = yes)
                                       %% if driver.rateset.ocb_set  == 'cal';
settings.obs_corr_matrix = -1;         %% just use nc_error (-1) or try to be fancy and use full cov matrix (+1)

settings.tie_sst_lowestlayer = +1;     %% tie together SST with lowest T(z)
settings.invtype         = 1;          %% pinv, see /home/sergio/MATLABCODE/oem_pkg/rodgers.m
settings.iNlays_retrieve = 97;         %% do all 97 layers
settings.descORasc = +1;               %% descending default
settings.iXJac = 0;                    %% const geo jacs, replace as needed CO2/CH4/N20;  
                                       %% +2 uses kCARTA varying geo/trace, +1 uses SARTA varying geo/trace jacs, 0 = constant kcarta jacs
settings.iDoStrowFiniteJac = -1;       %% -1 : do not change the time varying anomaly jacs                                done for all anomaly timesteps
                                       %% +1 stick to Sergio tracegas jacs = BT(1.001 X(t,latbin)) - BT(1.00 X(t,latbin)) interp in time
                                       %% +2 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin)))      interp in time .. 
                                       %% +3 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin)))      done for all anomaly timesteps
settings.iChSet = 1;                   %% +1 default, old chans (about 500)
                                       %% +2, new chans (about 400) with CFC11,CFC12      and weak WV, bad chans gone
                                       %% +3, new chans (about 400) w/o  CFC11 with CFC12 and weak WV, bad chans gone

allowedparams = [{'ocb_set'},{'numchan'},{'chan_LW_SW'},{'iChSet'},{'set_tracegas'},{'offsetrates'},...
			    {'addco2jacs'},{'obs_corr_matrix'},{'invtype'},{'tie_sst_lowestlayer'},{'iNlays_retrieve'},...
                            {'descORasc'},{'dataset'},{'iXJac'},{'co2lays'},{'iDoStrowFiniteJac'}];


%disp('settings before')
%settings

if nargin == 3
  optvar = fieldnames(topts);
  for i = 1 : length(optvar)
   if (length(intersect(allowedparams,optvar{i})) == 1)
     eval(sprintf('settings.%s = topts.%s;', optvar{i}, optvar{i}));
   else
     fprintf(1,'topts param not in allowed list ... %s \n',optvar{i});
     error('quitting ');
   end
 end
end

%disp('settings after')
settings

aux.invtype = settings.invtype;

driver.topts = topts;

%---------------------------------------------------------------------------
% Which latitude bin
ix = driver.iibin;
%---------------------------------------------------------------------------
% Fitting [obs][cal][bias], pick one
if settings.ocb_set == -1
  driver.rateset.ocb_set  = 'bias';
  disp('bias')
elseif settings.ocb_set == 0
  driver.rateset.ocb_set  = 'obs';
  disp('obs')
elseif settings.ocb_set == +1
  driver.rateset.ocb_set  = 'cal';
  disp('cal')
elseif abs(settings.ocb_set) > 1
  settings.ocb_set
  error('incorrect settings.ocb_set')
end
%---------------------------------------------------------------------------
% Raw rate data file        
if settings.dataset == 1
  disp('AIRS 16 year rates or anomalies')
  if settings.descORasc == +1 & driver.i16daytimestep < 0
    disp('doing descending latbin rates')
    driver.rateset.datafile  = 'convert_strowrates2oemrates_random_16_year_v32_clear_nucal.mat';
    driver.rateset.datafile  = 'convert_strowrates2oemrates_random_16_year_v32_clear_nucal_obs_cal_bias.mat';
    if settings.ocb_set == +1  & driver.i16daytimestep < 0
      %% convert_strowrates2oemrates_random_16_year_v32_clear_nucal.mat does not have cal,bias so have to do this
      %% convert_strowrates2oemrates_random_16_year_v32_clear_nucal_obs_cal_bias.mat is complete so really no need for this
      driver.rateset.datafile  = 'convert_strowrates2oemrates_random_16_year_v32_clear.mat';
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      %% convert_strowrates2oemrates_random_16_year_v32_clear_nucal.mat does not have cal,bias so have to do this
      %% convert_strowrates2oemrates_random_16_year_v32_clear_nucal_obs_cal_bias.mat is complete so really no need for this
      driver.rateset.datafile  = 'convert_strowrates2oemrates_random_16_year_v32_clear_nucal_obs_cal_bias.mat';
    end
  elseif settings.descORasc == -1 & driver.i16daytimestep < 0
    disp('doing ascending latbin rates')
    driver.rateset.datafile  = 'convert_strowrates2oemrates_random_16_year_v32_clearasc_nucal.mat'; %% note bias/cal are from "usual" not nucal
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 0
    disp('doing descending OBS ANOMALY')
    driver.rateset.datafile = ['ANOM_16dayavg/latbin_180dayavg_' num2str(driver.iibin) '.mat'];  
    driver.rateset.datafile = ['ANOM_16dayavg/latbin_0dayavg_' num2str(driver.iibin) '.mat'];  
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 1
    disp('doing descending CAL ANOMALY')
    driver.rateset.datafile = ['ANOM_16dayavg/latbin_180dayavg_' num2str(driver.iibin) '_cal.mat'];  
    driver.rateset.datafile = ['ANOM_16dayavg/latbin_0dayavg_' num2str(driver.iibin) '_cal.mat'];  
  end

elseif settings.dataset == 2
  disp('AIRS 11 year rates or anomalies, overlaopping from 2007 with IASI')
  if settings.descORasc == +1 & driver.i16daytimestep < 0
    disp('doing descending latbin rates')
    driver.rateset.datafile  = 'convert_strowrates2oemrates_clear_11year_iasitimespan_obs_cal_bias.mat';
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 0 & settings.descORasc == +1
    disp('doing descending OBS ANOMALY')
    driver.rateset.datafile = ['IASI_ANOM_16dayavg/latbin_180dayavg_' num2str(driver.iibin) '.mat'];  
    driver.rateset.datafile = ['IASI_ANOM_16dayavg/latbin_0dayavg_' num2str(driver.iibin) '.mat'];  
  end

elseif settings.dataset == 3
  disp('IASI2AIRS 11 year rates or anomalies, overlaopping from 2007 with IASI')
  if settings.descORasc == +1 & driver.i16daytimestep < 0
    disp('doing descending latbin rates')
    driver.rateset.datafile  = 'convert_strowrates2oemrates_clear_11_year_iasi2airs_obs.mat';
  end
end

% Lag-1 correlation file; if using rate least-squares errors
driver.rateset.ncfile   = '../oem_pkg/Test/all_lagcor.mat';
driver.rateset.ncfile   = driver.rateset.datafile;

% Get rate data, do Q/A elsewhere
driver = get_rates(driver);
%---------------------------------------------------------------------------
% Jacobian file: f = 2378x1 and M_TS_jac_all = 36x2378x200
% driver.jacobian.filename = '../../oem_pkg/Test/M_TS_jac_all.mat';
%% clear sky

iXJac = settings.iXJac;
%if driver.i16daytimestep > 0
%  iXJac = 0; %% const geo kcarta jcs
%  iXJac = 1; %% varying geo sarta jacs
%  iXJac = 2; %% varying geo kcarta jacs
%end

if driver.i16daytimestep < 0
  if settings.descORasc == +1
    driver.jacobian.filename = '../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/JUNK/kcarta_M_TS_jac_all_5_97_97_97_2645.mat';
    fprintf(1,'reading in constant kcarta jac file %s \n',driver.jacobian.filename)
  else
    %% for now assume same jacs
    driver.jacobian.filename = '../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/JUNK/kcarta_M_TS_jac_all_5_97_97_97_2645.mat';
    fprintf(1,'reading in constant kcarta jac file %s \n',driver.jacobian.filename)
  end
elseif driver.i16daytimestep > 0
  junk = num2str(driver.i16daytimestep,'%03d');
  if iXJac == 1
    %% sarta time vary jacs
    %asarta  = load('../MakeJacsSARTA/SARTA_AIRSL1c_Anomaly365_16/RESULTS/sarta_182_fixCFC_M_TS_jac_all_5_97_97_97_2645.mat');
    driver.jacobian.filename = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_Anomaly365_16_with_seasonal_OldSarta_largepert//RESULTS/sarta_' junk '_fixCFC_M_TS_jac_all_5_97_97_97_2645.mat']; %% old sarta
    driver.jacobian.filename = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_Anomaly365_16_no_seasonal_OldSarta_smallpert//RESULTS/sarta_' junk '_fixCFC_M_TS_jac_all_5_97_97_97_2645.mat'];   %% old sarta
    driver.jacobian.filename = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_Anomaly365_16/RESULTS/sarta_' junk '_fixCFC_M_TS_jac_all_5_97_97_97_2645.mat'];                         %% new sarta
    fprintf(1,'iXJac == 1 reading in timestep sarta jac file %s \n',driver.jacobian.filename)

  elseif iXJac == 2  
    %% kcarta time vary jac
    %akcarta = load('../MakeJacskCARTA_CLR/Anomaly365_16/RESULTS/kcarta_182_M_TS_jac_all_5_97_97_97_2645.mat');

    %% until June 26, 2019, the finite diff tracegas jacs used dQ = 0.1, dT = 1.0
    %% together with iXJac = 2,iDoStrowFiniteJac = -1 gives anomaly_0dayavg_results_strow_latbin_unc.mat
    driver.jacobian.filename = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8_dQpert0.1_dTpert1.0/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2645.mat'];  

    %% after June 27, 2019, the finite diff tracegas jacs used dQ = 0.001, dT = 0.01 but BAD tracegas profiles
    %% use with  together with iXJac = 2,iDoStrowFiniteJac = -1
    driver.jacobian.filename = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8_tillJuly01_2019/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2645.mat']; 

    %% after July 3 to July 16, 2019, the finite diff tracegas jacs used dQ = 0.001,dT = 0.01 GOOD tracegas profiles (basically the glatm.dat tracegas profiles for CO2/N2O/CH4 adjusted in time)
    %% the profiles still have seasonal
    %% use with  together with iXJac = 2,iDoStrowFiniteJac = -1
    driver.jacobian.filename = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8_July12_2019_Great_But_with_Seasonal/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2645.mat']; 

    %% after July 16, 2019, the finite diff tracegas jacs used dQ = 0.001,dT = 0.01 GOOD tracegas profiles (basically the glatm.dat tracegas profiles for CO2/N2O/CH4 adjusted in time)
    %% the profiles do not have seasonal
    %% use with  together with iXJac = 2,iDoStrowFiniteJac = -1
    driver.jacobian.filename = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2645.mat']; 

    fprintf(1,'iXJac == 2 reading in timestep kcarta jac file %s \n',driver.jacobian.filename)
  elseif iXJac == 0
    %% constant kcarta jacs
    driver.jacobian.filename = '../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/JUNK/kcarta_M_TS_jac_all_5_97_97_97_2645.mat';
    fprintf(1,'iXJac == 0 reading in constant kcarta jac file %s \n',driver.jacobian.filename)
  end
end

%% THIS IS DEFAULT -- 4 column trace gas (CO2/N2O/CH4/CFC11/CFC12), 1 stemp, (97x3) geo
driver.jacobian.varname  = 'M_TS_jac_all';
driver.jacobian.scalar_i = 1:6;
driver.jacobian.water_i  = 7:103;
driver.jacobian.temp_i   = 104:200;
driver.jacobian.ozone_i  = 201:297;
driver.jacobian.numlays  = 97;

% Get jacobians, and combine the 97 layer T(z)/WV(z)/O3(z) into N layers
jac               = load(driver.jacobian.filename);

%% add in extra column for CFC12 >>>>>>
if iXJac == 0 | iXJac == 1
  m_ts_jac0_noCFC12 = squeeze(jac.M_TS_jac_all(ix,:,:));
  m_ts_jac0 = zeros(2645,297);
  m_ts_jac0(:,1:4)   = m_ts_jac0_noCFC12(:,1:4);
  m_ts_jac0(:,6:297) = m_ts_jac0_noCFC12(:,5:296);

  xm_ts_jac_coljac   = m_ts_jac0(:,1:6);
  if driver.i16daytimestep > 0  
    xm_ts_jac_coljac = replace_time_cfc12jac(xm_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3);
  else
    xm_ts_jac_coljac = replace_time_cfc12jac(xm_ts_jac_coljac,driver.iibin,floor(365/2),3);
  end
  m_ts_jac0(:,5) = xm_ts_jac_coljac(:,5);
  clear xm_ts_jac_coljac xm_ts_jac_coljac

  qrenormjunk = zeros(1,length(jac.qrenorm)+1);
  qrenormjunk(1:4)   = jac.qrenorm(1:4);
  qrenormjunk(5)     = jac.qrenorm(4);
  qrenormjunk(6:297) = jac.qrenorm(5:296);
  jac.qrenorm = qrenormjunk;
else
  m_ts_jac0 = squeeze(jac.M_TS_jac_all(ix,:,:));
  %% oops forgot to fix qrenorm, will do later
  %% qrenormjunk = zeros(1,length(jac.qrenorm)+1);
  %% qrenormjunk(1:4)   = jac.qrenorm(1:4);
  %% qrenormjunk(5)     = jac.qrenorm(4);
  %% qrenormjunk(6:297) = jac.qrenorm(5:296);
  %% jac.qrenorm = qrenormjunk;
end

m_ts_jac_coljac   = m_ts_jac0(:,1:6);
driver.qrenorm  = jac.qrenorm;       %% set this default

if iNlays_retrieve <= 60
  [m_ts_jac_wv,qWV,layWV]  = combinejaclays(m_ts_jac0,driver.jacobian.water_i,jac.qrenorm,iNlays_retrieve);
  [m_ts_jac_t,qT,layT]     = combinejaclays(m_ts_jac0,driver.jacobian.temp_i, jac.qrenorm,iNlays_retrieve);
  [m_ts_jac_o3,qO3,layO3]  = combinejaclays(m_ts_jac0,driver.jacobian.ozone_i,jac.qrenorm,iNlays_retrieve);
else
  fprintf(1,'setting iNlays_retrieve ( > 60) from %2i to 97 \n',iNlays_retrieve);
  iNlays_retrieve = 97;
end

%% replace CO2,N2O,CH4 jacs
if driver.i16daytimestep > 0  
  %% this is for 365 anomaly time steps
  %% put in time varying Jacobian, err no more need to do this??? well sarta has older CO2/CH4 but let's comment this for now
  iDoStrowFiniteJac = +2; %% testing Strows finite difference jacs CO2(t)-CO2(370) ...    6/24-27/2019 interp in time, used to be +1
  iDoStrowFiniteJac = +3; %% testing Strows finite difference jacs CO2(t)-CO2(370) ...    6/24-27/2019 at all anom timesteps
  iDoStrowFiniteJac = +1; %% testing new Sergio finite diff jacs                          interp in time
  iDoStrowFiniteJac = -1; %% default, rely on time varying CO2/N20/CH4 jacs from kcarta,  done for all anom timsteps

  iDoStrowFiniteJac = settings.iDoStrowFiniteJac; %% from 6/29/2019

  if iXJac == 0 & iDoStrowFiniteJac == 1
    fprintf(1,'updating const kCARTA CO2/N2O/CH4 jacs with Sergio interpolated time varying jacs...\n');
    %% const kCARTA jacs, update the trace gases
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
  elseif iXJac == 2 & iDoStrowFiniteJac == 2
    fprintf(1,'updating time varying kCARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs -1,+2 for testing 6/24-27/2019...\n');
    fprintf(1,'  note before July2, the tracegas profile (CO2/N2O/CH4) was US Std shoehorned and multiplied, so ppm was quite wonky except at 500 mb');
    fprintf(1,'else turned off \n')
    %% const kCARTA jacs, update the trace gases
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% only used this on 6/26
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% added this on 6/27
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% added this on 6/27
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% added this on 8/21
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% added this on 8/29 ---> comment this out if iUgh == 4
  elseif iXJac == 2 & iDoStrowFiniteJac == 3
    fprintf(1,'updating time varying kCARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 6/24-27/2019...\n');
    fprintf(1,'  note before July2, the tracegas profile (CO2/N2O/CH4) was US Std shoehorned and multiplied, so ppm was quite wonky except at 500 mb');
    fprintf(1,'  note after July2, have improved the tracegas profile (CO2/N2O/CH4) so they are the same shape as glatm.dat');
    fprintf(1,'else turned off \n')
    %% const kCARTA jacs, update the trace gases
    %% kcarta strow finite jacs
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% only used this on 6/26
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% added this on 6/27
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% added this on 6/27
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% added this on 8/21
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% added this on 8/29 ---> comment this out if iUgh == 4 <<<>>><<<>>><<<>>>
  elseif iXJac == 1 & iDoStrowFiniteJac == 3
    fprintf(1,'updating time varying SARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 6/24-27/2019...\n');
    fprintf(1,'  note before July2, the tracegas profile (CO2/N2O/CH4) was US Std shoehorned and multiplied, so ppm was quite wonky except at 500 mb');
    fprintf(1,'  note after July2, have improved the tracegas profile (CO2/N2O/CH4) so they are the same shape as glatm.dat');
    fprintf(1,'else turned off \n')
    %% sarta strow finite jacs
    iVarType = -3; %% this uses SARTA  finitediff jacs, which I have shown are bad?? or good??
    iVarType = +3; %% this uses kCARTA finitediff jacs, which I have shown are good, just want to test
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/21
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/21 --> comment this out if iUgh == 4
  end
else    
  %% this is for 1 average rate
  iTimeStepUse = 1;            %% put in constant Jacobian, at timestep 1 (2002/09)
  iTimeStepUse = 365;          %% put in constant Jacobian, at timestep 365 (2018/08)
  iTimeStepUse = floor(365/2); %% put in constant Jacobian, half way through (365/2 ==> 2009/09)
  if iXJac == 0
    fprintf(1,'updating CO2/N2O/CH4 jacs ...\n');
    %% const kCARTA jacs, update the trace gases
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,iTimeStepUse);
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,iTimeStepUse);
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,iTimeStepUse);
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,iTimeStepUse); %% added this on 8/21
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,iTimeStepUse); %% added this on 8/21
  end
end

if settings.co2lays == 3
  m_ts_jac_coljac = replace_time_co2_3layjac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep);
end

if iNlays_retrieve <= 60
  m_ts_jac = [m_ts_jac_coljac m_ts_jac_wv m_ts_jac_t m_ts_jac_o3];
else
  m_ts_jac = m_ts_jac0;
  qWV = (1:iNlays_retrieve) + 6;
end

bad = find(isnan(m_ts_jac) | isinf(m_ts_jac));
if length(bad) > 0
  fprintf('oopsy foound %5i NaN or Inf in jacobian, resetting to 0 \n',length(bad))
  m_ts_jac(bad) = 0;
end

iNlays_retrieve0 = iNlays_retrieve;
iNlays_retrieve = length(qWV);
ixlays = 1:iNlays_retrieve;
if settings.co2lays == 1
  driver.jacobian.scalar_i = 1:6;
  driver.jacobian.wvjaclays_offset = 6;
  if iNlays_retrieve <= 60
    driver.qrenorm  = [jac.qrenorm(1:6) qWV qT qO3];
  end
elseif settings.co2lays == 3
  driver.jacobian.scalar_i = 1:8;
  driver.jacobian.wvjaclays_offset = 8;
  if iNlays_retrieve <= 60
    driver.qrenorm  = [jac.qrenorm(1) jac.qrenorm(1) jac.qrenorm(1) jac.qrenorm(2:6) qWV qT qO3];
  end
end

iBruteForce_co2adjustjac = +1;
iBruteForce_co2adjustjac = -1;
if iBruteForce_co2adjustjac > 0
  disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CO2 jac ---> co2 jac * 0.85 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
  if settings.co2lays == 1
    m_ts_jac(:,1)  = m_ts_jac(:,1) * 0.85;
  elseif settings.co2lays == 3
    m_ts_jac(:,1:3)  = m_ts_jac(:,1:3) * 0.85;
  end
  disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CO2 jac ---> co2 jac * 0.85 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
end

if iNlays_retrieve <= 60
  driver.jacobian.water_i  = max(driver.jacobian.scalar_i) + ixlays;  
  driver.jacobian.temp_i   = max(driver.jacobian.water_i) + ixlays;  
  driver.jacobian.ozone_i  = max(driver.jacobian.temp_i) + ixlays;  
  driver.jacobian.numlays  = iNlays_retrieve;
  driver.jacobian.wvjaclays_used   = layWV;
end

aux.m_ts_jac    = m_ts_jac;
f = jac.f;
clear jac
%---------------------------------------------------------------------------
% Good channel set
lenrates = length(driver.rateset.rates);

iChSet = 2; %% new chans
iChSet = 1; %% old chans (default)
iChSet = 3; %% new chans, but no CFC11
iChSet = topts.iChSet;

ch = find_the_oem_channels(f,lenrates,settings.numchan,settings.chan_LW_SW,iChSet);

driver.topts.iChSet = iChSet;
driver.jacobian.chanset = ch;
%---------------------------------------------------------------------------
% Apriori file
%driver.oem.apriori_filename = 'apriori_lls';

% Load in apriori
%xb = load(driver.oem.apriori_filename,'apriori');
%xb = xb.apriori;

%xb = zeros(296,1);
xb = zeros(driver.jacobian.wvjaclays_offset + iNlays_retrieve*3,1);

if abs(settings.set_tracegas) ~= 1
  settings.set_tracegas
  error('incorrect setting for overriding xb tracegas values');
end
if settings.set_tracegas == +1 & driver.i16daytimestep < 0
  disp('setting constant rates for tracegas apriori : CO2 = 2.2  CH4 = 4.5 N2O = 0.8 CFC = -1.25')
  if settings.co2lays == 1
    xb(1) = 2.2;  % Set CO2 apriori
    xb(2) = 1;
    xb(3) = 4.5;
    xb(4) = -1.0;
    xb(5) = -1.0;

    xb(1) = 2.2;        % Set CO2 apriori
    xb(2) = 0.8;        % set N2O 
    xb(3) = 4.5;        % set CH4
    xb(4) = -1.25 * 0;  % set CFC11, before Aug 23 the mult was 1
    xb(5) = -1.25 * 0;  % set CFC12, before Aug 23 the mult was 1

  elseif settings.co2lays == 3
    xb(1) = 2.2;        % Set CO2 apriori lower trop
    xb(2) = 2.2;        % Set CO2 apriori mid trop
    xb(3) = 2.2;        % Set CO2 apriori strat

    xb(4) = 0.8;        % set N2O 
    xb(5) = 4.5;        % set CH4
    xb(6) = -1.25 * 0;  % set CFC11, before Aug 23 the mult was 1
    xb(7) = -1.25 * 0;  % set CFC12, before Aug 23 the mult was 1
  end

elseif settings.set_tracegas == +1 & driver.i16daytimestep > 0
  junk = 365/16; %% days per timestep 
  junk = (driver.i16daytimestep-1)/junk;
  str = ['setting time varying rates for tracegas apriori : CO2 = 2.2  CH4 = 4.5 N2O = 0.8 CFC = -1.25 for ' num2str(junk) ' years']; 
  disp(str);
  if settings.co2lays == 1
    deltaT = 365/16; %% days per timestep

    %% default all this while getting good results
    xb(1) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori
    xb(2) = 0.8 * (driver.i16daytimestep-1)/deltaT * 1.0;    % set N2O 
    xb(3) = 4.5 * (driver.i16daytimestep-1)/deltaT * 1.0;    % set CH4
    xb(4) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;    % set CFC11, before Aug 23 the mult was 1
    xb(5) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;    % set CFC12, before Aug 23 the mult was 1

  elseif settings.co2lays == 3
    deltaT = 365/16; %% days per timestep

    xb(1) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori lower trop
    xb(2) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori mid trop
    xb(3) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori strat

    xb(4) = 0.8 * (driver.i16daytimestep-1)/deltaT * 1.0;        % set N2O 
    xb(5) = 4.5 * (driver.i16daytimestep-1)/deltaT * 1.0;        % set CH4
    xb(6) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;  % set CFC11, before Aug 23 the mult was 1
    xb(7) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;  % set CFC12, before Aug 23 the mult was 1
  end
end

[mm,nn] = size(xb);
if nn > 1
  xb = xb(:,driver.iibin);
end

% A Priori stored in aux.xb
aux.xb = xb./driver.qrenorm';

%---------------------------------------------------------------------------
% SARTA forward model and other "representation" errors
driver.oem.sarta_error = 0.0;
driver.oem.xb = xb;  %% note this is un-normalized xb

%---------------------------------------------------------------------------
if abs(settings.offsetrates) ~= 1
  settings.offsetrates
  error('incorrect setting for overriding obs rates by constant');
elseif settings.offsetrates > 0 & driver.i16daytimestep < 0
   disp('Offset linear rates (obs or cal or bias) by constant to get sensitivity to AIRS BT drift')
   driver.rateset.rates = driver.rateset.rates + 0.01;
elseif settings.offsetrates > 0 & driver.i16daytimestep > 0 & settings.ocb_set == 0 & sum(abs(driver.rateset.rates)) > 0
   disp('Offset anomaly (obs) by constant 0.01/year to get sensitivity to AIRS BT drift')
   oktimes = load('ok365times.mat');
   oktimes = oktimes.okdates(driver.i16daytimestep);   %% need to get the CO2 jac at this time!!!!
   fprintf(1,'  --> oktimes = %8.6f which is %8.6f years away from 2002.75 \n',oktimes,oktimes - 2002.75)
   driver.rateset.rates = driver.rateset.rates + 0.01 * (oktimes - 2002.75);
elseif settings.offsetrates > 0 & driver.i16daytimestep > 0 & settings.ocb_set == 1 & sum(abs(driver.rateset.rates)) > 0
   disp('Offset anomaly (cal) by constant 0.01/year to get sensitivity to AIRS BT drift')
   oktimes = load('ok365times.mat');
   oktimes = oktimes.okdates(driver.i16daytimestep);   %% need to get the CO2 jac at this time!!!!
   fprintf(1,'  --> oktimes = %8.6f which is %8.6f years away from 2002.75 \n',oktimes,oktimes - 2002.75)
   driver.rateset.rates = driver.rateset.rates + 0.01 * (oktimes - 2002.75);
end

%------------------------------------------------------------------------
if abs(settings.addco2jacs) ~= 1
  settings.addco2jacs
  error('incorrect setting for adding co2jacs to ERA calcrates');
end
if strfind(driver.rateset.ocb_set,'cal') & settings.addco2jacs > 0
  disp('Offset rates to add in CO2 jacs to ERA rates, to pretend there is CO2')
  jac = load('co2_kcartajac.mat');
  %  haha = driver.rateset.rates;
  %  baba = jac.co2jac(ix,:)';
  %  whos haha baba
  driver.rateset.rates = driver.rateset.rates + jac.co2jac(ix,:)';
end

%---------------------------------------------------------------------------
% Modify rates with lag-1 correlation errors or add to above
if driver.i16daytimestep < 0
  nc_cor = nc_rates(driver);
  driver.rateset.unc_rates = driver.rateset.unc_rates.*nc_cor;   %% THIS IS AN ARRAY
end

%driver.rateset.unc_rates = 0.001 * ones(size(driver.rateset.unc_rates));   %% THIS IS AN ARRAY

%%%%%%%%%%%%%%%%%%%%%%%%% >>>>>>
if settings.obs_corr_matrix > 0
  addpath /home/sergio/MATLABCODE
  thecov = load('/home/sergio/MATLABCODE/oem_pkg_run/Simulate_Calcs/thecov_clear');
  [f2645,i2645] = map_2834_to_2645;
  junk = thecov.thecov; junk = junk(i2645,i2645);

  %% see https://www.mathworks.com/help/finance/corr2cov.html
  %%junk = corr2cov(driver.rateset.unc_rates,junk);
  %% new, fast
  junk0 = junk;
  junk0 = diag(driver.rateset.unc_rates) * junk0 * diag(driver.rateset.unc_rates);

  %% new,slow
  %[mm,nn] = size(junk);
  %for ii = 1 : mm
  %  for jj = 1 : nn
  %    junk(ii,jj) = junk(ii,jj) * driver.rateset.unc_rates(ii) * driver.rateset.unc_rates(jj);
  %  end
  %end
  %sum(sum(junk-junk0))/nn/nn
  junk = junk0;  
  fprintf(1,'numchans, rank, mean rank, condition number of obs cov matrix = %6i %6i %8.6e %8.6e \n',nn,rank(junk),rank(junk)/length(junk),cond(junk))

  %% old
  %plot(f2645,driver.rateset.unc_rates.^2,'r',thecov.fairs(i2645),diag(junk))
  %%keyboard_nowindow
  %slope = (2645+1);
  %xind = 1:2645; diagind=slope*(xind-1)+1;
  %junk(diagind) = driver.rateset.unc_rates.^2;

  driver.rateset.unc_rates = junk;  %% THIS IS A SQUARE MATRTIX
end
%%%%%%%%%%%%%%%%%%%%%%%%% >>>>>>

% Modify with estimated error in freq + regress errors 
%driver.rateset.unc_rates = ones(2378,1)*0.001 +driver.rateset.unc_rates.*nc_cor;
%load corr_errors
%load 15yr_corr_errors_2314chans

%---------------------------------------------------------------------------
% Do rate Q/A (empty for now)
%---------------------------------------------------------------------------

build_cov_matrices

