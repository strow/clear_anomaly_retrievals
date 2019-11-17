function [driver,aux] = strow_override_defaults_latbins_AIRS(driver,topts);

settings.dataset = 1;       % (1) AIRS 16 year dataset (2) AIRS 11 year (IASI) (3) IASI2AIRS 11 year
settings.co2lays = 1;       % assume column jac for CO2, or 3 lays (gnd-500,500-trop,trop-TOA

settings.ocb_set = 0;       % 0 = obs, 1 = cal, -1 = bias
settings.numchan = 2645;    % L1b = 2378, L1c = 2645
settings.chan_LW_SW = 0;    % 0 is 640 to 1640, 1 is 700 to 1640, -1 is 640 to 2740

settings.set_tracegas = -1;    %% do we leave apriori as 0 or set CO2/N2o/CH4/CFC to be 2.2, 1, 4.5, -1
                               %% (-1 = no, +1 = yes)
settings.offsetrates  = -1;    %% do we add a constant offset to the spectral rates (-1 = no, +1 = yes)
                               %% if driver.rateset.ocb_set  == 'obs';
settings.addco2jacs   = -1;    %% do we add co2 jacs to the spectral rates (-1 = no, +1 = yes)
                               %% if driver.rateset.ocb_set  == 'cal';
settings.obs_corr_matrix = -1; %% just use nc_error (-1) or try to be fancy and use full cov matrix (+1)

settings.tie_sst_lowestlayer = +1;     %% tie together SST with lowest T(z)
settings.invtype         = 1;          %% pinv, see /home/sergio/MATLABCODE/oem_pkg/rodgers.m
settings.iNlays_retrieve = 97;         %% do all 97 layers
settings.descORasc = +1;               %% descending default
settings.iXJac = 0;                    %% const geo jacs, replace as needed CO2/CH4/N20;  
                                       %%   +1 uses kCARTA varying geo/trace, +2 uses SARTA varying geo/trace jacs
settings.iDoStrowFiniteJac = -1;       %% -1 : do not change the time varying anomaly jacs                                done for all anomaly timesteps
                                       %% +1 stick to Sergio tracegas jacs = BT(1.001 X(t,latbin)) - BT(1.00 X(t,latbin)) interp in time
                                       %% +2 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin)))      interp in time .. 
                                       %% +3 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin)))      done for all anomaly timesteps

allowedparams = [{'ocb_set'},{'numchan'},{'chan_LW_SW'},{'set_tracegas'},{'offsetrates'},...
			    {'addco2jacs'},{'obs_corr_matrix'},{'invtype'},{'tie_sst_lowestlayer'},{'iNlays_retrieve'},...
                            {'descORasc'},{'dataset'},{'iXJac'},{'co2lays'},{'iDoStrowFiniteJac'}];

if nargin == 2
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
elseif settings.ocb_set == 0
  driver.rateset.ocb_set  = 'obs';
elseif settings.ocb_set == +1
  driver.rateset.ocb_set  = 'cal';
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
driver.jacobian.filename = '../../oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_Oct2018_CLR/sarta_fixCFC_M_TS_jac_all_5_97_97_97.mat';
driver.jacobian.filename = '../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/JUNK/kcarta_M_TS_jac_all_5_97_97_97_2645.mat';

error('whoa need so many JAC options set here see eg strow_override_defaults_latbins_AIRS_fewlays.m')

driver.jacobian.varname  = 'M_TS_jac_all';
driver.jacobian.scalar_i = 1:5;
driver.jacobian.water_i  = 6:102;
driver.jacobian.temp_i   = 103:199;
driver.jacobian.ozone_i  = 200:296;
driver.jacobian.numlays  = 97;

% Get jacobians
jac             = load(driver.jacobian.filename);
aux.m_ts_jac    = squeeze(jac.M_TS_jac_all(ix,:,:));
driver.qrenorm  = jac.qrenorm;
f = jac.f;
clear jac
%---------------------------------------------------------------------------
% Good channel set
iNum = 2378;
iNum = 2645;
iNum = settings.numchan;

if iNum ~= 2378 & iNum ~= 2645
  settings.numchan
  error('incorrect number of chans');
end

lenrates = length(driver.rateset.rates);
if iNum ~= lenrates
  fprintf(1,'you want to use good chans for %4i chans but your spectral rates are length %4i \n',iNum,lenrates);
  error('please fix');
end
if iNum == 2378
  ch = choose_goodchans_from_2378;                      %% Old 2378 chans
else
  ch = choose_goodchans_from_2645(settings.chan_LW_SW); %% New 2645 chans
end

driver.jacobian.chanset = ch;
%driver.jacobian.chanset = ch(1:100:length(ch));

chans38 = [        41          54         181         273         317         359         445         449 ...
                  532         758         903         904        1000        1020        1034        1055 ...
                 1075        1103        1249        1282        1291        1447        1475        1557 ...
                 1604        1614        1618        1660        1790        1866        1867        1868 ...
                 1878        1888        2112        2140        2321        2333];
chans41 = [chans38 2325        2339        2353]; chans41 = sort(chans41);
driver.jacobian.chanset = chans41(chans41 < 2200);

%---------------------------------------------------------------------------
% Apriori file
driver.oem.apriori_filename = 'apriori_lls';

% Load in apriori
xb = load(driver.oem.apriori_filename,'apriori');
xb = xb.apriori;

xb = zeros(296,1);

if abs(settings.set_tracegas) ~= 1
  settings.set_tracegas
  error('incorrect setting for overriding xb tracegas values');
end
if settings.set_tracegas == +1
  xb(1) = 2.2;  % Set CO2 apriori
  xb(2) = 1;
  xb(3) = 4.5;
  xb(4) = -1.0;

  xb(1) = 2.2;    % Set CO2 apriori
  xb(2) = 0.8;    % set N2O 
  xb(3) = 4.5;    % set CH4
  xb(4) = -1.25;  % set CFC
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

% Convert radiance rates into bt rate
% deriv = drdbt(f,rad2bt(f,driver.rateset.r));
% driver.rateset.rates = driver.rateset.rates./(1000*deriv);
% driver.rateset.unc_rates = driver.rateset.unc_rates./(1000*deriv);
%---------------------------------------------------------------------------
% Load in freq corrections
%load Data/dbt_10year  % alldbt
%driver.rateset.rates = driver.rateset.rates-alldbt(ix,:)'/10;

%------------------------------------------------------------------------
if abs(settings.offsetrates) ~= 1
  settings.offsetrates
  error('incorrect setting for overriding obs rates by constant');
end
if driver.rateset.ocb_set  == 'obs' & settings.offsetrates > 0
   disp('Offset rates to get sensitivity to AIRS BT drift')
   driver.rateset.rates = driver.rateset.rates + 0.01;
end
%------------------------------------------------------------------------
%{
load ../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/JUNK/kcarta_M_TS_jac_all_5_97_97_97_2645.mat
co2jac = squeeze(M_TS_jac_all(:,:,1));
save co2_kcartajac.mat co2jac
%}

if abs(settings.addco2jacs) ~= 1
  settings.addco2jacs
  error('incorrect setting for adding co2jacs to ERA calcrates');
end
if driver.rateset.ocb_set  == 'cal' & settings.addco2jacs > 0
  disp('Offset rates to add in CO2 jacs to ERA rates, to pretend there is CO2')
  jac = load('co2_kcartajac.mat');
  %  haha = driver.rateset.rates;
  %  baba = jac.co2jac(ix,:)';
  %  whos haha baba
  driver.rateset.rates = driver.rateset.rates + jac.co2jac(ix,:)';
end

%---------------------------------------------------------------------------
% Modify rates with lag-1 correlation errors or add to above
nc_cor = nc_rates(driver);
driver.rateset.unc_rates = driver.rateset.unc_rates.*nc_cor;   %% THIS IS AN ARRAY

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
fprintf(1,'numchans, rank, mean rank, condiition number of obs cov matrix = %6i %6i %8.6e %8.6e \n',nn,rank(junk),rank(junk)/length(junk),cond(junk))

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
