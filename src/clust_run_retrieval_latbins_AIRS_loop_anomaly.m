%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS_loop_anomaly.m
%---------------------------------------------------------------------------

%% this is the timestep : 1: 365 (coincidecne : there are 365 days/year and
%% I did 16 day averages .... so 365/16 steps per year ... and 2002-2018 is
%% 16 years so total of 365/16 * 16 = 365 steps

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 20
%JOB = 116
%JOB = JUNK    %% for use with loop_clust_run_retrieval_latbins_AIRS_loop_anomaly.m

%---------------------------------------------------------------------------
addpath /home/sergio/MATLABCODE/oem_pkg
addpath Plotutils
%---------------------------------------------------------------------------
% Doing debug?
 driver.debug = false;
 driver.debug_dir = '../Debug';

% Open debug file if desired
 if driver.debug
    writelog('open');
 end;
%---------------------------------------------------------------------------
% Loop over anomaly time setps, 365 of them
%---------------------------------------------------------------------------
for iLat = 1 : 40
   driver.iibin     = iLat;
   driver.i16daytimestep = JOB;
   ix = iLat;

   driver.oem.dofit = true;
   driver.lls.dofit = false;
   driver.oem.nloop = 1;
   driver.oem.doplots = false;
%---------------------------------------------------------------------------
   % Override many settings and add covariance matrix
  topts = [];

  topts = struct;
  topts.numchan = 2645;
  topts.chan_LW_SW = 0;    %% just LW/MW
  %topts.chan_LW_SW = -1;  %% all chans

  topts.descORasc = +1;   %% descending default
  %topts.descORasc = -1;  %% ascending  new; note have not really re-done jacs

  %topts.offsetrates = +1;  %% add in offset of 0.01 K/yr; note for anomaly, if topts.ocb_set = +1;topts.offsetrates = +1
                           %% then we LINEARLY add in 0,01 K/year, depending on timestep
  topts.ocb_set = 0;        %% obs
  %topts.ocb_set = +1;       %% cal
  %topts.addco2jacs = +1;   %% add co2 = 2.2 ppmv/yr to cal

  topts.obs_corr_matrix = -1; %% use diagnol obs uncertainty, gives nice 2.0 ppmv/yr CO2 for all latbins
  %topts.obs_corr_matrix = +1; %% instead of diagnol obs cov, use full cov matrix

  topts.tie_sst_lowestlayer = -1;  %% no  cross cov between SST and lowest layer
  topts.tie_sst_lowestlayer = +1;  %% yes cross cov between SST and lowest layer

  topts.invtype = 1; %% pinv, DEFAULT
  %topts.invtype = 3; %% inverse, Texas A&M
  %topts.invtype = 4; %% Se ridge rgression
  %topts.invtype = 5; %% inverse, minimum eigenvalue

  topts.iNlays_retrieve = 97;
  topts.iNlays_retrieve = 10;
  topts.iNlays_retrieve = 20;

  if topts.ocb_set == 0
    driver.filename = ['OutputAnomaly_OBS/' num2str(iLat,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif topts.ocb_set == 1
    driver.filename = ['OutputAnomaly_CAL/' num2str(iLat,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  end

  if topts.iNlays_retrieve >= 97 & ~exist(driver.filename)
    [driver,aux] = strow_override_defaults_latbins_AIRS(driver,topts);
  elseif topts.iNlays_retrieve < 97 & ~exist(driver.filename)
    [driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,topts.iNlays_retrieve,topts);
  end

%---------------------------------------------------------------------------
  % Do the retrieval
  if ~exist(driver.filename)
     driver = retrieval(driver,aux);
  else
    driver.rateset.rates = zeros(2645,1);
  end
%---------------------------------------------------------------------------
  % Save retrieval output from this loop
  
  if sum(abs(driver.rateset.rates)) > 0 & ~exist(driver.filename)
    save(driver.filename,'-struct','driver');
  elseif sum(abs(driver.rateset.rates)) < eps & ~exist(driver.filename)
    fprintf(1,'not saving %s since sum(abs(driver.rateset.rates)) = 0 \n',driver.filename);
  elseif sum(abs(driver.rateset.rates)) < eps & exist(driver.filename)
    fprintf(1,'not saving %s since it already exists\n',driver.filename);
   end

%   alld(JOB) = driver;

%---------------------------------------------------------------------------
% Some simple output
   if sum(abs(driver.rateset.rates)) > 0
     fprintf('Scalar Retrievals from OEM latbin %2i timestep %3i \n',iLat,JOB)
     fprintf(1,'CO2   (ppm)   %5.3f  +- %5.3f \n',driver.oem.finalrates(1),driver.oem.finalsigs(1));
     fprintf(1,'N2O   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(2),driver.oem.finalsigs(2));
     fprintf(1,'CH4   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(3),driver.oem.finalsigs(3));
     fprintf(1,'CFC11 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(4),driver.oem.finalsigs(4));
     fprintf(1,'SST   (K)    %5.3f  +- %5.3f \n',driver.oem.finalrates(5),driver.oem.finalsigs(5));
     %---------------------------------------------------------------------------
     % Pull interesting variable out for quick look
     co2(JOB) = driver.oem.finalrates(1);
     co2_sigs(JOB) = driver.oem.finalsigs(1); 
%     o3(JOB) = driver.oem.finalrates(2); 
%     o3_sigs(JOB) = driver.oem.finalsigs(2); 
     n2o(JOB) = driver.oem.finalrates(2); 
     n2o_sigs(JOB) = driver.oem.finalsigs(2); 
     ch4(JOB) = driver.oem.finalrates(3); 
     ch4_sigs(JOB) = driver.oem.finalsigs(3); 
     cfc11(JOB) = driver.oem.finalrates(4); 
     cfc11_sigs(JOB) = driver.oem.finalsigs(4); 
     sst(JOB) = driver.oem.finalrates(5); 
     sst_sigs(JOB) = driver.oem.finalsigs(5); 
   end

   % Plot Results
%{
  if topts.iNlays_retrieve >= 97
    plot_retrieval_latbins
  else
    plot_retrieval_latbins_fewlays
  end
   %disp('Hit return for next latitude'); pause
   pause(0.1)
%}

end % end of latbin loop  
%---------------------------------------------------------------------------
% Close debug file
if driver.debug
   writelog('close')
end
%--------------------------------------------------------------------------

%{
if topts.iNlays_retrieve >= 97
  plot_all_latbins_anom
else
  plot_all_latbins_fewlays_anom
end
%}
