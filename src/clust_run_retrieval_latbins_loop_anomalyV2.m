%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS_loop_anomaly.m
%---------------------------------------------------------------------------

addpath /home/sergio/MATLABCODE
system_slurm_stats

t1x = tic;

%% this is the timestep : 1: 365 (coincidecne : there are 365 days/year and
%% I did 16 day averages .... so 365/16 steps per year ... and 2002-2018 is
%% 16 years so total of 365/16 * 16 = 365 steps

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 20

%%%%%%%%%% ANOM or RATES %%%%%%%%%%
%JOB = 20   %%% uncomment this when trying to fit for linear rates!!! fix change_important_topts_settings, and set <<< driver.i16daytimestep = -1 >>>;  below
%%%%%%%%%% ANOM or RATES %%%%%%%%%%

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
% for this JOB latbin (1:40), loop over 365 anomaly time setps
%---------------------------------------------------------------------------
iTimeStep0 = 11; iTimeStepE = 31;
iTimeStep0 = 21; iTimeStepE = 21;
iTimeStep0 =  1; iTimeStepE = 365;

for iTimeStep = iTimeStep0 : iTimeStepE
  disp(' ')
  fprintf(1,'timestep = %3i latbin = %2i \n',iTimeStep,JOB);

%------------------------------------------------------------------------
%% <<<<<<<    no real need to touch any of this  >>>>>>>>
  driver.iibin     = JOB;

  %%%%%%%%%% ANOM or RATES %%%%%%%%%%
  driver.i16daytimestep = -1;         %% for the rates, not anomalies, RUN BY HAND BY UN-COMMENTING THIS LINE and 
                                      %% on top JOB = 20, in change_important_topts_settings.m also set topts.set_tracegas = -1;
  driver.i16daytimestep = iTimeStep;  %% this is when doing anomaly
  %%%%%%%%%% ANOM or RATES %%%%%%%%%%

  iLat = driver.iibin;
  ix = driver.iibin;

  driver.oem.dofit = true;
  driver.lls.dofit = false;
  driver.oem.nloop = 2;
  driver.oem.nloop = 3;
  driver.oem.nloop = 1;
  driver.oem.doplots = false;
%---------------------------------------------------------------------------
  change_important_topts_settings  % Override many settings and add covariance matrix

  if topts.ocb_set == 0 & driver.i16daytimestep > 0
    driver.outfilename = ['OutputAnomaly_OBS/' num2str(iLat,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif topts.ocb_set == 1 & driver.i16daytimestep > 0
    driver.outfilename = ['OutputAnomaly_CAL/' num2str(iLat,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif driver.i16daytimestep < 0
    driver.outfilename = ['Output/test' int2str(iLat) '.mat'];
  end

  if topts.iNlays_retrieve >= 97 & ~exist(driver.outfilename)
    [driver,aux] = strow_override_defaults_latbins_AIRS(driver,topts);
  elseif topts.iNlays_retrieve < 97 & ~exist(driver.outfilename)
    [driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,topts.iNlays_retrieve,topts);
  end
  %%[driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,topts.iNlays_retrieve,topts);

%---------------------------------------------------------------------------
  % Do the retrieval
  if ~exist(driver.outfilename)
     driver = retrieval(driver,aux);
  else
    driver.rateset.rates = zeros(2645,1);
  end
%---------------------------------------------------------------------------
  % Save retrieval output from this loop

  if isfield(topts,'iFixTz_NoFit')
    if topts.iFixTz_NoFit > 0
      junk = 1:length(driver.oem.finalrates)+length(aux.orig_temp_i);

      %% put ozone into correct (expected) spot
      driver.oem.finalrates(aux.orig_ozone_i) = driver.oem.finalrates(aux.orig_temp_i);
      driver.oem.finalsigs(aux.orig_ozone_i)  = driver.oem.finalsigs(aux.orig_temp_i);

      %% put fixed/unchanging T anomaly
      driver.oem.finalrates(aux.orig_temp_i) = aux.FixTz_NoFit;
      driver.oem.finalsigs(aux.orig_temp_i)  = 0;

      driver.jacobian.temp_i  = driver.jacobian.ozone_i;
      driver.jacobian.ozone_i = driver.jacobian.ozone_i + driver.jacobian.numlays;

    end
  end

  if isfield(topts,'iFixO3_NoFit')
    if topts.iFixO3_NoFit > 0

      %% put ozone into correct (expected) spot
      driver.oem.finalrates(aux.orig_ozone_i) = aux.FixO3_NoFit;
      driver.oem.finalsigs(aux.orig_ozone_i)  = 0;

      driver.jacobian.ozone_i = driver.jacobian.temp_i + driver.jacobian.numlays;
    end
  end

  if sum(abs(driver.rateset.rates)) > 0 & ~exist(driver.outfilename)
    save(driver.outfilename,'-struct','driver');
  elseif sum(abs(driver.rateset.rates)) < eps & ~exist(driver.outfilename)
    fprintf(1,'not saving %s since sum(abs(driver.rateset.rates)) = 0 \n',driver.outfilename);
  elseif sum(abs(driver.rateset.rates)) < eps & exist(driver.outfilename)
    fprintf(1,'not saving %s since it already exists\n',driver.outfilename);
   end

%   alld(JOB) = driver;

%---------------------------------------------------------------------------
% Some simple output
   if sum(abs(driver.rateset.rates)) > 0
     fprintf('Scalar Retrievals from OEM latbin %2i timestep %3i \n',JOB,iTimeStep)
     if topts.co2lays == 1
       fprintf(1,'CO2   (ppm)   %5.3f  +- %5.3f \n',driver.oem.finalrates(1),driver.oem.finalsigs(1));
       fprintf(1,'N2O   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(2),driver.oem.finalsigs(2));
       fprintf(1,'CH4   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(3),driver.oem.finalsigs(3));
       fprintf(1,'CFC11 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(4),driver.oem.finalsigs(4));
       fprintf(1,'CFC12 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(5),driver.oem.finalsigs(5));
       fprintf(1,'SST   (K)    %5.3f  +- %5.3f \n',driver.oem.finalrates(6),driver.oem.finalsigs(6));
     elseif topts.co2lays == 3
       fprintf(1,'CO2 lower trop  (ppm)   %5.3f  +- %5.3f \n',driver.oem.finalrates(1),driver.oem.finalsigs(1));
       fprintf(1,'N2O   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(4),driver.oem.finalsigs(4));
       fprintf(1,'CH4   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(5),driver.oem.finalsigs(5));
       fprintf(1,'CFC11 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(6),driver.oem.finalsigs(6));
       fprintf(1,'CFC12 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(7),driver.oem.finalsigs(7));
       fprintf(1,'SST   (K)    %5.3f  +- %5.3f \n',driver.oem.finalrates(8),driver.oem.finalsigs(8));
     end

     %---------------------------------------------------------------------------
     % Pull interesting variable out for quick look
     if topts.co2lays == 1
       co2(JOB) = driver.oem.finalrates(1);
       co2_sigs(JOB) = driver.oem.finalsigs(1); 
       n2o(JOB) = driver.oem.finalrates(2); 
       n2o_sigs(JOB) = driver.oem.finalsigs(2); 
       ch4(JOB) = driver.oem.finalrates(3); 
       ch4_sigs(JOB) = driver.oem.finalsigs(3); 
       cfc11(JOB) = driver.oem.finalrates(4); 
       cfc11_sigs(JOB) = driver.oem.finalsigs(4); 
       cfc12(JOB) = driver.oem.finalrates(5); 
       cfc12_sigs(JOB) = driver.oem.finalsigs(5); 
       sst(JOB) = driver.oem.finalrates(6); 
       sst_sigs(JOB) = driver.oem.finalsigs(6); 
     elseif topts.co2lays == 3
       co2(JOB) = driver.oem.finalrates(1);
       co2_sigs(JOB) = driver.oem.finalsigs(1); 
       n2o(JOB) = driver.oem.finalrates(4); 
       n2o_sigs(JOB) = driver.oem.finalsigs(4); 
       ch4(JOB) = driver.oem.finalrates(5); 
       ch4_sigs(JOB) = driver.oem.finalsigs(5); 
       cfc11(JOB) = driver.oem.finalrates(6); 
       cfc11_sigs(JOB) = driver.oem.finalsigs(6); 
       cfc12(JOB) = driver.oem.finalrates(7); 
       cfc12_sigs(JOB) = driver.oem.finalsigs(7); 
       sst(JOB) = driver.oem.finalrates(8); 
       sst_sigs(JOB) = driver.oem.finalsigs(8); 
     end
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

t2x = toc(t1x);
fprintf(1,'time taken (in seconds) for a %3i layer retrievals = %8.4f \n',topts.iNlays_retrieve,t2x)

if driver.i16daytimestep < 0
  plot_all_latbins_fewlays
end
