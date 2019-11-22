%---------------------------------------------------------------------------

addpath ../utility

%% this is the timestep : 1: 365 (coincidecne : there are 365 days/year and
%% I did 16 day averages .... so 365/16 steps per year ... and 2002-2018 is
%% 16 years so total of 365/16 * 16 = 365 steps

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));  %% JOB = 1 : 40 (fixed latbin, loop over time)
% JOB = 19;

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
  driver.i16daytimestep = iTimeStep;
  iLat = driver.iibin;
  ix = driver.iibin;

  driver.oem.dofit = true;
  driver.lls.dofit = false;
  driver.oem.nloop = 2;
  driver.oem.nloop = 3;
  driver.oem.nloop = 1;
  driver.oem.doplots = false;
%---------------------------------------------------------------------------
  % Override many settings and add covariance matrix
  change_important_topts_settings

  if topts.ocb_set == 0
    driver.outfilename = ['OutputAnomaly_OBS/' num2str(iLat,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif topts.ocb_set == 1
    driver.outfilename = ['OutputAnomaly_CAL/' num2str(iLat,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
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
  
  [status,ghash] = githash;
  driver.githash = ghash;
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
     fprintf('Scalar Retrievals from OEM latbin %2i timestep %3i \n',iLat,JOB)
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
