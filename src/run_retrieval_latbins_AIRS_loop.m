%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS_loop.m
%---------------------------------------------------------------------------
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
% Loop over latitude bins
%---------------------------------------------------------------------------
for JOB = 1 : 40
   ix = JOB;
   driver.iibin = JOB;
   driver.i16daytimestep = -9999;

   driver.oem.dofit = true;
   driver.lls.dofit = false;
   driver.oem.nloop = 1;
   driver.oem.doplots = false;
   driver.i16daytimestep = -1; %% this is NOT anomaly

%---------------------------------------------------------------------------
   % Override many settings and add covariance matrix
  topts = [];

  topts = struct;

  topts.dataset = 1;  %% AIRS 16 year
  topts.dataset = 3;  %% IASI2AIRS 11 year overlap with IASI
  topts.dataset = 2;  %% AIRS 11 year overlap with IASI

  topts.numchan = 2645;
  topts.chan_LW_SW = 0;    %% just LW/MW
  %topts.chan_LW_SW = -1;  %% all chans

  topts.descORasc = +1;   %% descending default
  %topts.descORasc = -1;   %% ascending  new; note have not really re-done jacs

  topts.offsetrates = +1;  %% add in offset of 0.01 K/yr
  topts.ocb_set = +1;      %% cal
  %topts.ocb_set = -1;      %% bias
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
  if topts.iNlays_retrieve >= 97
    [driver,aux] = strow_override_defaults_latbins_AIRS(driver,topts);
  else
    [driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,topts.iNlays_retrieve,topts);
  end

%---------------------------------------------------------------------------
  % Do the retrieval
   driver = retrieval(driver,aux);
%---------------------------------------------------------------------------
  % Save retrieval output from this loop
  driver.filename = ['Output/test' int2str(driver.iibin)];
  save(driver.filename,'-struct','driver');

   alld(JOB) = driver;
%---------------------------------------------------------------------------
% Some simple output
fprintf('Scalar Retrievals from OEM latbin %2i \n',JOB)
   fprintf(1,'CO2   (ppm)   %5.3f  +- %5.3f \n',driver.oem.finalrates(1),driver.oem.finalsigs(1));
   fprintf(1,'N2O   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(2),driver.oem.finalsigs(2));
   fprintf(1,'CH4   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(3),driver.oem.finalsigs(3));
   fprintf(1,'CFC11 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(4),driver.oem.finalsigs(4));
   fprintf(1,'SST   (K)    %5.3f  +- %5.3f \n',driver.oem.finalrates(5),driver.oem.finalsigs(5));
   %---------------------------------------------------------------------------
   % Pull interesting variable out for quick look
   co2(JOB) = driver.oem.finalrates(1);
   co2_sigs(JOB) = driver.oem.finalsigs(1); 
%   o3(JOB) = driver.oem.finalrates(2); 
%   o3_sigs(JOB) = driver.oem.finalsigs(2); 
   n2o(JOB) = driver.oem.finalrates(2); 
   n2o_sigs(JOB) = driver.oem.finalsigs(2); 
   ch4(JOB) = driver.oem.finalrates(3); 
   ch4_sigs(JOB) = driver.oem.finalsigs(3); 
   cfc11(JOB) = driver.oem.finalrates(4); 
   cfc11_sigs(JOB) = driver.oem.finalsigs(4); 
   sst(JOB) = driver.oem.finalrates(5); 
   sst_sigs(JOB) = driver.oem.finalsigs(5); 

   % Plot Results
  if topts.iNlays_retrieve >= 97
    plot_retrieval_latbins
  else
    plot_retrieval_latbins_fewlays
  end
   %disp('Hit return for next latitude'); pause
   pause(0.1)

end % end of latbin loop  
%---------------------------------------------------------------------------
% Close debug file
if driver.debug
   writelog('close')
end
%--------------------------------------------------------------------------
if topts.iNlays_retrieve >= 97
  plot_all_latbins
else
  plot_all_latbins_fewlays
end
%plot_all_latbins_chooseOutput; 

% To save everything
% save filename alld, or if include special vars
% save obs_no_ozone_chans  alld ...
%        co2 co2_sigs ...
%        o3 o3_sigs   ...
%        n2o n2o_sigs ...
%        ch4 ch4_sigs ...
%        cfc11 cfc11_sigs ...
%        sst sst_sigs
