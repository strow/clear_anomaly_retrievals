addpath /home/sergio/MATLABCODE/oem_pkg

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS.m
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Select the latitude bin

%JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
JOB = 03;  %% debug
JOB = 20;  %% debug
%JOB = 38;  %% debug

ix = JOB;
fprintf(1,'JOB = %2i \n',ix);

driver.iibin = JOB;
%---------------------------------------------------------------------------
% Need oem_pkg
%---------------------------------------------------------------------------
% Doing debug?
driver.debug = false;
driver.debug_dir = '../Debug';

% Open debug file if desired
if driver.debug
  writelog('open');
end;
%---------------------------------------------------------------------------
% Perform OEM fit?
driver.oem.dofit = true;
driver.lls.dofit = false;

% Oem loops?  Just one if linear.
driver.oem.nloop = 1;

% Debug plots inside rodgers?
driver.oem.doplots = false;

driver.i16daytimestep = -1; %% this is NOT anomaly

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Override many settings and add covariance matrix
topts = struct;


topts.numchan = 2645;
topts.chan_LW_SW = 0;   %% just LW/MW
%topts.chan_LW_SW = -1;  %% all chans

%topts.offsetrates = +1;  %% add in offset of 0.01 K/yr
%topts.ocb_set = +1;      %% cal
%topts.addco2jacs = +1;   %% add co2 = 2.2 ppmv/yr to cal

topts.obs_corr_matrix = -1; %% use diagnol obs uncertainty 
%topts.obs_corr_matrix = +1; %% use full obs cov matrix

topts.tie_sst_lowestlayer = -1;  %% no  cross cov between SST and lowest layer
topts.tie_sst_lowestlayer = +1;  %% yes cross cov between SST and lowest layer

topts.invtype = 1; %% pinv, DEFAULT
%topts.invtype = 3; %% inverse, Texas A&M
%topts.invtype = 4; %% inverse, ridge regression
%topts.invtype = 5; %% inverse, minimum eigenvalue

topts.iNlays_retrieve = 97;
topts.iNlays_retrieve = 30;
topts.iNlays_retrieve = 10;
topts.iNlays_retrieve = 5;
topts.iNlays_retrieve = 20;
if topts.iNlays_retrieve >= 97
  [driver,aux] = strow_override_defaults_latbins_AIRS(driver,topts);
else
  [driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,topts.iNlays_retrieve,topts);
end

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Do the retrieval
driver = retrieval(driver,aux);
%---------------------------------------------------------------------------
% Save retrieval output
driver.filename = ['Output/test' int2str(driver.iibin)];
save(driver.filename,'-struct','driver');
%---------------------------------------------------------------------------
% Close debug file
if driver.debug
  writelog('close')
end
%---------------------------------------------------------------------------
% Some simple output
fprintf('Scalar Retrievals from OEM\n')
fprintf(1,'CO2   (ppm)   %5.3f  +- %5.3f \n',driver.oem.finalrates(1),driver.oem.finalsigs(1));
fprintf(1,'N2O   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(2),driver.oem.finalsigs(2));
fprintf(1,'CH4   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(3),driver.oem.finalsigs(3));
fprintf(1,'CFC11 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(4),driver.oem.finalsigs(4));
fprintf(1,'SST   (K)    %5.3f  +- %5.3f \n',driver.oem.finalrates(5),driver.oem.finalsigs(5));
%---------------------------------------------------------------------------
co2(JOB) = driver.oem.finalrates(1);
%---------------------------------------------------------------------------
% Plot Results
addpath Plotutils
if topts.iNlays_retrieve >= 97
  plot_retrieval_latbins
else
  plot_retrieval_latbins_fewlays
end

% for i=1:4;clf(i);end
%---------------------------------------------------------------------------

