% Build covariance matrices; fixed, water, temperature
pmat_size = driver.jacobian.numlays;

% Normalization depends on parameter
fnorm = driver.qrenorm(driver.jacobian.scalar_i);
wnorm = driver.qrenorm(driver.jacobian.water_i);
tnorm = driver.qrenorm(driver.jacobian.temp_i);
oznorm = driver.qrenorm(driver.jacobian.ozone_i);

% Make sure re-scale to Jacobians before squaring cov matrix
for i=1:pmat_size
  for j=1:pmat_size
    mat_od(i,j) = abs(i-j);
  end
end

% Defines tropopause index trpi
trop_index
if iNlays_retrieve < 90
  trpi = floor(trpi/(100/iNlays_retrieve));
end

%%%%%%%%%%%%%%%%%%%%
%% see ../Strow/strow_override_defaults_latbins_AIRS.m

%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLABCODE/oem_pkg_run/AIRS_AllSky_97_O3jacs_Apr2017/strow_build_cov_matrix.m
%%         lc   ct.lev1  ct.lev2   ct_wide   cw.lev1  cw.lev2    cw_wide alpha_T  alpha_w  alpha_oz
%%               TROP    STRAT               TROP     STRAT
%% if lc is tiny, mat_od ~ diagnol; if lc is large, mat_od is very broad
%%
%% if cX_wide is tiny, transition is very smooth and gradual ... if cX_wide is large, almost step function transition from trop to strat
%% reg ~ Scov + S_tikonov ~ (sig_q)^-1 + alpha L1'L1       where sig_q is uncertainty in param (eg K^2 for temperature)
%%   if (sig_q)^-1 = 1/(uncertainty)^2 >> alpha ==> covariance regularization
%%      (sig_q)^-1 = 1/(uncertainty)^2 << alpha ==> tikonov    regularization
%%   if sig_q -> 0   then you say you are VERY sure about a-priori ==> do not change ==> delta(param) --> 0
%%      sig_q -> INF then you say you are DO NOT TRUST    a-priori ==>        change ==> delta(param) --> bigly wigly
%%   if alpha -> 0   then you say you are DO NOT TRUST    a-priori ==>        change ==> delta(param) --> bigly wigly
%%      alpha -> INF then you say you are VERY sure about a-priori ==> do not change ==> delta(param) --> 0
%%
%%               sigT_t    sigT_s                sigWV_t   sigWV_s             sigO3_t   sigO3_s
%%         lc   ct.lev1  ct.lev2   ct_wide     cw.lev1  cw.lev2    cw_wide  coz.lev1  coz.lev2    coz_wide  alpha_T  alpha_w  alpha_oz

%%         lc   ct.lev1   ct.lev2   ct_wide   cw.lev1  cw.lev2    cw_wide  coz.lev1  coz.lev2 coz_wide  alpha_T  alpha_w  alpha_oz

%%% topts.obs_corr_matrix = -1, 10,20 lays, topts.invtype = 1,3
cov_set = [1.0  0.005/2   0.005/2   1/2       0.005/25     0.005/25   1/2      0.001/2   0.001/0.75     1/2        1E-5     1E-5  1E-5]; %pretty good for obs  YEAH YEAH YEAH
cov_set = [1.0  0.05/2    0.05/2    1/2       0.05/25      0.05/25    1/2      0.005/2   0.005/0.75     1/2        1E-5     1E-5  1E-5]; %works pretty well for topts.obs_corr_matrix = -1;
cov_set = [1.0  0.05/2    0.05/2    1/2       0.01/25      0.01/25    1/2      0.05/2    0.05/0.75      1/2        1E-1     1E-1  1E-1]; %works pretty well for topts.obs_corr_matrix = -1  and 10 lays, water not very wiggly
cov_set = [1.0  0.05/2    0.05/2    1/2       0.1/25       0.1/25     1/2      0.05/2    0.05/0.75      1/2        1E-1     1E-1  1E-1]; %works pretty well for topts.obs_corr_matrix = -1  and 10 lays, OBS AND ERA CALCS good fits (great for ERA calcs T and O3 need to slightly improve WV make it slightly more wiggly)
cov_set = [1.0  0.05/2    0.05/2    1/2       0.5/25       0.5/25     1/2      0.05/2    0.05/0.75      1/2        1E-1     1E-1  1E-1]; %works pretty well for topts.obs_corr_matrix = -1  and 10 lays, OBS AND ERA CALCS good fits (great for ERA calcs T and O3 need to slightly improve WV make it slightly less wiggly)
cov_set = [1.0  0.05/2    0.05/2    1/2       0.15/25      0.15/25    1/2      0.05/2    0.05/0.75      1/2        1E-1     1E-1  1E-1]; %works pretty well for topts.obs_corr_matrix = -1  and 10 lays, OBS AND ERA CALCS good fits (great for ERA calcs T and O3 need to slightly improve WV make it slightly less wiggly)

%{
%%% topts.obs_corr_matrix = +1, 10 lays, topts.invtype = 4
%cov_set = [1.0  1e-6      1e-6      1/2       1e-11         1e-11       1/2      1e-6      1e-6           1/2        1E+16     1E+16  1E+16]; %works pretty well for topts.obs_corr_matrix = +1  and 10 lays, OBS AND ERA CALCS good fits (great for ERA calcs T and O3 need to slightly improve WV make it slightly more wiggly)
%cov_set = [1.0  5e-12     5e-12     1/2       5e-13         1e-12       1/2      1e-12     1e-12          1/2        2e+11     5e+16  2e+11]; %works pretty well for topts.obs_corr_matrix = +1  and 10 lays, OBS AND ERA CALCS good fits (great for ERA calcs T and O3 need to slightly improve WV make it slightly more wiggly)
%cov_set = [1.0  5e-12     5e-12     1/2       2e-12         2e-12       1/2      2e-12     2e-12        1/2        1e+1      1e+1   2e+1 ]; %works pretty well for topts.obs_corr_matrix = +1  and 10 lays, ERA CALCS fits 750 chans
cov_set = [1.0  5e-5      5e+1      1/2       1e-9          1e-9        1/2      1e-4      1e-4         1/2        1e+2      1e+12   2e+2 ]; %works pretty well for topts.obs_corr_matrix = +1  and 20 lays, ERA CALCS fits 700 chans

%%% topts.obs_corr_matrix = +1, 20 lays, topts.invtype = 5, CALCS, iVers = 1 in oem_pkg/inverse_minimum_eigenvalue_matrix.m
cov_set = [1.0  5e-6      5e-6      1/2       1e-6          1e-6        1/2      1e-6      1e-6         1/2        1e+2      1e+2   2e+2 ]; %works pretty well for topts.obs_corr_matrix = +1  and 20 lays, ERA CALCS fits 700 chans
cov_set = [1.0  4e-5      4e-5      1/2       5e-6          5e-6        1/2      5e-6      5e-6         1/2        1e+2      1e+0   1e+0 ]; %works really well for topts.obs_corr_matrix = +1  and 20 lays, ERA CALCS fits 700 chans
cov_set = [1.0  8e-5      8e-5      1/2       1e-5          1e-5        1/2      5e-5      5e-5         1/2        1e+2      1e+0   1e+0 ]; %works really well for topts.obs_corr_matrix = +1  and 20 lays, ERA CALCS fits 700 chans

%%% topts.obs_corr_matrix = +1, 20 lays, topts.invtype = 5, CALCS, iVers = 2 in oem_pkg/inverse_minimum_eigenvalue_matrix.m
cov_set = [1.0  1e-1     1e-1      1/2       1e-1          1e-1        1/2      1e-1      1e-1         1/2        1e+0      1e+0   1e+0 ]; %works really well for topts.obs_corr_matrix = +1  and 20 lays, ERA CALCS fits 700 chans
cov_set = [1.0  1e-1     1e-1      1/2       1e-1          1e-1        1/2      1e-1      1e-1         1/2        1e+2      1e+2   1e+2 ]; %works really well for topts.obs_corr_matrix = +1  and 20 lays, ERA CALCS fits 700 chans

%%%% topts.obs_corr_matrix = +1, 20 lays, topts.invtype = 5, OBS
%cov_set = [1.0  1e-2      1e-2      1/2       1e-4          1e-4        1/2      1e-4      1e-4         1/2        1e+4      1e+5   1e+2 ]; %works really well for topts.obs_corr_matrix = +1  and 20 lays, OBS fits 700 chans
%cov_set = [1.0  1e-2      1e-2      1/2       1e-3          1e-3        1/2      1e-2      1e-2         1/2        1e+4      1e+5   1e+2 ]; %works really well for topts.obs_corr_matrix = +1  and 20 lays, OBS fits 700 chans
%}

fprintf(1,'cov_set lc = %8.6f \n',cov_set(1));
junk = cov_set([2 3 4  11]); fprintf(1,'      T  : sig_trop  sig_strat cwide alpha = %8.6e %8.6e %8.6e %8.6e \n',junk)
junk = cov_set([5 6 7  12]); fprintf(1,'      WV : sig_trop  sig_strat cwide alpha = %8.6e %8.6e %8.6e %8.6e \n',junk);
junk = cov_set([8 9 10 13]); fprintf(1,'      O3 : sig_trop  sig_strat cwide alpha = %8.6e %8.6e %8.6e %8.6e \n',junk)

%%%%%%%%%%%%%%%%%%%%%%%%%

% Relative off-diagonal
l_c = cov_set(1);;
mat_od = exp(-mat_od.^2./(1*l_c^2));

   ct(ix).trans1 = trpi(ix);
   ct(ix).trans2 = trpi(ix)+10;
   ct(ix).lev1 = cov_set(2);
   ct(ix).lev2 = cov_set(3);
   ct(ix).lev3 = ct(ix).lev2;
   ct(ix).width1 = cov_set(4);
   ct(ix).width2 = ct(ix).width1;

   cw(ix).trans1 = trpi(ix);
   cw(ix).trans2 = trpi(ix)+10;
   cw(ix).trans1 = cw(ix).trans1 + 10;  %% new move towards groud
   cw(ix).trans2 = cw(ix).trans2 + 10;  %% new move towards ground
   cw(ix).lev1 = cov_set(5);
   cw(ix).lev2 = cov_set(6);
   cw(ix).lev3 = cw(ix).lev2;
   cw(ix).width1 = cov_set(7);
   cw(ix).width2 = cw(ix).width1;

   coz(ix).trans1 = trpi(ix);
   coz(ix).trans2 = trpi(ix)+10;
   coz(ix).trans1 = coz(ix).trans1 - 40; %% new move towards TOA
   coz(ix).trans2 = coz(ix).trans2 - 40; %% new move towards TOA
   coz(ix).lev1 = cov_set(8);
   coz(ix).lev2 = cov_set(9);
   coz(ix).lev3 = coz(ix).lev2;
   coz(ix).width1 = cov_set(10);
   coz(ix).width2 = coz(ix).width1;

% Temperature level uncertainties, then scaled and squared
if ~exist('iFixTz_NoFit')
  tunc = cov2lev(ct(ix),driver.jacobian.numlays);
  t_sigma = (tunc./tnorm);
  tmat = (t_sigma'*t_sigma).*mat_od;
  driver.oem.tunc = tunc;
elseif exist('iFixTz_NoFit','var') & ~strfind(driver.rateset.ocb_set,'obs')
  tunc = cov2lev(ct(ix),driver.jacobian.numlays);
  t_sigma = (tunc./tnorm);
  tmat = (t_sigma'*t_sigma).*mat_od;
  driver.oem.tunc = tunc;
elseif exist('iFixTz_NoFit','var') & strfind(driver.rateset.ocb_set,'obs')
  %% from strow_override_defaults_latbins_AIRS_fewlays.m
  if iFixTz_NoFit > 0
    disp('making sure T is not fitted!!!!')

    %driver.oem.tunc = ones(size(tunc)) * 1.0e-16;  cov_set(11) = 1.0e10;  %% THIS COMBO WORKS, CO2, CFC,N2O, CH4 good but still see residual dT/dt
    %driver.oem.tunc = ones(size(tunc)) * 1.0e-16;  cov_set(11) = 1.0e16;  %% THIS COMBO WORKS, much smaller dT/dt, CO2 good, CFC,O3,WV bad
    %driver.oem.tunc = ones(size(tunc)) * 1.0e-16;  cov_set(11) = 1.0e12;  %% 

    driver.oem.tunc = []; cov_set(11) = 1.0e10;
  end
end

% Water level uncertainties, then scaled and squared
wunc = cov2lev(cw(ix),driver.jacobian.numlays);
w_sigma = (wunc./wnorm);
wmat = (w_sigma'*w_sigma).*mat_od;
driver.oem.wunc = wunc;

% ozone level uncertainties, then scaled and squared
ozunc = cov2lev(coz(ix),driver.jacobian.numlays);
oz_sigma = (ozunc./oznorm);
ozmat = (oz_sigma'*oz_sigma).*mat_od;
driver.oem.ozunc = ozunc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scalar uncertainties
%fmat definitions below
%            CO2(ppm) N2O(ppb) CH4(ppb) CFC11(ppt) CFC12 (ppt) Tsurf(K)    cng1 cng2 cpsz1 cpsz2]   
if settings.co2lays == 1
  fmatd = [2     0.1       2      1     1            0.1/3]*3;   %% scanty test time ONLY done Tue Jul 16 15:45:49 EDT 2019
  fmatd = [2     0.1       2      0.25  0.1          0.1];       %% try this, CFC --> -1 ppm/yr, pretty good CFC!!      : files 8/30/19  *newchans_newkcartajacs_finitediffV3_5_tweakCFCunc.mat 
  fmatd = [2     1.0       5      2     2            0.1];       %% try loosening  so all gases have unc ~ 1.0 year growth : files 8/30/19  *newchans_newkcartajacs_finitediffV3_5_allTraceGasUnc_1yearrate.mat
  fmatd = [2     1.0       5      2     2            0.1]/10;    %% try tightening so all gases have unc ~ 0.1 year growth : files 9/03/19  *newchans_newkcartajacs_finitediffV3_5_allTraceGasUnc_0p1yearrate.mat
  fmatd = [2     0.1       2      1     1            0.1];       %% 99.99999999999999999999% of the time, works pretty well but CFC11,CFC12 rates are too high by x2  -3.4 ppm/yr (July/Aug 2019) >>> ORIG, BEST
elseif settings.co2lays == 3
  fmatd = [2 2 2 0.1       2      1     1            0.1];
end

if settings.set_tracegas == 1 & settings.co2lays == 1 & driver.i16daytimestep < 0
  fmatd(1:5) = fmatd(1:5)*0.0000001;   %% have put in xb(1:3) so we should not change those values
elseif settings.set_tracegas == 1 & settings.co2lays == 3 & driver.i16daytimestep < 0
  fmatd(1:7) = fmatd(1:7)*0.0000001;   %% have put in xb(1:3) so we should not change those values
end

fmat  = diag(fmatd./fnorm); 
if exist('tmat','var')
  driver.oem.cov = blkdiag(fmat,wmat,tmat,ozmat);
else
  driver.oem.cov = blkdiag(fmat,wmat,ozmat);
end

if topts.tie_sst_lowestlayer > 0 & exist('tmat','var')
  %% now tie together surface temp with lowest layers
  wah1 = driver.oem.cov(driver.jacobian.scalar_i(end),driver.jacobian.scalar_i(end));
  wah2 = driver.oem.cov(driver.jacobian.temp_i(end),driver.jacobian.temp_i(end));
  %driver.oem.cov(driver.jacobian.scalar_i(end),driver.jacobian.temp_i(end-1):driver.jacobian.temp_i(end)) = sqrt(abs(wah1*wah2));
  %driver.oem.cov(driver.jacobian.temp_i(end-1):driver.jacobian.temp_i(end),driver.jacobian.scalar_i(end)) = sqrt(abs(wah1*wah2));
  driver.oem.cov(driver.jacobian.scalar_i(end),driver.jacobian.temp_i(end):driver.jacobian.temp_i(end)) = sqrt(abs(wah1*wah2));
  driver.oem.cov(driver.jacobian.temp_i(end):driver.jacobian.temp_i(end),driver.jacobian.scalar_i(end)) = sqrt(abs(wah1*wah2));
end

%---------------------------------------------------------------------
% Empirical regularization parameters and switches
driver.oem.reg_type = 'reg_and_cov'; % 'reg_and_cov','cov','reg' are other choices

% Separate reg weights for water, temperature profiles
if ~exist('iFixTz_NoFit')
  driver.oem.alpha_temp =  cov_set(11);
elseif exist('iFixTz_NoFit','var') & ~strfind(driver.rateset.ocb_set,'obs')
  driver.oem.alpha_temp =  cov_set(11);
elseif exist('iFixTz_NoFit','var') & strfind(driver.rateset.ocb_set,'obs')
  disp('oops iFixTz_NoFit exits, for driver.rateset.ocb_set = obs so not settiong driver.oem.alpha_temp')
end
driver.oem.alpha_water = cov_set(12);
driver.oem.alpha_ozone = cov_set(13);
