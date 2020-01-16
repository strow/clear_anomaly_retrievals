function m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,iiBin,i16daytimestep,iSwitchType);

if nargin == 3
  iSwitchType = 1; %% default, my original finite diff jacs of 10% pert, done almost all the time
end

if i16daytimestep < 0
  %% disp('use constant jac')
  else
  fprintf(1,'replacing ch4 jac with that at timestep %4i of 365 \n',i16daytimestep)
  lstrow = load('sarta_chans_for_l1c.mat');

  iOrigOrV2 = 1;  %% this is my original finite diff jacs of 10% pert, done almost all the time
  iOrigOrV2 = 2;  %% this is using Larrabee's idea of CO2(t)-370 , test 6/24
  iOrigOrV2 = 3;  %% this is using Larrabee's idea of CO2(t)-370 , test 7/01 underlying geo = time varying
  iOrigOrV2 = 4;  %% this is using Larrabee's idea of CO2(t)-370 , test 12/08 underlying geo = time varying, has Age of Air

  iOrigOrV2 = iSwitchType;

  iUgh = 4;   %% latest, no seasonal, without CFC11, works very well
  iUgh = 5;   %% latest, no seasonal, with CFC12
  if iOrigOrV2 == 4  
    disp('replacing CH4 jac with time interpolation kcarta, iOrigOrV2 == 4, Age of AIR')
    if iUgh == 5
      newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS_FiniteDiff_Try4/'];       %% no seasonal, redone
      newjacname = [newjacname '/kcarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_5_2645_V5.mat']; %% Dec 08,         better CO2/CH4/N2O/CFC11/CFC12 prof, fixed 2002/09 yay plus Age of AIR
    end
  elseif iOrigOrV2 == 3 
    disp('replacing CH4 jac with time interpolation kcarta, iOrigOrV2 == 3')
    if iUgh == 1
      newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS_FiniteDiff/'];
      newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS_FiniteDiff_Try1/'];
      newjacname = [newjacname '/kcarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_4_2645.mat'];   %% June 30/July 1, bad CO2/CH4/N2O prof
    elseif iUgh == 2
      newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS_FiniteDiff_Try2/'];
      newjacname = [newjacname '/kcarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_4_2645_V2.mat']; %% July 2,         better CO2/CH4/N2O prof, screwed 2002/09 ugh
    elseif iUgh == 3
      newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8_July12_2019_Great_But_with_Seasonal/RESULTS_FiniteDiff_Try3/'];  %% have seasonal    
      newjacname = [newjacname '/kcarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_4_2645_V3.mat']; %% July 11,        better CO2/CH4/N2O/CFC11       prof, fixed 2002/09 yay
    elseif iUgh == 4
      newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8_tillAug25_2019/RESULTS_FiniteDiff_Try3_noCFC12/']; %% no seasonal
      newjacname = [newjacname '/kcarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_4_2645_V3.mat']; %% July 11,        better CO2/CH4/N2O/CFC11       prof, fixed 2002/09 yay
    elseif iUgh == 5
      newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS_FiniteDiff_Try3/']; %% no seasonal
      newjacname = [newjacname '/kcarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_5_2645_V4.mat']; %% Aug 21,         better CO2/CH4/N2O/CFC11/CFC12 prof, fixed 2002/09 yay
    end
  elseif iOrigOrV2 == -3 
    disp('replacing CH4 jac with time interpolation sarta, iOrigOrV2 == 3')
    newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_Anomaly365_16/RESULTS_FiniteDiff_Try3/']; %% no seasonal
    newjacname = [newjacname '/sarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_3gas_2645_V3.mat']; %% July 11,         better CO2/CH4/N2O prof, fixed 2002/09 yay
  elseif iOrigOrV2 == 1
    disp('replacing CH4 jac with time interpolation, iOrigOrV2 == 1')
    newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CH4_1700_2000/ch4_jac_2834_latbin' num2str(iiBin) '.mat'];
  elseif iOrigOrV2 == 2
    disp('replacing CH4 jac with time interpolation, iOrigOrV2 == 2')
    newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CH4_1700_2000/ch4_jac_2834_latbin' num2str(iiBin) 'v2.mat'];
  else
    error('huh????')
  end

  new = load(newjacname);
  if abs(iOrigOrV2) == 3 | iOrigOrV2 == 4
    new_ch4_jac = squeeze(new.tracegas(iiBin,:,3));
  elseif iOrigOrV2 == 1 | iOrigOrV2 == 2
    new_ch4_jac_alltime = new.qcx(lstrow.ichan,:);
    new_iaCH4 = new.iaCH4;

    ch4_esrl = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/ch4_esrl.txt');
    new_timeCH4 = interp1(ch4_esrl(:,2)/1000,ch4_esrl(:,1),new_iaCH4,[],'extrap');
    % plot(ch4_esrl(:,1),ch4_esrl(:,2)/1000,'.-',new_timeCH4,new_iaCH4,'o-'); grid
    % whos new_ch4_jac_alltime new_timeCH4

%{
>> load anomaly_0dayavg_results_spectra_cal_constantCO2jac.mat
>> whos
  Name              Size               Bytes  Class     Attributes

  chanset         532x1                 4256  double
  iaTropics         1x20                 160  double
  okdates           1x365               2920  double
  okrtime           1x365               2920  double
  raaCal         2645x365            7723400  double
  raaObs         2645x365            7723400  double

>> save ok365times okdates okrtime
%}

    oktimes = load('ok365times.mat');
    oktimes = oktimes.okdates(i16daytimestep);   %% need to get the CH4 jac at this time!!!!

    for ii = 1 : 2645
      new_ch4_jac(ii) = interp1(new_timeCH4,new_ch4_jac_alltime(ii,:),oktimes);
    end
  end

  %plot(1:2645,m_ts_jac_coljac(:,3),'k.-',1:2645,new_ch4_jac,'r',1:2645,new_ch4_jac_alltime(:,1),'b',1:2645,new_ch4_jac_alltime(:,13),'c')
  m_ts_jac_coljac(:,3) = new_ch4_jac;
end
