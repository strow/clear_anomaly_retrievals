function m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,iiBin,i16daytimestep,iSwitchType);

if nargin == 3
  iSwitchType = 1; %% default, my original finite diff jacs of 10% pert, done almost all the time
end

if i16daytimestep < 0
  %% disp('use constant jac')
else
  fprintf(1,'replacing cfc11 jac with that at timestep %4i of 365 \n',i16daytimestep)
  lstrow = load('sarta_chans_for_l1c.mat');

  iOrigOrV2 = 1;  %% this is my original finite diff jacs of 10% pert, done almost all the time
  iOrigOrV2 = 2;  %% this is using Larrabee's idea of CFC11(t)-370 , test 6/24 underlying geo = average
  iOrigOrV2 = 3;  %% this is using Larrabee's idea of CFC11(t)-370 , test 7/01 underlying geo = time varying

  iOrigOrV2 = iSwitchType;

  iUgh = 4;   %% latest, no seasonal, without CFC11, works very well
  iUgh = 5;   %% latest, no seasonal, with CFC12
  if iOrigOrV2 == 3 
    disp('replacing CFC11 jac with time interpolation kcarta, iOrigOrV2 == 3')
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
    disp('replacing CFC11 jac with time interpolation sarta, iOrigOrV2 == 3')
    newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_Anomaly365_16_with_seasonal_OldSarta/RESULTS_FiniteDiff_Try3/']; %% with seasonal, old sarta
    newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_Anomaly365_16_with_seasonal_OldSarta/RESULTS_FiniteDiff_Try3/']; %% no seasonal, old sarta
    newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_Anomaly365_16/RESULTS_FiniteDiff_Try3/'];                        %% no seasonal, new sarta
    newjacname = [newjacname '/sarta_' num2str(i16daytimestep,'%03d') '_tracegas_finitediff_3gas_2645_V3.mat']; %% July 11,         better CO2/CH4/N2O prof, fixed 2002/09 yay
  elseif iOrigOrV2 == 1
    disp('replacing CFC11 jac with time interpolation, iOrigOrV2 == 1')
    newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CO2_370_385_400_415_12p8/co2_jac_2834_latbin' num2str(iiBin) '.mat'];
  elseif iOrigOrV2 == 2
    disp('replacing CFC11 jac with time interpolation, iOrigOrV2 == 2')
    newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CO2_370_385_400_415_12p8/co2_jac_2834_latbin' num2str(iiBin) 'v2.mat'];
  else
    error('huh????')
  end

  fprintf(1,' external cfc11jac = %s \n',newjacname);
  new = load(newjacname);
  if abs(iOrigOrV2) == 3
    new_cfc11_jac = squeeze(new.tracegas(iiBin,:,4));
  elseif iOrigOrV2 == 1 | iOrigOrV2 == 2

    new_cfc11_jac_alltime = new.qcx(lstrow.ichan,:);
    new_iaCFC11 = new.iaCFC11;

    yyjunk = 2000:2020;
    slope = (226-260)/(2020-2000); %% see https://www.esrl.noaa.gov/gmd/hats/combined/CFC11.html
    cfc11junk = 260 + (yyjunk-2000)*slope;
    new_timeCFC11 = interp1(cfc11junk,yyjunk,new_iaCFC11);       %% these are the jac times
    %  plot(yyjunk,cfc11junk,'.-',new_timeCFC11,new_iaCFC11,'o-')
    %  whos new_cfc11_jac_alltime new_timeCFC11

    oktimes = load('ok365times.mat');
    oktimes = oktimes.okdates(i16daytimestep);   %% need to get the CFC11 jac at this time!!!!

    for ii = 1 : 2645
      new_cfc11_jac(ii) = interp1(new_timeCFC11,new_cfc11_jac_alltime(ii,:),oktimes);
    end
  end

  %plot(1:2645,m_ts_jac_coljac(:,1),'k.-',1:2645,new_cfc11_jac,'r',1:2645,new_cfc11_jac_alltime(:,1),'b',1:2645,new_cfc11_jac_alltime(:,10),'c')
  m_ts_jac_coljac(:,4) = new_cfc11_jac;

end
