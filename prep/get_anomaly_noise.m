%% /home/strow/Work/Airs/Stability/get_anomaly_noise.m
addpath /asl/matlib/aslutil

load nedt_in_l1c_array
load_fairs

deriv = drdbt(fairs,250);
nen = deriv.*nedt;

load ~/Matlab/Stats/mlat_equal_area
maxilat = 40;

fstart = 'Data/Desc/statlat';
for latid = 1:maxilat
   latid
   g = load([fstart int2str(latid)],'count','robs');
   allc = g.count(:,1);
   robs = g.robs;
   clear g;
   % 16-day averages
   
   for i=1:362
     ii = (((i-1)*16)+1):(((i-1)*16)+16);
     robs_m = nanmean(robs(ii,:));
     count_m = nansum(allc(ii));
     btm = rad2bt(fairs,robs_m')';
     deriv = drdbt(fairs,btm);
     btn = nen./(deriv);
     if count_m < 10
        btn_avg(latid,i,:) = NaN(2645,1);
     else
        btn_avg(latid,i,:) = btn./sqrt(count_m);
     end
   end   
end

% btn_avag is 40 x 362 by 2645
% save btn_avg.mat btn_avg
