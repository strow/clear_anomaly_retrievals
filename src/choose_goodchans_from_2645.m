function [ch,fairs2645] = choose_goodchans_from_2645(chan_LW_SW)

%%% topts.numchan = 2645;
%%% topts.chan_LW_SW =  0;  %% just LW/MW DEFAULT
%%% topts.chan_LW_SW = -1;  %% all chans,LW/MW and SW <<<<<<<
%%% topts.chan_LW_SW = +1;  %% LW/MW but avoid deep 15 um 
%%% topts.chan_LW_SW = -2;  %% SW only
%%% topts.chan_LW_SW = +2;  %% LW/MW/SW but avoid deep 15 um

% chan_LW_SW       0 is 640 to 1640 (default), 1 is 700 to 1640, -1 is 640 to 2740
if nargin == 0
  chan_LW_SW = 1;
end

load indices_of_l1b_in_l1c.mat
%  l1b_ind_in_l1c       2314x1             18512  double
%  l1c_ind_for_l1b      2314x1             18512  double

%% get list of good L1b chans
iAB = 1;   %% use kg_fixed_ab  %% most of the time
iAB = 0;   %% use kg_ab0       %% 3/27/19 test

bind = load('good_chans_2016.mat');
if iAB == 1
  g_l1b = bind.kg_fixed_ab;       %% most of the time
elseif iAB == 0
  g_l1b = bind.kg_ab0;            %% 3/27/19 test
end
clear bind

[c ia ib] = intersect(g_l1b,l1b_ind_in_l1c);
g_l1c = l1c_ind_for_l1b(ib);
ch    = g_l1c;   %%  l1c indices of the good l1b channels
%whos ch

%%% new
good = 1:length(ch);

if iAB == 1
  load kbad_chans   %% from 640 to 1640 cm-1
  good = setdiff(1:length(ch),kbad);
elseif iAB == 0
  load kbad_chans0   %% from 640 to 1640 cm-1
  good = setdiff(1:length(ch),kbad0);
end

ch = ch(good);
%whos ch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
fairs = vchan2834;
load sarta_chans_for_l1c.mat
fairs2645 = fairs(ichan);

%%%%%%%%%%%%%%%%%%%%%%%%%
if chan_LW_SW == 0
  disp('use only LW/MW channels in choose_goodchans_from_2645.m')
  ahaLW = find(fairs2645 <= 1640);
  ch = intersect(ch,ahaLW);

  newbad_n2o_list_dec1_2019 = [1621 1622 1624 1644];
    ch = setdiff(ch,newbad_n2o_list_dec1_2019);
  no_more_cfc11 = find(fairs2645 >= 843 & fairs2645 <= 849);
  no_more_cfc11 = find((fairs2645 >= 840 & fairs2645 <= 856) | (fairs2645 >= 1070 & fairs2645 <= 1090));
  no_more_cfc11 = [];  %% find_the_oem_channels.m --- use iChSet == 3
    ch = setdiff(ch,no_more_cfc11);

%%%%%%%%%%%%%%%%%%%%%%%%%
elseif chan_LW_SW == +1
  disp('use LW/MW away from 700 cm-1 in choose_goodchans_from_2645.m')
  ahaLW = find(fairs2645 >= 690 & fairs2645 <= 1640);
  ahaLW = find(fairs2645 >= 700 & fairs2645 <= 1640);
  ch = intersect(ch,ahaLW);

  newbad_n2o_list_dec1_2019 = [1621 1622 1624 1644];
    ch = setdiff(ch,newbad_n2o_list_dec1_2019);
  no_more_cfc11 = find(fairs2645 >= 843 & fairs2645 <= 849);
  no_more_cfc11 = find((fairs2645 >= 840 & fairs2645 <= 856) | (fairs2645 >= 1070 & fairs2645 <= 1090));
  no_more_cfc11 = [];  %% find_the_oem_channels.m --- use iChSet == 3
    ch = setdiff(ch,no_more_cfc11);

%%%%%%%%%%%%%%%%%%%%%%%%%
elseif chan_LW_SW == -1
  disp('use LW and SW chans in choose_goodchans_from_2645.m');

  ahaLW = find(fairs2645 <= 1640);
  chLW = intersect(ch,ahaLW);

  newbad_n2o_list_dec1_2019 = [1621 1622 1624 1644];
    ch = setdiff(ch,newbad_n2o_list_dec1_2019);
  no_more_cfc11 = find(fairs2645 >= 843 & fairs2645 <= 849);
  no_more_cfc11 = find((fairs2645 >= 840 & fairs2645 <= 856) | (fairs2645 >= 1070 & fairs2645 <= 1090));
  no_more_cfc11 = [];  %% find_the_oem_channels.m --- use iChSet == 3
    ch = setdiff(ch,no_more_cfc11);

%  %% ahaSWhigh = find(fairs2645 >= 2310 & fairs2645 <= 2370);
%  junk = load('stratSW.mat');
%  ahaSWhigh = junk.iStratSW;
%  %ch = setdiff(ch,ahaSWhigh);
  ahaSWhigh1 = find(fairs2645 >= 2310 & fairs2645 <= 2370);
  ahaSWhigh1 = find(fairs2645 >= 1910 & fairs2645 <= 2870);
  ahaSWhigh1 = find(fairs2645 >= 2230 & fairs2645 <= 2770);  %% avoid CO chans!!!
  ahaSWhigh1 = find((fairs2645 >= 2230 & fairs2645 <= 2260) | (fairs2645 >= 2385 & fairs2645 <= 2760));  %% avoid CO and high alt CO2 chans!!!
  junk = load('stratSW.mat');
  ahaSWhigh2 = junk.iStratSW;
  chSW = setdiff(ahaSWhigh1,ahaSWhigh2);
  chacha = find(fairs2645(chSW) > 1910);  %% somehow a few LW chans still get in
  chSW = chSW(chacha);

  ch = union(chLW,chSW)

%%%%%%%%%%%%%%%%%%%%%%%%%
elseif chan_LW_SW == -2
  disp('use SW chans only in choose_goodchans_from_2645.m');
  ahaSWhigh1 = find(fairs2645 >= 2310 & fairs2645 <= 2370);
  ahaSWhigh1 = find(fairs2645 >= 1910 & fairs2645 <= 2870);
  ahaSWhigh1 = find(fairs2645 >= 2230 & fairs2645 <= 2770);  %% avoid CO chans!!!
  ahaSWhigh1 = find((fairs2645 >= 2230 & fairs2645 <= 2260) | (fairs2645 >= 2385 & fairs2645 <= 2760));  %% avoid CO and high alt CO2 chans!!!
  junk = load('stratSW.mat');
  ahaSWhigh2 = junk.iStratSW;
  ch = setdiff(ahaSWhigh1,ahaSWhigh2);
  chacha = find(fairs2645(ch) > 1910);  %% somehow a few LW chans still get in
  %ch = ch(6:end);                      %% somehow a few LW chans still get in
  ch = ch(chacha);

%%%%%%%%%%%%%%%%%%%%%%%%%
elseif chan_LW_SW == +2
  disp('use LW/MW away from 700 cm-1, plus SW chans, in choose_goodchans_from_2645.m')

  ahaLW = find(fairs2645 <= 1640);
  ahaLW = find(fairs2645 >= 700 & fairs2645 <= 1640);

  chLW = intersect(ch,ahaLW);

  newbad_n2o_list_dec1_2019 = [1621 1622 1624 1644];
    ch = setdiff(ch,newbad_n2o_list_dec1_2019);
  no_more_cfc11 = find(fairs2645 >= 843 & fairs2645 <= 849);
  no_more_cfc11 = find((fairs2645 >= 840 & fairs2645 <= 856) | (fairs2645 >= 1070 & fairs2645 <= 1090));
  no_more_cfc11 = [];  %% find_the_oem_channels.m --- use iChSet == 3
    ch = setdiff(ch,no_more_cfc11);

%  %% ahaSWhigh = find(fairs2645 >= 2310 & fairs2645 <= 2370);
%  junk = load('stratSW.mat');
%  ahaSWhigh = junk.iStratSW;
%  %ch = setdiff(ch,ahaSWhigh);
  ahaSWhigh1 = find(fairs2645 >= 2310 & fairs2645 <= 2370);
  ahaSWhigh1 = find(fairs2645 >= 1910 & fairs2645 <= 2870);
  ahaSWhigh1 = find(fairs2645 >= 2230 & fairs2645 <= 2770);  %% avoid CO chans!!!
  ahaSWhigh1 = find((fairs2645 >= 2230 & fairs2645 <= 2260) | (fairs2645 >= 2385 & fairs2645 <= 2760));  %% avoid CO and high alt CO2 chans!!!
  junk = load('stratSW.mat');
  ahaSWhigh2 = junk.iStratSW;
  chSW = setdiff(ahaSWhigh1,ahaSWhigh2);

  chacha = find(fairs2645(chSW) > 1910);  %% somehow a few LW chans still get in
  chSW = chSW(chacha);

  ch = union(chLW,chSW)

else
  error('unknown option')  
end

%clf; plot(fairs2645(ch),'+-'); grid; title(['good channels = ' num2str(length(ch))]); disp('ret'); pause
%keyboard_nowindow
%error('lkgsjlkjgs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
1) Suppose you have an l1c file (freqs, rads), then the “real” channels are
load rads_l1c                            %% rads_l1c = 1x2645  = l1c chans, some of which are real, others are fake
rads_real = rads_l1c(l1c_ind_for_l1b);   %% rad_real  = 2314, not 2378 wow

2) if you have a rads file with l1b info rads_l1b = 1x2378  and you want them “IN” a l1c array
  rads_l1c = zeros(1,2645);
  rads_l1c(l1c_ind_for_l1b) = rads_l1b(l1b_ind_in_l1c)
  these indices are all 2314 long, the number of l1b channels in l1c

3) My most recent *l1B” good channels, use kg_fixed_ab in good_chans_2016.mat

4) To get good chans in L1c

g_l1b = channel id’s for good l1b chans
[c ia ib] = intersect(g_l1b,l1b_ind_in_l1c);
g_l1c = l1c_ind_for_l1b(ib)
g_l1c are the l1c indices of the good l1b channels
%}

