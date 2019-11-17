function [ch,fairs2645] = choose_goodchans_from_2645(chan_LW_SW)

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
  ahaLW = find(fairs2645 <= 1640);
  ch = intersect(ch,ahaLW);
elseif chan_LW_SW == 1
  ahaLW = find(fairs2645 >= 690 & fairs2645 <= 1640);
  ch = intersect(ch,ahaLW);
elseif chan_LW_SW == -1
  ch = ch;
else
  error('unknown option')  
end

%clf; plot(fairs2645(ch),'+-'); grid; title(['good channels = ' num2str(length(ch))]); disp('ret'); pause
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

