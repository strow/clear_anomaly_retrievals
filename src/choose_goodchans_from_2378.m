function [ch,fairs2378] = choose_goodchans_from_2378()

%load /asl/s1/rates/Clear/good_chanset.mat 
%driver.jacobian.chanset = chanset;
%aha = find(isfinite(driver.rateset.rates));
%driver.jacobian.chanset = intersect(aha,chanset);

xyz = load('good_chans.mat');
xyz.woo = xyz.woo([1:274 276:end]);  % remove channel 312, 740 qbranch

load fitchans_noch4
xyz.woo = intersect(xyz.woo,k);

% SW only
%xyz.woo = xyz.woo([1:274 276:1034]);  % remove channel 312, 740 qbranch

ch = xyz.woo;
% Switch to l1c channels
%load kc; load kb;
%ch = kb(xyz.woo);

% Get rid of shortwave
ch = ch(1:918);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fairs2378 = instr_chans;
