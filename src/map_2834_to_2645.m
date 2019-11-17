function [f ichan] = map_2834_to_2645

hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
f = vchan2834;

load sarta_chans_for_l1c.mat
f = f(ichan);
