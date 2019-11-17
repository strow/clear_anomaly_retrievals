function ch = find_the_oem_channels(f,lenrates,settings_numchan,settings_chan_LW_SW,iChSet);

%iChSet = 1; %% till  August 20, 2019
%iChSet = 2; %% after August 20, 2019, uses CFC11/CFC12 chans
%iChSet = 3; %% after August 20, 2019, does not use CFC11 chans

if nargin == 4
  iChSet = 1;
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iNum = 2378;
iNum = 2645;
iNum = settings_numchan;

if iNum ~= 2378 & iNum ~= 2645
  settings_numchan
  error('incorrect number of chans');
end

%lenrates = length(driver.rateset.rates);
if iNum ~= lenrates
  fprintf(1,'you want to use good chans for %4i chans but your spectral rates are length %4i \n',iNum,lenrates);
  error('please fix');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iNum == 2378
  ch = choose_goodchans_from_2378;                      %% Old 2378 chans
else
  ch = choose_goodchans_from_2645(settings_chan_LW_SW); %% New 2645 chans
end

if iChSet == 1
  disp('find_the_oem_channels.m iChSet == 1 (old chans)')
  iless700 = find(f(ch) < 700);
    iless700 = iless700(1:5:length(iless700));
  imore700 = find(f(ch) >= 700 & f(ch) < 730);
    imore700 = imore700(1:5:length(imore700));
  imore730 = find(f(ch) >= 730);
  
  ch = ch([iless700; imore700; imore730]);
  ch = setdiff(ch,76);   %% 668.0348 cm-1 is a bad chan
    
%{
  driver.jacobian.chanset = ch;
  
  %plot(f,driver.rateset.rates,f(ch),driver.rateset.rates(ch),'o')
  
  %driver.jacobian.chanset = ch(1:10:length(ch));
  %driver.jacobian.chanset = ch(1:5:length(ch));
  
  chans38 = [        41          54         181         273         317         359         445         449 ...
                    532         758         903         904        1000        1020        1034        1055 ...
                   1075        1103        1249        1282        1291        1447        1475        1557 ...
                   1604        1614        1618        1660        1790        1866        1867        1868 ...
                   1878        1888        2112        2140        2321        2333];
  chans41 = [chans38 2325        2339        2353]; chans41 = sort(chans41);
  %driver.jacobian.chanset = chans41(chans41 < 2200);
%}
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif iChSet == 2
  %% after August 21, 2019
  %% see Channel_changes_for_anomaly_fits.txt
  disp('find_the_oem_channels.m iChSet == 2 (new chans)')

  iless700 = find(f(ch) < 700);
    iless700 = iless700(1:5:length(iless700));
  imore700 = find(f(ch) >= 700 & f(ch) < 730);
    imore700 = imore700(1:5:length(imore700));
  imore730 = find(f(ch) >= 730);
  
  ch = ch([iless700; imore700; imore730]);
  ch = setdiff(ch,76);   %% 668.0348 cm-1 is a bad chan

  addCFC11_weakWV = [634 635 637 1513 1520];
  rmBadChan = [289 647 648 649 650 652 654 655 660 663 664 669 671 673 674 676 682 685 689 690 ...
               1125 1152 1184 1194 ...
               1563 1586 1587 1588 1591 1592 1596 ...
               1619 1627 1633 1635 1636 1642 1643 1659 1694 1701 1702 2035 2012];
  ch = [ch; addCFC11_weakWV'];
  ch = setdiff(ch,rmBadChan);

elseif iChSet == 3
  %% after August 21, 2019
  %% see Channel_changes_for_anomaly_fits.txt
  disp('find_the_oem_channels.m iChSet == 3 (new chans w/o CFC)')

  iless700 = find(f(ch) < 700);
    iless700 = iless700(1:5:length(iless700));
  imore700 = find(f(ch) >= 700 & f(ch) < 730);
    imore700 = imore700(1:5:length(imore700));
  imore730 = find(f(ch) >= 730);
  
  ch = ch([iless700; imore700; imore730]);
  ch = setdiff(ch,76);   %% 668.0348 cm-1 is a bad chan

  addCFC11_weakWV = [634 635 637 1513 1520];
  rmBadChan = [289 647 648 649 650 652 654 655 660 663 664 669 671 673 674 676 682 685 689 690 ...
               1125 1152 1184 1194 ...
               1563 1586 1587 1588 1591 1592 1596 ...
               1619 1627 1633 1635 1636 1642 1643 1659 1694 1701 1702 2035 2012];
  ch = [ch; addCFC11_weakWV'];
  ch = setdiff(ch,rmBadChan);

  rmCFC11 = [[600:680] [1189:1280]];
  ch = setdiff(ch,rmCFC11);  

else
  iChSet
  error('iChSet = 1,2,3 only')
end
