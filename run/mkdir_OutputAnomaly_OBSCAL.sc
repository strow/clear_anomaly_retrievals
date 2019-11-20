echo " "

echo "rming any slurm*.out files"
/bin/rm slurm*.out

echo " "
if [ -d "OutputAnomaly_OBS" ] 
then
  echo "Directory OutputAnomaly_OBS exists, deleting files"
  /bin/rm -R OutputAnomaly_OBS/*/*.mat
else
  mkdir OutputAnomaly_OBS
  mkdir OutputAnomaly_OBS/01
  mkdir OutputAnomaly_OBS/02
  mkdir OutputAnomaly_OBS/03
  mkdir OutputAnomaly_OBS/04
  mkdir OutputAnomaly_OBS/05
  mkdir OutputAnomaly_OBS/06
  mkdir OutputAnomaly_OBS/07
  mkdir OutputAnomaly_OBS/08
  mkdir OutputAnomaly_OBS/09
  mkdir OutputAnomaly_OBS/10
  mkdir OutputAnomaly_OBS/11
  mkdir OutputAnomaly_OBS/12
  mkdir OutputAnomaly_OBS/13
  mkdir OutputAnomaly_OBS/14
  mkdir OutputAnomaly_OBS/15
  mkdir OutputAnomaly_OBS/16
  mkdir OutputAnomaly_OBS/17
  mkdir OutputAnomaly_OBS/18
  mkdir OutputAnomaly_OBS/19
  mkdir OutputAnomaly_OBS/20
  mkdir OutputAnomaly_OBS/21
  mkdir OutputAnomaly_OBS/22
  mkdir OutputAnomaly_OBS/23
  mkdir OutputAnomaly_OBS/24
  mkdir OutputAnomaly_OBS/25
  mkdir OutputAnomaly_OBS/26
  mkdir OutputAnomaly_OBS/27
  mkdir OutputAnomaly_OBS/28
  mkdir OutputAnomaly_OBS/29
  mkdir OutputAnomaly_OBS/30
  mkdir OutputAnomaly_OBS/31
  mkdir OutputAnomaly_OBS/32
  mkdir OutputAnomaly_OBS/33
  mkdir OutputAnomaly_OBS/34
  mkdir OutputAnomaly_OBS/35
  mkdir OutputAnomaly_OBS/36
  mkdir OutputAnomaly_OBS/37
  mkdir OutputAnomaly_OBS/38
  mkdir OutputAnomaly_OBS/39
  mkdir OutputAnomaly_OBS/40
fi

echo " "
if [ -d "OutputAnomaly_CAL" ] 
then
  echo "Directory OutputAnomaly_CAL exists, deleting files"
  /bin/rm -R OutputAnomaly_CAL/*/*.mat
else
  mkdir OutputAnomaly_CAL
  mkdir OutputAnomaly_CAL/01
  mkdir OutputAnomaly_CAL/02
  mkdir OutputAnomaly_CAL/03
  mkdir OutputAnomaly_CAL/04
  mkdir OutputAnomaly_CAL/05
  mkdir OutputAnomaly_CAL/06
  mkdir OutputAnomaly_CAL/07
  mkdir OutputAnomaly_CAL/08
  mkdir OutputAnomaly_CAL/09
  mkdir OutputAnomaly_CAL/10
  mkdir OutputAnomaly_CAL/11
  mkdir OutputAnomaly_CAL/12
  mkdir OutputAnomaly_CAL/13
  mkdir OutputAnomaly_CAL/14
  mkdir OutputAnomaly_CAL/15
  mkdir OutputAnomaly_CAL/16
  mkdir OutputAnomaly_CAL/17
  mkdir OutputAnomaly_CAL/18
  mkdir OutputAnomaly_CAL/19
  mkdir OutputAnomaly_CAL/20
  mkdir OutputAnomaly_CAL/21
  mkdir OutputAnomaly_CAL/22
  mkdir OutputAnomaly_CAL/23
  mkdir OutputAnomaly_CAL/24
  mkdir OutputAnomaly_CAL/25
  mkdir OutputAnomaly_CAL/26
  mkdir OutputAnomaly_CAL/27
  mkdir OutputAnomaly_CAL/28
  mkdir OutputAnomaly_CAL/29
  mkdir OutputAnomaly_CAL/30
  mkdir OutputAnomaly_CAL/31
  mkdir OutputAnomaly_CAL/32
  mkdir OutputAnomaly_CAL/33
  mkdir OutputAnomaly_CAL/34
  mkdir OutputAnomaly_CAL/35
  mkdir OutputAnomaly_CAL/36
  mkdir OutputAnomaly_CAL/37
  mkdir OutputAnomaly_CAL/38
  mkdir OutputAnomaly_CAL/39
  mkdir OutputAnomaly_CAL/40
fi

echo " "