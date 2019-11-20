while true; do  printf `date "+%Y/%m/%d:::%H:%M:%S"`; printf ' count so far :  '; ls OutputAnomaly_OBS/*/*.mat | wc -l ; ls OutputAnomaly_CAL/*/*.mat | wc -l; sleep 15; done
