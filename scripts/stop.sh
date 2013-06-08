FOLDER=/nv/hp5/jriousset6/data/
NAME=$(notdir $(CURDIR))
for file in ~/log$NAME.e*; 
do 
  JOBID=`echo ${file} | sed s/^.*\.e//`;
  canceljob $JOBID
  rm -rvf $file ${file/.e$JOBID/.o$JOBID}
done;
