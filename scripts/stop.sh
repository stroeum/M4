FOLDER=/nv/hp5/jriousset6/data/
NAME=$(basename $PWD)
for file in ~/log$NAME.e*; 
do 
  JOBID=`echo ${file} | sed s/^.*\.e//`;
  echo $JOBID
  canceljob $JOBID
  rm -rvf $file ${file/.e$JOBID/.o$JOBID}
done;
