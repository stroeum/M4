FOLDER=/nv/hp5/jriousset6/data/
for file in ~/log*.e*; 
do 
  JOBID=`echo ${file} | sed s/^.*\.e//`;
  canceljob $JOBID
  rm -rvf $file ${file/.e$JOBID/.o$JOBID}
done;
