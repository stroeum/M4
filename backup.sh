DATE=`date +%Y%m%d-%H%M`
FOLDER=`pwd`
NAME=`basename $CURRENT`
cd $FOLDER
E_BADARGS=2
if [ ! -n "$1" ]
then
  echo "Saving visualizations in ${HOME}/data/BAK/ as M4_RK_${DATE}.tar.gz"
  tar cvzf "${HOME}/data/BAK/M4_RK_${DATE}.tar.gz" viz_dir/*.*if* viz_dir/*.png viz_dir/M4*
else
  echo "Saving visualizations in ${HOME}/data/BAK/ as M4_RK_${DATE}_$1.tar.gz"
  cd viz_dir
  rm -rfv *.gif
  for file in *.miff;
    do mv "$file" "${file/.miff/_$1.miff}";
  done;
  for file in *.gif;
    do mv "$file" "${file/.gif/_$1.gif}";
  done;
  for file in *.png;
    do mv "$file" "${file/.png/_$1.png}";
  done;
  cd $FOLDER
  tar cvzf "${HOME}/data/BAK/M4_RK_${DATE}_$1.tar.gz" viz_dir/*.*if* viz_dir/*.png viz_dir/M4*
fi

#
