FOLDER=`pwd`
NAME=`basename $FOLDER`

cd $FOLDER
echo "Convert miff to gif"
rm -rfv *.gif
for file in *.miff;
#do convert -compress none -antialias -contrast -render -enhance -density 72 -despeckle -quality 100 -transparent white -verbose -delay 50 "$file" "${file/.miff/.gif}";
do convert -compress none -enhance -render -antialias -contrast -quality 100 -delay 50 -verbose -monitor "$file" "${file/.miff/.gif}";
done;
