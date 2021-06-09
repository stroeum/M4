DATE=`date +%Y%m%d-%H%M`
MYDIR=`pwd`

find ${MYDIR}         \( -name ".*.swp" -o -name ".nfs*"     \) -exec rm -rfv {} \;
find ${MYDIR}         \( -name "*.?if*" -o -name "*.png"     \) -exec rm -rfv {} \;
find ${MYDIR}/viz_dir \( -name "*.dat"  -o -name "*.general" \) -exec rm -rfv {} \;
find ${MYDIR}/output  \( -name "*.dat" 		                 \) -exec rm -rfv {} \;