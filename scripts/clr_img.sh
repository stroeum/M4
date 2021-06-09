DATE=`date +%Y%m%d-%H%M`
MYDIR=`pwd`

find ${MYDIR} \( -name ".*.swp" -o -name ".nfs*" \) -exec rm -rfv {} \;
find ${MYDIR} \( -name "*.?if*" -o -name "*.png" \) -exec rm -rfv {} \;
find ${MYDIR} \( -name "*.pdf"  -o -name "*.ps"  \) -exec rm -rfv {} \;
