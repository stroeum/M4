DATE=`date +%Y%m%d-%H%M`
MYDIR=`pwd`

find ${MYDIR} \( -name ".*.swp" -o -name ".nfs*"     \) -exec rm -rfv {} \;
find ${MYDIR} \( -name "main"   -o -name "convert" -o -name "*.o" \) -exec rm -rfv {} \;