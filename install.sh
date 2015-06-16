BASEDIR=$(dirname $0)
clear
echo 'build packages from source...'
echo '****************************************'
R CMD INSTALL --no-multiarch --with-keep.source $BASEDIR/
echo '****************************************'
echo ''

cat <<EOF
Do you want to install topOnto.HDO.db package and run an example for topOnto?(y/n)
EOF
read ans
[ -z ans ] && ans = 'y'
if [ $ans = 'y' ]; then
echo 'installing topOnto.HDO.db ...'
#(cd $BASEDIR/inst/scripts/create.db/ && sh run.sh)
R CMD INSTALL --no-multiarch --with-keep.source $BASEDIR/inst/scripts/create.db/topOnto.HDO.db/
echo '****************************************'
	Rscript $BASEDIR/inst/scripts/exampleRun.R
echo '****************************************'
fi
echo 'done!'
echo ''
