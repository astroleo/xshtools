##
## download_sort_xshoo.sh
##
##    sort XSHOOTER raw data
##
## OPTIONS
##    give any number of parameters (any text does it) and it will **not** look for the download file and not perform the download, but just sort the data
##
## CHANGELOG
##
##    2015-07-06   created from rescue_sort.sh (SINFOTOOLS)
##

set -e

if [[ $# -eq 0 ]]; then
	dl=true
else
	dl=false
fi

if [[ "$dl" == true ]]; then
	## confirm, execute and remove ESO download scripts
	downloaddir="$HOME/Downloads"
	cd $downloaddir
	
	esodlfile=`find * | egrep "^downloadRequest[0-9]{6}script.sh$"`
	nfiles=`echo $esodlfile | wc -w`
	if [ ! $nfiles == 1 ]; then
		echo "Found $nfiles ESO download files in $downloaddir."
		exit
	fi
	
	nfiles=$(cat $esodlfile | wc -l)
	
	echo $esodlfile
	cat $esodlfile
	
	read -p "Execute (and then delete) this download file [y/n]?" -n 1
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then exit; fi
fi

cd $XSHOODATAIN

if [[ "$dl" == true ]]; then
	sh $downloaddir/$esodlfile
fi

##
## let's not remove any non-FITS files for now, some may be necessary for pipeline...? (or perhaps contain other useful info)
##
#nfiles=$(find . -name "XSHOO.*.txt" | wc -l)
#if [[ ! $nfiles -eq 0 ]]; then
#	rm $XSHOODATAIN/XSHOO.*.txt
#fi
#
#nfiles=$(find . -name "M.XSHOOTER.*.fits" | wc -l)
#if [[ ! $nfiles -eq 0 ]]; then
#	rm $XSHOODATAIN/M.XSHOOTER.*.fits
#fi
#nfiles=$(find . -name "*.xml" | wc -l)
#if [[ ! $nfiles -eq 0 ]]; then
#	rm $XSHOODATAIN/*.xml
#fi

## if there are subdirectories move all files to main directory and remove subdirectory
subdirs=$(find . -type d -d 1)
for d in $subdirs; do
	mv $d/* .
	rmdir $d
done

echo "Now unzipping files..."
for zipfile in `ls *.Z`; do gunzip $zipfile; done

echo "Now moving/copying files to correct places..."
for f in `ls XSHOO.*.fits`
do
	dateobs=`echo $f | awk -F "." '{print $2}'`
	night=`whichnight_date.sh $dateobs`
	dir=$XSHOODATA/$night
	if [[ ! -d $dir ]]; then mkdir $dir; fi
	mv $f $dir
done

for f in `ls M.XSHOOTER.*.fits`
do
	dateobs=`echo $f | awk -F "." '{print $3}'`
	night=`whichnight_date.sh $dateobs`
	dir=$XSHOODATA/$night
	if [[ ! -d $dir ]]; then mkdir $dir; fi
	mv $f $dir
done

for f in `ls XSHOO.*.NL.txt`
do
	dateobs=`echo $f | awk -F "." '{print $2}'`
	night=`whichnight_date.sh $dateobs`
	dir=$XSHOODATA/$night
	if [[ ! -d $dir ]]; then mkdir $dir; fi
	mv $f $dir
done

for f in `ls XSHOO.*.xml`
do
	dateobs=`echo $f | awk -F "." '{print $2}'`
	night=`whichnight_date.sh $dateobs`
	dir=$XSHOODATA/$night
	if [[ ! -d $dir ]]; then mkdir $dir; fi
	mv $f $dir
done

if [[ "$dl" == true ]]; then
	cd $downloaddir
	if [[ -e $esodlfile ]]; then
		rm $esodlfile
		echo "removed ESO download file $downloaddir/$esodlfile"
	else
		echo "could not remove esodlfile: $esodlfile."
	fi
fi