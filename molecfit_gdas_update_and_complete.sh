##
## molecfit looks for GDAS profiles in 3 hour bins, but sometimes, for example
##    in the range Oct-Dec 2013 (within 2004-2016), the data are only given in
##    six hour bins. In that case, molecfit (as of today, 2 Feb 2016) crashes.
##    The workaround is to look for such missing files and then copie a close-by file
##    to the expected 3-hour binned location.
##
##    This script downloads the latest tar ball from ESO, then applies the 
##    workaround and repacks the many small files into a new tar ball that is 
##    then symlinked to the correct position in the molecfit directory
##
##    It has only been tested on Mac OS X 10.11 "El Capitan".
##

## settings
gdas_dir=${HOME}/molecfit/data/profiles/gdas/

if ! [[ $(uname) == "Darwin" ]]; then
	echo "On non-Mac systems you need to adapt the date commands before using this script."
	exit
fi

cd $gdas_dir

if [ ! -d tmp ]; then mkdir tmp; fi
cd tmp

tarname="gdas_profiles_C-70.4-24.6.tar.gz"
wget "ftp://ftp.eso.org/pub/dfs/pipelines/skytools/molecfit/gdas/$tarname"
tar xzf $tarname
rm $tarname
##
## remove all unnecessary files
rm C-70.4-24.6D2004*
rm C-70.4-24.6D2005*
rm C-70.4-24.6D2006*
rm C-70.4-24.6D2007*
rm C-70.4-24.6D2008*
rm C-70.4-24.6D2009*
rm C-70.4-24.6D2010*
rm C-70.4-24.6D2011*
rm C-70.4-24.6D2012*

date_start="2013-12-01"
##
## find last file (for break condition)
lastf=$(ls -r | head -n 1)
echo "Most recent available file is ${lastf}."

d=$date_start

while true; do
	hours="03 06 09 12 15 18 21"
	f1="C-70.4-24.6D${d}T00.gdas"
	
	for h in $hours; do
		f2="C-70.4-24.6D${d}T${h}.gdas"

		if ! [[ -e $f1 ]]; then
			if [[ -e $f2 ]]; then
				cp $f2 $f1
				echo "copied $f2 to $f1"
			else
				echo "Oops. Neither $f1 nor $f2 exist."
				exit 1
			fi
		elif ! [[ -e $f2 ]]; then
			cp $f1 $f2
			echo "copied $f1 to $f2"
		fi		
		f1=$f2
	done

	if [[ "$f2" == "$lastf" ]]; then break; fi

	YYYY=$(echo $d | awk -F "-" '{print $1}')
	MM=$(echo $d | awk -F "-" '{print $2}')
	DD=$(echo $d | awk -F "-" '{print $3}')
	d=$(date -v${YYYY}y -v${MM}m -v${DD}d -v "+1d" "+%Y-%m-%d")
	echo $d
done

flist=filelist.txt
ftar_fake=gdas_profiles_C-70.4-24.6_fake_complete.tar
ftar_real=gdas_profiles_C-70.4-24.6.tar.gz

ls > $flist
tar -T $flist -czf $ftar_fake
mv $ftar_fake ..
cd ..
rm $ftar_real
ln -s $ftar_fake $ftar_real
rm -r tmp