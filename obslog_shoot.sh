##
## sourcenights
##
## PURPOSE
## Query SQLite database for all nights in which a specified object has been observed
##
##
## set locale to use point as decimal symbol
export LC_NUMERIC="C"

dbfile="$OBSDB"
nobj=`sqlite3 $dbfile "select count(*) from sources where id='$1';"`
if [ ! $nobj -eq 1 ]; then
#	echo "No such source (or too many)."
	exit
fi

obj_name=`sqlite3 $dbfile "select name from sources where id='$1';"`
ra=`sqlite3 $dbfile "select ra from sources where id='$1';"`
dec=`sqlite3 $dbfile "select dec from sources where id='$1';"`
##
## define 10 arcsec search box, i.e. 0.0028
searchrad=0.0028
ra_min=$(echo "$ra $searchrad" | awk '{printf "%f", $1 - $2}')
ra_max=$(echo "$ra $searchrad" | awk '{printf "%f", $1 + $2}')
dec_min=$(echo "$dec $searchrad" | awk '{printf "%f", $1 - $2}')
dec_max=$(echo "$dec $searchrad" | awk '{printf "%f", $1 + $2}')
##
## get all nights in which object has been observed

nights=`sqlite3 $dbfile "select distinct night from shoot where ra between $ra_min and $ra_max and dec between $dec_min and $dec_max;"`

for night in $nights; do
	##
	## select distinct OBs
	OBs=`sqlite3 $dbfile "select distinct ob_name from shoot where ra between $ra_min and $ra_max and dec between $dec_min and $dec_max and night='$night';"`
	for ob in $OBs; do
		airmass=`sqlite3 $dbfile "select avg(airm_start) from shoot where ra between $ra_min and $ra_max and dec between $dec_min and $dec_max and night='$night' and ob_name = '$ob';"`
		seeing=`sqlite3 $dbfile "select avg(FWHM_start) from shoot where ra between $ra_min and $ra_max and dec between $dec_min and $dec_max and night='$night' and ob_name = '$ob';"`				
		echo "$1,$obj_name,$ob,$night,$airmass,$seeing"
	done
done