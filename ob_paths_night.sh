## night in, relevant file paths out...
night=$1
arms="UVB VIS NIR"

dir_reduced="/Users/leo/Projekte/LP-BAT/XSHOOTER/data/reflex_end_products/2015-11-12_pipeline_2-6-8"
dir_raw="/Volumes/astrodata2/XSHOODATA"

for arm in $arms; do

	obs=$(sqlite3 $OBSDB "select distinct ob_name from shoot where night=\"$night\" and arm=\"$arm\" and opti2_name=\"IFU\";")

	for ob in $obs; do
		if [[ $ob == "Calibration" ]]; then
			continue
		fi
	#	echo $ob
		arcfile=$(sqlite3 $OBSDB "select arcfile from shoot where night=\"$night\" and ob_name=\"$ob\" and arm=\"$arm\" and opti2_name=\"IFU\" limit 1;")
		dataset_id_eso=$(echo $arcfile | awk -F ".fits" '{print $1}')
		dataset_id_mac=$(echo $dataset_id_eso | sed 's/:/_/g')
	
	#	echo $dataset_id_eso
	#	echo $dataset_id_mac

		dprtype=$(sqlite3 $OBSDB "select dprtype from shoot where night=\"$night\" and ob_name=\"$ob\" and arm=\"$arm\" and opti2_name=\"IFU\" limit 1;")
	#	echo $dprtype
		if [[ $dprtype == "STD,TELLURIC" ]]; then
			dprstr="TELL"
		elif [[ $dprtype == "STD,FLUX" ]]; then
			dprstr="FLUX"
		elif [[ $dprtype == "OBJECT" ]]; then
			dprstr="SCI"
		else
			#echo "${dprtype}: ignored"
			continue
		fi
	
		dataset_raw="${dir_raw}/${night}/${dataset_id_eso}.fits"
		dataset_red="${dir_reduced}/${dataset_id_mac}_tpl/${ob}_${dprstr}_IFU_MERGE3D_DATA_OBJ_${arm}.fits"
		echo "$night $arm $dprstr $dataset_id_eso $dataset_raw $dataset_red"
	
	done
done