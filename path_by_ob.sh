##
## associate OB and file names
##
ob_name="NGC7582_1"
dateobs_obj_vis=$(sqlite3 $OBSDB "select dateobs from shoot where ob_name=\"$ob_name\" and dprtech=\"ECHELLE,IFU,OFFSET\" and arm = \"VIS\" order by dateobs limit 1;")
night=$(sqlite3 $OBSDB "select night from shoot where ob_name=\"$ob_name\" and dprtech=\"ECHELLE,IFU,OFFSET\" and arm = \"VIS\" order by dateobs limit 1;")

##
## find all tellurics of the night; print ra,dec



