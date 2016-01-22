## LB, 2016-01-19
## Python script to get redshift for a galaxy
## adapted from $SINFOTOOLS/redshift_line.py

from astroquery.ned import Ned

def get_redshift(object_name):
	Q=Ned.query_object(object_name)
	z=Q['Redshift'][0]
	return(z)
