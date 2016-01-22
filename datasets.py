##
## FUNCTION dataset_by_obname
##
## PURPOSE
##    get a dataset defintion (as astropy Table object) including DPIDs for a given dataset identified by OB name
##
## INPUT
##    ob_name
##
import sqlite3
import os
from astropy.io import ascii
from astropy.table import Table
import pdb

def dataset_by_obname(ob_name):
	dataset_definition = os.getenv("XDIR")+'/dataset_definition/'+ob_name+'.txt'
	d=ascii.read(dataset_definition)
	night=d['night'][0]
	ob_name_list=[d['object'][0],d['telluric'][0],d['flux'][0]]
	arms=["NIR","VIS","UVB"]
	dprlist=["SCI","TELL","FLUX"]
	
	dataset = Table(names=('dprtype', 'ob_name', 'arm', 'dpid'), dtype=('S10', 'S20', 'S3', 'S29'), meta={'night': night})

	conn=sqlite3.connect(os.getenv("OBSDB"))
	c=conn.cursor()
	
	for arm in arms:
		f_arm=dataset_definition.split('.txt')[0]+'_'+arm+'.txt'
		if os.path.isfile(f_arm):
			print("not yet coded...")
		else:
			with open(f_arm,'w') as f:
				for ob,dpr in zip(ob_name_list,dprlist):
					query = "select arcfile from shoot where night=\"" + night + \
						"\" and ob_name=\""+ ob + \
						"\" and arm=\"" + arm + \
						"\" and opti2_name=\"IFU\" limit 1;"
					c.execute(query)
					dpid=c.fetchone()
					dataset.add_row([dpr,ob,arm,dpid])
#					pdb.set_trace()
					file_string=ob+" "+dpid[0]+"\n"
					f.write(file_string)

	return(dataset)

##
## build QC pages by dataset
##
def QC_by_obname(ob_name):
	pass