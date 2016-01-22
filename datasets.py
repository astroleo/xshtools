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
import numpy as np

import pdb

class dataset:
	##
	## build dataset from previously saved files or from queries to $OBSDB if files do not yet exist
	##
	def __init__(self, ob_name):
		self.ob_name = ob_name
		self.dataset_definition = os.getenv("XDIR")+'/dataset_definition/'+self.ob_name+'.txt'
		d=ascii.read(self.dataset_definition)
		self.night=d['night'][0]
		ob_name_list=[d['object'][0],d['telluric'][0],d['flux'][0]]
		self.arms=["NIR","VIS","UVB"]
		self.dprlist=["SCI","TELL","FLUX"]
	
#		self.dataset = Table(names=('dprtype', 'ob_name', 'arm', 'dpid'), dtype=('a10', 'a20', 'a3', 'a29'), meta={'night': self.night})
		self.dataset = Table(names=('dprtype', 'ob_name', 'arm', 'dpid'), dtype=(np.dtype((str,10)), np.dtype((str,20)), np.dtype((str,3)), np.dtype((str,29))), meta={'night': self.night})

		conn=sqlite3.connect(os.getenv("OBSDB"))
		c=conn.cursor()
	
		for arm in self.arms:
			f_arm=self.dataset_definition.split('.txt')[0]+'_'+arm+'.txt'
			if os.path.isfile(f_arm):
				d_arm=ascii.read(f_arm,data_start=0,names=["ob","dpid"])
				for ob,dpid,dpr in zip(d_arm["ob"],d_arm["dpid"],self.dprlist):
					self.dataset.add_row([dpr,ob,arm,dpid])
			else:
				with open(f_arm,'w') as f:
					for ob,dpr in zip(ob_name_list,self.dprlist):
						query = "select arcfile from shoot where night=\"" + self.night + \
							"\" and ob_name=\""+ ob + \
							"\" and arm=\"" + arm + \
							"\" and opti2_name=\"IFU\" limit 1;"
						c.execute(query)
						dpid=c.fetchone()
						self.dataset.add_row([dpr,ob,arm,dpid])
						file_string=ob+" "+dpid[0]+"\n"
						f.write(file_string)
	
	##
	## build QC pages by dataset
	##
	def make_QC_page(self):
		html_QC = os.getenv("XDIR")+'/QC/datasets/'+self.ob_name+'.html'
		if os.path.isfile(html_QC):
			print(html_QC + ' exists.')
		else:
			with open(html_QC,'w') as html:
				html.write("<html>\n")
				html.write("<header>\n")
				html.write("<title>"+self.ob_name+"</title>\n")
				html.write("</header>\n")
				html.write("<body>\n")
				
				for arm in self.arms:
					html.write("<h1>"+arm+"</h1>")
					for row in self.dataset[np.where(self.dataset['arm']==arm)]:
						dpr = row['dprtype']
						ob = row['ob_name']
						dpid=row['dpid']

						html.write("<h2>" + dpr + ": " + ob + "</h2>\n")
						img="../plot/" + self.night + "/" + dpid.replace(":","_") + "_" + dpr + "_" + arm + ".png"
						html.write("<img src=\"" + img + "\">\n")
						html.write("<br>\n")
						html.write("<hr>\n")

				html.write("</body>\n")
				html.write("</html>\n")