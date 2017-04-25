##
## determine combined mask that masks all regions defined in a list of input files
##
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii
import glob

## from http://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
def merge_intervals(intervals):
	sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
	merged=[]
	
	for higher in sorted_by_lower_bound:
		if not merged:
			merged.append(higher)
		else:
			lower = merged[-1]
			if higher[0] <= lower[1]:
				upper_bound = max(lower[1], higher[1])
#				merged[-1] = (lower[0], upper_bound)
				merged[-1] = [lower[0], upper_bound]
			else:
				merged.append(higher)
	return merged

masks=glob.glob("*.sm")
allmasks=[]
y=1

for ix,m in enumerate(masks):
	a = ascii.read(m,data_start=1)
	ms = a["col1"].data
	me = a["col2"].data
	for ims,ime in zip(ms,me):
		plt.plot([ims,ime],[y,y],'k',linewidth=1)
	plt.text(10000,y+0.1,m.split("_")[1],fontsize=9,ha="right")
	y+=1

	onemask = np.transpose([ms.tolist(),me.tolist()]).tolist()
	allmasks.extend(onemask)

merged_mask = merge_intervals(allmasks)
for interval in merged_mask:
	plt.plot(interval,[y,y],'k',linewidth=5)
plt.text(10000,y+0.1,"combined",fontsize=9,ha="right")

plt.plot([5007,5007],[0,10],'k-',linewidth=0.5)
plt.plot([6563,6563],[0,10],'k-',linewidth=0.5)
	
plt.ylim([0,y+1])
plt.xlim([3500,10500])
plt.savefig("masks.pdf")
plt.clf()