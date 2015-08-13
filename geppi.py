import sys
import re
import urllib2
import math as m
import numpy as np
from colormap import gradient_interpolate


#######################################
# define rest api data mining functions

def BioGRID(gene):
	'''Pull PPI for a gene from BioGRID using the PSICQUIC REST API'''
	try:
		tab = urllib2.urlopen('http://tyersrest.tyerslab.com:8805/psicquic/webservices/current/search/interactor/{0}?species=human&format=tab25'.format(gene)).read()
	except urllib2.URLError:
		return []
	if not tab:
		return []
	table = [t.split('\t') for t in tab.splitlines()]
	def getGene(entry):
		genes = []
		if ':' in entry:
			for db in entry.split('|'):
				genes.append(db.split(':')[1])
			return genes[0]
		else:
			return ''
	inters = []
	for row in table:
		gene1 = getGene(row[2])
		gene2 = getGene(row[3])
		if gene1 and gene2 and gene1 != gene:
			inters.append(gene1)
		else:
			inters.append(gene2)
	return [i for i in inters if i in glist]


def SPIKE(gene):
	'''Pull PPI for a gene from SPIKE using the PSICQUIC REST API'''
	try:
		tab = urllib2.urlopen('http://spike.cs.tau.ac.il/psicquic-ws/webservices/current/search/interactor/{0}?species=human&format=tab25'.format(gene)).read()
	except urllib2.URLError:
		return []
	if not tab:
		return []
	table = [t.split('\t') for t in tab.splitlines()]
	def getGene(entry):
		finder = re.search(':[A-Za-z0-9\-]+\(gene name\)', entry)
		if finder:
			finder = finder.group()
			fgene = finder[1:].replace('(gene name)', '').upper()
			return fgene
		else:
			return ''
	inters = []
	for row in table:
		gene1 = getGene(row[2])
		gene2 = getGene(row[3])
		if gene1 and gene2 and gene1 != gene:
			inters.append(gene1)
		else:
			inters.append(gene2)
	return [i for i in inters if i in glist]


def IntAct(gene):
	'''Pull PPI for a gene from IntAct using the PSICQUIC REST API'''
	try:
		tab = urllib2.urlopen('http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/{0}?species=human&format=tab25'.format(gene)).read()
	except urllib2.URLError:
		return []
	if not tab:
		return []
	table = [t.split('\t') for t in tab.splitlines()]
	def getGene(entry):
		finder = re.search(':[A-Za-z0-9\-]+\(gene name\)', entry)
		if finder:
			finder = finder.group()
			fgene = finder[1:].replace('(gene name)', '').upper()
			return fgene
		else:
			return ''
	inters = []
	for row in table:
		gene1 = getGene(row[4])
		gene2 = getGene(row[5])
		if gene1 and gene2 and gene1 != gene:
			inters.append(gene1)
		else:
			inters.append(gene2)
	return [i for i in inters if i in glist]


def APID(gene):
	'''Pull PPI for a gene from APID using the PSICQUIC REST API'''
	try:
		tab = urllib2.urlopen('http://cicblade.dep.usal.es/psicquic-ws/webservices/current/search/interactor/{0}?species=human&format=tab25'.format(gene)).read()
	except urllib2.URLError:
		return []
	if not tab:
		return []
	table = [t.split('\t') for t in tab.splitlines()]
	def getGene(entry):
		genes = []
		if ':' in entry:
			for db in entry.split('|'):
				genes.append(db.split(':')[1])
			return genes[0]
		else:
			return ''
	inters = []
	for row in table:
		gene1 = getGene(row[4])
		gene2 = getGene(row[5])
		if gene1 and gene2 and gene1 != gene:
			inters.append(gene1)
		else:
			inters.append(gene2)
		return [i for i in inters if i in glist]


fname, ucolor, dcolor = sys.argv[1:4]
#fname, ucolor, dcolor = 'hedgehog.csv', '#CC3300', '#006699'

with open(fname, 'r') as f:
	csv = [l.split('\t') for l in f.read().splitlines()]
	cols = [list(r) for r in zip(*csv)]
	glist = cols[0]
	expv = map(float, cols[1])

emin, emax = min(expv), max(expv)

gdict = {}
for n, gene, val in zip(range(len(glist)), glist, expv):
	interactors = list(set(BioGRID(gene)+SPIKE(gene)+IntAct(gene)+APID(gene)))
	if interactors:
		gdict[gene] = {'id': n, 'interactors': interactors, 'expval': val}


###########
# write gml

# create grid for positions
lg = len(glist)/10
xa = np.arange(0, 2000, 200)
ya = np.arange(0, (lg+1)*200, 200)
xx, yy = np.meshgrid(xa, ya)
xx = xx.ravel()
yy = yy.ravel()

# begin write string
stream = '''graph
[
'''

# write each node
c = 0
for gene in gdict.keys():
	gd = gdict[gene]
	x = xx[c]
	y = yy[c]
	c += 1
	t = (gd['expval']-emin)/(emax-emin)
	color = gradient_interpolate(dcolor, finish_hex=ucolor, t=t)
	stream += '''node
[
	id	{0}
	label	\"{1}\"
	graphics
	[
		x	{2}
		y	{3}
		w	70
		h	30
		type	\"ellipse\"
		fill	\"{4}\"
		outline	\"#111111\"
	]
]
'''.format(gd['id'], gene, x, y, color)

# write edges
for gene in gdict.keys():
	gd = gdict[gene]
	for intg in gd['interactors']:
		source = gd['id']
		target = gdict[intg]['id']
		stream += '''edge
[
	source	{0}
	target	{1}
	label	\"pp\"
	graphics
	[
		width	2
		type	\"line\"
		fill	\"#333333\"
	]
]
'''.format(source, target)

# finish string
stream += ']'

with open(fname.split('.')[0]+'.gml', 'w') as f:
	f.write(stream)

