import sys
import re
import urllib2
import math as m
import numpy as np
from colormap import bezier_interpolate


class PPIQuery(object):
    """Pull PPI for a gene from DB using the PSICQUIC REST API"""
    def __init__(self, url=None, gene=None, glist=[], 
                       regex=None, replace=None, gene_locs=[0, 0]):
        self.url = url
        self.gene = gene
        self.glist = glist
        self.regex = regex
        self.replace = replace
        self.gene_locs = gene_locs
    
    def try_url(self):
        url = self.url
        try:
            url_data = urllib2.urlopen(url.format(self.gene)).read()
            return url_data
        except urllib2.URLError:
            return []
        if not url_data:
            return []
    
    def make_table(self, url_data):
        return [t.split('\t') for t in url_data.splitlines()]
    
    def get_gene(self, entry=[]):
        finder = re.search(self.regex, entry)
        if finder:
            finder = finder.group()
            fgene = finder.replace(self.replace, '').upper()
            return fgene
        else:
            return ''
    
    def get_interactors(self):
        url_data = self.try_url()
        if not url_data:
            return []
        table = self.make_table(url_data)
        interactors = []
        for row in table:
            gene1 = self.get_gene(entry=row[self.gene_locs[0]])
            gene2 = self.get_gene(entry=row[self.gene_locs[1]])
            if gene1 and gene2 and gene1 != self.gene:
                interactors.append(gene1)
            else:
                interactors.append(gene2)
        return [i for i in interactors if i in set(self.glist)]


def make_gml(fname, (dcolor, mcolor, ucolor)):
    '''geppi.py writes a GML file using a tab-delimited
gene list file with expression values and two hex colors
that can form a bezier gradient with white in the middle.

Usage: python geppi.py <tab-delim filename> [hex dcolor] [hex ucolor]
Example: python geppi.py hedgehog.csv "#003366" "#AA0000"

The gene list file should look like this:
BRCA1    1.23
BRCA2    0.56

It's best to use complementary colors for "ucolor" and "dcolor". 
If you want to specify colors, please specify both colors.'''
    
    # open the file
    with open(fname, 'r') as f:
        csv = [l.split('\t') for l in f.read().splitlines()]
        cols = [list(r) for r in zip(*csv)]
        glist = cols[0]
        expv = map(float, cols[1])
    
    # center values at zero
    emax = max(abs(min(expv)), abs(max(expv)))
    
    print('Populating interactions dictionary...')
    gdict = {}
    for n, gene, val in zip(range(len(glist)), glist, expv):
        biogrid = PPIQuery(url='http://tyersrest.tyerslab.com:8805/psicquic/webservices/current/search/interactor/{0}?species=human&format=tab25', 
                    gene=gene, glist=glist, 
                    regex='gene/locuslink:[A-Z]+[A-Z0-9\-]+', 
                    replace='gene/locuslink:', gene_locs=[2, 3]
                    )
        spike = PPIQuery(url='http://spike.cs.tau.ac.il/psicquic-ws/webservices/current/search/interactor/{0}?species=human&format=tab25', 
                    gene=gene, glist=glist, 
                    regex='gene/locuslink:[A-Z]+[A-Z0-9\-]+', 
                    replace='gene/locuslink:', gene_locs=[4, 5]
                    )
        intact = PPIQuery(url='http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/{0}?species=human&format=tab25', 
                    gene=gene, glist=glist, 
                    regex='[A-Z]+[A-Za-z0-9\-]+\(gene name\)', 
                    replace='(gene name)', gene_locs=[4, 5]
                    )
        apid = PPIQuery(url='http://cicblade.dep.usal.es/psicquic-ws/webservices/current/search/interactor/{0}?species=human&format=tab25', 
                    gene=gene, glist=glist, 
                    regex='gene/locuslink:[A-Z]+[A-Z0-9\-]+', 
                    replace='gene/locuslink:', gene_locs=[4, 5]
                    )
        
        interactors = list(set(
                        biogrid.get_interactors() + 
                        spike.get_interactors() + 
                        intact.get_interactors() + 
                        apid.get_interactors()
                        ))
        
        if interactors:
            gdict[gene] = {'id': n, 'interactors': interactors, 'expval': val}
    
    # create grid for positions
    lg = len(glist)/10
    xa = np.arange(0, 2000, 200)
    ya = np.arange(0, (lg+1)*200, 200)
    xx, yy = np.meshgrid(xa, ya)
    xx = xx.ravel()
    yy = yy.ravel()
    
    print('Begin writing GML file...')
    
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
        t = (gd['expval'] + emax) / 2 * emax
        color = bezier_interpolate((dcolor, mcolor, ucolor), t=t)
        stream += '''node
[
    id    {0}
    label    \"{1}\"
    graphics
    [
        x    {2}
        y    {3}
        w    70
        h    30
        type    \"ellipse\"
        fill    \"{4}\"
        outline    \"#111111\"
    ]
]
'''.format(gd['id'], gene, x, y, color)
    
    # write each edge
    edges = set([])
    for gene in gdict.keys():
        gd = gdict[gene]
        for intg in gd['interactors']:
            try:
                source = gd['id']
                target = gdict[intg]['id']
            except KeyError:
                continue
            if (source, target) not in edges and (target, source) not in edges:
                edges.update([(source, target)])
                stream += '''edge
[
    source    {0}
    target    {1}
    label    \"pp\"
    graphics
    [
        width    2
        type    \"line\"
        fill    \"#333333\"
    ]
]
'''.format(source, target)
    
    # finish string
    stream += ']'
    
    # write to file
    with open(fname.split('.')[0]+'.gml', 'w') as f:
        f.write(stream)
    
    print('Finished.')
    return


if __name__ == '__main__':
    try:
        fname = sys.argv[1]
    except IndexError:
        print(make_gml.__doc__)
        sys.exit()
    
    mcolor = '#FFFFFF'
    
    try:
        dcolor, ucolor = sys.argv[2:4]
    except ValueError:
        dcolor, ucolor = '#003366', 'AA0000'
    
    make_gml(fname, (dcolor, mcolor, ucolor))





