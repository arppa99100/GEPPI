# GEPPI
GEPPI stands for Genelist Extraction of Protein-Protein Interactions. GEPPI allows you to input a gene list with their corresponding expression levels and get out a GLM file that contains the specifications for all protein-protein interactions within the gene list along with formatting that allows for colorizing nodes based on log2 expression fold change. This GLM file can be directly imported into [Cytoscape](http://www.cytoscape.org/). The script `geppi.py` is a stand-alone that only uses modules from the Python Standard Library and a slightly modified [color gradient](http://bsou.io/posts/color-gradients-with-python) script that's included. It has been tested in Python 2.7.X, but I'll test on other versions if there's any interest.

## Installation and Usage
Put `geppi.py` and `colormap.py` into your active folder. To see the usage docstring without coming to this README all the time, just write `python geppi.py` and this will print to screen:


```
Usage: python geppi.py <tab-delim filename> [hex ucolor] [hex dcolor]
Example: python geppi.py hedgehog.csv "#CC3300" "#006699"

The gene list file should look like this:
BRCA1	1.23
BRCA2	0.56

"ucolor" and "dcolor" are going to be blended in the
middle, so it's best to use complementary colors. But
if you want to specify colors, please specify both colors.
```

To make it work, you only need to input the file you want to process. The two colors are optional but must be enclosed by quotation marks in order to be processed correctly. I've picked orange for upregulation and blue for downregulation as default colors.

## Importing into Cytoscape
What you end up with is a `filename.glm` file that can be imported and further customized by [Cytoscape](http://www.cytoscape.org/). Please follow the link to download and install if you don't have it already. Open Cytoscape and you'll notice that the first thing it does is ask if you would like to create an empty network, import from file, etc. Click the import option and load your GLM file. Click OK on the pop-up.

![Hedgehog Network](https://raw.githubusercontent.com/pam-bot/GEPPI/master/hedgehogGrid.png)
**Fig. 1** The [Sonic Hedgehog pathway](https://en.wikipedia.org/wiki/Sonic_hedgehog) genes and their interactions. Gene expression fold changes were randomly generated for this example. 

What you should see are colorized ellipses connected by lines indicating interactions. The shade of the colors indicate degrees of differential expression. For now, they're arranged in a grid. That's because it takes a lot more graph theory to arrange nodes by their connectivity than by a grid in the order they occur, and it's a step that is much easier to do in Cytoscape.

Here's a [tutorial](http://opentutorials.cgl.ucsf.edu/index.php/Tutorial:Introduction_to_Cytoscape#Laying_Out_Your_Network) on the many kinds of layouts available in Cytoscape. Whichever you choose depends on your needs, but I'll demonstrate the Edge-Weighted Spring Embedded Layout (it may be called something else in another version, but they're all Force-directed Layouts in general). You can see what each layout on default settings give by selecting them in the `Layout` tab.

Go to `Layout > Settings` and select the Edge-Weighted Spring Embedded Layout. A bunch of settings will pop up. Increasing the spring constant makes the graph more spread out, and so does increasing the spring resting length. I doubled both the defaults. I found I had a lot of collisions, so I increased the strength to avoid collisions to 100 (originally set at 0). The result looked pretty good to me:

![Hedgehog Network](https://raw.githubusercontent.com/pam-bot/GEPPI/master/hedgehogGraph.png)
**Fig. 2** The network on the adjusted Edge-Weighted Spring Embedded Layout.


