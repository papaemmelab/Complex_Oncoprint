This is to create oncoprints of variants (mutations [required], copy number variants [optional], structural variants [optional]).
Look at the example file for an example of how to call the function. 
Feb 2019: features added to show RESPONSE, and SAMPLE.SOURCE (for example you can choose B-cell, T-cell, or you can even pass on the disease type such as PV, ET, etc.). Also, for timepoint data, you have the option to show a heatbar of samples that belong to the same patient. 
Sep 2019: new features added to allow defining the second annotation ribbon as wish (e.g., RACE or DISEASE). Other visualization added and will be explained in an analysis session.
You can install/upgrade using:

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
