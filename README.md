# RetinaChoroid_SarsCoV-2_Kuenzeletal2021
To the scientific community:
This repository contains the relevant R codes for the analysis of human retina and choroid scRNA seq datasets. The original publication can be found here ### (to be filled in, when article goes online).

We thank Madhvi Menon, Shahin Mohammadi, Jose Davila-Velderrain et al. from Broad Institute of MIT and Harvard, as well as Andrew P. Voigt et al. from the Insitute for Vision Research from the University of Iowa, for providing the rich transcriptomic datasets to the scientific community. Raw and processed data of the human retina and choroid have previously been deposited through GEO (gene expression omnibus, NCBI) under the accession number GSE137537 (license: CC-BY 4.0) and GSE135922 (license: CC-BY-NC-ND 4.0, with permission of corresponding author Prof. Mullins). The raw gene counts for each dataset were downloaded from GEO database (https://www.ncbi.nlm.nih.gov/geo/). 

We reprocessed the data from raw quantification matrix following standard Seurat (v.3.1) clustering procedure. We included cells which expressed at least 100 features, and features, that could be detected in at least 3 cells. To take advantage of the improved Seurat pre-processing and normalization workflow, we used the “SCTransform” function, thereby we “corrected” log-normalized expression values across datasets. To reannotate the cells, multiple clusterings of different resolutions were generated among which the one best matching published clustering was picked and manual annotation was undertaken using marker genes (compare Suppl. Figure 1).

If you have any questions or require assistance with similar projects, please do not hesitate to contact the corresponding author Steffen E. Künzel (steffen-emil.kuenzel@charite.de). 
