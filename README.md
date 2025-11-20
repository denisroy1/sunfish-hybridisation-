# Sunfish-hybridisation-
Hybridising Sisters: Hybridisation between two sunfish species in urbanised areas.

Scripts repository used to analyse bluegill and pumpkinseed sunfish and their possible hybridisation in more modified habitats near Canada's second largest urban area. This repository contains the relevant scripts needed to analyse the data files for the manuscript titled:

### _Hybridizing sisters: Evidence of bluegill-pumpkinseed hybridization in modified habitats at the confluence of two major Canadian waterways_ ##
*Canadian Journal of Fisheries and Aquatic Sciences (submitted)*

Authored by: Anthony Gagliano, Éléonore Villain, and Denis Roy

### Scripts:

To see the main scripts for this research, go to the 'Scripts' folder. Script description:

Included herein, are the scripts used in the above study that link the hybridisation between bluegill (*Lepomis macrochirus*) and pumpkinseed (*Lepomis gibbosus*) sunfish around the Island of Montreal to habitat modification. The study compares the hybridisation reported here to that documented in Lake Opinicon Ontario. The data set includes raw genotypes that were analysed using software available online (cited in the paper). It also includes data on geometric morphometrics (TPS-formatted), diet analyses, stable isotope analyses, and modeling linking hybrid presence to environmental parameters. The datasets and files can be found [here](https://borealisdata.ca/dataset.xhtml?persistentId=doi%3A10.5683%2FSP3%2FYYANBO&version=DRAFT) on the Borealis McGill dataverse repository. Datasets 2-5 were analysed using custom scripts available in the scripts folder above.


* **Genetics**: 
Called Genotypes from the sequencing machine. Analysed using various online software, including ARLEQUIN version 3.5.1.2 (Excoffier and Lischer 2010), STRUCTURE (Pritchard et al. 2000), CLUMPAK (Kopelman et al. 2015), Hybrid Lab version 1.1 (Nielsen et al. 2006), and GeneClass2.0 (Rannala and Mountain 1997; Piry et al. 2004). Once the data were analysed by the softwares above - the *INDIVQ* file was used in the `STRUCTURERES.R` to create plots recovering the most likely number of genetic clusters. Other results from the Geneclass/Hybrilab assignments were noted and used directly in the `hybplotsMTL-OPI.R`, which plotted the results of parental and hybrid assignments. 

1a-sunfish-genotypes2025.csv

1b-Evannoplotdata.csv

1c-sunfish_indivq.csv


* **Geometric Morphometrics**: 
Fish images were digitised using the tpsUtil64 (version 1.81) and tpsDIG2w64 (version 2.32; Rohlf 2015) software available at https://www.sbmorphometrics.org/. Once digitised, the TPS produced files (2a and 2c) and their metadata (2b and 2d) were analysed using the script `sunfishGM.R`.

2a-sunfishmtl.TPS

2b-sunfishmtl_meta.csv

2c-sunfishopi.TPS

2d-sunfishopi_meta.csv


* **Diet - stomach content**:
Diet data for each site (Montreal/Opnincon) were accumulated into a spreadsheet workbook file that contained separate worksheets for each location. Each was then converted to a .csv file (files 3a and 3b) and analysed using the `sunfishDiet.R script`. 

3a-dietsdr_mtl.csv

3b-dietsdr_opi.csv


* **Stable Isotope data**:
Stable Isotope data from each individual was run through the Ecological Change and Environmental Stressors Laboratory at McGill University, which then provided the corrected and calibrated data. These are the data provided herein (files 4a and 4b). The Montreal (mtl) and the Opinicon (opi) data were then analysed using the `sunfishSI_mtl.R` and the `sunfishSI_opi.R` scripts, respectively.   

4a-sunfishSI_mtl.csv

4b-sunfishSI_opi.csv


* **Environmental Modeling**:
The Environmental modeling data was collected as outlined in the paper and is provided here in .csv format (file 5). It was analysed using the custom script `hyb-env.R`.

5-envcomp.csv

## Note: 
All data files are available without cost at the Borealis (McGill Dataverse) data repository:  
