# grasslandmicrobialdroughtlegacies
Data from the paper: Legacy effects of intensified drought on the soil microbiome in a mesic grassland

Authors: Leena Vilonen, Shabana Hoosein, Melinda Smith, Pankaj Trivedi

Corresponding author:leena.vilonen@colostate.edu

Journal: pending

Methods: Study site/climate conditions
	This study was conducted in 2018 and 2019 at the Konza Prairie Biological Station. Konza is a native, tallgrass prairie in the Flint Hills of northeastern Kansas (39.09° N, 96.48° W) with warm, wet summers and dry, cold winters. The climate is mesic with ~835 mm of mean annual precipitation (MAP) and ~75% of that rainfall occurring in the growing season (April – September). The climate was differed between 2018 and 2019: MAP was 811 mm in 2018 and 971 mm in 2019, with ~64% and 75% of rainfall occurring during the growing season in each year, respectively. The site was annually burned, ungrazed, and located on a relatively flat, upland with ~1 m or more of well-drained clay loam soils, characterized as silty clay, also referred to as Mollisols.
Experimental Design
	This study utilized a large-scale, well-replicated drought experiment (EDGE) established in 2013. EDGE imposed drought from 2014-2017 using large rainfall exclusion shelters (n=20 total), each 6 x 6 m in size and hydrologically isolated to a depth of ~ 1 m. For the chronic treatment, 10 shelters were covered with strips of clear corrugated polycarbonate spaced as to reduce each growing season rainfall event by 66% (April – September). For the intense treatment, 10 shelters completely excluded rainfall with clear corrugated polycarbonate panels until a similar amount of rainfall would be reduced as the chronic treatment reduces throughout the entire growing season (May – July). This resulted in an annual, average ~45% reduction in rainfall. Shelters were erected in May and typically taken down in July for the intense treatment and September for the chronic treatment. No shelters were erected for the control treatment (n = 10), but the plots were hydrologically isolated from one another and received ambient rainfall throughout the growing season. At the end of the 2017 growing season, the shelters were removed, allowing ambient precipitation to fall onto all of the plots in subsequent years. The plots each had four subplots with one subplot designated as a “destructive plot” as detailed in Griffin-Nolan et al. (2019), in which all soil samples were collected.
Soil sampling
	In 2018 and 2019, we collected bulk soil mid growing season (July), since we expected soil microbial activity to be the most active at the middle of the growing season (Carson and Zeglin, 2018). We homogenized four random soil core samples (15 cm depth x 5.7 cm diameter) collected from each “destructive plot”. Soils were immediately placed in a cooler, were sieved to 2mm within 24 hours, and frozen at -20° C. Ethanol was used to clean the soil-corer and sieve in between sample collection and sieving. DNA was extracted within a month of collection.
Soil moisture
We measured soil moisture in both the field and the lab. We used a hand-held TDR to measure in-situ soil moisture to a depth of 15 cm at the time of soil sampling. We additionally dried field-collected soil (the same soil used to measure nutrients) for 48 hours at 60°C to calculate moisture and soil wet soil/dry conversion factors
DNA extraction
We extracted DNA using the DNeasy PowerSoil Kit (Qiagen, Hilden, Germany) following the manufacturer’s instructions, with the exception of increasing soil weight from 0.25g to 0.5 g and using nuclease free water. The DNA sample was tested for quality and quantity using a NanoDrop Lite Spechtrophotometer (Fischer Scientific, MA, USA). Extractions were repeated if 260/280 absorbance ratio was below 1.6. The DNA ranged from 40-85 ng/L and all samples were diluted to 10 ng/L and stored at -80 C.
qPCR
We quantified the numbers of copies of genes present of the extracted, diluted DNA using quantitative polymerase chain reaction (qPCR). We used 96-well plates using 2 L of diluted DNA with forward and reverse primers and Bioline 2x SensiFAST SYBR No-ROX Mix. The forward and reverse primes we used were Eub338- Eub 518 and ITS 1-5.8S, respectively, to amplify bacterial (16S) and fungal (ITS) genes. The plates with DNA, forward and reverse primers, and Bioline mix were then run on a CFX96 Touch Deep Well Real-Time PCR System (Bio-Rad, CA, USA). Resulting melt curves were visualized and any un-amplified samples were re-run. Number of copies of genes were calculated by dividing by initial DNA concentration and weight of soil added for the initial DNA extraction. The data was then ln-transformed for statistical analysis. We used R statistical program to run one-way ANOVAs in the car package and then ran post-hoc Tukey adjusted tests using the emmeans package.
Amplicon sequencing 
	Part of the DNA extracted from above was set aside for DNA sequencing to measure the diversity and community structure of bacteria and fungi. This DNA was amplified using a T100 PCR Thermal Cycler (Bio-Rad, CA, USA) using platinum and DNA marker bar codes specific to each sample to tag each sample to isolate samples after pooling and primer sets 515F/806R (Caporaso et al., 2012) and ITS1/ITS2R (Caporaso et al. 2012) to amplify a portion of the bacterial 16S rRNA gene and fungal ITS1 region. Bacteria (16S) and fungi (ITS) were then sequenced at the University of Colorado Anschutz Medical Campus Genomics Shared Resource. 16S (region v4; 515f-806r) and ITS (ITS1f-ITS2) paired-end 250-read sequencing was performed through amplicon sequencing using Illumina MiSeq Sequencing (Illumina Inc., CA, USA) to measure bacterial and fungal-associated sequences. Resulting data was returned as multiplexed FASTQ files for downstream analysis. 
Bioinformatics 
	The sequencing results were received in the FASTQ format and processed using CUTADAPT to remove adapters from sequences. The USEARCH v.11 pipeline was used for demultiplexing, denoising (UNOISE; Edgar 2016), quality filtering (UCHIME; Edgar 2011) and 97% operational taxonomic unit (OUT) generation (UPARSE; Edgar 2013). Specifically, we used cutadapt to remove adapters and primers (Martin 2011). Quality filtering was assessed using fastQC (Andrews 2010) and sequences were discarded if they had a low quality score (Q < 20), were short in length (<100 bp), or if they contained ambiguous nucleotides. We then merged paired-end reads. OTUs were clustered using DADA2 and denoised using uNoise 3 (Kang et al. 2021) and counted at the sample level. We assigned taxonomy to the OTUs using USEARCH and UCLUST against the SILVA (Quast et al. 2013) database for 16S sequences and UNITE (Nilsson et al. 2018) for ITS. We also removed sequences matching mitochondrial or chloroplast samples, using protocols per Edgar 2016. The OTU table was then exported as a text file and used for further downstream analyses in R statistical software.
Network analysis 
	To visualize combined network analyses with both bacteria and fungi, we calculated coocurrence network metrics using CoNet (Faust and Raes 2016) within Cytoscape v 3.7.2 software (Shannon et al. 2003). We created the visual representations of these coocurrence networks using the platform Gephi 0.9.2 (Bastian et al. 2009). OTUs with zero abundance and less than 60% occupancy were removed to limit noise in the visualizations. Interactions with a degree (number of linkages per node) less than 10 were also removed (Liu et al. 2011). 
Statistical analyses
	Due to the differences in climatic conditions between 2018 and 2019, we analyzed each year separately. All statistical analyses were performed in R using the OUT table created. We used mctoolsr (https://github.com/leffj/mctoolsr/) to upload our meta data and OUT table. We first measured alpha and beta diversity to test for differences in diversity between the treatments. We rarefied data to 3,000 sequences per sample for both alpha and beta diversity using the single_rarefy function in mctoolrs. We measured alpha diversity by calculating richness (number of OTUs) and Shannon’s Diversity using the diversity function from the vegan package. We then created generalized linear models with drought treatment as a fixed factor to measure significance between treatments for each alpha diversity measurement. We ran one-way ANOVAs and post-hoc Tukey tests using the emmeans package to test for significance between treatments. We then ran Permanovas in the vegan package to test for significant differences between the treatments and used pairwise Permanovas to test for significance between treatments using the package pairwiseadonis. We then measured beta diversity using Canonical Analysis of Principal Coordinates (CAP) using the CAPdiscrim function in the BiodiversityR package using treatment as a fixed factor. 
	To determine if there were differences from the OUT level to the phylum level in bacteria and fungi, we used fixed-effect negative binomial generalized linear models (GLM) from the MASS package in R. First, we transformed the data to relative abundance by dividing count per sample by the total per sample. We then normalized the data using TMM normalization and used the summarize_taxonomy function in mctoolsr to measure differences across the treatments for different taxonomic levels. We used the emmeans package to estimate abundance and standard error of each OTU. Emmeans calculates ln(counts), thus we back-transformed to counts. We were then able to run one-way ANOVAS across treatments depending on the taxonomic level and performed FDR adjustments. We extracted Tukey adjusted post-hoc comparisons across the treatments for phylum, order, class, and family  easur and looked for differences within these groups using the emmeans package. 


Data files: OTU_16S_withtax_2018.xlsx, OTU_16S_withtax_2019.xlsx,OTU_ITS_withtax_2018.xlsx,OTU_ITS_withtax_2019.xlsx, MappingFile_16S_2018.xlsx,MappingFile_16S_2019.xlsx, MappingFile_ITS_2018.xlsx, MappingFile_ITS_2019.xlsx,qpcr.csv, EDGEmicrobiomedateexploration.R

OTU_16S_withtax_2018.xlsx:
This file has the OTU table for bacteria (16S sequences) with phylogeny in 2018.
File aspects:
#OTU ID - The OTU that is associated with the count
Leena1_16S to Leena 30_16S - Sample is the header and numbers in the column are the counts of each OTU in that sample
taxonomy - phlogeny for each OTU

OTU_16S_withtax_2019.xlsx:
This file has the OTU table for bacteria (16S sequences) with phylogeny in 2019.
File aspects:
#OTU ID - The OTU that is associated with the count
Leena1_16S to Leena 30_16S - Sample is the header and numbers in the column are the counts of each OTU in that sample
taxonomy - phlogeny for each OTU

OTU_ITS_withtax_2018.xlsx:
This file has the OTU table for fungi (ITS sequences) with phylogeny in 2018.
File aspects:
#OTU ID - The OTU that is associated with the count
Leena1_ITS to Leena 30_ITS - Sample is the header and numbers in the column are the counts of each OTU in that sample
taxonomy - phlogeny for each OTU

OTU_ITS_withtax_2019.xlsx:
This file has the OTU table for fungi (ITS sequences) with phylogeny in 2019.
#OTU ID - The OTU that is associated with the count
Leena1_ITS to Leena 30_ITS - Sample is the header and numbers in the column are the counts of each OTU in that sample
taxonomy - phlogeny for each OTU

MappingFile_16S_2018.xlsx:
This file has the accompanying meta data for each sample for the bacteria (16S sequences) in 2018.
File aspect:
SampleID - The sample that matches the OTU table
Locus - 16S or ITS
Sample Type - PT1 = envrironmental
Sample Location - Kansas
Treatment - INT = intense, CHR = chronic, or CON = control. See methods for descriptions of the treatments

MappingFile_16S_2019.xlsx:
This file has the accompanying meta data for each sample for the bacteria (16S sequences) in 2019.
File aspect:
SampleID - The sample that matches the OTU table
Locus - 16S or ITS
Sample Type - PT1 = envrironmental
Sample Location - Kansas
Treatment - INT = intense, CHR = chronic, or CON = control. See methods for descriptions of the treatments

MappingFile_ITS_2018.xlsx:

This file has the accompanying meta data for each sample for the fungi (ITS sequences) in 2018.

File aspect:
SampleID - The sample that matches the OTU table
Locus - 16S or ITS
Sample Type - PT1 = envrironmental
Sample Location - Kansas
Treatment - INT = intense, CHR = chronic, or CON = control. See methods for descriptions of the treatments

MappingFile_ITS_2019.xlsx:

This file has the accompanying meta data for each sample for the fungi (ITS sequences) in 2019.

File aspect:
SampleID - The sample that matches the OTU table
Locus - 16S or ITS
Sample Type - PT1 = envrironmental
Sample Location - Kansas
Treatment - INT = intense, CHR = chronic, or CON = control. See methods for descriptions of the treatments

qpcr.csv:

This file contains the qPCR data for bacteria and fungi. 

File columns:
Sample - Which sample does this belong to. Samples are from 1-30.
Block - Which block the sample belongs to. There were 10 blocks and three samples per block.
Treatment - Which treatment the sample belongs to. INT is intense, CHR is chronic, and CON is control. See methods for what these means.
Bacteria - Copy number of bacteria in the sample
Fungi - Copy number of Fungi in the sample
Total - Copy number of bacteria and fungi.

EDGEmicrobiomedateexploration.R:

This file has all the data analyses run. Email for any questions.
