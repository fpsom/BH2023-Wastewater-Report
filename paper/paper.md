---
title: 'Standards for a multimodal data wastewater surveillance process'
title_short: 'Standards for a multimodal data wastewater surveillance process'
tags:
  - Wastewater surveillance
authors:
  - name: Fotis Psomopoulos
    orcid: 0000-0002-8021-9162
    affiliation: 1
  - name: Konstantinos Kyritsis
    orcid: 0000-0001-8035-341X
    affiliation: 1
  - name: Ivan Topolsky
    orcid: 0000-0002-7561-0810
    affiliation: 2
  - name: Amy Heather Fitzpatrick
    orcid: 0000-0002-1883-0489
    affiliation: 3
  - name: Gabriele Leoni
    orcid: 0000-0002-4899-5284
    affiliation: 4
affiliations:
  - name: Institute of Applied Biosciences, Centre for Research and Technology Hellas, Thessaloniki, Greece
    index: 1
  - name: Computational Biology Group, SIB Swiss Institute of Bioinformatics, Basel,Switzerland
    index: 2
  - name: University College Dublin
    index: 3
  - name: European Commission, Joint Research Centre (JRC), Ispra, Italy
    index: 4
date: 3 November 2023
cito-bibliography: paper.bib
event: ELIXIR BioHackathon 2023
biohackathon_name: "ELIXIR BioHackathon 2023"
biohackathon_url:   "https://biohackathon-europe.org/"
biohackathon_location: "Barcelona, Spain, 2023"
group: Wastewater Surveillance
# URL to project git repo --- should contain the actual paper.md:
git_url: https://github.com/fpsom/BH2022-Wastewater-Report
# This is the short authors description that is used at the
# bottom of the generated paper (typically the first two authors):
authors_short: Fotis Psomopoulos & Ivan Topolsky \emph{et al.}
---


<!--

The paper.md, bibtex and figure file can be found in this repo:

  https://github.com/journal-of-research-objects/Example-BioHackrXiv-Paper

To modify, please clone the repo. You can generate PDF of the paper by
pasting above link (or yours) in

  http://biohackrxiv.genenetwork.org/

-->

# Introduction

Wastewater-Based Epidemiology (WBE) holds significant promise as an early-warning system for various pathogens, with the potential to contribute to public health ([@. While considerable wet-lab efforts have been dedicated to this field, its limits and actual potential remain undefined when viewed from a dry-lab perspective.

Wastewater influent into a Wastewater Treatment Plant (WWTP) is a complex mixture of sewage from the catchment area, encompassing a diverse array of nucleic acid originating from Eukaryotic, Prokaryotic, and Viral organisms. Although the precise composition of wastewater remains poorly documented, it has been established that these microbial communities are highly fragmented. Additionally, Viral organisms are found in lower relative abundance compared to their Prokaryotic counterparts, primarily due to their smaller genomic size. These characteristics pose unique challenges for dry-lab and bioinformatic analysis.

Furthermore, the variability of wastewater samples across seasons, geographical locations, and days of the week adds another layer of complexity. The choice of wet-lab preparation methods and High Throughput Sequencing (HTS) techniques has been shown to impact the detection of specific microorganisms. Taking all of these components into account, it is clear that there are likley challenges in detecting and characterising the ground-truth in this sample type. 

As a response to these challenges, our biohackathon project is dedicated to exploring the potential of WBE in addressing specific biological questions. Our approach involves reviewing metadata associated with various read archieve wastewater sequencing data uploads and developing a prototype workflow to assess the accuracy of Antimicrobial Resistance (AMR) detection in wastewater. This research aims to shed light on the full capabilities and limitations of Wastewater-Based Epidemiology from a dry-lab perspective. The end-goal is a workflow that permits novice users to establish the probability of detecting their specific ground-truth in a shotugn metagenomics wastewater sequencing dataset.





# Metadata standards for wastewater HTS data upload

To collectively address the challenge of harmonizing metadata for NGS data generated from sequencing wastewater samples across major deposition databases (such as SRA and ENA), we conducted a community survey. The objective of this survey was to identify the essential metadata that should be included during the deposition of NGS data and to establish a common alignment point. The [survey]((https://ec.europa.eu/eusurvey/runner/ELIXIR-BH2023-Project31)) featured questions adapted from the ENA checklist (please note that these questions are a slight modification of the ENA checklist available at [link](https://www.ebi.ac.uk/ena/browser/view/ERC000036)). The questions can be viewed in the table below. 


Number |  Metadata name |  Metadata Description |  Metadata Type 
-----------|------------------------|--------------------------------|---------------------
1 | Name of the sampling site |  Refers to the name of the site/station where data/sample collection is performed | Free text
2 | Nucleic acid extraction |  A link to a literature reference, electronic resource or a standard operating procedure (SOP), that describes the material separation to recover the nucleic acid fraction from a sample. | Free text
3 | Nucleic acid amplification |  A link to a literature reference, electronic resource or a standard operating procedure (SOP), that describes the enzymatic amplification (PCR, TMA, NASBA) of specific nucleic acids |  Free text
4 | Investigation type | Nucleic Acid Sequence Report is the root element of all MIxS compliant reports as standardised by Genomic Standards Consortium | one of the: bacteria-archaea, eukaryote, metagenome, metagenome-assembled genome, metatranscriptome, mimarks-specimen, mimarks-survey, organelle, plasmid, single amplified genome, uncultivated viral genomes, virus
5 | Surveillance target | Valid species level NCBI taxon id of the organism being surveyed for (if any) e.g. for Escherichia coli enter 562. |Ontology choice (e.g. ncbitaxon)
6 | Collection date | The date the sample was collected with the intention of sequencing, either as an instance (single point in time) or interval. In case no exact time is available, the date/time can be right truncated i.e. all of these are valid ISO8601 compliant times: 2008-01-23T19:23:10+00:00; 2008-01-23T19:23:10; 2008-01-23; 2008-01; 2008. | Date
7 | Geographic location (country and/or sea) | The location the sample was collected from with the intention of sequencing, as defined by the country or sea. Country or sea names should be chosen from the INSDC country list (http://insdc.org/country.html). | Country code/list
8 | Geographic location (latitude) |  The geographical origin of the sample as defined by latitude. The values should be reported in decimal degrees and in WGS84 system | Number
9 | Geographic location (longitude) | The geographical origin of the sample as defined by longitude. The values should be reported in decimal degrees and in WGS84 system | Number
10 | Geographic location (region and locality) | The geographical origin of the sample as defined by the specific region name followed by the locality name. | Free text
11 | Amount or size of sample collected | The total amount or size (volume (ml), mass (g) or area (m2) ) of sample collected | Number
12 | Sample storage duration | duration for which sample was stored | Number
13 | Sample storage temperature | temperature at which sample was stored, e.g. -80 | Number
14 | Sample storage location | location at which sample was stored, usually name of a specific freezer/room | Free text
15 | Sampling time point | Sampling time point | Free text
16 | Sample transportation temperature | transportation temperature from sample site to storage | Number
17 | Sample transportation date | transportation/shipping date of the sample. Format:YYYY-MM-DD | Regular expression
18 | Sample transportation time | transportation time from sample site to storage | Number
19 | Receipt date | Date on which the sample was received. Format:YYYY-MM-DD. Please provide the highest precision possible. If the sample was received by the institution and not collected, the 'receipt date' must be provided instead. Either the 'collection date' or 'receipt date' must be provided. If available, provide both dates. | Regular expression
20 | Sewage type | Type of sewage based on origin. | one of the: wastewater treatment plant (municipal or industrial), open sewer line, river, stream, stagnant pool, or other
21 | Temperature | temperature of the sample at time of sampling | Number
22 | Area of sampling site | Please indicate if there are specific facilities in the area covered by the sewage sample. For example: farming, slaughterhouse(s), industry, hospital(s) or any other facility. | Free text
23 | Size of the catchment area | Refers to the size of the area that is drained by the sampled sewage system in square km. | Number
24 | Population size of the catchment area | Refers to the number of people living in the area covered by the sewage system | Number



We have compiled a list of critical information that we consider essential for facilitating a thorough comparison of wastewater analysis results and evaluating the robustness and reproducibility of wastewater workflows. These key metadata points include:

- Amplicon probe map (version) (e.g., inserts or primers bedfile)
- Quality control measures performed (e.g., filters applied to raw reads, inclusion of read deduplication steps, which can impact variant calling)
- Sequencing coverage
- Variant caller used (e.g., GATK4 HaplotypeCaller, which cannot call variants with Allele Frequency < 0.5, versus LowFreq, which can, but has higher false positive rates, and ShoRAH/VILOCA, known for its resilience at low frequencies)
- Version of the lineages database used (relevant if Pangolin, Freyja, or other tools were employed)

# LIMBO - Low Input Microbial Biological Observations workflow
One of the greater challenges working with wastewater sequencing data is determining at what frequency/coverage of detection can the user/analyst have confidence in results. In this case, we propose a workflow LIMBO, which permits the user to generate scenarios of low input microbial presence and incorporate end-points for classification accuracy such as taxonomic classification, detection of Single Nucleotide Variants (SNVs) or identification of novel taxa/variants/lineages. The output is a measure of the ground-truth as a probability. 

For this prototype we have taken shotgun metagenomics sequencing of wastewater for Antimicrobial Resistence (AMR) detection as a case-study. 

Proposed workflow: 
User data input with varying parameters
Sequencing method  -”amplicon”, “shotgun_metagenomics”, “metatranscriptomics”
Read length (short / long)
Sequencing strategy (Paired-End / Single-End )
Distribution of taxa within a multi fasta e.g “normal”, “beta” etc
Taxa relative abundance (Provided in reads percentages, coverage values, relative proportions or read counts) 
Path to input fasta file(s)
Simulation of HTS data using MeSS (able to take into account all user data 
parameter mentioned above)
Possibility to generate multiple replicates from the same input set using a normal distribution and a standard deviation parameter
Post sequencing processing determined by sequencing input 
Shotgun metagenomics - Alignment approach
Shotgun metagenomics - SPades assembly
End-point detection - AMR genes/plasmids with AMRFinderPlus, Abricate, CARD
Output - probability of obtaining ground-truth

## Acknowledgements

Some of the authors were funded by ELIXIR, the research infrastructure for life-science data, to join the BioHackathon Europe.

## References
