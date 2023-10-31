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
git_url: https://github.com/fpsom/BH2023-Wastewater-Report
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


Wastewater-Based Epidemiology (WBE) began to flourish in 2001 and was initially used for the detection of illegal drugs. Historically though, around the 1950's WBE was implemented for the detection of infectious diseases, when the first studies to detect poliovirus and enteroviruses in sewage samples were conducted. These studies were based on the fact that since viruses cannot grow outside of host cells, they lose the ability to evolve in wastewater, and thus the concentration values that were found in sewages are related to those that were initially egested by the population.

During the COVID-19 pandemic WBE became a community-wide monitoring tool which provided real-time data concerning both public health status, environmental monitoring as well as food safety. Moreover, in addition to the traditional wastewater surveillance metrics that rely on biochemical indicators, the pandemic led to the increased use of high-throughput sequencing data on wastewater samples. This led to a number of key challenges and questions in the field, ranging from the data capture efforts through biosensors, until the data management and data analysis efforts

The main goal of this project was to define a framework for addressing specific health policy-related questions based on multi-modal wastewater surveillance, including a critical review on the applicable standards for (meta)data. The effort continued on the outputs of the relevant BH2022 project on wastewater surveillance [@ref? ], as well as on the collective expertise of the [ELIXIR Wastewater Surveillance Working Group](https://www.covid19dataportal.org/partners?activeTab=Working%20groups). Ultimately, the proposed framework will be demonstrated using existing publicly available data.


# Considerations around defining the biological composition of wastewater




# Metadata for wastewater sequence data

We identified a list of relevant information that, to us, are necessary to properly compare the results of wastewater analyses and/or to assess wastewater workflows robustness and reproducibility. These are:

Amplicon probe map (version) (e.g.: inserts or primers bedfile)
Sequencing coverage
QC performed (i.e., filter used on raw reads, whether reads deduplication steps were included or not -> it affects variant calling)
Variant caller used (e.g., GATK4 HaplotypeCaller cannot call variants with Allele Freq < 0.5, LowFreq on the opposite can, but has higher false positive rates, ShoRAH/VILOCA is more resilient at low-frequencies)
Version of the lineages database used (in case Pangolin/Freyja/Others were used)

In order to do a community-driven review of the metadata that are 

In order to identify the key metadata that should be used during the deposition of NGS data generated after sequencing wastewater samples, with a goal of finding an optimal point of alignment across the various major deposition databases (such as SRA and ENA), we designed and run a short community survey (https://ec.europa.eu/eusurvey/runner/ELIXIR-BH2023-Project31) with the following questions (please note that these are a minor adaptation of the ENA checklist https://www.ebi.ac.uk/ena/browser/view/ERC000036):


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


## Acknowledgements

Some of the authors were funded by ELIXIR, the research infrastructure for life-science data, to join the BioHackathon Europe.

## References
