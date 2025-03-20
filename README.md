# Bitter crab disease monitoring in Bering Sea snow crab
This dataset combines the former AFSC Pathobiology Program bitter crab disease detection efforts (2014-2017; conventional PCR) and Virginia Institute of Marine Science bitter crab disease detection efforts (2018-2023; conventional and digital PCR) to 1) estimate annual disease prevalence in index monitoring sites, 2) investigate potential drivers of bitter crab disease, and 3) assess the accuracy and sensitivity of visual diagnostic methods. Many analyses using this combined dataset are leveraged from Fedewa et al. 2025, [Bitter crab disease dynamics in eastern Bering Sea Tanner and snow crab: An underestimated and emergent stressor](https://www.sciencedirect.com/science/article/pii/S016578362500044X). 

# Project PIs: 
Erin Fedewa, Hamish Small, Maya Groner, Reyn Yoshioka

# Funding: 
NOAA National Cooperative Research Program

# METADATA:
================================================

Bibliographic
-------------

| Published     | 03/13/2025   |
| ------------- | ------------- |
| Keywords      | snow crab |
|               |   diagnostic methods |
|               |  Bering Sea |
|               | bitter crab disease |
|               | Hematodinium |

 

Coverage
--------

### Temporal

| Begin    | 2014-06-01 |
| ------------- | ------|
| End   | 2023-08-31 |

 
### Spatial

| LME     |                     |
| ------------- | ------------- |
|                | Eastern Bering Sea |
|                | Northern Bering Sea |

![Rplot](https://github.com/Erin-Fedewa-NOAA/BCS_2014-2017/blob/master/figs/n_year.png)


Attributes
----------
One master dataset has been produced for further modeling via the "append_haul" script. "pcr_haul_master.csv" includes all crab biometric, disease diagnosis and haul level data, with data attributes listed below. 

| Name    |    Description   |   Unit    |
| ------- | ---------------- | ---------- |
|  `cruise` | Cruise ID for Bering Sea bottom trawl surveys. YEAR-01 designates EBS surveys, YEAR-02 designates NBS surveys  |  numeric
| `gis_station`   | Alpha-numeric designation for the station established in the design of AFSC standardized surveys | numeric/text
|  `area_swept`   |   Unit of effort for AFSC bottom trawl surveys: computed by distance towed*mean net width   |   numeric, in ha
|  `cpue`   |   Station-level snow crab density   |   numeric, crab/nmi^2
|  `spno` |   Unique specimen ID. First four numbers correspond to year of collection   |  numeric
|  `species_name` | Species name of specimen sampled   |  text
|  `sex` | Sex of specimen sampled. 1=Male, 2-Female   |  numeric
|  `size`  |  Carapace width of specimen sampled | numeric, in mm
|  `chela`   |   Chela height (males only) used to determine maturity status   | numeric, in mm
 | `shell_cond`  |   Shell condition of specimen sampled. 0=premolt/molting, 1=softshell, 2=newshell, 3=oldshell, 4=very oldshell, 5=graveyard, 9=unk | numeric
 |  `egg_color` | Egg color of clutch (females only). 0=No eggs, 2=Purple, 3=Brown, 4=Orange, 5=Purple-Brown, 6=Pink  |  numeric
 |  `egg_cond` | Egg condition of clutch (females only). 0=No Eggs, 1=Uneyed eggs, 2=Eyed eggs, 3=Dead eggs, 4=Empty eggs cases, 5=Hatching eggs  |  numeric
 |  `clutch` | Size of clutch (females only). 0=Immature: no eggs, 1=Mature: no eggs, 2=Trace to 1/8 full, 3=1/4 full, 4=1/2 full, 5=3/4 full, 6=Full   |  numeric
|  `collection_comments`    |    Notes on datasheet from at-sea samplers   |   text
|  `pcr_result` | Bitter crab disease diagnosis via conventional PCR assay for detection of Hematodinium spp. DNA. 0=uninfected, 1=infected, 3=undetermined   |  numeric
| `year`     |        Year of specimen collection | numeric
|  `pcr_result` | Bitter crab disease diagnosis via PCR assay for detection of Hematodinium spp. DNA. 0=uninfected, 1=infected, 3=undetermined   |  numeric
|  `index_site` | Bering Sea sites established for bitter crab disease monitoring. Sites 4-6 = EBS snow crab index sites, Site 7 = NBS snow crab index site. In 2018, the NBS index site was designated as "NBS" because this was a rapid response survey and sampling was not conducted in standard survey grid   |  numeric
|  `general_location` | Large marine ecosystem. EBS=eastern Bering Sea, NBS=northern Bering Sea   |  text
|  `collected_by` | Agency taking hemolymph samples. SAP = NOAA AFSC Shellfish Assessment Program   |  text
|  `dna_plate_no` | Plate number containing hemolymph sample   |  numeric
|  `dna_well_no` | Individual well plate number containing hemolymph sample   |  alpha-numeric
|  `host_tissue` | Type of biological sample collected for PCR   |  text
|  `preservative` | Method of preservation for tissue/blood sample   |  text
|  `sample_status` | Disposition of tissue/blood sample   |  text
|  `c_v_h` | Combined string with cruise, vessel and haul for collections   |  numeric 
|  `vessel`  |     ID number of the vessel used to collect data for that haul associated with vessel name    |   numeric
|  `haul`      |  Uniquely identifies a sampling event (haul) within an AFSC cruise. It is a sequential number, in chronological order of occurrence |  numeric
|  `visual_positive` | Infection status via visual diagnosis. NA=not recorded, 0=not infected, 1=infected   |  numeric
|  `prioritization` | Dataset has been filtered for prioritization=1, as these were samples prioritized for VIMS PCR analysis   |  numeric
|  `dna_quant`  |  Concentration of DNA in hemolymph sample extracted (2018+ data only)  | numeric, in ng/microL  
|  `ratio_260_280`  |  Measurement of the absorbance of light at 260 and 280 nm of extracted DNA. Used to determine the quality/purity, with ~1.8-2 ideal (2018+ data only)  | numeric  
|  `nssu_pcr`  |  General metazoan PCR primer set (nSSU A and nSSU B) used in conventional PCR as a nominal control for the presence of amplifiable DNA in the samples, where 0=failure to produce an amplification product/no host DNA cells present.  (2018+ data only)  | numeric  
|  `faint_pos`  |  Subjective assessment of the intensity of the ethidium bromide stained PCR amplification product (2018+ data only). F=faint positive  | numeric  
|  `no_valid_partitions`  |  Number of partitions where the reference fluorescence signal is detected and therefore this partition can be used in subsequent digital PCR calculations. (2018+ data only)  | numeric  
|  `no_positive_partitions`  |  Number of partitions where fluorescent signal from amplification of the target DNA (Hematodinium) was detected during digital PCR. (2018+ data only)  | numeric  
|  `copies_ul`  |  Poisson statistics-based calculation of the copy number of the target DNA sequence. (2018+ data only)  | numeric 
|  `maturity`  |  Maturity of specimen sampled, as determined by chela (males) or clutch morphology (females). 0=Immature, 1=Mature  | numeric  
| `start_date`     |        Date of sampling | date, month/day/year
| `mid_latitude`       |   Latitude of specimen collection. Designates latitude at start of haul for AFSC standardized surveys    | numeric
|  `mid_longitude`    | Longitude of specimen collection. Designates longitude at start of haul for AFSC standardized surveys | decimal degree
 | `bottom_depth`    |    Bottom depth at station for AFSC standardized surveys  | numeric, in m
|  `gear_temperature`   |    Bottom temperature at sampling station | degree C


:::

Distribution
------------

  File                                              Format    
  ------------------------------------------------- -------- -------------------------------------------------
  `pcr_haul_master`   `csv`    [Download](https://github.com/Erin-Fedewa-NOAA/BCS_2014-2017/tree/master/data)
 ------------------------------------------------- -------- -------------------------------------------------

 # NOAA License
This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

Software code created by U.S. Government employees is not subject to copyright in the United States (17 U.S.C. §105). The United States/Department of Commerce reserve all rights to seek and obtain copyright protection in countries other than the United States for Software authored in its entirety by the Department of Commerce. To this end, the Department of Commerce hereby grants to Recipient a royalty-free, nonexclusive license to use, copy, and create derivative works of the Software outside of the United States.
