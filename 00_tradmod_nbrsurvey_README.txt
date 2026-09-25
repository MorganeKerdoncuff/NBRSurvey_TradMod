<Help text is included in angle brackets and should be deleted before saving. We recommend you to add "00_" in front of the ReadMe file name (e.g. "00_README.txt"), which will make the file appear on the top of the file overview.>

<DataverseNO README File Template --- General --- Version: 2.4 (2024-03-21)>

This README file was generated on 2026-04-23 by Morgane KERDONCUFF.
Last updated: [YYYY-MM-DD].


-------------------
GENERAL INFORMATION
-------------------
// Title of Dataset: 
// DOI: 
// Contact Information
<The person to be contacted for questions about the dataset>
     // Name: Morgane KERDONCUFF
     // Institution: University of Bergen
     // Email: morgane.kerdoncuff@pm.me
     // ORCID: 0000-0003-2223-1857

<Whenever applicable, the following information should be registered in the metadata schema of DataverseNO. In the text below, remove fields/lines that are not applicable, and leave the rest unchanged. >
// Contributors: See metadata field Contributor.
// Data Type: See metadata field Data Type.
// Date of Collection: See metadata field Date of Collection.
// Geographic location: See metadata section Geospatial Metadata.
// Funding sources: See metadata section Funding Information.

// Description of dataset:
<(Short) description of what the dataset is about, including reference to related project(s) and publication(s), if applicable. Should correspond to the information entered in the metadata fields Description and Related Publication.>
This dataset is associated to the TradMod project (2018-2023) "From traditional resource use to modern industrial production: Holistic management in western Norway", funded by the Research Council of Norway (grant number 280299).
It contains biological, environmental and management data collected in 2019 and 2020 in the Nordhordland UNESCO Biosphere Reserve, Western Norway.
This dataset was used for the following publications: 
Demeaux et al. (2024). Just graze it! Biodiversity, nectar and forage resources in cultural landscapes grazed by different livestock species. Ecosystems and People. 20. 10.1080/26395916.2024.2311176;
Kerdoncuf (2026). Marking the landscape : Ecological assessment of small-scale grazing systems in the fjord region of the Nordhordland UNESCO Biosphere Reserve (doctoral thesis). https://hdl.handle.net/11250/5332751
One manuscript using this dataset is currently under review at Basic and Applied Ecology, and a data paper describing the method of data collection for this dataset is currently in writing.

--------------------------
METHODOLOGICAL INFORMATION
--------------------------
<Note! If the documentation referred to is not openly available through a persistent URL, it must be added here or uploaded as file(s) to the dataset.>

<Note! It may generally be considered appropriate to have overlap in the methods section of a research data README file with citation of the original article. See Committee on Publication Ethics (COPE) guidance on text recycling: https://publicationethics.org/resources/guidelines-new/text-recycling-guidelines-editors-0.>

// Description of sources and methods used for collection/generation of data:
<Include links or references to sources, publications, reports or other documentation (e.g. survey questionnaires, interview protocols, Preregistrations or Registered Reports) containing (experimental) study design or protocols, or other collection techniques used, as well as personnel involved in data collection/generation.>
Ecological data was collected in 2019 and 2020, from June to August, in the Nordhordland UNESCO Biosphere Reserve (60°47'N, 5°15'E), Western Norway.
We visited a total of 45 sites of semi-natural grasslands and heathlands distributed along a west-east gradient, from the coastline to subalpine areas.
The sampling design is described in details in Kerdoncuf (2026). Marking the landscape : Ecological assessment of small-scale grazing systems in the fjord region of the Nordhordland UNESCO Biosphere Reserve (doctoral thesis). https://hdl.handle.net/11250/5332751.

In each site, we laid out a nested sampling design within a representative area of 400 m^2^, preferentially 20 m x 20 m but when necessary 40 m x 10 m, avoiding rocks, waterlogged ground, trees and site edges.
According to the same exclusion rules, we demarcated three plots of 3 m x 3 m in the sampling area, at least 6 m apart from each other. 
Each plot was divided into nine subplots of 1 m^2^, five dedicated for non-destructive data recording (e.g. plant assemblage) and four for destructive sampling (e.g. soil, beetle assemblage).

For each site, we collected GPS coordinates, using the coordinate reference system ETRS89/UTM zone 32N, and elevation.
We determined slope angle and aspect between the lowest and highest ends of the sampling area using an inclinometer.
In the sampling areas, we determined slope angle and aspect between the lowest and highest points using an inclinometer.
We estimated the average percent cover of rock, mud, trees/tall shrubs, low shrubs, forbs, monocotyledons, bryophytes, ferns, and lichens in the 400 m^2^, and counted and measured the total length of livestock paths.

In non-destructive subplots, we estimated the average percent cover of bare ground, rock, litter, dead wood, dung, vascular plants, bryophytes, lichens, and flower bossom.
We assessed average vegetation height by taking three random measurements within the subplot, and we also measured the tallest plant to record maximum vegeration height.
We identified all plant species, specifically recording blossoming ones, and estimated their percent covers within the subplot.
When necessary, specimens were collected, coded and brought for microscopic or expert identification.
All field identifications and percent cover assessments were made by one or the other of the same two observers across the two field seasons.

In destructive subplots, we collected soil, aboveground biomass and dung-associated beetles.
We performed several measurements of soil physicochemical properties.
On site, we tested soil resistance to penetration twice in each subplot by dropping a sharpened metal rod of 43.4 cm length (diam. 2.2 cm; weight 1.3 kg) from 1 m above ground, through a PVC tube to ensure vertical fall.
The visible standing part of the stick was measured and subtracted from the total stick length to obtain the length of penetration into the soil.
We also recorded whether the stick hit the bedrock.
Before dropping the stick, the surrounding vegetation around a 5 cm radius was cleared to prevent resistance from the aboveground vegetation. 
The tip was regularly sharpened and remeasured to account for the decrease in length, which was of 1.6 % (0.7 cm) by the end of the data collection.

In each subplot, we also sampled two series of soil cores.
We collected three cores for the determination of bulk density and water gravimetric content using PVC tubes of 57.23 cm^3^ (diam. 3.7 cm, height 10 cm).
To limit soil compression during collection, we sharpened the edges of the PVC tubes and cut the surrounding root mat with a knife while pushing down the cores. 
Samples were covered with plastic foil and stored at 6°C before processing in the laboratory. 
If necessary, the soil volume of each core was corrected for holes or unfilled space.
Cheesecloth was attached to the bottom of soil cores with an elastic band to prevent soil from falling during processing.

We took the following measurements: 
initial weight of fresh soil;
weight after water saturation for 24h;
weight after air drying for 24h;
weight after air drying for 48h;
weight after oven drying at 105°C for 48h.
During the processing of 2020 samples, we encountered an issue with the weighing scale which affected the accuracy of some measurements, especially when soil core volume was low.
We therefore recommend to filter the data and only keep weight measurements of cores containing at least 50 cm^3^ of soil.

// Methods for processing the data:
<If data other than raw data are provided, describe how the submitted data were processed from the raw or collected/generated data. The documentation of methods used for data processing should include (if applicable): details that may influence reuse or replication efforts; data cleaning and analysis syntax; code or algorithms, with commenting to explain steps taken, to reproduce all reported findings; de-identification procedures for sensitive human subjects or endangered species data. If applicable, code, algorithm or command files used to create derived data should be included in the dataset and referred to in this section.>
Raw data was anonymised and cleaned using the *janitor* package to ensure that variables are consistent and follow naming convention standards.
One site (bog), which did not belong to the same habitat type than other sites (grasslands, heathlands), was entirely removed from the data.
Latin names of plant species and beetle families were corrected and recoded when mispelled.
When applicable, records containing wrong data entry (e.g. for soil bulk density) were removed from the dataset.

<Remove sections below that are not applicable.>

// Facility-, instrument- or software-specific information needed to interpret the data: 
<If not covered above, include full name and version of software, and any necessary packages or libraries needed to read and interpret the data, e.g. to run scripts. For experimental data, specify and describe the facilities and instruments used in the experiment(s).>

// Standards and calibration information: 

// Environmental/experimental conditions: 

// Describe any quality-assurance procedures performed on the data: 


--------------------
DATA & FILE OVERVIEW
--------------------
// File List: 
<List all files (or folders, as appropriate for dataset organization) contained in the dataset. Where appropriate, a file overview may be provided by explaining file naming conventions, instead of listing individual files. For each file (or folder), provide a brief description of what data it contains, and of the file format (e.g. plain text) if not obvious from file extension (e.g. .txt). If necessary, also include system and hardware requirements needed to open and read the file.>

// Relationship between files, if important: 

// Is this an updated version of a dataset published on DataverseNO? yes/no
<This question is relevant if you make changes to the original dataset on file level. It does not relate to (minor) changes in the metadata. If the answer is yes, fill out the information below and repeat for each file that was updated. If the answer is no, leave the lines below unchanged. These might be used for updates to the dataset in the future> 
     Version number of dataset:     
     File name: 
     Why was the file updated? 
     When was the file updated (YYYY-MM-DD)?: 
     What was changed? 


-----------------------------------------
DATA-SPECIFIC INFORMATION FOR: [FILENAME]
-----------------------------------------
<If file-level documentation is provided elsewhere (e.g., as integrated file metadata e.g. ), it is sufficient to refer to this documentation below.>

<Repeat this section for each dataset, folder or file, as appropriate. Recurring items may also be explained in a common initial section.>

<For TABULAR data, provide a data dictionary/code book containing the following information:>
// Variable/Column List: 
<List variable/column name(s), description(s), unit(s) of measurement, decimal separator (comma or point), value labels, and source(s) as appropriate for each.>
// Missing data codes: 
<Define codes or symbols used to indicate missing data.>
// Specialized formats or other abbreviations used: 

<For QUALITATIVE data (e.g. interviews, images), provide a data list:>
// Data List: 
<The data list should contain the following information, as appropriate: interview identifier (ID); age; gender; occupation, organisation; location; place of interview; date of interview; transcript file name. Information included should be detailed enough to enable sub-setting and filtering, but not detailed enough to enable identification of participants where confidentiality has been promised.>

// Contextual Information:
<Where applicable, provide biographical or other contextual information relating to the interview or interviewee, either embedded at the beginning of the transcript files, as a summary page, or added below, as appropriate.>


--------------------------
SHARING/ACCESS INFORMATION
--------------------------
<Whenever applicable, the following information should be registered in the metadata schema of DataverseNO. In the text below, remove fields that are not applicable, and leave the rest unchanged. >
// Licenses/Restrictions: See Terms tab.
// Links to publications that cite or use the data: See metadata field Related Publication.
// Links/relationships to related data sets: See metadata field Related Dataset.
// Data sources: See metadata field Data Source.
// Recommended citation: See citation generated by repository.


<
Acknowledgements. This README file template is adapted from the following documents/resources:

AJPS. ‘American Journal of Political Science Qualitative Data Verification Checklist’. Wiley, 11 March 2016. https://ajps.org/wp-content/uploads/2019/01/ajps-qualdata-checklist-ver-1-0.pdf.

Berkeley Library, University of California. How to Write a Good Documentation. Available at: http://guides.lib.berkeley.edu/how-to-write-good-documentation

Cornell University. Guide to writing "readme" style metadata. Available at: https://data.research.cornell.edu/content/readme#bestpractices

Corti, Louise, Veerle Van den Eynden, Libby Bishop, Matthew Woollard, Maureen Haaker, and Scott Summers. Managing and Sharing Research Data: A Guide to Good Practice. 2nd edition. Los Angeles: SAGE, 2019.

Dryad. Best practices for creating reusable data publications. Available at: https://datadryad.org/stash/best_practices#describe

PurpleBooth. A template to make good README.md. Available at: https://gist.github.com/PurpleBooth/109311bb0361f32d87a2

University of Bath. Working with data: Data Documentation and Metadata. Available at: https://library.bath.ac.uk/research-data/working-with-data/data-documentation-metadata

Zalando. Zalando's README Template. Available at: https://github.com/zalando/zalando-howto-open-source/blob/master/READMEtemplate.md#readme
>
