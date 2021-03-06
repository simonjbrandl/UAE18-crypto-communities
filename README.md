-   [Summary](#summary)
-   [Details & Contents](#details-contents)
-   [Additional information and
    resources](#additional-information-and-resources)

Summary
=======

Extreme environmental conditions reduce coral reef fish biodiversity and productivity
-------------------------------------------------------------------------------------

### Simon J. Brandl, Jacob L. Johansen, Jordan M. Casey, Luke Tornabene, Renato A. Morais & John A. Burt

### Nature Communications, 2020, <a href="https://doi.org/10.1038/s41467-020-17731-2" class="uri">https://doi.org/10.1038/s41467-020-17731-2</a>

Tropical ectotherms are hypothesized to be vulnerable to environmental
changes, but cas- cading effects of organismal tolerances on the
assembly and functioning of reef fishes are largely unknown. Here, we
examine differences in organismal traits, assemblage structure, and
productivity of cryptobenthic reef fishes between the world’s hottest,
most extreme coral reefs in the southern Arabian Gulf and the nearby,
but more environmentally benign, Gulf of Oman. We show that assemblages
in the Arabian Gulf are half as diverse and less than 25% as abundant as
in the Gulf of Oman, despite comparable benthic composition and live
coral cover. This pattern appears to be driven by energetic deficiencies
caused by environmental extremes and distinct prey resource availability
rather than absolute thermal tolerances. As a consequence, production,
transfer, and replenishment of biomass through cryptobenthic fish
assemblages is greatly reduced on Earth’s hottest coral reefs. Extreme
environmental con- ditions, as predicted for the end of the 21st
century, could thus disrupt the community structure and productivity of
a critical functional group, independent of live coral loss.

Details & Contents
==================

The materials in this folder permit the reproduction of the results
presented in Brandl et al. “Extreme environmental conditions reduce
coral reef fish biodiversity and productivity,” published in Nature
Communications:
<a href="https://doi.org/10.1038/s41467-020-17731-2" class="uri">https://doi.org/10.1038/s41467-020-17731-2</a>.
The following elements are provided:

    * UAE18 crypto communities.Rproj: R-project file
    * UAE18_Rscript.Rmd: R-markdown file that contains all code necessary to reproduce the analyses. The script is divided into 11 chunks, that follow the flow of the manuscript:

      1. SETUP
      
      2. LOAD PACKAGES
      
      3. LOAD DATASETS AND SET AESTHETICS
      
      4. TEMPERATURE MAP OF AREA
      
      5. DETAILED TEMPERATURE DATA FROM IN SITU LOGGERS
      
      6. COMMUNITY STRUCTURE: DIVERSITY, ABUNDANCE, BIOMASS
      
      7. COMMUNITY COMPOSITION: FISH AND BENTHOS
      
      8. PHYSIOLOGY: CTmax and CTmin
      
      9. DIET METABARCODING: NETWORK ANALYSES
      
      10. DIET METABARCODING: RAREFACTION
      
      11. SIZE STRUCTURE AND ABUNDANCES
      
      12. GROWTH MODELING
      

Each chunk has subsections pertaining to data wrangling, data analysis,
and figures yielded by each analysis. Chunk 4 requires download of the
MODIS Aqua temperature data
(<a href="https://oceandata.sci.gsfc.nasa.gov/MODIS-Aqua/Mapped/Daily/4km/sst" class="uri">https://oceandata.sci.gsfc.nasa.gov/MODIS-Aqua/Mapped/Daily/4km/sst</a>)
from 2010 to 2018, which can be processed using the auxilliary script
provided and are provided on FigShare
(<a href="https://figshare.com/projects/Cryptobenthic_fish_assemblages_in_the_United_Arab_Emirates/81644" class="uri">https://figshare.com/projects/Cryptobenthic_fish_assemblages_in_the_United_Arab_Emirates/81644</a>).
Chunk 12 requires the functions provided in the rfishprod package, which
is included.

      * Data: Folder that includes all datasets required and produced by the analysis. All data are described in detail in Chunk 3 of the Rmd file.
      * Brandletal_UAE_NatComms_AuxilliaryScript.R: Script to process temperature data from the MODIS Aqua database.
      * rfishprod: folder that contains functions accessed during the modeling of growth and mortality in Chunk 12. Created by RA Morais. 
      

Additional information and resources
====================================

All code written by Simon J. Brandl
(<a href="mailto:simonjbrandl@gmail.com" class="email">simonjbrandl@gmail.com</a>
and
<a href="https://github.com/simonjbrandl" class="uri">https://github.com/simonjbrandl</a>).
All other data necessary to reproduce the paper are available on
FigShare:
<a href="https://figshare.com/projects/Cryptobenthic_fish_assemblages_in_the_United_Arab_Emirates/81644" class="uri">https://figshare.com/projects/Cryptobenthic_fish_assemblages_in_the_United_Arab_Emirates/81644</a>.
