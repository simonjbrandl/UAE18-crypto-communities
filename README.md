The materials in this folder permit the reproduction of the results
presented in Brandl et al. “Extreme environmental conditions reduce
coral reef fish biodiversity and productivity,” published in Nature
Communications:
<a href="https://doi.org/10.1038/s41467-020-17731-2" class="uri">https://doi.org/10.1038/s41467-020-17731-2</a>.
The following elements are provided: • UAE18 crypto communities.Rproj:
R-project file • UAE18\_Rscript.Rmd: R-markdown file that contains all
code necessary to reproduce the analyses. The script is divided into 11
chunks, that follow the flow of the manuscript: 1. SETUP 2. LOAD
PACKAGES 3. LOAD DATASETS AND SET AESTHETICS 4. TEMPERATURE MAP OF AREA
5. DETAILED TEMPERATURE DATA FROM IN SITU LOGGERS 6. COMMUNITY
STRUCTURE: DIVERSITY, ABUNDANCE, BIOMASS 7. COMMUNITY COMPOSITION: FISH
AND BENTHOS 8. PHYSIOLOGY: CTmax and CTmin 9. DIET METABARCODING:
NETWORK ANALYSES 10. DIET METABARCODING: RAREFACTION 11. SIZE STRUCTURE
AND ABUNDANCES 12. GROWTH MODELING

      Each chunk has subsections pertaining to data wrangling, data analysis, and figures yielded by each analysis. Chunk 4 requires download of the MODIS Aqua temperature data (https://oceandata.sci.gsfc.nasa.gov/MODIS-Aqua/Mapped/Daily/4km/sst) from 2010 to 2018, which can be processed using the auxilliary script provided and are provided on FigShare (https://figshare.com/projects/Cryptobenthic_fish_assemblages_in_the_United_Arab_Emirates/81644). Chunk 12 requires the functions provided in the rfishprod package, which is included.
      • Data: Folder that includes all datasets required and produced by the analysis. All data are described in detail in Chunk 3 of the Rmd file.
      • Brandletal_UAE_NatComms_AuxilliaryScript.R: Script to process temperature data from the MODIS Aqua database.
      • rfishprod: folder that contains functions accessed during the modeling of growth and mortality in Chunk 12. Created by RA Morais. 
      

All code written by Simon J. Brandl
(<a href="mailto:simonjbrandl@gmail.com" class="email">simonjbrandl@gmail.com</a>
and
<a href="https://github.com/simonjbrandl" class="uri">https://github.com/simonjbrandl</a>).
All other data necessary to reproduce the paper are available on
FigShare:
<a href="https://figshare.com/projects/Cryptobenthic_fish_assemblages_in_the_United_Arab_Emirates/81644" class="uri">https://figshare.com/projects/Cryptobenthic_fish_assemblages_in_the_United_Arab_Emirates/81644</a>.
