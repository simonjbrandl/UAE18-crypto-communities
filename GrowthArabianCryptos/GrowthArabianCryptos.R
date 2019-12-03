##############################################################################
####### Predicting growth and mortality and calculating productivity #########
############################################################## Dec/2019 ######

library (devtools)
library (tidyverse)

options (stringsAsFactors = F, strip.white = T, header = T)

#------------------------------------------#
## Setting file paths and loading package ##
#------------------------------------------#

pkgpath <- gsub ('Fish Growth UAE', 
			 	 'Documents/SJB1/Projects/UAE18 crypto communities/rfishprod', getwd ())
		## Package location in my computer

devtools::load_all(pkgpath)


###########################################
########### Processing Dataset ############
###########################################

#----------------#
## Loading data ##
#----------------#

## Initial traits ##
predata <- read.csv ('traitlist_allfish.csv')
	## Individual level dataset sent by Simon
	
lenwei <- read.csv ('famlevparamsSJB.csv')
	## Merging with length-weight regression parameters at the family level, needs refining 


#-------------------------#
## Merging length-weight ##
#-------------------------#

data <- left_join (predata, lenwei, by = 'Family') %>%
			mutate (Method = 'Otolth',
					Size = TL/10) %>%
				tidytrait (dataset = db)


#-----------------------------------------------------------------------------#
## Applying rfishprod functions to predict growth and mortality coefficients ##
#-----------------------------------------------------------------------------#

## Defining model formula ##
fmod <- formula (~ sstmean + MaxSizeTL + Diet + Movimentation + Method)
set.seed (31)

## Predicting VBGM parameters ##
gr <- predKmax (data, dataset = db, fmod = fmod, params = xgboostparams, niter = 1000, return = 'pred')

## Predicting individual size-specific mortality rates #
gr$pred$Md <- with (gr$pred, 
					predM (Lmeas = TL, Lmax = MaxSizeTL, Kmax = Kmax, 
						   Lr = 1, temp = sstmean, method = 'Function'))

### Final object ###
datafin <- gr$pred


#------------------------------------------------------#
## Checking if sizes and maximum sizes are compatible ##
#------------------------------------------------------#
			
if (any (datafin$Size > datafin$MaxSizeTL)) { "OhOhOh SOMETHING IS WRONG"} else { "ALL GOOD MATE!!"}
## They are not, so... ##

## Maximum sizes observed for these species ##
sppMaxSize <- datafin [which (datafin$Size > datafin$MaxSizeTL),] %>% 
					group_by (Tax) %>%
						summarise (MaxSize = max (Size) + 0.05) %>%
							select (Tax, MaxSize) %>%
								as.data.frame


## Now including in the final dataset ##
datafin [datafin$Tax %in% sppMaxSize$Tax, 'MaxSizeTL'] <- sppMaxSize[match (datafin [datafin$Tax %in% sppMaxSize$Tax, 'Tax'], sppMaxSize$Tax), 'MaxSize']


## Double-checking ##
if (any (datafin$Size > datafin$MaxSizeTL)) { "OhOhOh SOMETHING IS WRONG"} else { "ALL GOOD MATE!!"}


#-------------------------------------------------#
## Calculating growth and loss through mortality ##
#-------------------------------------------------#

datafin <- datafin %>% mutate (EstWeight = a * (Size ^ b),
							   EstGrowth = somaGain (a = a, b = b, Lmeas = Size, t = 1, 
							   						 Lmax = MaxSizeTL, Kmax = Kmax, silent = FALSE),
							   RelGr = EstGrowth / EstWeight,
							   Growth = (W * RelGr),
							   MortLoss = somaLoss (M = Md, Wei = W, t = 1)) %>%
							rename (Weight = 'W') %>%
								select (Location:TL, Full_ID:Genspe, Kmax, Md, Weight,Growth, MortLoss)
							   
#write.csv (datafin, 'dataproc.csv', row.names = F)

### STARTING FROM HERE IN CASE YOU DIDN'T WANT TO BOOTSTRAP GROWTH PREDICTIONS ###

datafin <- read.csv ('dataproc.csv')



#--------------------------------------------------------------------------#
## Tidying up and calculating Productivity, Consumed Biomass and Turnover ##
#--------------------------------------------------------------------------#

aggreg <- datafin %>% group_by (Location, Sitename, Sitenumber, Habitat) %>%
				summarise (Biom = sum (Weight),
						   Prod = sum (Growth),
						   Cons = sum (MortLoss),
						   Turn = (Prod / Biom) + (Cons / Biom)) %>%
					as.data.frame

## Standing biomass is in g/area
## Productivity is in g/area/day
## Consumed biomass is in g/area/day
## Turnover is in 1/day						   

## Function for plotting ##

myplot <- function (df, x, y, leg, coltab = c('red', 'blue'), jit = 0.5, ...) {

xv <- as.numeric (as.factor (df[,x])) 
at <- range (xv)
xlim <- at + c(-0.5, 0.5)

plot (jitter (xv, jit), df [, y], pch = 16, xaxt = 'n', xlab = '',
		col = coltab [xv], xlim = xlim, ...)
	axis (1, at = 1:2, labels = levels (as.factor (df[, x])))
	legend ('topleft', leg, bty = 'n')

}


## Labels ##
vars <- c ('Biom', 'Prod', 'Cons', 'Turn')
lab <- setNames (c('Standing Biomass', 'Productivity', 'Consumed biomass', 'Turnover'), vars)

## Plotting ##
pdf ('ArabianCryptos.pdf', width = 8, height = 2.5, useDingbats = F)
par (mfrow = c(1,4), mar = c(2,3.2,0,1), mgp = c(1.9, 0.4,0), oma = c (0,0,1,0))

for (i in vars) {

myplot (aggreg, x = 'Location', y = i, leg = lab [i], ylab = i)

}

dev.off ()

