########################################################################
#NLA 2007 analyses - Comparison of water-column and surface sediment   #
#For manuscript, Winegardner et al.                                    #
#R version 3.0.2                                                       #
########################################################################

########################################################################
#Section 1 - Packages                                                  #
########################################################################
#Packages
library(vegan) #ordination, transformations etc. 
library(packfor) #for different type of forward selection
library(ggplot2) #plotting and visuals
library(moments) #skewness
library(PCNM) #PCNM
library(car)  #simple correlations
library(spdep) #spatial
library(ape) #spatial
library(ade4) #spatial
library(AEM) #spatial
library(spacemakeR) #spatial
library(geoR) #BoxCox transformations
library(gridExtra) #Arranging visuals in ggplot2

########################################################################
#Section 2 - Workflow                                                  #
########################################################################

#a) Import site info, environmental data and species data. 

#b) Hellinger transformation of diatom species relative abundance data. 

#c) Test environmental variables for normality, Box-Cox 
#transformations. 

#d) Forward selection of environmental variables and created reduced
# environmental dataset. 

#e) Generate PCNM axes and 3rd order polynomials using coordinates. 
#Forward selection to reduce spatial variable dataset. 

#f) PCAs and OLS regression of PC scores (whole dataset). 

#g) RDAs and OLS regression of RDA scores (whole dataset). 

#h) Procrustes analyses. 

#i) Repeat analyses with planktonic species only. 

#j) Repeat analyses using visit data. 

#k) Bray-Curtis dissimilarities. 

##########################################################################

##########################################################################
#Section 3 - Data, transformations and normality                         #
##########################################################################
#a) Import site info, environmental data and species data. 

#b) Hellinger transformation of diatom species relative abundance data. 

#c) Test environmental variables for normality, Box-Cox 
#transformations. 

#Site and physical data
sites.phys<- read.csv(file.choose()) #nla2007_siteinfo.csv
rownames(site.phys)<- as.character(site.phys[,1])

#Env
env.data<- read.csv(file.choose()) #nla2007_env_468.csv
rownames(env.data)<- as.character(env.data[,1]) 
#Can be standardized using decostand(  , "standardize") 

#Diatom species
wc.diat<- read.csv(file.choose()) #nla2007_diat_468wc.csv 
rownames(wc.diat)<- as.character(wc.diat[,1])

ssed.diat<- read.csv(file.choose()) #nla2007_diat_468sed.csv
rownames(ssed.diat)<- as.character(ssed.diat[,1]) 
#Can be Hellinger transformed using decostand( , "hell") 

#Normality of environmental variables and Box-Cox 
ph<-env.data$pH_field #PH
cond<-envd.ata$Cond #Conductivity
turb<-env.data$Turb #Turbidity
DOC<-env.data$DOC #Dissolved organic carbon
NH4<-env.data$NH4 #Ammonium
NO3_NO2<-env.data$NO3_NO2 #Nitrate + Nitrite by flow injection analysis
NTL<-env.data$NTL #Total nitrogen
PTL<-env.data$PTL #Total phosphorous 
SO4<-env.data$SO4 #Sulphate
Ca<-env.data$CA #Calcium
Mg<-env.data$MG #Magnesium
col<-env.data$Colour #Colour
SIO2<-env.data$SIO2 #Silica
H<-env.data$H #H+
OH<-env.data$OH #Hydroxide 
ion.bal<-env.data$BALANCE2 #Ion balance 
secchi<-env.data$SECMEAN #Mean secchi disk depth
Cl<-env.data$CL #Chlorine
chla<-env.data$CHLA #Chlorophyll a

#Transforming variables

#ph - stays as raw data

#cond
shapiro.test(cond)
hist(cond) #non-normal
bcoxcond = boxcoxfit(cond)
bcoxcond #record parameters
cond.tr = (cond^bcoxcond$lambda - 1)/bcoxcond$lambda #transformed conductivity
shapiro.test(cond.tr)
hist(cond.tr)
#use cond.tr

#turb
shapiro.test(turb)
hist(turb) #non-normal
bcoxturb = boxcoxfit(turb)
bcoxturb
turb.tr = (turb^bcoxturb$lambda - 1)/bcoxturb$lambda #transformed turbidity
shapiro.test(turb.tr)
hist(turb.tr)
#use turb.tr

#DOC
shapiro.test(DOC)
hist(DOC) #non-normal
bcoxDOC = boxcoxfit(DOC)
bcoxDOC
DOC.tr = (cond^bcoxDOC$lambda - 1)/bcoxDOC$lambda #transformed DOC
shapiro.test(DOC.tr)
hist(DOC.tr)
#use DOC.tr

#NH4
shapiro.test(NH4)
hist(NH4) #non-normal
bcoxNH4 = boxcoxfit(NH4)
bcoxNH4
NH4.tr = (NH4^bcoxcond$lambda - 1)/bcoxNH4$lambda #transformed NH4
shapiro.test(NH4.tr)
hist(NH4.tr)
NH4.log<-log(NH4)
hist(NH4.log)
#use NH4.tr

#NO3_NO2
shapiro.test(NO3_NO2)
hist(NO3_NO2) #non-normal
#requires positive data
NO3_NO2.pos<-NO3_NO2-min(NO3_NO2)+0.1
bcoxNO3_NO2.pos = boxcoxfit(NO3_NO2.pos)
bcoxNO3_NO2.pos
NO3_NO2.pos.tr = (NO3_NO2.pos^bcoxNO3_NO2.pos$lambda - 1)/bcoxNO3_NO2.pos$lambda #transformed NO3_NO2 (positive)
hist(NO3_NO2.pos.tr)
#can't normalize --> REMOVE from analysis. 

#NTL
shapiro.test(NTL)
hist(NTL) #non-normal
bcoxNTL = boxcoxfit(NTL)
bcoxNTL
NTL.tr = (NTL^bcoxNTL$lambda - 1)/bcoxNTL$lambda #transformed NTL
shapiro.test(NTL.tr)
hist(NTL.tr)
#use NTL.tr

#PTL
shapiro.test(PTL)
hist(PTL) #non-normal
bcoxPTL = boxcoxfit(PTL)
bcoxPTL
PTL.tr = (PTL^bcoxPTL$lambda - 1)/bcoxPTL$lambda #transformed PTL
shapiro.test(PTL.tr)
hist(PTL.tr)
#use PTL.tr

#SO4
shapiro.test(SO4)
hist(SO4) #non-normal
bcoxSO4 = boxcoxfit(SO4)
bcoxSO4
SO4.tr = (SO4^bcoxSO4$lambda - 1)/bcoxSO4$lambda #transformed SO4
shapiro.test(SO4.tr)
hist(SO4.tr)
#use SO4.tr

#Ca
shapiro.test(Ca)
hist(Ca) #non-normal
bcoxCa = boxcoxfit(Ca)
bcoxCa
Ca.tr = (Ca^bcoxCa$lambda - 1)/bcoxCa$lambda #transformed Ca
shapiro.test(Ca.tr)
hist(Ca.tr)
#use Ca.tr

#Mg
shapiro.test(Mg)
hist(Mg) #non-normal
bcoxMg = boxcoxfit(Mg)
bcoxMg
Mg.tr = (Mg^bcoxMg$lambda - 1)/bcoxMg$lambda #transformed Mg
shapiro.test(Mg.tr)
hist(Mg.tr)
#use Mg.tr

#col
shapiro.test(col)
hist(col) #non-normal
#requires positive data
col.pos<-col-min(col)+0.1
bcoxcol = boxcoxfit(col.pos)
bcoxcol
col.pos.tr = (col.pos^bcoxcol$lambda - 1)/bcoxcol$lambda #transformed colour
shapiro.test(col.pos.tr)
hist(col.pos.tr)
#use col.pos.tr

#SIO2
shapiro.test(SIO2)
hist(SIO2) #non-normal
bcoxSIO2 = boxcoxfit(SIO2)
bcoxSIO2
SIO2.tr = (SIO2^bcoxSIO2$lambda - 1)/bcoxSIO2$lambda #transformed SIO2
shapiro.test(SIO2.tr)
hist(SIO2.tr)
#use SIO2.tr

#H
shapiro.test(H)
hist(H) #non-normal
#requires positive data
H.pos<-H-min(H)+0.1
bcoxH = boxcoxfit(H.pos)
bcoxH
H.pos.tr = (H.pos^bcoxH$lambda - 1)/bcoxH$lambda #transformed positive H
shapiro.test(H.pos.tr)
hist(H.pos.tr)
#Can't normalize, REMOVE from dataset

#OH
shapiro.test(OH)
hist(OH) #non-normal
bcoxOH = boxcoxfit(OH)
bcoxOH
OH.tr = (OH^bcoxOH$lambda - 1)/bcoxOH$lambda #transformed OH
shapiro.test(OH.tr)
hist(OH.tr)
#Use OH.tr

#ion.bal
shapiro.test(ion.bal)
hist(ion.bal) #non-normal
#requires positive data
ion.bal.pos<-ion.bal-min(ion.bal)+0.1
bcoxionbal = boxcoxfit(ion.bal.pos)
bcoxionbal
ion.bal.pos.tr = (ion.bal.pos^bcoxionbal$lambda - 1)/bcoxionbal$lambda #transformed Ion balance
hist(ion.bal.pos.tr)
#Difficult to normalize, REMOVE from analyses

#secchi
shapiro.test(secchi)
hist(secchi) #non-normal
bcoxsecchi = boxcoxfit(secchi)
bcoxsecchi
secchi.tr = (secchi^bcoxsecchi$lambda - 1)/bcoxsecchi$lambda #transformed mean Secchi disk depth
shapiro.test(secchi.tr)
hist(secchi.tr) #and confirmed with Shapiro-Wilkes
#Use secchi.tr

#Cl
shapiro.test(Cl)
hist(Cl) #non-normal
bcoxCl = boxcoxfit(Cl)
bcoxCl
Cl.tr = (Cl^bcoxCl$lambda - 1)/bcoxCl$lambda #transformed Cl
shapiro.test(Cl.tr)
hist(Cl.tr)
#Use Cl.tr

#chla
shapiro.test(chla)
hist(chla) #non-normal
bcoxchla = boxcoxfit(chla)
bcoxchla
chla.tr = (chla^bcoxchla$lambda - 1)/bcoxchla$lambda #transformed Chl a
shapiro.test(chla.tr)
hist(chla.tr)
#use chla.tr

#zmax
shapiro.test(zmax)
hist(zmax) #non-normal
bcoxzmax = boxcoxfit(zmax)
bcoxzmax
zmax.tr = (zmax^bcoxzmax$lambda - 1)/bcoxzmax$lambda #transformed Zmax
shapiro.test(zmax.tr)
hist(zmax.tr)
#use zmax.tr

#Variables to retain and combine into a env/morphometric dataset
#Raw variables to standardize using decostand
#ph
#cond
#turb
#DOC
#NH4
#PTL
#SO4
#Ca
#Mg
#col
#OH
#zmax
#Use with Box-Cox transformation
#SIO2.tr
#NTL.tr
#secchi.tr
#chla.tr
#Cl.tr

##########################################################################

##########################################################################
#Section 4 - Forward selection of environmental and spatial variables    #
##########################################################################
#d) Forward selection of environmental variables and created reduced
# environmental dataset. 

#e) Generate PCNM axes and 3rd order polynomials using coordinates. 
#Forward selection to reduce spatial variable dataset. 

#Code for this section is reproduced in part from 
#Borcard et al. (2011) "Numerical Ecology with R" 

#Use cbind and as.data.frame to create matrix of transformed env
#variables from Section 3. --> env.data2

#Fwd selection of env variables using water-column diatom data
#Packfor
wc.rda.full<-rda(wc.diat[,#:#], env.data2[,#:#]) #Add Hellinger if not done
(wc.full.R2a<- RsquareAdj(wc.rda.full)$adj.r.squared) 
wc.fwd<-forward.sel(wc.diat[,#:#], env.data2[,#:#])
wc.sign.fwd<- sort(wc.fwd$order) #significant variables and order
env.data.red<- env.data2[,c(wc.sign.fwd)] #New reduced env matrix
colnames(env.data.red)

#Ordistep
wc.rda.full <- rda(wc.diat[,#:#] ~ ., env.data2[,#:#]) #Add Hellinger if not done
ordistep(wc.rda.full) 
wc.rda.mod<-ordistep(rda(wc.diat[,#:#] ~ 1, env.data2[,#:#]), scope = formula(wc.rda.full))

#Generate spatial variables and selection
sites.xy<-data.frame(sites.phys[,2:3]) #use Albers as coordinate projections
xy.dl<-dist(sites.xy)

sites.PCNM.auto<-PCNM(xy.dl)
sites.PCNM.auto$expected_Moran 
sites.PCNM.auto$Moran_I
(select<- which(sites.PCNM.auto$Moran_I$Positive == TRUE))
length(select)
sites.PCNM.pos<- as.data.frame(sites.PCNM.auto$vectors)[,select]

#Detrend water-column diatom data for generating spatial variables
wcdiat.det<-resid(lm(as.matrix(wc.diat) ~ ., data = sites.xy))

wcdiat.PCNM.rda<-rda(wcdiat.det, sites.PCNM.pos) ##Now using positive only. 
anova.cca(diat.PCNM.rda) #If analysis is significant, compute adj R2
# and run fwd selection on PCNM variables. 

(wc.space.R2a<- RsquareAdj(wcdiat.PCNM.rda)$adj.r.squared) 
(wcdiat.PCNM.fwd<- forward.sel(wcdiat.det, as.matrix(sites.PCNM.pos), adjR2thresh=wc.space.R2a))
(no.sig.PCNM.wc<- nrow(wcdiat.PCNM.fwd))
#Can then continue on with reduced set of PCNM variables as spatial dataset. 

#Variation partitioning with varpart() as per Borcard et al. (2011)
#Third order polynomials for coordinates can be created using:
variable.poly<- poly(variable, degree=3) 

###########################################################################

###########################################################################
#Section 5 - Whole dataset analyses                                       #
###########################################################################
#f) PCAs and OLS regression of PC scores (whole dataset). 

#g) RDAs and OLS regression of RDA scores (whole dataset). 

#h) Procrustes analyses. 

#PCAs
#PCA 1 - water-column diatoms
wc.pca<- rda(wc.diat[,#:#])
wc.pcscores<- data.frame(scores(wc.pca, choices = 1, display="sites")) 
#can extract PC2 as well, choices = 1:2

#PCA 2 - surface sediment diatoms
ssed.pca<- rda(ssed.diat[,#:#])
ssed.pcscores<- data.frame(scores(ssed.pca, choices = 1, display="sites"))

pcdiat.scores<-data.frame(cbind(wc.pcscores, ssed.pcscores))
colnames(pcdiat.scores)<-c("WC_PC1", "SSed_PC1")
#Correlation water-column PC1 with surface sediment PC1
model<-lm(pcdiat.scores$WC_PC1~pcdiat.scores$SSed_PC1)
summary(model)

#Example plot#
plotpca<- qplot(y=pcdiat.scores$SSed_PC1, x=pcdiat.scores$WC_PC1, data=pcdiat.scores)
model_coef<- coef(model)
#plotpca<- plotpca + geom_abline(intercept=model_coef[1], slope=model_coef[2])
plotpca<- plotpca + labs(x="WC diatoms PC1 scores", y="SSed diatoms PC1 scores")
plotpca<- plotpca + theme_bw()
plotpca<- plotpca + theme(axis.text.x = element_text(colour="black",size=15))
plotpca<- plotpca + theme(axis.text.y = element_text(colour="black",size=15))
plotpca<- plotpca + theme(axis.title.x = element_text(size = rel(1.8), angle=00))
plotpca<- plotpca + theme(axis.title.y = element_text(size = rel(1.8), angle=90))
plotpca<- plotpca + annotate("text", x = 0.4, y = -0.5, label = "Adj R2 = 0.23", size=6)
##

#PCA 3 - environmental variables (water-column)
env.pca<- rda(env.data.red[,#:#]) 
env.pcscores<- data.frame(scores(env.pca, choices = 1, display="sites") 

#PCA 4 - spatial variables (water-column) 
space.pca<- rda(pcnm.red[,#:#])
space.pcscores<- data.frame(scores(space.pca, choices = 1, display="sites") 

#Regression as for PCA 1 and 2. 

#RDAs
#RDA 1 - environmental variables
#Active ordination with wc diatoms and env variables
wc.env.rda<- rda(wc.diat[,#:#], env.data.red[,#:#]) #Add Hellinger if not done
#Extract wc site scores
wc.env.scores<- data.frame(scores(wc.env.rda, choices = 1:2, display="sites")

#Passivley add (weighted-averaging) surface sediment diatoms and extract scores
ssed.envpred.rda<- predict(wc.env.rda, type="wa", newdata=(ssed.diat[,#:#]), scaling=2)
#this produces a set of predicted scores

#RDA 2 - spatial variables
#Active ordination with wc diatoms and spatial variables
wc.space.rda<- rda(wc.diat[,#:#], pcnm.red[,#:#]) #Add Hellinger if not done
#Extract wc site scores
wc.space.scores<- data.frame(scores(wc.space.rda, choices = 1:2, display="sites")

#Passivley add (weighted-averaging) surface sediment diatoms and extract scores
ssed.spacepred.rda<- predict(wc.space.rda, type="wa", newdata=(ssed.diat[,#:#]), scaling=2)
#this produces a set of predicted scores

#OLS regression
#Use extracted scores and extracted predicted scores for lm()
#lm1: wc.env.scores ~ ssed.envpred.rda (RDA1 and 2) 
#lm2: wc.space.scores ~ ssed.spacepred.rda (RDA1 and 2) 
#lm3: wc.env.scores ~ wc.space.scores (RDA1 and 2) 

#Procrustes
#Can compare output from 2 ordinations, either PCA or RDA (if both active)
#Generally extract 4 axes, choices = 1:4 for scores

protest(wc.pcscores, ssed.pcscores) 
residuals(procrustes(wc.pcscores, ssed.pcscores)) 
#number of sites and species in each ordination must match
#can extract residuals and plot separately

###########################################################################

###########################################################################
#Section 6 - Planktonic species only                                      #
###########################################################################
#i) Repeat analyses with planktonic species only. 

#All of the analyses from section 4 and 5 can be repeated using matrix
#of only planktonic species (as per nla2007_diat_classification.csv)

###########################################################################

###########################################################################
#Section 7 - Analyses with time visit data                                #
###########################################################################
#j) Repeat analyses using visit data. 

#k) Bray-Curtis dissimilarities. 

#Time steps 
#use nla2007_env_51.csv, nla2007_diat_51wc.csv and nla2007_diat_51sed.csv
#created average between visits 1 and 2

#Bray-Curtis dissimilarities
#Create dissimilarity matrices using vegdist( , "bray")

############################################################################

