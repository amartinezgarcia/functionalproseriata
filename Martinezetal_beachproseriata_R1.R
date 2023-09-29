######################################
# Martinez, Fontaneto, Curini-Galletti
# Ecological filters in reflective beaches swash zone
# version 1: 16.01.2023 in Verbania
# final version: 18.03.2023 in Arrecife
######################################


library(dplyr)
library(ggplot2)
library(BAT)


# Visualization of spatial objects
library(lme4)
library(car)
library(nlme)
library(lme4)
library(MuMIn)
library(car)
library(vegan)


setwd("~/Dropbox/_Papers/_READY/submitted - Martinez et al - Proseriata western mediterranean/(5) Ecology (sept2023)/Script")


"%ni%" <- Negate("%in%")

############ Part 0: ARRANGING DATA FILES -------------------------------


######## 0.1. Preparation of the community matrix

comm <- read.csv2("community.csv",
                  sep = ";", dec = ".")
      rownames(comm) <- comm$X
      comm <- t(comm[,-1])



######## 0.2. Preparation file of stations data -------------------------------

      station <- read.csv2("stations2.csv",
                         sep = ";", dec = ".")


      ## Set the correct variable types
      station[,c(4:7,16:18)]<-lapply(station[,c(4:7,16:18)],factor)
      str(station)

      # Number of beaches per country for the text
      a <- table(station$Country)
      (a/3) # check console
      rm(a)




######## 0.3. Analyes of climatic and oceanographic variables -------------------------------

        climatic <- read.csv2("Stations_climatic.csv",
                             sep = ";", dec = ".")

        str(climatic)

        #psych::pairs.panels(climatic[2:ncol(climatic)], cex.cor = 5)

        climatic00 <- climatic[c("BIO1_Mean_Temp",
                                 "BIO19_Precipitation.of.Coldest.Quarter",
                                 "BIO4_Temp_Seasonality",
                                 "BIO7_Temp_Range")]

        #psych::pairs.panels(climatic00)


        oceanographic <- read.csv2("Stations_copernicus.csv",
                              sep = ";", dec = ".")

        oceanographic00 <- oceanographic[,c(1,2,5,8,11,14,17,20,23,26,
                                            29,32,35,38,41,44,47)]

        #psych::pairs.panels(oceanographic00[2:ncol(oceanographic00)], cex.cor = 3)


        oceanographic00 <- oceanographic00[c("salinity_avg",
                                      "temperature_avg",
                                      "pprod_avg",
                                      "nitrate_avg",
                                      "northward_avg",
                                      "easteward_avg")]



        macroecological <- cbind(oceanographic00,climatic00)

        write.table(cbind(oceanographic,climatic),
                    "TableS_predictors.csv",
                    sep = ";", dec = ".")

        rm(oceanographic,oceanographic00,climatic,climatic00) #clean the house



######## 0.4. Preparation of the dataset with ecological covariates -------------------------------

        comm.env <- station[c("code","name","level","mean","Length","latitude","longitude")]

        richness <- cbind(comm.env,macroecological)
        nspec <- as.data.frame(rowSums(comm))
        nspec <- cbind(row.names(comm),nspec); colnames(nspec) <- c("code","Taxon")

        richness <- merge(nspec,richness, by="code")
        rm(nspec)

        richness[c(5,6,9:18)] <- lapply(richness[c(5,6,9:18)], function(x) c(scale(x)))

        richness$latitude <- as.numeric(richness$latitude)
        richness$longitude <- as.numeric(richness$longitude)

        richness <- richness %>% replace(is.na(.), 0)

        ### Calculation of endemic species

        comm.beach <- rowsum(comm, (seq_len(nrow(comm)) + 2) %/% 3)
        endemic <- as.data.frame(colSums(ifelse(comm.beach>1,1,comm.beach)))
        colnames(endemic) <- "nrecords"
        endemic$names <- rownames(endemic)
        endemic2 <- as.vector(endemic[ which (endemic$nrecords == 1),2])


        comm.beach.endemic <- comm.beach[ , which (colnames(comm.beach) %in% endemic2)]
        comm.beach.endemic <- ifelse(comm.beach.endemic>1,1,comm.beach.endemic)
        beach.endemic.richness <- rowSums(comm.beach.endemic)
        beach.endemic.richness <- rep(beach.endemic.richness, each  = 3)

        richness$nendemic <- beach.endemic.richness

        rm(comm.beach.endemic,endemic,endemic2)

        # Few data for the manuscript
        print(paste0("Range of the species per beach = ",
                     range(richness$Taxon)))


        print(paste0("Number of beach with 0-1 species =",
                     nrow( richness[ which (richness$Taxon < 4),])/length(richness$Taxon)))


        # write.table(richness,"Table_S1.csv", sep = ";", dec = ".", row.names = F) ## Supplementary table


           rm(macroecological, comm.env)

######## 0.5. Preparation of the trait matrix -------------------------------


traits <- read.csv2("traits.csv", dec = ".", sep = ";")

    # setting species names as row names
    row.names(traits) <- traits$species
    traits <- traits[,-c(1:2)]


    # delete the pair columns containing explanation
    traits <- traits[,2*(-c(1:14))]


    # set categorical and continuous traits
    traits[,c(3:14)]<-lapply(traits[,c(3:14)],factor)


    traits00 <- traits[c("length","width", "adhesiveness",
                     "pharynx", "eyes", "body_pigmentation",
                     "statocyst", "ciliation", "penis", "accessory_spiny_organs")]


    traits00$length <- scale(traits00$length)
    traits00$width <- scale(traits00$width)


    str(traits00)


    # Calculation of gower distance matrix

      gower.mat <- gower(traits00)
      euc.pco <- labdsv::pco(gower.mat,k=4)

      SUM <- sum(euc.pco$eig[euc.pco$eig > 0])

      barplot(euc.pco$eig)

      (euc.pco$eig[1]+euc.pco$eig[2]+euc.pco$eig[3]+
          euc.pco$eig[4])/sum(euc.pco$eig[euc.pco$eig > 0])  # 0.792


    # Generation of the gower.based trait matrix
      traits02 <- euc.pco$points # <<-- Trait matrix based on Gower distance
      colnames(traits02) <- c(paste(rep("axes_",4),1:4,sep=''))



    # match community and traits
      traits02 <- traits02[match(colnames(comm),rownames(traits02)),]


rm(gower.mat, traits00) ## keep it tidy
rm(euc.pco)


######## 0.6. Calculation omega parameter -------------------------------

## Omega = Hb (T/(0.0165646*mean grain)-0.825744 )

## T = 5 seg in the Mediterranean (https://doi.org/10.1016/j.ecss.2018.09.027)
## Hb = 20 - 100 cm


morphodynamic <- station %>%
  group_by(name) %>%
  summarise(Granulometry_mean = mean(mean))



morphodynamic$veldec <- (0.0165646*morphodynamic$Granulometry_mean)-0.825744
Ts <- 5
Hbmin <- 20 ## normal high
Hbmax <- 100 ## storm weather

morphodynamic$omega_min <- Hbmin/(Ts*morphodynamic$veldec)
morphodynamic$omega_max <- Hbmax/(Ts*morphodynamic$veldec)

range(morphodynamic$omega_min)
count(morphodynamic[morphodynamic$omega_min < 1,])/length(morphodynamic$omega_min)
count(morphodynamic[morphodynamic$omega_min > 1,])

cor(morphodynamic$Granulometry_mean,morphodynamic$omega_min)
cor(morphodynamic$Granulometry_mean,morphodynamic$omega_max)

write.table(morphodynamic, "Table_omega.csv", sep = ";", dec=".", row.names = F)

rm(morphodynamic)




##########################################################################################
##### Hypothesis 1: Differences in species richness --------------------------------------
##########################################################################################

##### 1.1. Modelling species richness version environmental parameters ----------

f.tax <- formula(Taxon ~ level + mean + Length +
                   salinity_avg + temperature_avg + pprod_avg +
                   nitrate_avg + northward_avg + easteward_avg +
                   BIO1_Mean_Temp + BIO19_Precipitation.of.Coldest.Quarter +
                   BIO7_Temp_Range + BIO4_Temp_Seasonality)



Tax.B0 <- nlme::gls(f.tax, data = richness)
Tax.B0A <- nlme::gls(f.tax, correlation = corSpher(form = ~ longitude + latitude, nugget = TRUE), data = richness)
Tax.B0B <- nlme::gls(f.tax, correlation = corLin(form = ~ longitude + latitude, nugget = TRUE), data = richness)
Tax.B0C <- nlme::gls(f.tax, correlation = corRatio(form = ~ longitude + latitude, nugget = TRUE), data = richness)
Tax.B0D <- nlme::gls(f.tax, correlation = corGaus(form = ~ longitude + latitude, nugget = TRUE), data = richness)
Tax.B0E <- nlme::gls(f.tax, correlation = corExp(form = ~ longitude + latitude, nugget = TRUE), data = richness)

AIC(Tax.B0, Tax.B0A, Tax.B0B, Tax.B0C, Tax.B0D, Tax.B0E)

    anova(Tax.B0,Tax.B0E)
    summary(Tax.B0E)
    Anova(Tax.B0E)


rm(Tax.B0, Tax.B0A, Tax.B0B, Tax.B0C, Tax.B0D, Tax.B0E)
rm(f.tax, f.tax00)



## Figure 1c: Boxplot showing the differences between beach levels


(Fig1c <- ggplot(richness,
             aes(x=level, y=Taxon, fill = level)) +
              geom_boxplot(alpha=0.6) +
              scale_fill_manual(values = c("#1f2a4d","#fe6e02", "#ffcc77")) +
              labs(title="a",x="Habitat", y = "Number species") +
              geom_jitter(aes(stroke=0.5,
                              alpha = 30), width=0.2, height = 0.2) +
              coord_flip() +
              theme_classic() + theme(legend.position = "none"))



  ## relationship between grain size and number of species
  # (Fig1S <- ggplot(richness,
  #           aes(x=mean, y=Taxon)) +
  #           geom_point(alpha=0.6, aes(color=level)) +
  #           scale_colour_manual(values = c("#1f2a4d","#fe6e02", "#ffcc77")) +
  #           labs(title="a",x="Mean grain size (um)", y = "Number species") +
  #           geom_smooth(method=lm , color="black", se=FALSE) +
  #           theme_classic() + theme(legend.position = "none"))

    ### boxplot keeping the beach identity
    (FigS1 <- ggplot(richness,
              aes(x=level, y=Taxon, fill = level)) +
              geom_boxplot(alpha=0.6) +
              scale_fill_manual(values = c("#1f2a4d","#fe6e02", "#ffcc77")) +
              labs(title="a",x="Habitat", y = "Number species") +
              geom_jitter(aes(stroke=0.5,
                         alpha = 30), width=0.2, height = 0.2) +
              theme_classic() + theme(legend.position = "none") +

              geom_line(aes(group = name), alpha = 0.6, colour = "black", data = richness))


rm(Fig1c,Fig1S)

##### 1.2 Modelling species richness at each beach level (Supplementary) ----------

f.tax00 <- formula(Taxon ~ mean + Length +
                     salinity_avg + temperature_avg + pprod_avg +
                     nitrate_avg + northward_avg + easteward_avg +
                     BIO1_Mean_Temp + BIO19_Precipitation.of.Coldest.Quarter +
                     BIO7_Temp_Range + BIO4_Temp_Seasonality)


####  Swash level

uni.swash <- richness[ which (richness$level == "swash"),]

Tax.B0 <- nlme::gls(f.tax00,data = uni.swash) # the model without spatial structure
Tax.B0A <- nlme::gls(f.tax00, correlation = corSpher(form = ~ longitude + latitude, nugget = TRUE), data = uni.swash)
#Tax.B0B <- nlme::gls(f.tax00, correlation = corLin(form = ~ longitude + latitude, nugget = TRUE), data = uni.swash)
Tax.B0C <- nlme::gls(f.tax00, correlation = corRatio(form = ~ longitude + latitude, nugget = TRUE), data = uni.swash)
Tax.B0D <- nlme::gls(f.tax00, correlation = corGaus(form = ~ longitude + latitude, nugget = TRUE), data = uni.swash)
Tax.B0E <- nlme::gls(f.tax00, correlation = corExp(form = ~ longitude + latitude, nugget = TRUE), data = uni.swash)

AIC(Tax.B0, Tax.B0A, Tax.B0C, Tax.B0D,Tax.B0E)
summary(Tax.B0)


rm(Tax.B0, Tax.B0A, Tax.B0D, Tax.B0C, Tax.B0E, uni.swash)


##### Shoaling level

uni.shoaling <- richness[ which (richness$level == "mid"),]

Tax.B0 <- nlme::gls(f.tax00,data = uni.shoaling) # the model without spatial structure
Tax.B0A <- nlme::gls(f.tax00, correlation = corSpher(form = ~ longitude + latitude, nugget = TRUE), data = uni.shoaling)
#Tax.B0B <- nlme::gls(f.tax00, correlation = corLin(form = ~ longitude + latitude, nugget = TRUE), data = uni.shoaling)
Tax.B0C <- nlme::gls(f.tax00, correlation = corRatio(form = ~ longitude + latitude, nugget = TRUE), data = uni.shoaling)
Tax.B0D <- nlme::gls(f.tax00, correlation = corGaus(form = ~ longitude + latitude, nugget = TRUE), data = uni.shoaling)
Tax.B0E <- nlme::gls(f.tax00, correlation = corExp(form = ~ longitude + latitude, nugget = TRUE), data = uni.shoaling)

AIC(Tax.B0, Tax.B0A, Tax.B0C, Tax.B0D,Tax.B0E)
anova(Tax.B0,Tax.B0A)
summary(Tax.B0A)


rm(Tax.B0, Tax.B0A, Tax.B0C, Tax.B0D, Tax.B0E, uni.shoaling)


##### Subtidal level

uni.subtidal <- richness[ which (richness$level == "deep"),]

Tax.B0 <- nlme::gls(f.tax00,data = uni.subtidal) # the model without spatial structure
Tax.B0A <- nlme::gls(f.tax00, correlation = corSpher(form = ~ longitude + latitude, nugget = TRUE), data = uni.subtidal)
#Tax.B0B <- nlme::gls(f.tax00, correlation = corLin(form = ~ longitude + latitude, nugget = TRUE), data = uni.subtidal)
Tax.B0C <- nlme::gls(f.tax00, correlation = corRatio(form = ~ longitude + latitude, nugget = TRUE), data = uni.subtidal)
Tax.B0D <- nlme::gls(f.tax00, correlation = corGaus(form = ~ longitude + latitude, nugget = TRUE), data = uni.subtidal)
Tax.B0E <- nlme::gls(f.tax00, correlation = corExp(form = ~ longitude + latitude, nugget = TRUE), data = uni.subtidal)

AIC(Tax.B0, Tax.B0A, Tax.B0C,Tax.B0D,Tax.B0E)
anova(Tax.B0,Tax.B0A)

summary(Tax.B0A)

rm(Tax.B0, Tax.B0A, Tax.B0C,Tax.B0D,Tax.B0E)

rm(uni.subtidal,uni.swash, uni.shoaling)


######## 1.3 Modelling endemism against environmental parameters ----------------

f.end <- formula(nendemic ~ level + mean + Length +
                   salinity_avg + temperature_avg + pprod_avg +
                   nitrate_avg + northward_avg + easteward_avg +
                   BIO1_Mean_Temp + BIO19_Precipitation.of.Coldest.Quarter +
                   BIO7_Temp_Range + BIO4_Temp_Seasonality)



End.B0 <- nlme::gls(f.end, data = richness)
#End.B0A <- nlme::gls(f.end, correlation = corSpher(form = ~ longitude + latitude, nugget = TRUE), data = richness)
End.B0B <- nlme::gls(f.end, correlation = corLin(form = ~ longitude + latitude, nugget = TRUE), data = richness)
End.B0C <- nlme::gls(f.end, correlation = corRatio(form = ~ longitude + latitude, nugget = TRUE), data = richness)
End.B0D <- nlme::gls(f.end, correlation = corGaus(form = ~ longitude + latitude, nugget = TRUE), data = richness)
End.B0E <- nlme::gls(f.end, correlation = corExp(form = ~ longitude + latitude, nugget = TRUE), data = richness)

AIC(End.B0, End.B0B, End.B0C, End.B0D, End.B0E)

anova(End.B0,End.B0E)
summary(End.B0E)
Anova(End.B0E)


rm(End.B0, End.B0B, End.B0C, End.B0D, End.B0E)
rm(f.end)

##### 1.4 Figure 1A: Plotting richness and endemism in the map--------------------------------

comm.beach <- rowsum(comm, (seq_len(nrow(comm)) + 2) %/% 3)
comm.beach <- ifelse(comm.beach>1,1,comm.beach)
  beach.richness <- rowSums(comm.beach)

endemic <- as.data.frame(colSums(comm.beach))
    colnames(endemic) <- "nrecords"
    endemic$names <- rownames(endemic)
    endemic2 <- as.vector(endemic[ which (endemic$nrecords == 1),2])


comm.beach.endemic <- comm.beach[ , which (colnames(comm.beach) %in% endemic2)]
    comm.beach.endemic <- ifelse(comm.beach.endemic>1,1,comm.beach.endemic)
    beach.endemic.richness <- rowSums(comm.beach.endemic)


map.richness <-   data.frame(beach = unique(richness$name),
                            species = beach.richness,
                            nonendemic = (beach.richness - beach.endemic.richness)/beach.richness,
                            endemic = beach.endemic.richness/beach.richness,
                            latitude = richness[seq(1, nrow(richness), 3), 7],
                            longitude = richness[seq(1, nrow(richness), 3), 8],
                            col1 = rep("black",length(beach.richness)),
                            col2 =rep("white",length(beach.richness)))


pies <- read.csv2("pies.csv")

pies<- pies[c("beach","pie.x","pie.y")]

map.richness2 <- merge(map.richness, pies, by = "beach")
map.richness2$pie.x <- as.numeric(map.richness2$pie.x)
map.richness2$pie.y <- as.numeric(map.richness2$pie.y)
map.richness2 <- map.richness2[complete.cases(map.richness2),]

rm(endemic, comm.beach.endemic,beach.richness)


WestMed<-marmap::getNOAA.bathy(lon1=-6,lon2=18, lat1=33,lat2=46,resolution=4)

      blues<-c("lightsteelblue4","lightsteelblue3", "lightsteelblue2","lightsteelblue1")
      greys<-c(grey(0.6),grey(0.93),grey(0.99))



plot(WestMed, image = TRUE,
     land = TRUE, n=1,
     bpal = list(c(0, max(WestMed), greys), c(min(WestMed), 0, blues)))

for (i in 1:nrow(map.richness2)){
  space.pies(map.richness2$longitude[i], map.richness2$latitude[i],
             pie.slices = map.richness2[i,c(3:4)],
             pie.colors=map.richness2[i,c(7:8)],
             pie.radius=map.richness2$species[i]/25,
             pie.space =10,
             seg.lwd = 0.4,
             coord = map.richness2[i,c(9:10)])
}

## Requieres editting outside R
points(beaches$good.lon, beaches$good.lat,
       pch = 21, col = "black", bg = "yellow", cex = 0.5, lwd = 0.2)


rm(blues,endemic2,grys,SUM)
rm(map.richness,map.richness2,pies,beach.endemic.richness)

############# ############# ############# ############# ############# ############# #############
############# Hypothesis 2: Relation between traits and levels  --------------------
############# ############# ############# ############# ############# ############# #############


###### 2.1 Main test: Differences in trait composition across levels --------------------------------

richness01 <- richness[c("level","name","mean","easteward_avg")]
rownames(richness01) <- richness$code


Trait.multi <- mvabund::traitglm(comm, richness01, traits02, method="glm1path") ## Go for a drink
Trait.multi$fourth.corner
anova_trait.multi <- anova(Trait.multi)



######## 2.2 Variance partition between geology and geography ----------------------------------------

library(ade4)
library(spdep)
library(adespatial)
library(adegraphics)

comm.xy <- richness[c("longitude", "latitude")]; colnames(comm.xy) <- c("x","y")


## Neightbour matrix

nbgab <- graph2nb(gabrielneigh(comm.xy), sym = TRUE)
distgab <- nbdists(nbgab, comm.xy)
fdist <- lapply(distgab, function(x) 1 - x/max(dist(comm.xy)))
listwgab <- nb2listw(nbgab, glist = fdist) ## weighted neighbour matrix


mem.gab <- mem(listwgab) ### moran eigenvector matrix


pca.hell <- dudi.pca(comm, scale = FALSE, scannf = FALSE, nf = 2)
moran.randtest(pca.hell$li, listw = listwgab)

# test.rda <- randtest(rda.hell)
# test.rda


## select of the MEMs
mem.gab.sel <- mem.select(pca.hell$tab, listw = listwgab)
mem.gab.sel$global.test
mem.gab.sel$summary


## Selection amongst different spatial weighting matrix
cand.lw <- listw.candidates(comm.xy, nb = c("gab", "rel"), weights = c("bin", "flin"))
sel.lw <- listw.select(pca.hell$tab, candidates = cand.lw, nperm = 99)
sel.lw$candidates
sel.lw$best.id
lw.best <- cand.lw[[sel.lw$best.id]]
sel.lw$best


rda.hell <- pcaiv(pca.hell, sel.lw$best$MEM.select, scannf = FALSE)
test.rda <- randtest(rda.hell)
test.rda
plot(test.rda)

nrow(sel.lw$best$MEM.selec)
nrow(richness[,c(4,5,8:18)])

vp1 <- varpart(pca.hell$tab, richness[,c(4,5,8:18)], sel.lw$best$MEM.select)
vp1
plot(vp1, bg = c(3, 5), Xnames = c("environment", "spatial"))



######### 2.3 Effect of ecology ------------------------------------------------------------

library(mvabund)

comm <- mvabund(comm)
meanvar.plot(comm)


mod.multi <- manyglm(comm ~ level + name + mean + easteward_avg,
                  family = "negative_binomial", data = richness) ## Takes 8 hrs

plot(mod.multi)

anova_mod.multi_adj <-  anova(mod.multi, p.uni = "adjusted")

        saveRDS(anova_mod.multi_adj,"anova_mod.multi_adj.rds")
        anova <- readRDS("anova_mod.multi_adj.rds")


  uni.p <- as.data.frame(t(anova$uni.p))
  uni.p$species <- row.names(uni.p)

  uni.test <- as.data.frame(t(anova$uni.test))
  uni.test$species <- row.names(uni.test)

  uni.analyses <- merge(uni.test,uni.p, by = "species")

  colnames(uni.analyses) <- c("species","Dev.intercept","Dev.level",
                              "Dev.beach","Dev.mean","Dev.easteward",
                              "p.intercept","p.level","p.name","p.mean",
                              "p.easteward")


write.table(uni.analyses,"anovamulti_species.csv",row.names = F, sep=";", dec=".")
write.table(anova_mod.multi_adj$table,"anovamulti_summary.csv",row.names = F, sep=";", dec=".")






############ 2.4 Effect of geography  -------------------------------



## Geographical distances using the centroid

  WestMed<-marmap::getNOAA.bathy(lon1=-6,lon2=18, lat1=33,lat2=46,resolution=1)

  beaches <- read.csv2("beaches_coordinatesall.csv", dec = ".") ## coordinates corrected for Marmap
  beaches = beaches[seq(1, nrow(beaches), 3), ]
  colnames(beaches) <- c("name","good.lon","good.lat")


  richness <- merge(richness,beaches, by="name", all.x = T)
  richness <- richness[order(richness$code),]

  depths <- get.depth(WestMed,richness$good.lon, richness$good.lat, locator = F)
  depths$depth > -1 ## all good

  rm(depths)

## Euclidean distances between beaches

  d.euclidean <- as.matrix(geodist::geodist(richness[c("good.lon", "good.lat")],
                                             measure = "geodesic"))


  station.level = richness[ which (richness$level == "swash"),]
  d.euclidean.level <- as.matrix(geodist::geodist(station.level[c("good.lon", "good.lat")],
                                                   measure = "geodesic"))


  ## Least cost distance between beaches, considering 1,100 m depth

  dist.lc <- marmap::lc.dist(trans.mat(WestMed, min.depth= -1,  max.depth =-200),
                             richness[c("good.lon", "good.lat")], res = "dist")


  station.level = richness[ which (richness$level == "swash"),]
  dist.lc.level <-  marmap::lc.dist(trans.mat(WestMed, min.depth= -1,  max.depth =-200),
                                    station.level[c("good.lon", "good.lat")], res = "dist")

### Beta diversity


  st.swash <- richness[ which (richness$level == "swash"),]
  st.shoaling <- richness[ which (richness$level == "mid"),]
  st.subtidal <- richness[ which (richness$level == "deep"),]

    beta.tax <-   BAT::beta(comm, abund=F)
    beta.swash <- BAT::beta(comm[which(rownames(comm) %in% st.swash$code),], abund=F)
    beta.shoaling <-   BAT::beta(comm[which(rownames(comm) %in% st.shoaling$code),], abund=F)
    beta.subtidal <-  BAT::beta(comm[which(rownames(comm) %in% st.subtidal$code),], abund=F)


###### Mantel tests

     vegan::mantel(beta.tax$Btotal, dist.lc,
                   method = "spearman", permutations = 999)

     vegan::mantel(beta.swash$Btotal, dist.lc.level,
                   method = "spearman", permutations = 999)

     vegan::mantel(beta.shoaling$Btotal, dist.lc.level,
                   method = "spearman", permutations = 999)

     vegan::mantel(beta.subtidal$Btotal, dist.lc.level,
                   method = "spearman", permutations = 999)



##### 2.5. Figure 3: Mantel tests relationships --------------


dist.lc.level.m <- as.matrix(dist.lc.level)
levels <- c("st.swash", "st.shoaling", "st.subtidal")



distance.swash <- data.frame()
for (i in 1:ncol(dist.lc.level.m)){

  pairs <- data.frame(distance = dist.lc.level.m[,i],
                      station1 = rep(st.swash$code[i], ncol(dist.lc.level.m)),
                      station2 = st.swash$code)
  distance.swash <- rbind(distance.swash,pairs)
  }



betasw <- as.matrix(beta.swash$Btotal)
beta.swash <- data.frame()
for (i in 1:ncol(betasw)){
  pairs <- data.frame(beta = betasw[,i],
                      station1 = rep(colnames(betasw)[i],length(betasw[,i])),
                      station2 = rownames(betasw))
  beta.swash <- rbind(beta.swash,pairs)}

distance.swash <- cbind(beta.swash,distance.swash)
distance.swash <- distance.swash[,c(1:4)]

dist1 <- ggplot(distance.swash,
             aes(x=distance, y=beta)) +
              geom_point(alpha=0.4, color = "#ffcc77") +
              labs(title="a",x="Beta swash", y = "Distance swash") +
              geom_smooth(method=lm , color="black", se=FALSE) +
              theme_classic() + theme(legend.position = "none")




distance.mid <- data.frame()
  for (i in 1:ncol(dist.lc.level.m)){

    pairs <- data.frame(distance = dist.lc.level.m[,i],
                        station1 = rep(st.shoaling$code[i], ncol(dist.lc.level.m)),
                        station2 = st.shoaling$code)
    distance.mid <- rbind(distance.mid,pairs)
  }



betash <- as.matrix(beta.shoaling$Btotal)
beta.sh <- data.frame()
for (i in 1:ncol(betash)){
  pairs <- data.frame(beta = betash[,i],
                      station1 = rep(colnames(betash)[i],length(betash[,i])),
                      station2 = rownames(betash))
  beta.sh <- rbind(beta.sh,pairs)}


distance.mid <- cbind(beta.sh,distance.mid)
distance.mid <- distance.mid[,c(1:4)]


dist2 <- ggplot(distance.mid,
                aes(x=distance, y=beta)) +
  geom_point(alpha=0.4, color = "#fe6e02") +
  labs(title="a",x="Beta breaking", y = "Distance breaking") +
  geom_smooth(method=lm , color="black", se=FALSE) +
  theme_classic() + theme(legend.position = "none")


distance.deep <- data.frame()
  for (i in 1:ncol(dist.lc.level.m)){
    pairs <- data.frame(distance = dist.lc.level.m[,i],
                        station1 = rep(st.subtidal$code[i], ncol(dist.lc.level.m)),
                        station2 = st.subtidal$code)
    distance.deep <- rbind(distance.deep,pairs)
  }



betasub <- as.matrix(beta.subtidal$Btotal)
beta.sub <- data.frame()
for (i in 1:ncol(betasub)){
  pairs <- data.frame(beta = betasub[,i],
                      station1 = rep(colnames(betasub)[i],length(betash[,i])),
                      station2 = rownames(betasub))
  beta.sub <- rbind(beta.sub,pairs)}

distance.deep <- cbind(beta.sub,distance.deep)
distance.deep <- distance.deep[,c(1:4)]



dist3 <- ggplot(distance.deep,
                aes(x=distance, y=beta)) +
  geom_point(alpha=0.4, color = "#1f2a4d") +
  labs(title="c", x="Beta Shoaling", y = "Distance Shoaling") +
  geom_smooth(method=lm , color="black", se=FALSE) +
  theme_classic() + theme(legend.position = "none")



gridExtra::grid.arrange(dist1,dist2,dist3,ncol=3)




rm(dist1,dist2,dist3,d.euclidean,d.euclidean.level,
   dist.lc,dist.lc.level, dist.lc.level.m, dist.lc1)
rm(distance.deep,distance.k, distance.mid,distance.swash, distances.centroid,level.k)
rm(st.shoaling,st.swash,st.subtidal,station.level,station2, pairs,code.swash)
rm(beatsh,betashoaling,betasub,betasw,beta.tax,beta.sh,
   beta.shoaling,beta.sub,beta.subtidal, beta.swash,betash,beta.mid)
rm(i,k,greys,blues,beach.endemic.richness,endemic2,SUM,beach.level,distance,levels)

############ ############ ############ ############ ############ ############ ############ ############
############ Hypothesis 3: Differences in trait space ---------------------------------------------
############ ############ ############ ############ ############ ############ ############ ############


##### 3.1. Hypervolume calculation ------------------------------------------------------------------------

comm.hab <- rbind(swash = colSums(comm[which (row.names(comm) %in%
                                    as.vector(station[ which (station$level == "swash"),1])),]),
                  mid = colSums(comm[which (row.names(comm) %in%
                                    as.vector(station[ which (station$level == "mid"),1])),]),
                  deep = colSums(comm[which (row.names(comm) %in%
                                    as.vector(station[ which (station$level == "deep"),1])),]))



  kernel.hab.abund <- kernel.build(comm=comm.hab,trait=traits02, abund=TRUE,cores=3,
                              method="gaussian")


  rich.hab.abund <- kernel.alpha(kernel.hab.abund)

          Alpha.hab.abund <- data.frame(Richness=rich.hab.abund,
                                        habitat=c("swash","mid","deep"))




######## 3.2. Null Modelling ----------------------------------


source("abspre.R")

#
# comm.hab1 <- comm.hab %>% replace(comm.hab > 0, 1)
#
# species.freq <- as.data.frame(t(comm.hab)/348)
# species <- as.vector(unique(rownames(species.freq)))
#
# habitat.rich <-  rowSums(comm.hab1)
# str(habitat.rich)
#
# nreplicates = 10
#
# for (k in 1:2){
#
#   random.communities.all <- data.frame()
#
#   for (i in 1:length(habitat.rich)){
#
#     richness.i <- habitat.rich[i]
#     habitat.i <- names(habitat.rich)[i]
#
#     random.comm.habitat.i <-  data.frame()
#
#     for (p in 1:nreplicates){
#       random.community <- as.vector(sample(species, richness.i))
#       random.community01 <- data.frame(species=random.community,
#                                        replicate = rep(paste0(habitat.i,p),richness.i),
#                                        island = rep(habitat.i,richness.i))
#       random.comm.habitat.i <- rbind(random.community01,random.comm.habitat.i)
#     }
#     random.communities.all <- rbind(random.comm.habitat.i,random.communities.all)
#   }
#
#   comm.random <- abspres(data=random.communities.all, sites.col="replicate", sp.col = "species", keep.n=F)
#
#   row.names(comm.random) <- comm.random[,1]; comm.random <- comm.random[,-1]
#
#   kernelFD.random <- kernel.build(comm=comm.random,trait=traits02,abund=FALSE,cores=3,method="gaussian")
#
#   saveRDS(kernelFD.random, file = paste0("kernelFD.random.proseriata",k,".rds"))
#
#
#
#   # Extracting hypervolumes for each habitat
#   HVr.swash <- hypervolume::hypervolume_join(kernelFD.random@HVList[1:nreplicates])
#   HVr.mid <- hypervolume::hypervolume_join(kernelFD.random@HVList[(nreplicates+1):(2*nreplicates)])
#   HVr.deep <- hypervolume::hypervolume_join(kernelFD.random@HVList[(2*nreplicates+1):(3*nreplicates)])
#
#
#   rich.random <- kernel.alpha(kernelFD.random)
#   even.random <- kernel.evenness(comm=kernelFD.random)
#   habitat.random <- c(rep("swash",nreplicates),rep("mid",nreplicates),rep("deep",nreplicates))
#
#   Alpha.null <- data.frame(island=habitat.random,Richness=rich.random,
#                            Even=even.random)
#   write.csv2(Alpha.null,paste0("Alpha.null.new",k,".csv"),sep=";", dec="." )
#   print(paste(k, "of", 10))
# }


Alpha.null <- read.csv2("Nullmodelling2.csv", dec = ".")



## Swash zone
null.swash <- Alpha.null[which(Alpha.null$island=="swash"),]
mean(null.swash$Richness)
ses(Alpha.hab.abund$Richness[1], as.vector(null.swash$Richness))

Null.swash.R <- ggplot(null.swash, aes(x=Richness)) +
  geom_histogram(fill="#ffcc77") +
  geom_vline(aes(xintercept=Alpha.hab.abund$Richness[1]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=Alpha.hab.abund$Richness.weighted[1]),color="black", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean(Richness)), color="red", linetype="dashed", size=0.5)+
  xlab("Functional richness")+ylab("Density")+labs(title="Swash zone")+
  theme_bw()+ theme(legend.position = "none",plot.title = element_text(face="bold"))


## Shoaling zone

null.mid <- Alpha.null[which(Alpha.null$island=="mid"),]
ses(Alpha.hab.abund$Richness[2], as.vector(null.mid$Richness))

mean(null.mid$Richness)

Null.mid.R <- ggplot(null.mid, aes(x=Richness)) +
  geom_histogram(fill="#fe6e02") +
  geom_vline(aes(xintercept=Alpha.hab.abund$Richness[2]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=Alpha.hab.abund$Richness.weighted[2]),color="black", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean(Richness)), color="red", linetype="dashed", size=0.5)+
  xlab("Functional richness")+ylab("Density")+labs(title="Slope")+
  theme_bw()+ theme(legend.position = "none",plot.title = element_text(face="bold"))



## Subtidal zone
null.deep <- Alpha.null[which(Alpha.null$island=="deep"),]
ses(Alpha.hab.abund$Richness[3], as.vector(null.deep$Richness))

mean(null.deep$Richness)

Null.deep.R <- ggplot(null.deep, aes(x=Richness)) +
  geom_histogram(fill="#1f2a4d") +
  geom_vline(aes(xintercept=Alpha.hab.abund$Richness[3]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=Alpha.hab.abund$Richness.weighted[3]),color="black", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean(Richness)), color="red", linetype="dashed", size=0.5)+
  xlab("Functional richness")+ylab("Density")+labs(title="Deep")+
  theme_bw()+ theme(legend.position = "none",plot.title = element_text(face="bold"))


gridExtra::grid.arrange(Null.swash.R, Null.mid.R, Null.deep.R, ncol = 3, nrow=1)

rm(null.deep,Null.deep.R,null.mid,Null.mid.R,null.swash,Null.swash.R)
rm(kernel.hab.abund,this.station)

###### 3.3. Species functional contribution -------------------------------------


#contri.pro.abund2 <- kernel.contribution(kernel.hab.abund, func="one out")
# write.csv2(contri.pro.abund2,"(0) Data/contribution_proseriata_hab_abund_oneout.csv")

contribution <- read.csv2("contribution_proseriata_hab_abund_oneout.csv",
                          dec = ",", row.names = 1)


contribution.t <- as.data.frame(t(contribution))
contribution.t[is.na(contribution.t)] <- 0


contribution.v <- data.frame()
for (i in 1:ncol(contribution.t)){
  this.col <- as.numeric(contribution.t[,i])
  this.col.name <- rep((colnames(contribution.t)[i]),nrow(contribution.t))
  this.station <- data.frame(row.names(contribution.t),this.col.name,this.col)
  contribution.v <- rbind(contribution.v,this.station)}

colnames(contribution.v) <- c("ID","habitat","contribution")

contribution.v$habitat = factor(contribution.v$habitat)
contribution.v <- contribution.v[!is.na(contribution.v$contribution),]

contribution.v <- contribution.v[ which (contribution.v$contribution != 0), ]


con4 <- (ggplot(data=contribution.v, aes(x=contribution, group=habitat, fill=habitat)) +
           geom_density(alpha=.4) +
           xlab("species contribution") +
           ylab("Density")+
           scale_fill_manual(values = c("#ffcc77","#fe6e02", "#1f2a4d"))+
           labs(title="a") +
           theme_bw() + theme(legend.position = "none"))

std <- function(x) sd(x)/sqrt(length(x))

contribution.sum <- contribution.v %>%
  group_by(habitat) %>%
  dplyr::summarise(average = mean(contribution),
                   se = std(contribution),
                   minimun = min(contribution),
                   maximum = max(contribution))


contribution.model <- aov(contribution ~ habitat, data = contribution.v)
performance::check_model(contribution.model)
TukeyHSD(contribution.model)


rm(contribution.model,contribution.sum,contribution.t,contribution.v)
rm(Alpha.hab.abund,Alpha.null,con4,contribution)



############ ############ ############ ############ ############ ############ ############ ############
############ Hypothesis 4: Functional traits analyses  --------- ------------ ------------ ------------
############ ############ ############ ############ ############ ############ ############ ############


traits.discrete <- traits[,c(3:14)]

      probabilities <- list()
      n_trait <- ncol(traits.discrete)
      n_communities  <- nrow(comm)

      for (i in 1 : n_trait){

        n_states <- as.character(unique(traits.discrete[,i]))
        sum_states <- length(n_states)

        name <- colnames(traits.discrete)[i]

        traits.value <- matrix(NA,n_communities,(2*sum_states)+1)
        #traits.value <- matrix(NA,n_communities,sum_states)
        n_total <- append(paste("perc", n_states, sep="_"),"total")
        colnames(traits.value) <- append(n_states, n_total)
        rownames(traits.value) <- rownames(comm)
        #colnames(traits.value) <- n_states


        for (j in 1 : sum_states){

          state.j <- as.data.frame(traits.discrete[traits.discrete[,i] == n_states[j],])
          state.names <- row.names(state.j)
          states <- c()
          total <- c()
          perc <- c()

          for (k in 1:n_communities) {

            this.community <- comm[k,]
            state.values <- this.community[ which (names(this.community) %in% state.names)]
            perc    <- append(perc, (sum(state.values)/sum(this.community))) # calculating the probability of a given state by abundance
            states <- append(states, sum(state.values))
            total <- append(total, sum(this.community))
          }

          traits.value[,j] <- states
          traits.value[,sum_states+j] <- perc
          traits.value[,(2*sum_states)+1] <- total
          traits.value[is.na(traits.value)] <- 0
        }

        probabilities[[i]] <- traits.value #storing the result

      }

      names(probabilities) <- colnames(traits.discrete)


      rm(traits.value, state.j)

library(MASS)
library(emmeans)


###### 4.1. Stickinesss -------------------------

### We merge c(0, 1), poorly adhesive; c(2,3) adhesive

stickiness <- probabilities[[1]]

      stickiness <- data.frame(non_adhesive = stickiness[,2] + stickiness[,3],
                                 adhesive = stickiness[,1] + stickiness[,4],
                                 station = rownames(comm),
                                 habitat = station$level)
      
      stickiness <- merge(stickiness, station[,c(1,3)], by.x = "station", by.y="code")

      stickiness$habitat <- factor(stickiness$habitat,
                                   levels=c("swash","mid","deep"))


      stickiness <- stickiness %>% replace(is.na(.), 0)


        t1 <- ggplot(stickiness[ which (stickiness$non_adhesive != 0 &
                                          stickiness$adhesive != 0 ),],
                           aes(x=habitat,
                           y=adhesive/(non_adhesive + adhesive),
                           fill = habitat)) +
                  geom_boxplot(alpha=0.6, cex = 0.2) +
                  scale_fill_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                  labs(x="", y = "") +
                  geom_jitter(aes(colour = habitat), width=0.2) +
                  scale_colour_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                  coord_flip() +
                  scale_y_continuous(position = "right",
                                     breaks = seq(0,1, by = 0.2)) +
                  theme_classic() + theme(legend.position = "none",
                                          axis.text.y=element_blank())
        
   m.sticky <- glmer(cbind(adhesive, non_adhesive) ~ habitat + (1|name), family=binomial(logit), data=stickiness)
          Anova(m.sticky)
          summary(multcomp::glht(m.sticky, mcp(habitat="Tukey")))
          


      rm(m.sticky,m_means,stickiness)

###### 4.2. Body shapped -------------------------

      ## flattened yes or no

      shape <- probabilities[[4]]

      shape <- data.frame(flat = shape[,2],
                          no_flat = shape[,1] + shape[,3],
                          station = rownames(comm),
                          habitat = station$level)
      
      shape <- merge(shape, station[,c(1,3)], by.x = "station", by.y="code")

      shape <- shape %>% replace(is.na(.), 0)
      shape$habitat <- factor(shape$habitat,
                              levels=c("deep","mid","swash"))



      t2 <- ggplot(shape[ which (shape$flat != 0 &
                                   shape$no_flat != 0 ),],
                     aes(x=habitat, y=flat/(no_flat + flat), fill = habitat)) +
                    geom_boxplot(alpha=0.6, cex = 0.2) +
                    scale_fill_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                    labs(x="", y = "") +
                    geom_jitter(aes(colour = habitat), width=0.2) +
                    scale_colour_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                    coord_flip() +
                    scale_y_continuous(position = "right",
                                       breaks = seq(0,1, by = 0.2)) +
                    theme_classic() + theme(legend.position = "none",
                                            axis.text.y=element_blank())


      m.flat <- glmer(cbind(flat, no_flat) ~ habitat + (1|name), family=binomial(logit), data=shape)
      Anova(m.flat)
      summary(multcomp::glht(m.flat, mcp(habitat="Tukey")))

      rm(m.flat,m_means,shape)


###### 4.3. Cephalic sensory area  -------------------------

      sensory <- probabilities[[5]]; head(sensory)

      sensory <- data.frame(present = sensory[,2],
                          absent = sensory[,1],
                          station = rownames(comm),
                          habitat = station$level)
      
      sensory <- merge(sensory, station[,c(1,3)], by.x = "station", by.y="code")

      sensory <- sensory %>% replace(is.na(.), 0)

      sensory$habitat <- factor(sensory$habitat,
                                levels=c("deep","mid","swash"))


      t3 <- ggplot(sensory <- sensory[ which (sensory$present != 0 &
                                                sensory$absent != 0 ),],
                     aes(x=habitat, y=present/(absent+present), fill = habitat)) +
                    geom_boxplot(alpha=0.6, cex = 0.2) +
                    scale_fill_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                    labs(x="", y = "") +
                    geom_jitter(aes(colour = habitat), width=0.2) +
                    scale_colour_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                    coord_flip() +
                    scale_y_continuous(position = "right",
                                       breaks = seq(0,1, by = 0.2)) +
                    theme_classic() + theme(legend.position = "none",
                                            axis.text.y=element_blank())

      m.sensory <- glmer(cbind(present, absent) ~ habitat + (1|name), family=binomial(logit), data=sensory)
      Anova(m.sensory)
      summary(multcomp::glht(m.sensory, mcp(habitat="Tukey")))

      rm(m.sensory,m_means,sensory)



###### 4.4. Brain capsule  -------------------------

      brain <- probabilities[[3]]; head(brain)

      brain <- data.frame(present = brain[,2],
                            absent = brain[,1] + brain[,3],
                            station = rownames(comm),
                            habitat = station$level)

      brain <- merge(brain, station[,c(1,3)], by.x = "station", by.y="code")
      
      brain$habitat <- factor(brain$habitat,
                              levels=c("deep","mid","swash"))

      brain <- brain %>% replace(is.na(.), 0)


      t4 <- ggplot(brain[ which (brain$present != 0 &
                                   brain$absent != 0 ),],
                   aes(x=habitat, y=present/(absent+present), fill = habitat)) +
                  geom_boxplot(alpha=0.6, cex = 0.2) +
                  scale_fill_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                  labs(x="", y = "") +
                  geom_jitter(aes(colour = habitat), width=0.2) +
                  scale_colour_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                  coord_flip() +
                  scale_y_continuous(position = "right",
                                     breaks = seq(0,1, by = 0.2)) +
                  theme_classic() + theme(legend.position = "none",
                                          axis.text.y=element_blank())


      m.brain <- glmer(cbind(present, absent) ~ habitat + (1|name), family=binomial(logit), data=brain)
      Anova(m.brain)
      summary(multcomp::glht(m.brain, mcp(habitat="Tukey")))

      rm(m.brain,m_means,brain)


###### 4.5. Ciliary pattern  -------------------------

      cilia <- probabilities[[9]]; head(cilia)

      cilia <- data.frame(complete = cilia[,1],
                          ventral = cilia[,2],
                          station = rownames(comm),
                          habitat = station$level)
      
      cilia <- merge(cilia, station[,c(1,3)], by.x = "station", by.y="code")

      cilia$habitat <- factor(cilia$habitat,
                              levels=c("deep","mid","swash"))

      cilia <- cilia %>% replace(is.na(.), 0)


      t5 <- ggplot(cilia[ which (cilia$complete != 0 &
                                   cilia$ventral != 0 ),],
                   aes(x=habitat, y=ventral/(complete+ventral), fill = habitat)) +
                  geom_boxplot(alpha=0.6, cex = 0.2) +
                  scale_fill_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                  labs(x="", y = "") +
                  geom_jitter(aes(colour = habitat), width=0.2) +
                  scale_colour_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                  coord_flip() +
                  scale_y_continuous(position = "right",
                                     breaks = seq(0,1, by = 0.2)) +
                  theme_classic() + theme(legend.position = "none",
                                          axis.text.y=element_blank())


      m.cilia <- glmer(cbind(ventral, complete) ~ habitat + (1|name), family=binomial(logit), data=cilia)
      Anova(m.cilia)

      summary(multcomp::glht(m.cilia, mcp(habitat="Tukey")))

      rm(m.cilia,m_means,cilia)



###### 4.6. Feeding  -------------------------

      feeding <- probabilities[[2]]; head(feeding)

      feeding <- data.frame(scavenger = feeding[,2],
                          predator = feeding[,1] + feeding[,3],
                          station = rownames(comm),
                          habitat = station$level)
      
      feeding <- merge(feeding, station[,c(1,3)], by.x = "station", by.y="code")
      
      feeding$habitat <- factor(feeding$habitat,
                                levels=c("deep","mid","swash"))


      feeding <- feeding %>% replace(is.na(.), 0)


      t6 <- ggplot(feeding[ which (feeding$scavenger != 0 &
                                   feeding$predator != 0 ),],
                   aes(x=habitat, y=scavenger/(predator+scavenger), fill = habitat)) +
                  geom_boxplot(alpha=0.6, cex = 0.2) +
                  scale_fill_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                  labs(x="", y = "") +
                  geom_jitter(aes(colour = habitat), width=0.2) +
                  scale_colour_manual(values = c("#1f2a4d","#fe6e02","#ffcc77")) +
                  coord_flip() +
                  scale_y_continuous(position = "right",
                                     breaks = seq(0,1, by = 0.2)) +
                  theme_classic() + theme(legend.position = "none",
                                          axis.text.y=element_blank())

     m.feed <- glmer(cbind(scavenger, predator) ~ habitat + (1|name), family=binomial(logit), data=feeding)
      Anova(m.feed)

      summary(multcomp::glht(m.feed, mcp(habitat="Tukey")))

      rm(m.feed,m_means,feeding)



gridExtra::grid.arrange(t3,t4,t5,t2,t6,t1, ncol = 1, nrow=6)

 rm(t1,t2,t3,t4,t5,t6)

 rm(probabilities,traits.discrete)

