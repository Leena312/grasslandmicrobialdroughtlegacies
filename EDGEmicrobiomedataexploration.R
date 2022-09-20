##alpha diversity

library(mctoolsr)
library(vegan)
library(car)
library(metacoder)
library(emmeans)
library(multcompView)



##2018 16S Shannon
input <- load_taxa_table("~/Desktop/R/OTU_16S_withtax_2018.txt", "~/Desktop/R/MappingFile_16S_2018.txt")

##Decision for rarefication

colSums(input$data_loaded)
inputr = single_rarefy(input, 2500)

##By Treatment
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- diversity(t(INT$data_loaded), index = "shannon")
CHRshan <- diversity(t(CHR$data_loaded), index = "shannon")
CONshan <- diversity(t(CON$data_loaded), index = "shannon")

colors<-c("orange","dark gray","blue")

pdf(file="Shannon Diversity 16S treatment 2018",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Shannon's Diversity", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##2018 GLM 16S

allshan <- diversity(t(inputr$data_loaded), index = "shannon")

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")



emmeans(allshan_glm,pairwise~Treatment)



##Richness 16S 2018
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- specnumber(t(INT$data_loaded))
CHRshan <- specnumber(t(CHR$data_loaded))
CONshan <- specnumber(t(CON$data_loaded))

colors<-c("orange","dark gray","blue")

pdf(file="Richness 16S treatment 2018",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Richness", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##2018 GLM 16S species richness

allshan <- specnumber(t(inputr$data_loaded))

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")

emmeans(allshan_glm,pairwise~Treatment)

##Simpson
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- diversity(t(INT$data_loaded), index = "simpson")
CHRshan <- diversity(t(CHR$data_loaded), index = "simpson")
CONshan <- diversity(t(CON$data_loaded), index = "simpson")

colors<-c("orange","dark gray","blue")

pdf(file="Simpson Diversity 16S treatment 2018",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Simpson Diversity", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##2018 GLM 16S

allshan <- diversity(t(inputr$data_loaded), index = "simpson")

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")



emmeans(allshan_glm,pairwise~Treatment)


##2018 ITS Shannon
input <- load_taxa_table("~/Desktop/R/OTU_ITS_with tax_2018.txt", "~/Desktop/R/MappingFile_ITS_2018.txt")

##Decision for rarefication

colSums(input$data_loaded)
inputr = single_rarefy(input, 2600)

##By Treatment
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- diversity(t(INT$data_loaded), index = "shannon")
CHRshan <- diversity(t(CHR$data_loaded), index = "shannon")
CONshan <- diversity(t(CON$data_loaded), index = "shannon")

colors<-c("orange","dark gray","blue")

pdf(file="Shannon Diversity 16S treatment 2018",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Shannon's Diversity", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##2018 GLM 16S

allshan <- diversity(t(inputr$data_loaded), index = "shannon")

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")



emmeans(allshan_glm,pairwise~Treatment)



##Richness ITS 2018
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- specnumber(t(INT$data_loaded))
CHRshan <- specnumber(t(CHR$data_loaded))
CONshan <- specnumber(t(CON$data_loaded))

colors<-c("orange","dark gray","blue")

pdf(file="Richness ITS treatment 2018",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Richness", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##2018 GLM ITS species richness

allshan <- specnumber(t(inputr$data_loaded))

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")

emmeans(allshan_glm,pairwise~Treatment)

##Simpson
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- diversity(t(INT$data_loaded), index = "simpson")
CHRshan <- diversity(t(CHR$data_loaded), index = "simpson")
CONshan <- diversity(t(CON$data_loaded), index = "simpson")

colors<-c("orange","dark gray","blue")

pdf(file="Simpson Diversity ITS treatment 2018",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Simpson's Diversity", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##2018 GLM 16S

allshan <- diversity(t(inputr$data_loaded), index = "simpson")

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")



emmeans(allshan_glm,pairwise~Treatment)

##2019 16S Shannon
input <- load_taxa_table("~/Desktop/R/OTU_16S_withtax_2019.txt", "~/Desktop/R/MappingFile_16S_2019.txt")

##Decision for rarefication

colSums(input$data_loaded)
inputr = single_rarefy(input, 4000)

##By Treatment
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- diversity(t(INT$data_loaded), index = "shannon")
CHRshan <- diversity(t(CHR$data_loaded), index = "shannon")
CONshan <- diversity(t(CON$data_loaded), index = "shannon")

colors<-c("orange","dark gray","blue")

pdf(file="Shannon Diversity 16S treatment 2019",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Shannon's Diversity", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##GLM 16S

allshan <- diversity(t(inputr$data_loaded), index = "shannon")

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")



emmeans(allshan_glm,pairwise~Treatment)



##Richness 16S 2018
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- specnumber(t(INT$data_loaded))
CHRshan <- specnumber(t(CHR$data_loaded))
CONshan <- specnumber(t(CON$data_loaded))

colors<-c("orange","dark gray","blue")

pdf(file="Richness 16S treatment 2019",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Richness", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##2018 GLM 16S species richness

allshan <- specnumber(t(inputr$data_loaded))

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")

emmeans(allshan_glm,pairwise~Treatment)

##Simpson
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- diversity(t(INT$data_loaded), index = "simpson")
CHRshan <- diversity(t(CHR$data_loaded), index = "simpson")
CONshan <- diversity(t(CON$data_loaded), index = "simpson")

colors<-c("orange","dark gray","blue")

pdf(file="Simpson Diversity 16S treatment 2019",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Simpson Diversity", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##2018 GLM 16S

allshan <- diversity(t(inputr$data_loaded), index = "simpson")

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")



emmeans(allshan_glm,pairwise~Treatment)


##2019 ITS Shannon
input <- load_taxa_table("~/Desktop/R/OTU_ITS_withtax_2019.txt", "~/Desktop/R/MappingFile_ITS_2019.txt")

##Decision for rarefication

colSums(input$data_loaded)
inputr = single_rarefy(input, 3900)

##By Treatment
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- diversity(t(INT$data_loaded), index = "shannon")
CHRshan <- diversity(t(CHR$data_loaded), index = "shannon")
CONshan <- diversity(t(CON$data_loaded), index = "shannon")

colors<-c("orange","dark gray","blue")

pdf(file="Shannon Diversity 16S treatment 2019",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Shannon's Diversity", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##2018 GLM 16S

allshan <- diversity(t(inputr$data_loaded), index = "shannon")

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")



emmeans(allshan_glm,pairwise~Treatment)



##Richness ITS 2019
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- specnumber(t(INT$data_loaded))
CHRshan <- specnumber(t(CHR$data_loaded))
CONshan <- specnumber(t(CON$data_loaded))

colors<-c("orange","dark gray","blue")

pdf(file="Richness ITS treatment 2019",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Richness", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##2018 GLM ITS species richness

allshan <- specnumber(t(inputr$data_loaded))

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")

emmeans(allshan_glm,pairwise~Treatment)

##Simpson
INT = filter_data(inputr, 'Treatment', keep_vals = 'INT')
CHR = filter_data(inputr, 'Treatment', keep_vals = 'CHR')
CON = filter_data(inputr, 'Treatment', keep_vals = 'CON')

INTshan <- diversity(t(INT$data_loaded), index = "simpson")
CHRshan <- diversity(t(CHR$data_loaded), index = "simpson")
CONshan <- diversity(t(CON$data_loaded), index = "simpson")

colors<-c("orange","dark gray","blue")

pdf(file="Simpson Diversity ITS treatment 2019",
    width=4,
    height=5)

a <- boxplot(INTshan, CHRshan, CONshan,col=colors, main = "Treatment", ylab = "Simpson's Diversity", xlab = "Treament", names = c("Intense", "Chronic", "Control"))

dev.off()


##2018 GLM 16S

allshan <- diversity(t(inputr$data_loaded), index = "simpson")

allshan <- allshan[match(names(allshan), row.names(inputr$map_loaded))]
allshan_df <- cbind(inputr$map_loaded, allshan)

inputr$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")



emmeans(allshan_glm,pairwise~Treatment)

##CAPS


library(dplyr)
library(BiodiversityR)
library (vegan)
library (MASS)

#2019 16S
input <- load_taxa_table("~/Desktop/R/OTU_16S_withtax_2018.txt", "~/Desktop/R/MappingFile_16S_2018.txt")

inputr = single_rarefy(input, 3000)

otuTAB<-inputr$data_loaded

otuTAB <- as.data.frame (otuTAB)
otuTABtranspose<-t(otuTAB)
mta<-input$map_loaded
mtafinal<-mta[-c(7,27),]
mtafinal$Treatment <- as.factor (mtafinal$Treatment)

o.dist <-vegdist (otuTABtranspose)
Ordination.model1 <- CAPdiscrim (o.dist ~ Treatment, data = mtafinal,
                                 dist="bray", axes=2, m=0, add=FALSE)

pdf(file="16S 2019",
    width=4,
    height=5)

plot1 <- ordiplot(Ordination.model1, type="none")
ordisymbol(plot1, mtafinal, "Treatment", legend=TRUE,legend.x="topright")

ordiellipse(Ordination.model1, groups = mtafinal$Treatment, draw = "polygon", col=c("dark gray","blue","orange"))
dev.off()

##Decision for rarefication
#2018 16S

input <- load_taxa_table("~/Desktop/R/OTU_16S_withtax_2018.txt", "~/Desktop/R/MappingFile_16S_2018.txt")
colSums(input$data_loaded)
inputr = single_rarefy(input, 3000)

otuTAB<-inputr$data_loaded

otuTAB <- as.data.frame (otuTAB)
otuTABtranspose<-t(otuTAB)
mta<-input$map_loaded
mtafinal<-mta[-c(7,14,15,16,25,26,27,30),]
mtafinal$Treatment <- as.factor (mtafinal$Treatment)

o.dist <-vegdist (otuTABtranspose)
Ordination.model1 <- CAPdiscrim (o.dist ~ Treatment, data = mtafinal,
                                 dist="bray", axes=2, m=0, add=FALSE)

pdf(file="16S 2018",
    width=4,
    height=5)


plot1 <- ordiplot(Ordination.model1, type="none")
ordisymbol(plot1, mtafinal, "Treatment", legend=TRUE,legend.x="topleft")
ordiellipse(Ordination.model1, groups = mtafinal$Treatment, draw = "polygon", col=c("dark gray","blue","orange"))
dev.off()

inputr$map_loaded$Treatment <- as.factor (inputr$map_loaded$Treatment)

dm <- calc_dm(t(inputr$data_loaded))
calc_pairwise_permanovas(dm, inputr$map_loaded, "Treatment")


install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") 
library(pairwiseAdonis)
pair.mod<-pairwise.adonis(t(otuTAB),factors=mtafinal$Treatment)
pair.mod


##ITS 2018
input <- load_taxa_table("~/Desktop/R/OTU_ITS_with tax_2018.txt", "~/Desktop/R/MappingFile_ITS_2018.txt")
inputr = single_rarefy(input, 2600)

otuTAB<-inputr$data_loaded

otuTAB <- as.data.frame (otuTAB)
otuTABtranspose<-t(otuTAB)
mta<-input$map_loaded
mtafinal<-mta[-c(7,14,18,21,28,29),]
mtafinal$Treatment <- as.factor (mtafinal$Treatment)

o.dist <-vegdist (otuTABtranspose)
Ordination.model1 <- CAPdiscrim (o.dist ~ Treatment, data = mtafinal,
                                 dist="bray", axes=2, m=0, add=FALSE)

pdf(file="ITS 2018",
    width=4,
    height=5)
v<-c("orange","blue","dark gray")

plot1 <- ordiplot(Ordination.model1, type="none")
ordisymbol(plot1, mtafinal, "Treatment", legend=TRUE,legend.x="topleft",colors=v)
ordiellipse(Ordination.model1, groups = mtafinal$Treatment, draw = "polygon", col=c("dark gray","blue","orange"))
dev.off()

##ITS 2019
input <- load_taxa_table("~/Desktop/R/OTU_ITS_withtax_2019.txt", "~/Desktop/R/MappingFile_ITS_2019.txt")
colSums(input$data_loaded)
inputr = single_rarefy(input, 3200)

otuTAB<-inputr$data_loaded

otuTAB <- as.data.frame (otuTAB)
otuTABtranspose<-t(otuTAB)
mta<-input$map_loaded
mtafinal<-mta[-c(2),]
mtafinal$Treatment <- as.factor (mtafinal$Treatment)
pair.mod<-pairwise.adonis(t(otuTAB),factors=mtafinal$Treatment)
pair.mod
o.dist <-vegdist (otuTABtranspose)
Ordination.model1 <- CAPdiscrim (o.dist ~ Treatment, data = mtafinal,
                                 dist="bray", axes=2, m=0, add=FALSE)

pdf(file="ITS 2019",
    width=4,
    height=5)


plot1 <- ordiplot(Ordination.model1, type="none")
ordisymbol(plot1, mtafinal, "Treatment", legend=TRUE,legend.x="topleft")
ordiellipse(Ordination.model1, groups = mtafinal$Treatment, draw = "polygon", col=c("dark gray","blue","orange"))
dev.off()

##Permanova

library(vegan)
library(devtools)


model1<-adonis(otuTABtranspose~Treatment,data=mtafinal,permutations=999,method="bray")
model1





##nb.glms

library(mctoolsr)
library(edgeR)
library(vegan)
library(MASS)
library(lme4)
library(colorspace)
library(car)
library(doBy)
library(lmerTest)
library(lattice)
library(ggplot2)
library(ggthemes)
library(scales)
library(viridis)
library(RColorBrewer)
library(colorRamps)
library(reshape)
library(Hmisc)
library(emmeans)


zotu_ITS_input<-load_taxa_table("~/Desktop/R/OTU_ITS_with tax_2018.txt", "~/Desktop/R/MappingFile_ITS_2018.txt")
zotu_16S_input <- load_taxa_table("~/Desktop/R/OTU_16S_withtax_2018.txt", "~/Desktop/R/MappingFile_16S_2018.txt")

##for relative abundance
zotu_16S_input$data_loaded = apply(zotu_16S_input$data_loaded,2,function(x){100*x/sum(x)})

zotu_ITS_input$data_loaded = apply(zotu_ITS_input$data_loaded,2,function(x){100*x/sum(x)})



sort(colSums(zotu_16S_input$data_loaded))
sort(colSums(zotu_ITS_input$data_loaded))

##TMM Normalization
#16S
y <- zotu_16S_input$data_loaded
# y <- zotu_16S_input$data_loaded[rowSums(zotu_16S_input$data_loaded,na.rm = TRUE)>500,]
z <- zotu_16S_input$taxonomy_loaded[rownames(zotu_16S_input$taxonomy_loaded) %in% rownames(y),]
x <- DGEList(y)
x <- calcNormFactors(x)
x <- estimateCommonDisp(x)
zotu_16S_TMM <- list(x$pseudo.counts,zotu_16S_input$map_loaded, z)
names(zotu_16S_TMM) <- names(zotu_16S_input)

#ITS
y <- zotu_ITS_input$data_loaded
# y <- zotu_ITS_input$data_loaded[rowSums(zotu_ITS_input$data_loaded,na.rm = TRUE)>100,]
z <- zotu_ITS_input$taxonomy_loaded[rownames(zotu_ITS_input$taxonomy_loaded) %in% rownames(y),]
x <- DGEList(y)
x <- calcNormFactors(x)
x <- estimateCommonDisp(x)
zotu_ITS_TMM <- list(x$pseudo.counts,zotu_ITS_input$map_loaded, z)
names(zotu_ITS_TMM) <- names(zotu_ITS_input)

zotu_16S_TMM$map_loaded$Treatment <- as.factor(zotu_16S_TMM$map_loaded$Treatment)


zotu_ITS_TMM$map_loaded$Treatment <- as.factor(zotu_ITS_TMM$map_loaded$Treatment)


##Taxonomy

tax_NBGLMs <- function(data.set) {
    results <- vector()
    for(i in 1:2){           #splits dataset for each taxonomic group
        print(paste(i, "of 2"))
        input_tax_level = summarize_taxonomy(data.set, level = i, report_higher_tax = FALSE, relative = TRUE)
        test.data.raw <-  cbind(data.set$map_loaded,t(input_tax_level))
        start.col <- dim(data.set$map_loaded)[2] + 1
        
        for(j in start.col:ncol(test.data.raw)){         #Does the glm and lsmean for each taxonomic group
            #print(j)
            
            test.glm <- tryCatch({
                suppressWarnings(glm.nb(test.data.raw[,j] ~ Treatment, data = test.data.raw))
            }, error = function(e){
                return(rep(NA, 8))
            })
            if(is(test.glm) != "negbin") {
                results <- rbind(results, test.glm);
                rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
                next
            }
            
            test.aov <- tryCatch({
                Anova(test.glm, test.statistic = "F")
            }, error = function(e){
                return(rep(NA, 8))
            })#, warning = function(w){
            #return(rep(NA, 20))
            #})
            if(is(test.aov)[1] != "anova") {
                results <- rbind(results, test.aov);
                rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
                next
            }
            
            test.lsm <- emmeans(test.glm, specs = "Treatment")
            pwc<-pairs(test.lsm)
            a <- c(test.aov$'Sum Sq', test.aov$Df, test.aov$F, test.aov$`Pr(>F)`,summary(test.lsm)$emmean, summary(test.lsm)$SE,summary(pwc)$p.value)
            results <- rbind(results, a)
            rownames(results)[nrow(results)] <- colnames(test.data.raw)[j]
            
        }
        
    }
    colnames(results) <- c(apply(expand.grid(rownames(test.aov), colnames(test.aov)), 1, paste, collapse="."), 
                           apply(expand.grid(summary(test.lsm)$Treatment, colnames(summary(test.lsm))[2:3]), 1, paste,collapse="."),
                           summary(pwc)$contrast)
    return(results)
}



#16S
zotu_NBGLM_tax_16S <- tax_NBGLMs(zotu_16S_TMM)

#ITS
zotu_NBGLM_tax_ITS <- tax_NBGLMs(zotu_ITS_TMM)

###Individual OTUs
otu_NBGLMs <- function(data.set) {
    
    test.data.raw <-  cbind(data.set$map_loaded,t(data.set$data_loaded))
    start.col <- dim(data.set$map_loaded)[2] + 1
    
    results <- vector()
    
    for(j in start.col:ncol(test.data.raw)){         #Does the glm and lsmean for each taxonomic group
        
        if(j%%100 == 0) { print(paste(j, "out of", ncol(test.data.raw))) }
        
        if(sum(test.data.raw[,j], na.rm = TRUE) < 100) { next }
        
        #model <- formula(test.data.raw[,j] ~ RepNumber + Accession)
        
        test.glm <- tryCatch({
            suppressWarnings(glm.nb(test.data.raw[,j] ~ Treatment, data = test.data.raw))
        }, error = function(e){
            return(rep(NA, 8))
        })
        if(is(test.glm) != "negbin") {
            #results <- rbind(results, test.glm);
            #rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
            next
        }
        
        test.aov <- tryCatch({
            Anova(test.glm, test.statistic = "F")
        }, error = function(e){
            return(rep(NA, 8))
        })#, warning = function(w){
        #return(rep(NA, 20))
        #}
        if(is(test.aov)[1] != "anova") {
            #results <- rbind(results, test.aov);
            #rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
            next
        }
        
        test.lsm <- LSmeans(test.glm, effect = "Treatment")
        a <- c(test.aov$'Sum Sq', test.aov$Df, test.aov$F, test.aov$`Pr(>F)`,test.lsm$coef[,1], test.lsm$coef[,2])
        results <- rbind(results, a)
        rownames(results)[nrow(results)] <- colnames(test.data.raw)[j]
        
    }
    
    
    colnames(results) <- c(apply(expand.grid(rownames(test.aov), colnames(test.aov)), 1, paste, collapse="."), 
                           apply(expand.grid(test.lsm$grid$Treatment, colnames(test.lsm$coef)[1:2]), 1, paste,
                                 collapse="."))
    return(results)
}

#16S
zotu_NBGLM_raw_16S <- otu_NBGLMs(zotu_16S_TMM)

#ITS
zotu_NBGLM_raw_ITS <- otu_NBGLMs(zotu_ITS_TMM)

##CHANGE NAMES
#Transform LSMean estimates back to counts instead of ln(counts)
#16S
zotu_NBGLM_raw_16S[,grep(".emmean|.SE", colnames(zotu_NBGLM_raw_16S))] <- exp(zotu_NBGLM_raw_16S[,grep(".emmean|.SE", colnames(zotu_NBGLM_raw_16S))])
zotu_NBGLM_tax_16S[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_16S))] <- exp(zotu_NBGLM_tax_16S[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_16S))])

#ITS
zotu_NBGLM_raw_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_raw_ITS))] <- exp(zotu_NBGLM_raw_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_raw_ITS))])
zotu_NBGLM_tax_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_ITS))] <- exp(zotu_NBGLM_tax_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_ITS))])

# FDR Adjustments
#16S
p.cols <- grep(".Pr[:(:][:>:]F[:):]", colnames(zotu_NBGLM_raw_16S))
w <- apply(zotu_NBGLM_raw_16S[,p.cols], 2, function(x){ p.adjust(x, method = 'fdr') })
colnames(w) <- paste(colnames(w), "FDR", sep = ".")
zotu_NBGLM_raw_16S <- cbind(zotu_NBGLM_raw_16S, w)
p.cols <- grep(".Pr[:(:][:>:]F[:):]", colnames(zotu_NBGLM_tax_16S))
w <- apply(zotu_NBGLM_tax_16S[,p.cols], 2, function(x){ p.adjust(x, method = 'fdr') })
colnames(w) <- paste(colnames(w), "FDR", sep = ".")
zotu_NBGLM_tax_16S <- cbind(zotu_NBGLM_tax_16S, w)

#ITS
p.cols <- grep(".Pr[:(:][:>:]F[:):]", colnames(zotu_NBGLM_raw_ITS))
w <- apply(zotu_NBGLM_raw_ITS[,p.cols], 2, function(x){ p.adjust(x, method = 'fdr') })
colnames(w) <- paste(colnames(w), "FDR", sep = ".")
zotu_NBGLM_raw_ITS <- cbind(zotu_NBGLM_raw_ITS, w)
p.cols <- grep(".Pr[:(:][:>:]F[:):]", colnames(zotu_NBGLM_tax_ITS))
w <- apply(zotu_NBGLM_tax_ITS[,p.cols], 2, function(x){ p.adjust(x, method = 'fdr') })
colnames(w) <- paste(colnames(w), "FDR", sep = ".")
zotu_NBGLM_tax_ITS <- cbind(zotu_NBGLM_tax_ITS, w)

rm(p.cols, w)




# Eta-Squared Calcs
####May want to combine this process with the above functions.
etaSq.calc <- function(nbglm.data){
    results <- vector()
    SS.cols <- grep("Sum Sq$", colnames(nbglm.data))
    for(i in 1:nrow(nbglm.data)){
        eta.sq <- nbglm.data[i,SS.cols] / sum(nbglm.data[i,SS.cols], na.rm = T)
        results <- rbind(results, eta.sq)
    }
    colnames(results) <- gsub("Sum Sq", "etaSq", colnames(results))
    rownames(results) <- rownames(nbglm.data)
    return(results)
}

#16S
zotu_NBGLM_raw_16S <- cbind(zotu_NBGLM_raw_16S,  etaSq.calc(zotu_NBGLM_raw_16S))
zotu_NBGLM_tax_16S <- cbind(zotu_NBGLM_tax_16S,  etaSq.calc(zotu_NBGLM_tax_16S))

#ITS
zotu_NBGLM_raw_ITS <- cbind(zotu_NBGLM_raw_ITS,  etaSq.calc(zotu_NBGLM_raw_ITS))
zotu_NBGLM_tax_ITS <- cbind(zotu_NBGLM_tax_ITS,  etaSq.calc(zotu_NBGLM_tax_ITS))



#Write Results to CSV
write.csv(zotu_NBGLM_raw_16S, "zotu_NBGLM_raw_16S_2018.csv")
write.csv(zotu_NBGLM_tax_16S, "zotu_NBGLM_tax_16S_2018.csv")

write.csv(zotu_NBGLM_raw_ITS, "zotu_NBGLM_raw_ITS_2018.csv")
write.csv(zotu_NBGLM_tax_ITS, "zotu_NBGLM_tax_ITS_2018.csv")

##2019

zotu_ITS_input<-load_taxa_table("~/Desktop/R/OTU_ITS_withtax_2019.txt", "~/Desktop/R/MappingFile_ITS_2019.txt")
zotu_16S_input <- load_taxa_table("~/Desktop/R/OTU_16S_withtax_2019.txt", "~/Desktop/R/MappingFile_16S_2019.txt")

zotu_16S_input$data_loaded = apply(zotu_16S_input$data_loaded,2,function(x){100*x/sum(x)})

zotu_ITS_input$data_loaded = apply(zotu_ITS_input$data_loaded,2,function(x){100*x/sum(x)})

sort(colSums(zotu_16S_input$data_loaded))
sort(colSums(zotu_ITS_input$data_loaded))

##TMM Normalization
#16S
y <- zotu_16S_input$data_loaded
# y <- zotu_16S_input$data_loaded[rowSums(zotu_16S_input$data_loaded,na.rm = TRUE)>500,]
z <- zotu_16S_input$taxonomy_loaded[rownames(zotu_16S_input$taxonomy_loaded) %in% rownames(y),]
x <- DGEList(y)
x <- calcNormFactors(x)
x <- estimateCommonDisp(x)
zotu_16S_TMM <- list(x$pseudo.counts,zotu_16S_input$map_loaded, z)
names(zotu_16S_TMM) <- names(zotu_16S_input)

#ITS
y <- zotu_ITS_input$data_loaded
# y <- zotu_ITS_input$data_loaded[rowSums(zotu_ITS_input$data_loaded,na.rm = TRUE)>100,]
z <- zotu_ITS_input$taxonomy_loaded[rownames(zotu_ITS_input$taxonomy_loaded) %in% rownames(y),]
x <- DGEList(y)
x <- calcNormFactors(x)
x <- estimateCommonDisp(x)
zotu_ITS_TMM <- list(x$pseudo.counts,zotu_ITS_input$map_loaded, z)
names(zotu_ITS_TMM) <- names(zotu_ITS_input)

zotu_16S_TMM$map_loaded$Treatment <- as.factor(zotu_16S_TMM$map_loaded$Treatment)


zotu_ITS_TMM$map_loaded$Treatment <- as.factor(zotu_ITS_TMM$map_loaded$Treatment)


##Taxonomy

tax_NBGLMs <- function(data.set) {
    results <- vector()
    for(i in 1:2){           #splits dataset for each taxonomic group
        print(paste(i, "of 2"))
        input_tax_level = summarize_taxonomy(data.set, level = i, report_higher_tax = FALSE, relative = TRUE)
        test.data.raw <-  cbind(data.set$map_loaded,t(input_tax_level))
        start.col <- dim(data.set$map_loaded)[2] + 1
        
        for(j in start.col:ncol(test.data.raw)){         #Does the glm and lsmean for each taxonomic group
            #print(j)
            
            test.glm <- tryCatch({
                suppressWarnings(glm.nb(test.data.raw[,j] ~ Treatment, data = test.data.raw))
            }, error = function(e){
                return(rep(NA, 8))
            })
            if(is(test.glm) != "negbin") {
                results <- rbind(results, test.glm);
                rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
                next
            }
            
            test.aov <- tryCatch({
                Anova(test.glm, test.statistic = "F")
            }, error = function(e){
                return(rep(NA, 8))
            })#, warning = function(w){
            #return(rep(NA, 20))
            #})
            if(is(test.aov)[1] != "anova") {
                results <- rbind(results, test.aov);
                rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
                next
            }
            
            test.lsm <- emmeans(test.glm, specs = "Treatment")
            pwc<-pairs(test.lsm)
            a <- c(test.aov$'Sum Sq', test.aov$Df, test.aov$F, test.aov$`Pr(>F)`,summary(test.lsm)$emmean, summary(test.lsm)$SE,summary(pwc)$p.value)
            results <- rbind(results, a)
            rownames(results)[nrow(results)] <- colnames(test.data.raw)[j]
            
        }
        
    }
    colnames(results) <- c(apply(expand.grid(rownames(test.aov), colnames(test.aov)), 1, paste, collapse="."), 
                           apply(expand.grid(summary(test.lsm)$Treatment, colnames(summary(test.lsm))[2:3]), 1, paste,collapse="."),
                           summary(pwc)$contrast)
    return(results)
}


#16S
zotu_NBGLM_tax_16S <- tax_NBGLMs(zotu_16S_TMM)

#ITS
zotu_NBGLM_tax_ITS <- tax_NBGLMs(zotu_ITS_TMM)

###Individual OTUs
otu_NBGLMs <- function(data.set) {
    
    test.data.raw <-  cbind(data.set$map_loaded,t(data.set$data_loaded))
    start.col <- dim(data.set$map_loaded)[2] + 1
    
    results <- vector()
    
    for(j in start.col:ncol(test.data.raw)){         #Does the glm and lsmean for each taxonomic group
        
        if(j%%100 == 0) { print(paste(j, "out of", ncol(test.data.raw))) }
        
        if(sum(test.data.raw[,j], na.rm = TRUE) < 100) { next }
        
        #model <- formula(test.data.raw[,j] ~ RepNumber + Accession)
        
        test.glm <- tryCatch({
            suppressWarnings(glm.nb(test.data.raw[,j] ~ Treatment, data = test.data.raw))
        }, error = function(e){
            return(rep(NA, 8))
        })
        if(is(test.glm) != "negbin") {
            #results <- rbind(results, test.glm);
            #rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
            next
        }
        
        test.aov <- tryCatch({
            Anova(test.glm, test.statistic = "F")
        }, error = function(e){
            return(rep(NA, 8))
        })#, warning = function(w){
        #return(rep(NA, 20))
        #}
        if(is(test.aov)[1] != "anova") {
            #results <- rbind(results, test.aov);
            #rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
            next
        }
        
        test.lsm <- LSmeans(test.glm, effect = "Treatment")
        a <- c(test.aov$'Sum Sq', test.aov$Df, test.aov$F, test.aov$`Pr(>F)`,test.lsm$coef[,1], test.lsm$coef[,2])
        results <- rbind(results, a)
        rownames(results)[nrow(results)] <- colnames(test.data.raw)[j]
        
    }
    
    
    colnames(results) <- c(apply(expand.grid(rownames(test.aov), colnames(test.aov)), 1, paste, collapse="."), 
                           apply(expand.grid(test.lsm$grid$Treatment, colnames(test.lsm$coef)[1:2]), 1, paste,
                                 collapse="."))
    return(results)
}

#16S
zotu_NBGLM_raw_16S <- otu_NBGLMs(zotu_16S_TMM)

#ITS
zotu_NBGLM_raw_ITS <- otu_NBGLMs(zotu_ITS_TMM)


#Transform LSMean estimates back to counts instead of ln(counts)
#16S
zotu_NBGLM_raw_16S[,grep(".estimate|.se", colnames(zotu_NBGLM_raw_16S))] <- exp(zotu_NBGLM_raw_16S[,grep(".estimate|.se", colnames(zotu_NBGLM_raw_16S))])
zotu_NBGLM_tax_16S[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_16S))] <- exp(zotu_NBGLM_tax_16S[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_16S))])

#ITS
zotu_NBGLM_raw_ITS[,grep(".estimate|.se", colnames(zotu_NBGLM_raw_ITS))] <- exp(zotu_NBGLM_raw_ITS[,grep(".estimate|.se", colnames(zotu_NBGLM_raw_ITS))])
zotu_NBGLM_tax_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_ITS))] <- exp(zotu_NBGLM_tax_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_ITS))])

# FDR Adjustments
#16S
p.cols <- grep(".Pr[:(:][:>:]F[:):]", colnames(zotu_NBGLM_raw_16S))
w <- apply(zotu_NBGLM_raw_16S[,p.cols], 2, function(x){ p.adjust(x, method = 'fdr') })
colnames(w) <- paste(colnames(w), "FDR", sep = ".")
zotu_NBGLM_raw_16S <- cbind(zotu_NBGLM_raw_16S, w)
p.cols <- grep(".Pr[:(:][:>:]F[:):]", colnames(zotu_NBGLM_tax_16S))
w <- apply(zotu_NBGLM_tax_16S[,p.cols], 2, function(x){ p.adjust(x, method = 'fdr') })
colnames(w) <- paste(colnames(w), "FDR", sep = ".")
zotu_NBGLM_tax_16S <- cbind(zotu_NBGLM_tax_16S, w)

#ITS
p.cols <- grep(".Pr[:(:][:>:]F[:):]", colnames(zotu_NBGLM_raw_ITS))
w <- apply(zotu_NBGLM_raw_ITS[,p.cols], 2, function(x){ p.adjust(x, method = 'fdr') })
colnames(w) <- paste(colnames(w), "FDR", sep = ".")
zotu_NBGLM_raw_ITS <- cbind(zotu_NBGLM_raw_ITS, w)
p.cols <- grep(".Pr[:(:][:>:]F[:):]", colnames(zotu_NBGLM_tax_ITS))
w <- apply(zotu_NBGLM_tax_ITS[,p.cols], 2, function(x){ p.adjust(x, method = 'fdr') })
colnames(w) <- paste(colnames(w), "FDR", sep = ".")
zotu_NBGLM_tax_ITS <- cbind(zotu_NBGLM_tax_ITS, w)

rm(p.cols, w)




# Eta-Squared Calcs
####May want to combine this process with the above functions.
etaSq.calc <- function(nbglm.data){
    results <- vector()
    SS.cols <- grep("Sum Sq$", colnames(nbglm.data))
    for(i in 1:nrow(nbglm.data)){
        eta.sq <- nbglm.data[i,SS.cols] / sum(nbglm.data[i,SS.cols], na.rm = T)
        results <- rbind(results, eta.sq)
    }
    colnames(results) <- gsub("Sum Sq", "etaSq", colnames(results))
    rownames(results) <- rownames(nbglm.data)
    return(results)
}

#16S
zotu_NBGLM_raw_16S <- cbind(zotu_NBGLM_raw_16S,  etaSq.calc(zotu_NBGLM_raw_16S))
zotu_NBGLM_tax_16S <- cbind(zotu_NBGLM_tax_16S,  etaSq.calc(zotu_NBGLM_tax_16S))

#ITS
zotu_NBGLM_raw_ITS <- cbind(zotu_NBGLM_raw_ITS,  etaSq.calc(zotu_NBGLM_raw_ITS))
zotu_NBGLM_tax_ITS <- cbind(zotu_NBGLM_tax_ITS,  etaSq.calc(zotu_NBGLM_tax_ITS))



#Write Results to CSV
write.csv(zotu_NBGLM_raw_16S, "zotu_NBGLM_raw_16S_2019.csv")
write.csv(zotu_NBGLM_tax_16S, "zotu_NBGLM_phylum_16S_2018.csv")

write.csv(zotu_NBGLM_raw_ITS, "zotu_NBGLM_raw_ITS_2019.csv")
write.csv(zotu_NBGLM_tax_ITS, "zotu_NBGLM_phylum_ITS_2018.csv")



###############################################
######   Plotting
###############################################




####Stacked barcharts for taxa
#16S 2018


LSmeans.cols <- grep(".estimate", colnames(zotu_NBGLM_tax_16S))
tax.rows <- grep("f_", rownames(zotu_NBGLM_tax_16S))
temp.df <- zotu_NBGLM_tax_16S[tax.rows, LSmeans.cols]

temp.df.m <- melt(temp.df, id.var="Phylum")

bacteria2018phylum<-read.csv(file.choose(),header=TRUE)
bacteria2018phylum<-bacteria2018phylum[-c(5:10)]

bacteria2018phylum.m <- melt(bacteria2018phylum, id.var="Phylum")

a<-ggplot(bacteria2018phylum.m, aes(x = variable, y = value, fill = Phylum)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = colorRampPalette(solarized_pal()(8))(11)) + 

    
    theme(
        axis.text = element_text(color='black',size=10),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_text(size=18),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

bacteria2019phylum<-read.csv(file.choose(),header=TRUE)
bacteria2019phylum<-bacteria2019phylum[-c(5:10)]
bacteria2019phylum.m <- melt(bacteria2019phylum, id.var="Phylum")

b<-ggplot(bacteria2019phylum.m, aes(x = variable, y = value, fill = Phylum)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = colorRampPalette(solarized_pal()(8))(10)) + 
    
    
    theme(
        axis.text = element_text(color='black',size=10),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_text(size=18),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

fungi2018phylum<-read.csv(file.choose(),header=TRUE)
fungi2018phylum<-fungi2018phylum[-c(5:10)]
fungi2018phylum.m <- melt(fungi2018phylum, id.var="Phylum")


c<-ggplot(fungi2018phylum.m, aes(x = variable, y = value, fill = Phylum)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = colorRampPalette(solarized_pal()(8))(10)) + 
    
    
    theme(
        axis.text = element_text(color='black',size=10),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_text(size=18),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

fungi2019phylum<-read.csv(file.choose(),header=TRUE)
fungi2019phylum<-fungi2019phylum[-c(5:10)]
fungi2019phylum.m <- melt(fungi2018phylum, id.var="Phylum")

d<-ggplot(fungi2019phylum.m, aes(x = variable, y = value, fill = Phylum)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = colorRampPalette(solarized_pal()(8))(10)) + 
    
    
    theme(
        axis.text = element_text(color='black',size=10),
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_text(size=18),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

library(grid)
library(gridExtra)
library(patchwork)

bargraphs<-(a+b)/(c+d)

ggsave(filename = "phylumgraphs.pdf",
       plot = bargraphs,
       bg = "transparent",
       width = 12, height = 8, units = "in",
       dpi = 600)


##ITS 2018

LSmeans.cols <- grep(".estimate", colnames(zotu_NBGLM_tax_ITS))
tax.rows <- grep("f_", rownames(zotu_NBGLM_tax_ITS))
temp.df <- zotu_NBGLM_tax_ITS[tax.rows, LSmeans.cols]
write.csv(temp.df, "ITS2019family.csv")

fungi2018phylum<-read.csv(file.choose(),header=TRUE)
fungi2018phylum.m <- melt(fungi2018phylum, id.var="Family")

fungi2018<-ggplot(fungi2018phylum.m, aes(x = variable, y = value, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_y_continuous(labels = percent_format()) +
  
    theme_fivethirtyeight() + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(filename = "2019fungifamily.pdf",
       plot = fungi2018,
       bg = "transparent",
       width = 36, height = 24, units = "in",
       dpi = 600)

#16S 2019


LSmeans.cols <- grep(".estimate", colnames(zotu_NBGLM_tax_16S))
tax.rows <- grep("c_", rownames(zotu_NBGLM_tax_16S))
temp.df <- zotu_NBGLM_tax_16S[tax.rows, LSmeans.cols]
write.csv(temp.df, "16S2019phylum.csv")

bacteria2019phylum<-read.csv(file.choose(),header=TRUE)
bacteria2019phylum.m <- melt(bacteria2019phylum, id.var="Phylum")



bacteria2019<-ggplot(bacteria2019phylum.m, aes(x = variable, y = value, fill = Phylum)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = colorRampPalette(solarized_pal()(8))(29)) + 
    theme_fivethirtyeight() + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(filename = "2019bacteriaphylum.pdf",
       plot = bacteria2019,
       bg = "transparent",
       width = 14, height = 18, units = "in",
       dpi = 600)

##ITS 2019

LSmeans.cols <- grep(".estimate", colnames(zotu_NBGLM_tax_ITS))
tax.rows <- grep("p_", rownames(zotu_NBGLM_tax_ITS))
temp.df <- zotu_NBGLM_tax_ITS[tax.rows, LSmeans.cols]
write.csv(temp.df, "ITS2019phylum.csv")

fungi2019phylum<-read.csv(file.choose(),header=TRUE)
fungi2019phylum.m <- melt(fungi2019phylum, id.var="Phylum")

fungi2019<-ggplot(fungi2019phylum.m, aes(x = variable, y = value, fill = Phylum)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = colorRampPalette(solarized_pal()(8))(25)) + 
    theme_fivethirtyeight() + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(filename = "2019fungiphylum.pdf",
       plot = fungi2019,
       bg = "transparent",
       width = 14, height = 18, units = "in",
       dpi = 600)


##barcharts of significance

phylum2018<-read.csv(file.choose(),header=TRUE)
library(ggplot2)
library(dplyr)
library(tidyverse)

fungi2019phylumrozello<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota"&Phylum!="d__Fungi; p__Basidiomycota"&Phylum!="d__Fungi; p__Basidiobolomycota"&Phylum!="d__Fungi; p__Calcarisporiellomycota")
fungi2019ascomycota<-subset(bacteriasig,Phylum!="d__Fungi; p__Rozellomycota"&Phylum!="d__Fungi; p__Basidiomycota"&Phylum!="d__Fungi; p__Basidiobolomycota"&Phylum!="d__Fungi; p__Calcarisporiellomycota")
fungi2019basidiomycota<-subset(bacteriasig,Phylum!="d__Fungi; p__Rozellomycota"&Phylum!="d__Fungi; p__Ascomycota"&Phylum!="d__Fungi; p__Basidiobolomycota"&Phylum!="d__Fungi; p__Calcarisporiellomycota")
fungi2019basidioblymcoaa<-subset(bacteriasig,Phylum!="d__Fungi; p__Rozellomycota"&Phylum!="d__Fungi; p__Ascomycota"&Phylum!="d__Fungi; p__Basidiomycota"&Phylum!="d__Fungi; p__Calcarisporiellomycota")
fungi2019calcar<-subset(bacteriasig,Phylum!="d__Fungi; p__Rozellomycota"&Phylum!="d__Fungi; p__Ascomycota"&Phylum!="d__Fungi; p__Basidiomycota"&Phylum!="d__Fungi; p__Basidiobolomycota")

fungi2019phylumrozello<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota"&Phylum!="d__Fungi; p__Basidiomycota"&Phylum!="d__Fungi; p__Basidiobolomycota"&Phylum!="d__Fungi; p__Calcarisporiellomycota")
fungi2019ascomycota<-subset(bacteriasig,Phylum!="d__Fungi; p__Rozellomycota"&Phylum!="d__Fungi; p__Basidiomycota"&Phylum!="d__Fungi; p__Basidiobolomycota"&Phylum!="d__Fungi; p__Calcarisporiellomycota")
fungi2019basidiomycota<-subset(bacteriasig,Phylum!="d__Fungi; p__Rozellomycota"&Phylum!="d__Fungi; p__Ascomycota"&Phylum!="d__Fungi; p__Basidiobolomycota"&Phylum!="d__Fungi; p__Calcarisporiellomycota")
fungi2019basidioblymcoaa<-subset(bacteriasig,Phylum!="d__Fungi; p__Rozellomycota"&Phylum!="d__Fungi; p__Ascomycota"&Phylum!="d__Fungi; p__Basidiomycota"&Phylum!="d__Fungi; p__Calcarisporiellomycota")
fungi2019calcar<-subset(bacteriasig,Phylum!="d__Fungi; p__Rozellomycota"&Phylum!="d__Fungi; p__Ascomycota"&Phylum!="d__Fungi; p__Basidiomycota"&Phylum!="d__Fungi; p__Basidiobolomycota")

fungiphyluma<-phylum2018[1:3,]
fungiphylumb<-phylum2018[4:6,]
fungiphylumc<-phylum2018[7:9,]
fungiphylumd<-phylum2018[10:12,]
fungiphylume<-phylum2018[13:15,]
fungiphylumf<-phylum2018[16:18,]

a<-ggplot(fungiphyluma,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
  
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

b<-ggplot(fungiphylumb,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
       
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


c<-ggplot(fungiphylumc,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

d<-ggplot(fungiphylumd,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        legend.position = "none",
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
       
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

e<-ggplot(fungiphylume,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

f<-ggplot(fungiphylumf,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

library(grid)
library(gridExtra)
library(patchwork)

phylum2018<-a+b
ggsave(filename = "family2019bacteriarelative.pdf",
       plot = a,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

a<-ggplot(phylum2018,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

bacteriasig<-read.csv(file.choose(),header=TRUE)

fungi2018a<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota; c__Xylonomycetes"&Phylum!="d__Fungi; p__Basidiobolomycota; c__Basidiobolomycetes"&Phylum!="d__Fungi; p__Basidiomycota; c__Agaricomycetes"&Phylum!="d__Fungi; p__Calcarisporiellomycota; c__Calcarisporiellomycetes")
fungi2018b<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota; c__Taphrinomycetes"&Phylum!="d__Fungi; p__Basidiobolomycota; c__Basidiobolomycetes"&Phylum!="d__Fungi; p__Basidiomycota; c__Agaricomycetes"&Phylum!="d__Fungi; p__Calcarisporiellomycota; c__Calcarisporiellomycetes")
fungi2018c<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota; c__Taphrinomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Xylonomycetes"&Phylum!="d__Fungi; p__Basidiomycota; c__Agaricomycetes"&Phylum!="d__Fungi; p__Calcarisporiellomycota; c__Calcarisporiellomycetes")
fungi2018d<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota; c__Taphrinomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Xylonomycetes"&Phylum!="d__Fungi; p__Basidiobolomycota; c__Basidiobolomycetes"&Phylum!="d__Fungi; p__Calcarisporiellomycota; c__Calcarisporiellomycetes")
fungi2018e<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota; c__Taphrinomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Xylonomycetes"&Phylum!="d__Fungi; p__Basidiobolomycota; c__Basidiobolomycetes"&Phylum!="d__Fungi; p__Basidiomycota; c__Agaricomycetes")

fungi2019a<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota; c__Orbiliomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Sordariomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Xylonomycetes"&Phylum!="d__Fungi; p__Rozellomycota; c__Rozellomycotina_cls_Incertae_sedis")
fungi2019b<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota; c__Leotiomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Sordariomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Xylonomycetes"&Phylum!="d__Fungi; p__Rozellomycota; c__Rozellomycotina_cls_Incertae_sedis")
fungi2019c<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota; c__Leotiomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Orbiliomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Xylonomycetes"&Phylum!="d__Fungi; p__Rozellomycota; c__Rozellomycotina_cls_Incertae_sedis")
fungi2019d<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota; c__Leotiomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Orbiliomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Sordariomycetes"&Phylum!="d__Fungi; p__Rozellomycota; c__Rozellomycotina_cls_Incertae_sedis")
fungi2019e<-subset(bacteriasig,Phylum!="d__Fungi; p__Ascomycota; c__Leotiomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Orbiliomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Sordariomycetes"&Phylum!="d__Fungi; p__Ascomycota; c__Xylonomycetes")


a<-ggplot(fungi2019a,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

b<-ggplot(fungi2019b,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


c<-ggplot(fungi2019c,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

d<-ggplot(fungi2019d,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

e<-ggplot(fungi2019e,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

class2019<-(a+b+c)/(d+e)

ggsave(filename = "class2019.pdf",
       plot = class2019,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##order. dropping rows should be easeir like this

class2018<-read.csv(file.choose(),header=TRUE)

fungiclassa<-class2018[1:3,]
fungiclassb<-class2018[4:6,]
fungiclassc<-class2018[7:9,]
fungiclassd<-class2018[10:12,]
fungiclasse<-class2018[13:15,]
fungiclassf<-class2018[16:18,]
fungiclassg<-class2018[19:21,]
fungiclassh<-class2018[22:24,]
fungiclassi<-class2018[25:27,]
fungiclassj<-class2018[28:30,]
fungiclassk<-class2018[31:33,]
fungiclassl<-class2018[34:36,]
fungiclassm<-class2018[37:39,]
fungiclassn<-class2018[40:42,]
fungiclasso<-class2018[43:45,]
fungiclassp<-class2018[46:48,]

a<-ggplot(fungiclassa,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

b<-ggplot(fungiclassb,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


c<-ggplot(fungiclassc,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

d<-ggplot(fungiclassd,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

e<-ggplot(fungiclasse,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

f<-ggplot(fungiclassf,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


g<-ggplot(fungiclassg,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

h<-ggplot(fungiclassh,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

i<-ggplot(fungiclassi,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

j<-ggplot(fungiclassj,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


k<-ggplot(fungiclassk,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

l<-ggplot(fungiclassl,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

m<-ggplot(fungiclassm,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

n<-ggplot(fungiclassn,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


o<-ggplot(fungiclasso,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

p<-ggplot(fungiclassp,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

order2018<-(a+b+c+d)/(e+f+g+h)/(i+j+k+l)/(m+n+o+p)+plot_layout(widths=c(4,4))

ggsave(filename = "order2018.pdf",
       plot = order2018,
       bg = "transparent",
       width = 8, height = 20, units = "in",
       dpi = 600)

class2018<-read.csv(file.choose(),header=TRUE)

fungiclassa<-class2018[1:3,]
fungiclassb<-class2018[4:6,]
fungiclassc<-class2018[7:9,]
fungiclassd<-class2018[10:12,]
fungiclasse<-class2018[13:15,]
fungiclassf<-class2018[16:18,]
fungiclassg<-class2018[19:21,]
fungiclassh<-class2018[22:24,]
fungiclassi<-class2018[25:27,]
fungiclassj<-class2018[28:30,]


a<-ggplot(fungiclassa,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

b<-ggplot(fungiclassb,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


c<-ggplot(fungiclassc,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        legend.position = "none",
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

d<-ggplot(fungiclassd,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

e<-ggplot(fungiclasse,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
      
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

f<-ggplot(fungiclassf,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
       
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


g<-ggplot(fungiclassg,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        legend.position="none",
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

h<-ggplot(fungiclassh,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

i<-ggplot(fungiclassi,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

j<-ggplot(fungiclassj,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))




order2019<-(a+b+c+d+e)/(f+g+h+i+j)

ggsave(filename = "order2019.pdf",
       plot = order2019,
       bg = "transparent",
       width = 8, height = 10, units = "in",
       dpi = 600)


phylum2018<-read.csv(file.choose(),header=TRUE)
library(ggplot2)
library(dplyr)
library(tidyverse)

bacteriaclassa<-phylum2018[1:3,]
bacteriaclassb<-phylum2018[4:6,]
bacteriaclassc<-phylum2018[7:9,]



a<-ggplot(bacteriaclassa,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

b<-ggplot(bacteriaclassb,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

c<-ggplot(bacteriaclassc,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))



library(grid)
library(gridExtra)
library(patchwork)

phylumbac2019<-a+b+c

ggsave(filename = "order2019ITSrelativeabundance.pdf",
       plot = phylumbac2019,
       bg = "transparent",
       width = 14, height = 8, units = "in",
       dpi = 600)

phylum2018<-read.csv(file.choose(),header=TRUE)

fungiordera<-phylum2018[1:3,]
fungiorderb<-phylum2018[4:6,]
fungiorderc<-phylum2018[7:9,]
fungiorderd<-phylum2018[10:12,]
fungiordere<-phylum2018[13:15,]
fungiorderf<-phylum2018[16:18,]



a<-ggplot(fungiordera,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

b<-ggplot(fungiorderb,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


c<-ggplot(fungiorderc,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        legend.position = "none",
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

d<-ggplot(fungiorderd,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

e<-ggplot(fungiordere,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

f<-ggplot(fungiorderf,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))




order2018<-(a+b+c)/(d+e+f)

ggsave(filename = "order201916Srelativeabundance.pdf",
       plot = order2018,
       bg = "transparent",
       width = 14, height = 8, units = "in",
       dpi = 600)

phylum2018<-read.csv(file.choose(),header=TRUE)

fungifamilya<-phylum2018[1:3,]
fungifamilyb<-phylum2018[4:6,]
fungifamilyc<-phylum2018[7:9,]
fungifamilyd<-phylum2018[10:12,]
fungifamilye<-phylum2018[13:15,]
fungifamilyf<-phylum2018[16:18,]
fungifamilyg<-phylum2018[19:21,]
fungifamilyh<-phylum2018[22:24,]
fungifamilyi<-phylum2018[25:27,]
fungifamilyj<-phylum2018[28:30,]


a<-ggplot(fungifamilya,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

b<-ggplot(fungifamilyb,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


c<-ggplot(fungifamilyc,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
     
        axis.title.y= element_blank(),
        legend.position = "none",
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

d<-ggplot(fungifamilyd,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

e<-ggplot(fungifamilye,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.y= element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

f<-ggplot(fungifamilyf,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=14),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


g<-ggplot(fungifamilyg,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        legend.position="none",
        
        axis.title.x= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

h<-ggplot(fungifamilyh,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position="none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

i<-ggplot(fungifamilyi,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        legend.position="none",
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

j<-ggplot(fungifamilyj,aes(Phylum,Estimate,fill=Treatment)) +
    geom_bar(stat="identity",width=0.5,color="black",aes(color="Treatment"),position=position_dodge())+
    scale_fill_manual("",breaks=c("CHR","CON","INT"),values=c("dark gray","blue","orange")) +
    
    geom_errorbar(aes(ymin=Estimate-std, ymax=Estimate+std), width=.2,
                  position=position_dodge(.5))+
    ylab("Abundance") +
    
    theme(
        axis.title.x= element_blank(),
        
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))


family2018<-(a+b+c)/(d+e+f)/(g+h+i)/(j)

ggsave(filename = "family201916Srelativeabundance.pdf",
       plot = family2018,
       bg = "transparent",
       width = 14, height = 8, units = "in",
       dpi = 600)



library(ggplot2)

Opporunistic<-read.csv(file.choose(),header=TRUE)

b<-ggplot(Opporunistic,aes(Treatment,fill=factor(Type,levels=c("Sensitive","Opportunistic","Resistant")))) +
    geom_bar(width=0.5,color="gray")+
    
 
    ylab("Count") +
    
    theme(
        
        axis.text = element_text(color='black',size=14),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "BarChartOpportunistic16S2019.pdf",
       plot = b,
       bg = "transparent",
       width = 8, height = 12, units = "in",
       dpi = 600)

library(dplyr)
library(tidyverse)

Ascomycotagraph<-subset(Opporunistic,Phylum!="Fungi"&Phylum!="Basidiobolomycota"&Phylum!="Basidiomycota"&Phylum!="Calcarisporiellomycota"
                        &Phylum!="Entorrhizomycota"&Phylum!="Glomeromycota"&Phylum!="Mortierellomycota"&Phylum!="Mucoromycota"&Phylum!="Rozellomycota")

AscomycotagraphCHR<-subset(Ascomycotagraph,Treatment!="INT")
AscomycotagraphINT<-subset(Ascomycotagraph,Treatment!="CHR")


a<-table(AscomycotagraphCHR$Type)
b<-table(AscomycotagraphINT$Type)

Ascomycotapiegraph<-read.csv(file.choose(),header=TRUE)

d<-ggplot(Ascomycotapiegraph, aes(x="", y=Amount, group=Type, color=Type, fill=Type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + facet_wrap(~ Treatment) +
    theme(
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        
        plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "Ascomycota2019.pdf",
       plot = d,
       bg = "transparent",
       width = 12, height = 6, units = "in",
       dpi = 600)

Basidiobolomycotagraph<-subset(Opporunistic,Phylum!="Fungi"&Phylum!="Ascomycota"&Phylum!="Basidiomycota"&Phylum!="Calcarisporiellomycota"
                        &Phylum!="Entorrhizomycota"&Phylum!="Glomeromycota"&Phylum!="Mortierellomycota"&Phylum!="Mucoromycota"&Phylum!="Rozellomycota")

BasidiobolomycotaCHR<-subset(Basidiobolomycotagraph,Treatment!="INT")
BasidiobolomycotaINT<-subset(Basidiobolomycotagraph,Treatment!="CHR")

a<-table(BasidiobolomycotaCHR$Type)
b<-table(BasidiobolomycotaINT$Type)

Basidiobolomycotagraph<-read.csv(file.choose(),header=TRUE)

e<-ggplot(Basidiobolomycotagraph, aes(x="", y=Amount, group=Type, color=Type, fill=Type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + facet_wrap(~ Treatment) +
    theme(
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        
        plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "Basidiobolomycota.pdf",
       plot = e,
       bg = "transparent",
       width = 12, height = 6, units = "in",
       dpi = 600)


Basidiomycotagraph<-subset(Opporunistic,Phylum!="Fungi"&Phylum!="Ascomycota"&Phylum!="Basidiobolomycota"&Phylum!="Calcarisporiellomycota"
                               &Phylum!="Entorrhizomycota"&Phylum!="Glomeromycota"&Phylum!="Mortierellomycota"&Phylum!="Mucoromycota"&Phylum!="Rozellomycota")

BasidiomycotaCHR<-subset(Basidiomycotagraph,Treatment!="INT")
BasidiomycotaINT<-subset(Basidiomycotagraph,Treatment!="CHR")

a<-table(BasidiomycotaCHR$Type)
b<-table(BasidiomycotaINT$Type)

Basidiomycotagraph<-read.csv(file.choose(),header=TRUE)

f<-ggplot(Basidiomycotagraph, aes(x="", y=Amount, group=Type, color=Type, fill=Type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + facet_wrap(~ Treatment) +
    theme(
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        
        plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "Basidiomycota.pdf",
       plot = f,
       bg = "transparent",
       width = 12, height = 6, units = "in",
       dpi = 600)

Calcarisporiellomycota<-subset(Opporunistic,Phylum!="Fungi"&Phylum!="Ascomycota"&Phylum!="Basidiobolomycota"&Phylum!="Basidiomycota"
                           &Phylum!="Entorrhizomycota"&Phylum!="Glomeromycota"&Phylum!="Mortierellomycota"&Phylum!="Mucoromycota"&Phylum!="Rozellomycota")

CalcarisporiellomycotaCHR<-subset(Calcarisporiellomycota,Treatment!="INT")
CalcarisporiellomycotaINT<-subset(Calcarisporiellomycota,Treatment!="CHR")

a<-table(CalcarisporiellomycotaCHR$Type)
b<-table(CalcarisporiellomycotaINT$Type)

Calcarisporiellomycotagraph<-read.csv(file.choose(),header=TRUE)

g<-ggplot(Calcarisporiellomycotagraph, aes(x="", y=Amount, group=Type, color=Type, fill=Type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + facet_wrap(~ Treatment) +
    theme(
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        
        plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "Calcarisporiellomycota.pdf",
       plot = g,
       bg = "transparent",
       width = 12, height = 6, units = "in",
       dpi = 600)

Entorrhizomycota<-subset(Opporunistic,Phylum!="Fungi"&Phylum!="Ascomycota"&Phylum!="Basidiobolomycota"&Phylum!="Basidiomycota"
                               &Phylum!="Calcarisporiellomycota"&Phylum!="Glomeromycota"&Phylum!="Mortierellomycota"&Phylum!="Mucoromycota"&Phylum!="Rozellomycota")

EntorrhizomycotaCHR<-subset(Entorrhizomycota,Treatment!="INT")
EntorrhizomycotaINT<-subset(Entorrhizomycota,Treatment!="CHR")

a<-table(EntorrhizomycotaCHR$Type)
b<-table(EntorrhizomycotaINT$Type)

Entorrhizomycotagraph<-read.csv(file.choose(),header=TRUE)

h<-ggplot(Entorrhizomycotagraph, aes(x="", y=Amount, group=Type, color=Type, fill=Type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + facet_wrap(~ Treatment) +
    theme(
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        
        plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "Entorrhizomycota.pdf",
       plot = h,
       bg = "transparent",
       width = 12, height = 6, units = "in",
       dpi = 600)

Glomeromycota<-subset(Opporunistic,Phylum!="Fungi"&Phylum!="Ascomycota"&Phylum!="Basidiobolomycota"&Phylum!="Basidiomycota"
                         &Phylum!="Calcarisporiellomycota"&Phylum!="Entorrhizomycota"&Phylum!="Mortierellomycota"&Phylum!="Mucoromycota"&Phylum!="Rozellomycota")

GlomeromycotaCHR<-subset(Glomeromycota,Treatment!="INT")
GlomeromycotaINT<-subset(Glomeromycota,Treatment!="CHR")

a<-table(GlomeromycotaCHR$Type)
b<-table(GlomeromycotaINT$Type)

Glomeromycotagraph<-read.csv(file.choose(),header=TRUE)

i<-ggplot(Glomeromycotagraph, aes(x="", y=Amount, group=Type, color=Type, fill=Type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + facet_wrap(~ Treatment) +
    theme(
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        
        plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "Glomeromycota.pdf",
       plot = i,
       bg = "transparent",
       width = 12, height = 6, units = "in",
       dpi = 600)

Mortierellomycota<-subset(Opporunistic,Phylum!="Fungi"&Phylum!="Ascomycota"&Phylum!="Basidiobolomycota"&Phylum!="Basidiomycota"
                      &Phylum!="Calcarisporiellomycota"&Phylum!="Entorrhizomycota"&Phylum!="Glomeromycota"&Phylum!="Mucoromycota"&Phylum!="Rozellomycota")

MortierellomycotaCHR<-subset(Mortierellomycota,Treatment!="INT")
MortierellomycotaINT<-subset(Mortierellomycota,Treatment!="CHR")

a<-table(MortierellomycotaCHR$Type)
b<-table(MortierellomycotaINT$Type)

Mortierellomycotagraph<-read.csv(file.choose(),header=TRUE)

j<-ggplot(Mortierellomycotagraph, aes(x="", y=Amount, group=Type, color=Type, fill=Type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + facet_wrap(~ Treatment) +
    theme(
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        
        plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "Mortierellomycota.pdf",
       plot = j,
       bg = "transparent",
       width = 12, height = 6, units = "in",
       dpi = 600)

Mucoromycota<-subset(Opporunistic,Phylum!="Fungi"&Phylum!="Ascomycota"&Phylum!="Basidiobolomycota"&Phylum!="Basidiomycota"
                          &Phylum!="Calcarisporiellomycota"&Phylum!="Entorrhizomycota"&Phylum!="Glomeromycota"&Phylum!="Mortierellomycota"&Phylum!="Rozellomycota")

MucoromycotaCHR<-subset(Mucoromycota,Treatment!="INT")
MucoromycotaINT<-subset(Mucoromycota,Treatment!="CHR")

a<-table(MucoromycotaCHR$Type)
b<-table(MucoromycotaINT$Type)

Mucoromycotagraph<-read.csv(file.choose(),header=TRUE)

k<-ggplot(Mucoromycotagraph, aes(x="", y=Amount, group=Type, color=Type, fill=Type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + facet_wrap(~ Treatment) +
    theme(
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        
        plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "Mucoromycota.pdf",
       plot = k,
       bg = "transparent",
       width = 12, height = 6, units = "in",
       dpi = 600)


Rozellomycota<-subset(Opporunistic,Phylum!="Fungi"&Phylum!="Ascomycota"&Phylum!="Basidiobolomycota"&Phylum!="Basidiomycota"
                     &Phylum!="Calcarisporiellomycota"&Phylum!="Entorrhizomycota"&Phylum!="Glomeromycota"&Phylum!="Mortierellomycota"&Phylum!="Mucoromycota")

RozellomycotaCHR<-subset(Rozellomycota,Treatment!="INT")
RozellomycotaINT<-subset(Rozellomycota,Treatment!="CHR")

a<-table(RozellomycotaCHR$Type)
b<-table(RozellomycotaINT$Type)

Rozellomycotagraph<-read.csv(file.choose(),header=TRUE)

k<-ggplot(Rozellomycotagraph, aes(x="", y=Amount, group=Type, color=Type, fill=Type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + facet_wrap(~ Treatment) +
    theme(
        
        axis.text = element_text(color='black',size=8),
        axis.title = element_text(color='black',size=12),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank(),
        panel.background = element_rect(fill=NA,color='black'),
        panel.grid = element_blank(),
        
        plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "Rozellomycota.pdf",
       plot = k,
       bg = "transparent",
       width = 12, height = 6, units = "in",
       dpi = 600)

##Number per phylum significant

numsig<-read.csv(file.choose(),header=TRUE)
library(ggplot2)

a<-ggplot(numsig, aes(Phylum,Numbersig, fill=Phylum)) + 
    ylab("Total Significant ASV per phyla (count)") +
    geom_bar(stat="identity", color="black")+
    scale_x_discrete(limits=c("Basidiomycota","Ascomycota","Rozellomycota","Mortierellomycota","Basidiobolomycota","Calcarisporiellomycota","Entorrhizomycota","Glomeromycota","Mucoromycota")) +
    theme_classic()+
    theme(
        legend.position="none", 
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size= 8.5),
        axis.ticks = element_line(color='black'),
        panel.background = element_rect(fill=NA),panel.grid = element_blank())

b<-ggplot(numsig, aes(Phylum,ratio, fill=Phylum)) + 
    ylab("Significant ASV per phyla/Total ASV per phyla (percent)") +
    geom_bar(stat="identity", color="black")+
    scale_x_discrete(limits=c("Basidiomycota","Ascomycota","Rozellomycota","Mortierellomycota","Basidiobolomycota","Calcarisporiellomycota","Entorrhizomycota","Glomeromycota","Mucoromycota")) +
    theme_classic()+
    theme(
        legend.position="none", 
        axis.text = element_text(color='black',size=6),
        axis.title = element_text(color='black',size= 8.5),
        axis.ticks = element_line(color='black'),
        panel.background = element_rect(fill=NA),panel.grid = element_blank())

library(grid)
library(gridExtra)
library(patchwork)

plot1<-a/b
ggsave(filename = "numberofsignificant2019.pdf",
       plot = plot1,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)


