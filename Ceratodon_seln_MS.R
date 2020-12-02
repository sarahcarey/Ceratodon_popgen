
### R Script for Carey et al., 202X on selection on Ceratodon UV sex chromosomes

# using R version 3.5.3 

library(ggplot2)
#v3.3.1

################### dS plot ###################

dn_vs_ds <- read.csv("Data/ceratodon_ds_dn_july2019.csv",header=TRUE)


png("Figures/ds_Oct2020.png", width = 8, height = 8, units = 'in', res = 300)

stat_box_data <- function(x, upper_limit = 1.05) {
  return( 
    data.frame(
      y = 0.9 * upper_limit,
      label = paste(format(round(mean(x), 3), big.mark = ",", decimal.mark = ".", scientific = FALSE))
    )
  )
}

b <- ggplot(dn_vs_ds, (aes(x = Sex, y = dS, fill=Sex))) +
  geom_boxplot(color=c("gray0", "gray50"), fill="white", alpha=0, lwd=2)
b  + scale_fill_manual(values=c("peru", "midnightblue")) +
  xlab("") +
  ylab("dS") +
  theme(axis.title.x = element_text(size=30)) +
  theme(axis.title.y = element_text(size=30))+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.75, jitter.height = 0), 
             aes(group=Sex), alpha=0.75, shape=19, size=2, colour="gray25") +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.000001, cex=10) + 
  guides(fill=FALSE) + 
  theme(axis.text.x= element_text(size=30))

dev.off()

# test for signifcant differences
wilcox.test(dn_vs_ds$dS ~ dn_vs_ds$Sex)

################### dN plot ###################

dn_vs_ds <- read.csv("Data/ceratodon_ds_dn_july2019.csv",header=TRUE)

png("Figures/dn_Oct2020.png", width = 8, height = 8, units = 'in', res = 300)

stat_box_data <- function(x, upper_limit = 0.4) {
  return( 
    data.frame(
      y = 1 * upper_limit,
      label = paste(format(round(mean(x), 3), big.mark = ",", decimal.mark = ".", scientific = FALSE))
    )
  )
}

b <- ggplot(dn_vs_ds, (aes(x = Sex, y = dN, fill=Sex))) +
  geom_boxplot(color=c("gray0", "gray50"), fill="white", alpha=0, lwd=2)
b  + scale_fill_manual(values=c("peru", "midnightblue")) +
  xlab("") +
  ylab("dN") +
  theme(axis.title.x = element_text(size=30)) +
  theme(axis.title.y = element_text(size=30))+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.75, jitter.height = 0), 
             aes(group=Sex), alpha=0.75, shape=19, size=2, colour="gray25") +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.000001, cex=10) + 
  guides(fill=FALSE) + 
  theme(axis.text.x= element_text(size=30))

dev.off()

# test for signifcant differences
wilcox.test(dn_vs_ds$dS ~ dn_vs_ds$Sex)

################### Population genetics analyses ###################

################### Nucleotide diversity ###################

library(PopGenome)
#v2.7.5

## split VCF files into chromosomes
VCF_split_into_scaffolds("Ceratodon_bwa.autochloro.filtered.md.vcf", "Ceratodon_bwa_AC_md")
VCF_split_into_scaffolds("Ceratodon_bwa.U.filtered.md.vcf", "Ceratodon_bwa_U_md")
VCF_split_into_scaffolds("Ceratodon_bwa.V.filtered.md.vcf", "Ceratodon_bwa_V_md")
VCF_split_into_scaffolds("Ceratodon_ngm.autochloro.filtered.md.vcf", "Ceratodon_ngm_AC_md")
VCF_split_into_scaffolds("Ceratodon_ngm.U.filtered.md.vcf", "Ceratodon_ngm_U_md")
VCF_split_into_scaffolds("Ceratodon_ngm.V.filtered.md.vcf", "Ceratodon_ngm_V_md")

# have to move each chromosome into its own directory for PopGenome to run

## various populations to run depending on the analysis
male <- c("Alaska_21.1.1_M1",
          "Dur_19.2.8_M1",
          "Dur_4.13.5_M2",
          "Port_4_M2",
          "Uconn_12.2.4_F2",
          "Port_8_M1",
          "Dur_13.8.1_M3",
          "Equ_E13.E3.14_M1")

female <- c("Alaska_21.1.8_F1",
            "Dur_16.6.2_F1",
            "Dur_4.13.7_F2",
            "Port_11_F2",
            "Uconn_15.12.12_F1",
            "Dur_13.8.10_F3",
            "Equ_E13.E3.1_F1",
            "Port_7_F1")

eastF<- c("Dur_16.6.2_F1",
          "Dur_4.13.7_F2",
          "Uconn_15.12.12_F1",
          "Dur_13.8.10_F3")

westF<- c("Alaska_21.1.8_F1",
          "Port_11_F2",
          "Port_7_F1")

eastM <- c("Dur_19.2.8_M1",
           "Dur_4.13.5_M2",
           "Uconn_12.2.4_F2",
           "Dur_13.8.1_M3")

westM <- c("Alaska_21.1.1_M1",
           "Port_4_M2",
           "Port_8_M1")

east <- c("Dur_19.2.8_M1",
           "Dur_4.13.5_M2",
           "Uconn_12.2.4_F2",
           "Dur_13.8.1_M3",
          "Dur_16.6.2_F1",
          "Dur_4.13.7_F2",
          "Uconn_15.12.12_F1",
          "Dur_13.8.10_F3")

west <- c("Alaska_21.1.1_M1",
           "Port_4_M2",
           "Port_8_M1",
          "Alaska_21.1.8_F1",
          "Port_11_F2",
          "Port_7_F1")

keep <- c("Alaska_21.1.1_M1",
          "Dur_19.2.8_M1",
          "Dur_4.13.5_M2",
          "Port_4_M2",
          "Uconn_12.2.4_F2",
          "Port_8_M1",
          "Dur_13.8.1_M3",
          "Equ_E13.E3.14_M1",
          "Alaska_21.1.8_F1",
          "Dur_16.6.2_F1",
          "Dur_4.13.7_F2",
          "Port_11_F2",
          "Uconn_15.12.12_F1",
          "Dur_13.8.10_F3",
          "Equ_E13.E3.1_F1",
          "Port_7_F1")

outgroup <- c("Chile_1.1_M1",
              "Chile_2.12_F1")



populations <- list(male,female)
populations <- list(keep,outgroup)
populations <- list(eastF,westF)
populations <- list(eastM,westM)
populations <- list(east,west)
populations <- list(male,female,outgroup)


# call in chromosome to analyze
GENOME.class <- readData("Ceratodon_bwa_AC_md/1", format = "VCF", 
                         include.unknown = TRUE, FAST = FALSE)

# set populations
GENOME.class <- set.populations(GENOME.class, populations,
                                diploid = FALSE, triploid = FALSE, 
                                tetraploid = FALSE)


# Analyses by chromosome
# Populations to use:
# use 'keep' for autos
# use 'male' or 'female' for sex chroms
# use 'east' and 'west' for auto and chloroplast Fst
# use 'eastM' and 'westM' or 'eastF' and 'westF' for sex chroms

GENOME.class <- neutrality.stats(GENOME.class) 
GENOME.class <- F_ST.stats(GENOME.class) 
get.sum.data(GENOME.class)

# Wu and Watterson's theta
thetaW <- GENOME.class@theta_Watterson
thetaW <- thetaW/GENOME.class@n.sites
thetaW

#Tajima's theta
thetaT <- GENOME.class@theta_Tajima
thetaT <- thetaT/GENOME.class@n.sites
thetaT

# Pi
Pi <- GENOME.class@nuc.diversity.within
Pi <- Pi/GENOME.class@n.sites
Pi

#Tajima's D 
TajimaD <- GENOME.class@Tajima.D
TajimaD

# get segregating sites
segsites <- GENOME.class@n.segregating.sites
segsites

## Fst
FST <- GENOME.class@nuc.F_ST.pairwise
FST


################### Sliding window analyses ###################

################### Neutrality tests ###################

###################  Autosomes ################### 
for (i in 1:12) {
  
GENOME.class <- readData(paste("Ceratodon_bwa_AC_md/",i, sep=""), format = "VCF", include.unknown=TRUE, FAST=FALSE)

keep <- c("Alaska_21.1.1_M1",
          "Dur_19.2.8_M1",
          "Dur_4.13.5_M2",
          "Port_4_M2",
          "Uconn_12.2.4_F2",
          "Port_8_M1",
          "Dur_13.8.1_M3",
          "Equ_E13.E3.14_M1",
          "Alaska_21.1.8_F1",
          "Dur_16.6.2_F1",
          "Dur_4.13.7_F2",
          "Port_11_F2",
          "Uconn_15.12.12_F1",
          "Dur_13.8.10_F3",
          "Equ_E13.E3.1_F1",
          "Port_7_F1")

populations <- list(keep)  

slide.GENOME.class <- sliding.window.transform(GENOME.class,width=100000,jump=10000,type=2)

slide.GENOME.class <- set.populations(slide.GENOME.class, populations,
                                diploid = FALSE, triploid = FALSE, 
                                tetraploid = FALSE)

slide.GENOME.class <- neutrality.stats(slide.GENOME.class) 
slide.GENOME.class <- F_ST.stats(slide.GENOME.class) 

# Wu and Watterson's theta
thetaW <- slide.GENOME.class@theta_Watterson
thetaW <- thetaW/100000

# Tajima's D
TajimaD <- slide.GENOME.class@Tajima.D

# Pi
Pi <- slide.GENOME.class@nuc.diversity.within
Pi <- Pi/100000


regions <- as.data.frame(slide.GENOME.class@region.names)
regions_split <- data.frame(do.call('rbind', strsplit(as.character(regions[,1]),' ',fixed=TRUE)))
chr <- rep(GENOME.class@region.names,each=length(thetaW[,1]))
regions_chr <- cbind(regions_split,chr,thetaW,TajimaD,Pi)
regions_cleaned <- subset(regions_chr, select=-c(X2,X4))
write.table(regions_cleaned, "neutrality.autos.txt", append=TRUE, col.names=FALSE, sep=",")

}

################### U chromosome ################### 

GENOME.class <- readData("Ceratodon_bwa_U_md/female", format = "VCF", include.unknown=TRUE, FAST=FALSE)

female <- c("Alaska_21.1.8_F1",
            "Dur_16.6.2_F1",
            "Dur_4.13.7_F2",
            "Port_11_F2",
            "Uconn_15.12.12_F1",
            "Dur_13.8.10_F3",
            "Equ_E13.E3.1_F1",
            "Port_7_F1")

populations <- list(female)

slide.GENOME.class <- sliding.window.transform(GENOME.class,width=100000,jump=10000,type=2)

slide.GENOME.class <- set.populations(slide.GENOME.class, populations,
                                      diploid = FALSE, triploid = FALSE, 
                                      tetraploid = FALSE)

slide.GENOME.class <- neutrality.stats(slide.GENOME.class)
slide.GENOME.class <- F_ST.stats(slide.GENOME.class) 


# Wu and Watterson's theta
thetaW <- slide.GENOME.class@theta_Watterson
thetaW <- thetaW/100000


# Tajima's D
TajimaD <- slide.GENOME.class@Tajima.D

# Pi
Pi <- slide.GENOME.class@nuc.diversity.within
Pi <- Pi/100000


regions <- as.data.frame(slide.GENOME.class@region.names)
regions_split <- data.frame(do.call('rbind', strsplit(as.character(regions[,1]),' ',fixed=TRUE)))
chr <- rep(GENOME.class@region.names,each=length(thetaW[,1]))
regions_chr <- cbind(regions_split,chr,thetaW,TajimaD,Pi)
regions_cleaned <- subset(regions_chr, select=-c(X2,X4))
write.table(regions_cleaned, "neutrality.U.txt", col.names=FALSE, sep=",")


################### V chromosome ###################

GENOME.class <- readData("Ceratodon_bwa_V_md/male", format = "VCF", include.unknown=TRUE, FAST=FALSE)

male <- c("Alaska_21.1.1_M1",
          "Dur_19.2.8_M1",
          "Dur_4.13.5_M2",
          "Port_4_M2",
          "Uconn_12.2.4_F2",
          "Port_8_M1",
          "Dur_13.8.1_M3",
          "Equ_E13.E3.14_M1")

populations <- list(male)

slide.GENOME.class <- sliding.window.transform(GENOME.class,width=100000,jump=10000,type=2)

slide.GENOME.class <- set.populations(slide.GENOME.class, populations,
                                      diploid = FALSE, triploid = FALSE, 
                                      tetraploid = FALSE)

slide.GENOME.class <- neutrality.stats(slide.GENOME.class)
slide.GENOME.class <- F_ST.stats(slide.GENOME.class) 


# Wu and Watterson's theta
thetaW <- slide.GENOME.class@theta_Watterson
thetaW <- thetaW/100000

# Tajima's D
TajimaD <- slide.GENOME.class@Tajima.D

Pi <- slide.GENOME.class@nuc.diversity.within
Pi <- Pi/100000


regions <- as.data.frame(slide.GENOME.class@region.names)
regions_split <- data.frame(do.call('rbind', strsplit(as.character(regions[,1]),' ',fixed=TRUE)))
chr <- rep(GENOME.class@region.names,each=length(thetaW[,1]))
regions_chr <- cbind(regions_split,chr,thetaW,TajimaD,Pi)
regions_cleaned <- subset(regions_chr, select=-c(X2,X4))
write.table(regions_cleaned, "neutrality.V.txt", col.names=FALSE, sep=",")


###################  Fst autosomes ################### 
for (i in 1:12) {
  
  GENOME.class <- readData(paste("Ceratodon_bwa_AC_md/",i, sep=""), format = "VCF", include.unknown=TRUE, FAST=FALSE)
  
  east <- c("Dur_19.2.8_M1",
            "Dur_4.13.5_M2",
            "Uconn_12.2.4_F2",
            "Dur_13.8.1_M3",
            "Dur_16.6.2_F1",
            "Dur_4.13.7_F2",
            "Uconn_15.12.12_F1",
            "Dur_13.8.10_F3")
  
  west <- c("Alaska_21.1.1_M1",
            "Port_4_M2",
            "Port_8_M1",
            "Alaska_21.1.8_F1",
            "Port_11_F2",
            "Port_7_F1")
  
  populations <- list(east,west)
  
  
  slide.GENOME.class <- sliding.window.transform(GENOME.class,width=100000,jump=10000,type=2)
  
  slide.GENOME.class <- set.populations(slide.GENOME.class, populations,
                                        diploid = FALSE, triploid = FALSE, 
                                        tetraploid = FALSE)
  
  
  slide.GENOME.class <- F_ST.stats(slide.GENOME.class)
  
  FST <- t(slide.GENOME.class@nuc.F_ST.pairwise)

  
  regions <- as.data.frame(slide.GENOME.class@region.names)
  regions_split <- data.frame(do.call('rbind', strsplit(as.character(regions[,1]),' ',fixed=TRUE)))
  chr <- rep(GENOME.class@region.names,each=length(FST[,1]))
  regions_chr <- cbind(regions_split,chr,FST)
  regions_cleaned <- subset(regions_chr, select=-c(X2,X4))
  write.table(regions_cleaned, "Fst.autos.txt", append=TRUE, col.names=FALSE, sep=",")
  
}

################### Fst U chromosome ################### 

GENOME.class <- readData("Ceratodon_bwa_U_md/female", format = "VCF", include.unknown=TRUE, FAST=FALSE)

eastF<- c("Dur_16.6.2_F1",
          "Dur_4.13.7_F2",
          "Uconn_15.12.12_F1",
          "Dur_13.8.10_F3")

westF<- c("Alaska_21.1.8_F1",
          "Port_11_F2",
          "Port_7_F1")

populations <- list(eastF,westF)

slide.GENOME.class <- sliding.window.transform(GENOME.class,width=100000,jump=10000,type=2)

slide.GENOME.class <- set.populations(slide.GENOME.class, populations,
                                      diploid = FALSE, triploid = FALSE, 
                                      tetraploid = FALSE)

slide.GENOME.class <- F_ST.stats(slide.GENOME.class) 

FST <- t(slide.GENOME.class@nuc.F_ST.pairwise)

regions <- as.data.frame(slide.GENOME.class@region.names)
regions_split <- data.frame(do.call('rbind', strsplit(as.character(regions[,1]),' ',fixed=TRUE)))
chr <- rep(GENOME.class@region.names,each=length(FST[,1]))
regions_chr <- cbind(regions_split,chr,FST)
regions_cleaned <- subset(regions_chr, select=-c(X2,X4))
write.table(regions_cleaned, "Fst.female.U.txt", col.names=FALSE, sep=",")



################### Fst V chromosome ################### 

GENOME.class <- readData("Ceratodon_bwa_V_md/male", format = "VCF", include.unknown=TRUE, FAST=FALSE)

eastM <- c("Dur_19.2.8_M1",
           "Dur_4.13.5_M2",
           "Uconn_12.2.4_F2",
           "Dur_13.8.1_M3")

westM <- c("Alaska_21.1.1_M1",
           "Port_4_M2",
           "Port_8_M1")

populations <- list(eastM,westM)

slide.GENOME.class <- sliding.window.transform(GENOME.class,width=100000,jump=10000,type=2)

slide.GENOME.class <- set.populations(slide.GENOME.class, populations,
                                      diploid = FALSE, triploid = FALSE, 
                                      tetraploid = FALSE)

slide.GENOME.class <- F_ST.stats(slide.GENOME.class) 

FST <- t(slide.GENOME.class@nuc.F_ST.pairwise)

regions <- as.data.frame(slide.GENOME.class@region.names)
regions_split <- data.frame(do.call('rbind', strsplit(as.character(regions[,1]),' ',fixed=TRUE)))
chr <- rep(GENOME.class@region.names,each=length(FST[,1]))
regions_chr <- cbind(regions_split,chr,FST)
regions_cleaned <- subset(regions_chr, select=-c(X2,X4))
write.table(regions_cleaned, "Fst.male.V.txt", col.names=FALSE, sep=",")


#################### Sliding window plot ####################

library(karyoploteR)
#v1.8.8

r40 <- read.table("Data/R40_genome_lengths_sex.txt", header=TRUE)
popgen_data <- read.csv("Data/ceratodon_popgen_Oct29.csv", header=TRUE)

popgen.GR <- toGRanges(data.frame(chr=popgen_data$Chromosome, start=popgen_data$Start, end=popgen_data$Stop))
custom.genome <- toGRanges(data.frame(chr=r40$Chromosome, start=r40$ChromStart, end=r40$ChromEnd))

png("Figures/neutrality_windows_Oct30.png", width = 9, height = 6, units = 'in', res = 300)

kp <- plotKaryotype(genome=custom.genome,pin=8, plot.type = 4,labels.plotter = NULL)

kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 5, tick.col="black", cex=0.5,
                 minor.tick.dist = 5000000, minor.tick.len = 5, minor.tick.col = "black", add.units = F)
kp <- kpDataBackground(kp, data.panel = 1, color="gray95")

kp <-kpPlotLoess(kp, data=popgen.GR, data.panel = 1, col=("gray0"), y=popgen_data$thetaW,
                 ymin=0, ymax=0.02, span=0.03, r0=0.75, r1=0.99)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.5, numticks=3, 
       labels=c("0","0.01","0.02"), r0=0.75, r1=0.99, side=1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.5, numticks=3, 
       labels=c("0","0.01","0.02"), r0=0.75, r1=0.99, side=2)
kp <-kpPlotLoess(kp, data=popgen.GR, data.panel = 1, col=("gray20"), y=popgen_data$thetaPi,
                 ymin=0, ymax=0.02, span=0.03, r0=0.5, r1=0.7)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.5, numticks=3, 
       labels=c("0","0.01","0.02"), r0=0.5, r1=0.7, side=1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.5, numticks=3, 
       labels=c("0","0.01","0.02"), r0=0.5, r1=0.7, side=2)
kp <-kpPlotLoess(kp, data=popgen.GR, data.panel = 1, col=("gray40"), y=popgen_data$tajimaD,
                 ymin=-2, ymax=2, span=0.03, r0=0.25, r1=0.45)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.5, numticks=3, 
       labels=c("-2","0","2"), r0=0.25, r1=0.45, side=1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.5, numticks=3, 
       labels=c("-2","0","2"), r0=0.25, r1=0.45, side=2)
kp <-kpPlotLoess(kp, data=popgen.GR, data.panel = 1, col=("gray60"), y=popgen_data$Fst,
                 ymin=0, ymax=0.75, span=0.03, r0=0, r1=0.2)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.5, numticks=3, 
       labels=c("0","0.375","0.75"), r0=0, r1=0.2, side=1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.5, numticks=3, 
       labels=c("0","0.375","0.75"), r0=0, r1=0.2, side=2)

dev.off()

# some additional asthesics done in Inkscape (adding chromosome numbers, metrics to y-axis)

####################### Testing significant differences between U, V, ####################### 
### and autos by random sampling of sliding windows (i.e., bootstrapping)

neutrality_data <- read.csv("Data/ceratodon_popgen_Oct29.csv", header=TRUE)


theta_U <- subset(neutrality_data, subset=neutrality_data$Chromosome == "U")
theta_V <- subset(neutrality_data , subset=neutrality_data$Chromosome == "V")
theta_auto <- subset(neutrality_data, subset=neutrality_data$Chromosome == "Chr01" | neutrality_data$Chromosome == "Chr02" 
                     | neutrality_data$Chromosome == "Chr03" | neutrality_data$Chromosome == "Chr04"
| neutrality_data$Chromosome == "Chr05" | neutrality_data$Chromosome =="Chr06" 
| neutrality_data$Chromosome == "Chr07"| neutrality_data$Chromosome == "Chr08" 
| neutrality_data$Chromosome == "Chr09" | neutrality_data$Chromosome == "Chr10" 
| neutrality_data$Chromosome == "Chr11"| neutrality_data$Chromosome == "Chr12")

# random sampling of 1000 windows with replacement
theta_U_rand <- theta_U[sample(nrow(theta_U), replace=TRUE, size=1000), ]
theta_V_rand <- theta_V[sample(nrow(theta_V), replace=TRUE, size=1000), ]
theta_auto_rand <- theta_auto[sample(nrow(theta_auto), replace=TRUE, size=1000), ]

rand_data <- rbind(theta_U_rand,theta_V_rand,theta_auto_rand)


# test significance
pairwise.wilcox.test(rand_data$thetaW, rand_data$chr_bin, p.adjust.method = "BH")

pairwise.wilcox.test(rand_data$thetaPi, rand_data$chr_bin, p.adjust.method = "BH")

pairwise.wilcox.test(rand_data$tajimaD, rand_data$chr_bin, p.adjust.method = "BH")

pairwise.wilcox.test(rand_data$Fst, rand_data$chr_bin, p.adjust.method = "BH")


####### Generating 99% confidence intervals using randomly selected windows from above ####### 

### theta
a <- mean(theta_auto_rand$thetaW)
s <- sd(theta_auto_rand$thetaW)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right


a <- mean(theta_U_rand$thetaW)
s <- sd(theta_U_rand$thetaW)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right

a <- mean(theta_V$thetaW)
s <- sd(theta_V$thetaW)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right


### Pi
a <- mean(theta_auto_rand$thetaPi)
s <- sd(theta_auto_rand$thetaPi)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right


a <- mean(theta_U_rand$thetaPi)
s <- sd(theta_U_rand$thetaPi)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right

a <- mean(theta_V$thetaPi)
s <- sd(theta_V$thetaPi)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right

### Tajima's D
a <- mean(theta_auto_rand$tajimaD)
s <- sd(theta_auto_rand$tajimaD)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right


a <- mean(theta_U_rand$tajimaD)
s <- sd(theta_U_rand$tajimaD)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right

a <- mean(theta_V$tajimaD)
s <- sd(theta_V$tajimaD)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right

### Fst
a <- mean(theta_auto_rand$Fst)
s <- sd(theta_auto_rand$Fst)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right


a <- mean(theta_U_rand$Fst)
s <- sd(theta_U_rand$Fst)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right

a <- mean(theta_V$Fst)
s <- sd(theta_V$Fst)
n <- 1000

error <- qnorm(0.995)*s/sqrt(n)
left <- a-error
right <- a+error
left
right


################## McDonald Krietman Test ##################

################## Autosomes MKT ##################

## autosomes 1-9
for (i in 1:9) {
  
  GENOME.class <- readData(paste("Ceratodon_bwa_AC_md/",i, sep=""), format="VCF", include.unknown=TRUE, FAST=FALSE, 
                           gffpath="../GFF_split")
  
  keep <- c("Alaska_21.1.1_M1",
            "Dur_19.2.8_M1",
            "Dur_4.13.5_M2",
            "Port_4_M2",
            "Uconn_12.2.4_F2",
            "Port_8_M1",
            "Dur_13.8.1_M3",
            "Equ_E13.E3.14_M1",
            "Alaska_21.1.8_F1",
            "Dur_16.6.2_F1",
            "Dur_4.13.7_F2",
            "Port_11_F2",
            "Uconn_15.12.12_F1",
            "Dur_13.8.10_F3",
            "Equ_E13.E3.1_F1",
            "Port_7_F1")
  
  outgroup <- c("Chile_1.1_M1",
                "Chile_2.12_F1")
  
  populations <- list(keep,outgroup)
  
  # set synonymous sites
  GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="../Ceratodon_purpureus_R40_plusU.plusChloro.fasta")
  
  # set populations
  GENOME.class <- set.populations(GENOME.class, populations,
                                  diploid = FALSE, triploid = FALSE, 
                                  tetraploid = FALSE)
  
  # by gene
  genes <- splitting.data(GENOME.class, subsites="gene")
  genes <- MKT(genes,do.fisher.test=TRUE)
  gene_names <- get.feature.names(genes, chr=paste("Chr0",i, sep=""), gff.file=paste("../GFF_split/Chr0",i, sep=""))
  MKT_gene_names <- cbind(get.MKT(genes),gene_names)
  
  write.csv(MKT_gene_names, "MKT.auto.csv", append=TRUE)
}

## autosomes 10-12
for (i in 10:12) {
  
  GENOME.class <- readData(paste("Ceratodon_bwa_AC_md/",i, sep=""), format="VCF", include.unknown=TRUE, FAST=FALSE, 
                           gffpath="../GFF_split")
  
  keep <- c("Alaska_21.1.1_M1",
            "Dur_19.2.8_M1",
            "Dur_4.13.5_M2",
            "Port_4_M2",
            "Uconn_12.2.4_F2",
            "Port_8_M1",
            "Dur_13.8.1_M3",
            "Equ_E13.E3.14_M1",
            "Alaska_21.1.8_F1",
            "Dur_16.6.2_F1",
            "Dur_4.13.7_F2",
            "Port_11_F2",
            "Uconn_15.12.12_F1",
            "Dur_13.8.10_F3",
            "Equ_E13.E3.1_F1",
            "Port_7_F1")
  
  outgroup <- c("Chile_1.1_M1",
                "Chile_2.12_F1")
  
  populations <- list(keep,outgroup)
  
  # set synonymous sites
  GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="../Ceratodon_purpureus_R40_plusU.plusChloro.fasta")
  
  # set populations
  GENOME.class <- set.populations(GENOME.class, populations,
                                  diploid = FALSE, triploid = FALSE, 
                                  tetraploid = FALSE)
  
  # by gene
  genes <- splitting.data(GENOME.class, subsites="gene")
  genes <- MKT(genes,do.fisher.test=TRUE)
  gene_names <- get.feature.names(genes, chr=paste("Chr",i, sep=""), gff.file=paste("../GFF_split/Chr",i, sep=""))
  MKT_gene_names <- cbind(get.MKT(genes),gene_names)
  
  write.csv(MKT_gene_names, "MKT.auto.csv", append=TRUE)
}


################## U MKT ##################

GENOME.class <- readData("Ceratodon_bwa_U_md/female", format = "VCF", include.unknown=TRUE, FAST=FALSE, 
                         gffpath="../GFF_split")

female <- c("Alaska_21.1.8_F1",
            "Dur_16.6.2_F1",
            "Dur_4.13.7_F2",
            "Port_11_F2",
            "Uconn_15.12.12_F1",
            "Dur_13.8.10_F3",
            "Equ_E13.E3.1_F1",
            "Port_7_F1")

outgroup <- c("Chile_2.12_F1")

populations <- list(female,outgroup)

# set synonymous sites
GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="../Ceratodon_purpureus_R40_plusU.plusChloro.fasta")

# set populations
GENOME.class <- set.populations(GENOME.class, populations,
                                diploid = FALSE, triploid = FALSE, 
                                tetraploid = FALSE)

# by gene
genes <- splitting.data(GENOME.class, subsites="gene")
genes <- MKT(genes,do.fisher.test=TRUE)
gene_names <- get.feature.names(genes, chr="U", gff.file="../GFF_split/U")
MKT_gene_names <- cbind(get.MKT(genes),gene_names)

write.csv(MKT_gene_names, "MKT.U.csv", append=TRUE)

################## V MKT ##################

GENOME.class <- readData("Ceratodon_bwa_V_md/male", format = "VCF", include.unknown=TRUE, FAST=FALSE, 
                         gffpath="../GFF_split")

male <- c("Alaska_21.1.1_M1",
          "Dur_19.2.8_M1",
          "Dur_4.13.5_M2",
          "Port_4_M2",
          "Uconn_12.2.4_F2",
          "Port_8_M1",
          "Dur_13.8.1_M3",
          "Equ_E13.E3.14_M1")

outgroup <- c("Chile_1.1_M1")



populations <- list(male,outgroup)

# set synonymous sites
GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="../Ceratodon_purpureus_R40_plusU.plusChloro.fasta")

# set populations
GENOME.class <- set.populations(GENOME.class, populations,
                                diploid = FALSE, triploid = FALSE, 
                                tetraploid = FALSE)

# by gene
genes <- splitting.data(GENOME.class, subsites="gene")
genes <- MKT(genes,do.fisher.test=TRUE)
gene_names <- get.feature.names(genes, chr="V", gff.file="../GFF_split/V")
MKT_gene_names <- cbind(get.MKT(genes),gene_names)

write.csv(MKT_gene_names, "MKTtest.V.csv", append=TRUE)


################## DoS plot ##################

library(ggplot2)

dos_data <- read.csv("Data/ceratodon_MKT_sig.csv", header=TRUE)

png("Figures/DoS_Oct2020.png", width = 8, height = 8, units = 'in', res = 300)

b <- ggplot(dos_data, (aes(x = Chromosome, y = DoS, fill=Chromosome))) +
  geom_boxplot(color=c("gray65","gray0", "gray50"), fill="white", alpha=0, lwd=2)
b  + scale_fill_manual(values=c("peru", "midnightblue","black")) +
  geom_abline(intercept = 0, slope = 0, lwd=2) +
  xlab("") +
  ylab("DoS") +
  theme(axis.title.x = element_text(size=30)) +
  theme(axis.title.y = element_text(size=30))+
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=20)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.75, jitter.height = 0), 
             aes(group=Chromosome), alpha=0.75, shape=19, size=2, colour="gray25") + 
  theme(legend.position="none")
 
dev.off()

