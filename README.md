---
title: Analyzing Human Genome Variation
author: Cassidy Enloe
date: March 2, 2021
---

## 1. Importing Allele Frequency Data

```
> options(width=150)
> freq <- read.table("freq.df", header=TRUE)
> head(freq)
  pop  dist superpop     lat      long                  popname CHROM      POS N_ALLELES N_CHR  Allele_A Allele_G
1 ACB 13.19      AFR 13.1776  -59.5412       African_Carib_BBDS    15 48426484         2   192 0.1041670 0.895833
2 ASW -8.78      AFR 36.1070 -112.1130  African_Ancestry_SW_USA    15 48426484         2   122 0.1885250 0.811475
3 BEB 23.68      SAS 23.6850   90.3563    Bengali_in_Bangladesh    15 48426484         2   172 0.5348840 0.465116
4 CDX 22.01      EAS 22.0088 -100.7971              Chinese_Dai    15 48426484         2   186 0.0000000 1.000000
5 CEU 62.28      EUR 39.3210 -111.0937 Utah_Resid_from_NWEurope    15 48426484         2   198 1.0000000 0.000000
6 CHB 23.13      EAS 39.9042  116.4074              Han_Chinese    15 48426484         2   206 0.0291262 0.970874
> names(freq)
 [1] "pop"       "dist"      "superpop"  "lat"       "long"      "popname"   "CHROM"     "POS"       "N_ALLELES" "N_CHR"     "Allele_A"  "Allele_G" 
> dim(freq)
[1] 26 12
```

## 2. Accessing Rows and Columns

```
> freq[1:3,]
  pop  dist superpop     lat      long                 popname CHROM      POS N_ALLELES N_CHR Allele_A Allele_G
1 ACB 13.19      AFR 13.1776  -59.5412      African_Carib_BBDS    15 48426484         2   192 0.104167 0.895833
2 ASW -8.78      AFR 36.1070 -112.1130 African_Ancestry_SW_USA    15 48426484         2   122 0.188525 0.811475
3 BEB 23.68      SAS 23.6850   90.3563   Bengali_in_Bangladesh    15 48426484         2   172 0.534884 0.465116
> freq[,1:3]
   pop  dist superpop
1  ACB 13.19      AFR
2  ASW -8.78      AFR
3  BEB 23.68      SAS
4  CDX 22.01      EAS
5  CEU 62.28      EUR
6  CHB 23.13      EAS
7  CHS 24.48      EAS
8  CLM  6.24      AMR
9  ESN 10.22      AFR
10 FIN 61.92      EUR
11 GBR 56.49      EUR
12 GIH 22.26      SAS
13 GWD 13.44      AFR
14 IBS 40.46      EUR
15 ITU 11.13      SAS
16 JPT 35.69      EAS
17 KHV 14.06      EAS
18 LWK -0.02      AFR
19 MSL  8.46      AFR
20 MXL 23.63      AMR
21 PEL -9.19      AMR
22 PJL 31.55      SAS
23 PUR 18.22      AMR
24 STU  7.87      SAS
25 TSI 43.77      EUR
26 YRI 10.16      AFR
```

## 3. Summary Statistics

```
> summary(freq$Allele_A)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00000 0.04065 0.38382 0.43479 0.77825 1.00000 
> summary(freq$Allele_G)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0000  0.2218  0.6162  0.5652  0.9594  1.0000 
> range(freq$lat)
[1] -9.1900 61.9241
> range(freq$long)
[1] -112.1130  139.6917
```

## 4. Histograms

```
> pdf("hist_AlleleA.pdf", width=8, height=8)
> hist(freq$Allele_A)
> dev.off()
```
<center>
<img src="hist_AlleleA.pdf" width=500px></img>
</center>

Adjust the bin size.
```
> pdf("hist_AlleleA_bins.pdf", width=8, height=8)
> hist(freq$Allele_A,breaks=20)
> dev.off()
```
<center>
<img src="hist_AlleleA_bins.pdf" width=500px></img>
</center>

## 5. Scatter Plots

```
> pdf("AlleleA_Scatter_open.pdf", width=8, height=8)
> plot(freq$Allele_A, freq$lat)
> dev.off()
```

<center>
<img src="AlleleA_Scatter_open.pdf" width=500px></img>
</center>


Adjust the point style.
```
> pdf("AlleleA_Scatter_closed.pdf", width=8, height=8)
> plot(freq$Allele_A, freq$lat, xlab="f(A) rs1426654", ylab="Latitude", pch=16, cex=0.8, col="salmon")
> dev.off()
```
<center>
<img src="AlleleA_Scatter_closed.pdf" width=500px></img>
</center>

Assign point color based on superpopulation.
```
> myColors <- c(AFR="red", AMR="blue", EAS = "darkgreen", EUR = "salmon", SAS="black")
> pdf("rs1426654_freq.pdf", width=8, height=8) 
> plot(freq$Allele_A, freq$lat, xlab="Freq(A) rs1426654", ylab="Latitude", pch=16, cex=0.8, col=myColors[freq$superpop], xlim=c(0,1), main="Latitudinal Variation in rs1426654 among 26 human populations")
> legend("topleft", c("African", "Admixed American", "East Asian", "European", "South Asian"), cex=0.8, col=c("red", "blue", "darkgreen", "salmon", "black"), pch=16, inset=0.02)
> dev.off()
```

<center>
<img src="rs126654_freq.pdf" width=500px></img>
</center>


## 6. Plotting Data on Geological Maps

Packages to install.
```
install.packages("maps")
install.packages("mapdata")
install.packages("scales")
install.packages("mapplots")
library(maps)
library(mapdata)
library(scales)
library(mapplots)
```

Drawing the world map outline.

```
> pdf("WorldMap_outline.pdf", width=8, height=8)
> map('worldHires', xlim=c(-120,142), ylim=c(-12,72), col='gray', fill=FALSE)
> box()
> dev.off()
```
<center>
<img src="WorldMap_outline.pdf" width=500px></img>
</center>



Plotting the human population.

```
> pdf("WorldMap_pop.pdf", width=8, height=8)
> map('worldHires', xlim=c(-120,142), ylim=c(-12,72), col='gray', fill=FALSE)
> points(freq$long, freq$lat, pch=16, col="salmon")
> box()
> dev.off()
```
<center>
<img src="WorldMap_pop.pdf" width=500px></img>
</center>


Adjusting the point size to correlate to the allele frequency.
```
> pdf("WorldMap_adjpop.pdf", width=8, height=8)
> map('worldHires', xlim=c(-120,142), ylim=c(-12,72), col='gray', fill=FALSE)
> points(freq$long, freq$lat, pch=16, cex=freq$Allele_A*1.5, col="blue")
> box()
> dev.off()
```

<center>
<img src="WorldMap_adjpop.pdf" width=500px></img>
</center>


Using pie charts to show relative allele frequencies.

```
> pdf("WorldPie_Final.pdf", width=10, height=7)
> map('worldHires', xlim=c(-120,142), ylim=c(-12,72), col='gray', fill=FALSE)
> 
> for(m in 1:26){
+     add.pie(z=c(freq$Allele_A[m], freq$Allele_G[m]), 
+             x=freq$long[m], y=freq$lat[m], 
+             radius=freq$N_CHR[m]/100, 
+             col=c(alpha("orange", 0.6), alpha("purple", 0.6)), labels="")
+     
+     m=m+1
+ }
> text(freq$long, freq$lat, labels=freq$superpop, cex=0.5, pos=1)
> box():if expand("%") == ""|browse confirm w|else|confirm w|endif

> legend('topright', bty='1', c("Freq. Allele A", "Freq. Allele G"), 
+        pch=16, col=c(alpha("orange", 0.6), alpha("blue", 0.6)), pt.cex=1, cex=0.7)
> title(main="Global Distribution of rs1426654 Alleles", font.main=1, cex.main=0.9)
> dev.off()
```
<center>
<img src="WorldPie_Final.pdf" width=500px></img>
</center>

