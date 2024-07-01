# Associated with Ward et al. (2024) "Home on the range: a multi-isotope investigation of ungulate resource partitioning at Ashfall Fossil Beds, Nebraska, USA" _Palaeogeography, Palaeoclimatology, Paleoecology_

dir <- c("")
setwd(dir)

# install packages
{
  install.packages("PMCMRplus", dependencies = T)
  install.packages("MASS", dependencies = T)
  install.packages("shipunov", dependencies = T)
}

# load packages
{
  library("PMCMRplus")
  library("MASS")
  library("shipunov")
}

# load in data 
{
  AFB_data <- read.csv("AFB_data.csv", header = T, sep = ",", stringsAsFactors = F, fileEncoding = "UTF-8")
  Oxygen_data <- read.csv("Oxygen_data.csv", header = T, fileEncoding = "UTF-8", stringsAsFactors = F)
  MacFadden_data <- read.csv("MacFadden_data.csv", header = T, fileEncoding = "UTF-8", stringsAsFactors = F)
}

# objects to make plotting easier
{
  taxa.list <- c("Pliohippus pernix", "Cormohipparion occidentale","Pseudhipparion gratum", "Procamelus grandis", "Protolabis heterodontus", "Longirostromeryx wellsi", "Teleoceras major")
  C.lim <- c(-10,-6.9)
  O.lim.VPDB <- c(-8.1,-1.3)
  VPDB.lim <- c(-10, -0.8)
  O.lim.VSMOW <- c(22.6,29.6)
  Sr.lim <- c(0.7086, 0.709)
  color.list <- c("darkorange", "deepskyblue", "indianred3", "aquamarine3", "paleturquoise4", "firebrick4", "lightgoldenrod4", "black", "black")
  color.list.full <- c("darkorange", "deepskyblue", "indianred3", "aquamarine3", "paleturquoise4", "firebrick4", "lightgoldenrod4", "khaki3", "lightskyblue3", "olivedrab", "saddlebrown", "lavender", "slateblue1")
  Sr.iso.uncertainty <- 0.00003
  C.lab <- expression(bold(delta^13*"C"["VPDB"]* " (\u2030)"))
  O.lab <- expression(bold(delta^18*"O"["VPDB"]* " (\u2030)"))
  Sr.lab <- expression(bold(""^87*"Sr/"^86*"Sr"))
}

# carbon predictive model uncertainty
{
  # calculating enrichment values (epsilon) and 90CI
  {
    # equations from Tejada-Lara et al. (2018)
    foregut.epsilon <- function(BM){exp(2.34 + 0.05 * log(BM))}
    hindgut.epsilon <- function(BM){exp(2.42 + 0.032 * log(BM))}
    
    # calculate each species' epsilon
    Teleo.ep <- hindgut.epsilon(950) # middle point of Mead (2000)
    Cormohip.ep <- hindgut.epsilon(150) # MacFadden (1986)
    Pseudhip.ep <- hindgut.epsilon(100) # Parker et al. (2018)
    Pliohip.ep <- hindgut.epsilon(150) # MacFadden (1986)
    Procam.ep <- foregut.epsilon(310) # Jerison (1971)
    Protolab.ep <- foregut.epsilon(95) # Jerison (1971)
    Longi.ep <- foregut.epsilon(15) # Clementz et al. (2008)
    
    # make list of all values
    weighted.ep <- c(rep(Teleo.ep, times = 13), rep(Cormohip.ep, times = 3), rep(Pseudhip.ep, times = 3), rep(Pliohip.ep, times = 3), Procam.ep, Protolab.ep, rep(Longi.ep, times = 3))
    
    # calculations
    weighted.mean <- mean(weighted.ep)
    weighted.sd <- sd(weighted.ep)
    weighted.SEM <- weighted.sd / length(weighted.ep)^0.5
    weighted.90CI <- 1.64 * weighted.SEM
  }
  
  # composite uncertainty
  {
    # Tipple et al. (2010) 90CI for 11.8 Ma 3 Myr rolling average is 0.34
    (0.34^2 + weighted.90CI^2)^0.5
  }
}

# non-parametric analyses
{
  # K-W, post hoc, and F-K tests BY GENUS
    
    # carbon
    kruskal.test(formula = d13C ~ as.factor(taxon), data = AFB_data)
    summary(dscfAllPairsTest(formula = d13C ~ as.factor(taxon), data = AFB_data))
    fligner.test(formula = d13C ~ as.factor(taxon), data = AFB_data)
    
    # oxygen
    kruskal.test(formula = d18O.VPDB ~ as.factor(taxon), data = AFB_data)
    summary(dscfAllPairsTest(formula = d18O.VPDB ~ as.factor(taxon), data = AFB_data))
    fligner.test(formula = d18O.VPDB ~ as.factor(taxon), data = AFB_data)
    
    # strontium
    kruskal.test(formula = Sr.iso ~ as.factor(taxon), data = AFB_data)
    summary(dscfAllPairsTest(formula = Sr.iso ~ as.factor(taxon), data = AFB_data))
    fligner.test(formula = Sr.iso ~ as.factor(taxon), data = AFB_data)
    
  # K-W, post hoc, and F-K tests BY FAMILY
    
    # carbon
    kruskal.test(formula = d13C ~ as.factor(family), data = AFB_data)
    summary(dscfAllPairsTest(formula = d13C ~ as.factor(family), data = AFB_data))
    fligner.test(formula = d13C ~ as.factor(family), data = AFB_data)
    
    # oxygen
    kruskal.test(formula = d18O.VPDB ~ as.factor(family), data = AFB_data)
    summary(dscfAllPairsTest(formula = d18O.VPDB ~ as.factor(family), data = AFB_data))
    fligner.test(formula = d18O.VPDB ~ as.factor(family), data = AFB_data)
    
    # strontium
    kruskal.test(formula = Sr.iso ~ as.factor(family), data = AFB_data)
    summary(dscfAllPairsTest(formula = Sr.iso ~ as.factor(family), data = AFB_data))
    fligner.test(formula = Sr.iso ~ as.factor(family), data = AFB_data)
  }

# Figure 2
{
    # calculations for plotting
    {
      # means per isotope; # first by species, then add families
      Sr.mean.s <- tapply(AFB_data$Sr.iso, factor(AFB_data$taxon, levels = unique(AFB_data$taxon)), mean)
      Sr.means <- c(Sr.mean.s, tapply(AFB_data$Sr.iso, factor(AFB_data$family, levels = c("Equidae", "Camelidae")), mean))
      C.mean.s <- tapply(AFB_data$d13C, factor(AFB_data$taxon, levels = unique(AFB_data$taxon)), mean)
      C.means <- c(C.mean.s, tapply(AFB_data$d13C, factor(AFB_data$family, levels = c("Equidae", "Camelidae")), mean))
      O.mean.s <- tapply(AFB_data$d18O.VPDB, factor(AFB_data$taxon, levels = unique(AFB_data$taxon)), mean)
      O.means <- c(O.mean.s, tapply(AFB_data$d18O.VPDB, factor(AFB_data$family, levels = c("Equidae", "Camelidae")), mean))
      # standard deviations
      Sr.s.d <- tapply(AFB_data$Sr.iso, factor(AFB_data$taxon, levels = unique(AFB_data$taxon)), sd)
      Sr.sd <- c(Sr.s.d, tapply(AFB_data$Sr.iso, factor(AFB_data$family, levels = c("Equidae", "Camelidae")), sd))
      C.s.d <- tapply(AFB_data$d13C, factor(AFB_data$taxon, levels = unique(AFB_data$taxon)), sd)
      C.sd <- c(C.s.d, tapply(AFB_data$d13C, factor(AFB_data$family, levels = c("Equidae", "Camelidae")), sd))
      O.s.d <- tapply(AFB_data$d18O.VPDB, factor(AFB_data$taxon, levels = unique(AFB_data$taxon)), sd)
      O.sd <- c(O.s.d, tapply(AFB_data$d18O.VPDB, factor(AFB_data$family, levels = c("Equidae", "Camelidae")), sd))
    }
    
    # plotting
    {
      # cairo_pdf("figure 2.pdf", width = 4, height = 8)
      # cairo sometimes doesn't work on some computers
      # the dimensions used in the publication are width = 4, height = 8
      layout(matrix(c(1,2,3), ncol = 1, byrow = T), heights = c(5,5,1.5))
      # C v O
      par(mar = c(4,5,1,1), cex.lab = 1.5)
      plot(C.means, O.means, col = color.list, pch = c(17, 9, 17, 4, 4, 17, 3, 17, 4), xlab = C.lab, ylab = O.lab, ylim = O.lim.VPDB, xlim = C.lim, cex = 2, cex.axis = 1.25)
      points(x = AFB_data$d13C, y = AFB_data$d18O.VPDB, col = AFB_data$color, pch = AFB_data$pch)
      for (i in 1:10) {
        C0 <- C.means[i] - C.sd[i]
        C1 <- C.means[i] + C.sd[i]
        O0 <- O.means[i] - O.sd[i]
        O1 <- O.means[i] + O.sd[i]
        arrows(C0, O.means[i], C1, O.means[i], col = color.list[i], length = 0.05, angle = 90, code = 3)
        arrows(C.means[i], O0, C.means[i], O1, col = color.list[i], length = 0.05, angle = 90, code = 3) }
      # Sr v O
      par(mar = c(4,5,1,1), cex.lab = 1.5)
      plot(Sr.means, O.means, col = color.list, pch = c(17, 9, 17, 4, 4, 17, 3, 17, 4), ylab = O.lab, xlab = Sr.lab, xlim = Sr.lim, ylim = O.lim.VPDB, cex = 2, cex.axis = 1.25)
      points(y = AFB_data$d18O.VPDB, x = AFB_data$Sr.iso, col = AFB_data$color, pch = AFB_data$pch)
      for (i in 1:10) {
        O0 <- O.means[i] - O.sd[i]
        O1 <- O.means[i] + O.sd[i]
        Sr0 <- Sr.means[i] - Sr.sd[i]
        Sr1 <- Sr.means[i] + Sr.sd[i]
        arrows(y0 = O0, x0 = Sr.means[i], y1 = O1, x1 = Sr.means[i], col = color.list[i], length = 0.05, angle = 90, code = 3)
        arrows(y0 = O.means[i], x0 = Sr0, y1 = O.means[i], x1 = Sr1, col = color.list[i], length = 0.05, angle = 90, code = 3)}
      # legend
      par(mar = c(0,0.4,0,0.4))
      plot.new()
      legend("center", legend = c(expression(bold("Rhinocerotidae")), expression(italic("Teleoceras major")), "male M3 (N=5)", "female M3 (N=8)", "",expression(bold("Camelidae")), "Procamelus (N=1)", "Protolabis (N=1)", expression(bold("Blastomerycinae")), "Longirostromeryx (N=3)", expression(bold("Equidae")), "Cormohipparion (N=3)", "Pliohippus (N=3)", "Pseudhipparion (N=3)"), col = c("black", color.list[7], color.list[7], color.list[7], NA, "black", color.list[4:5], "black", color.list[2], "black", color.list[c(1,3,6)]), pch = c(3, 15, 15, 16, NA, 4, 4, 4, 9, 9, 17, 17, 17, 17), cex = 0.9, ncol = 3, pt.cex = 1.4, xpd = T, bty = "n")
      # dev.off()
    }
  }

# Linear Discriminant Analysis
{
  lda.data <- AFB_data[,c("family", "Sr.iso", "d13C", "d18O.VPDB")]
  # CV = cross-validation
  CV.lda <- lda(family ~ Sr.iso + d13C + d18O.VPDB, AFB_data, tol = 7.0e-06, CV = T)
  lda <- lda(family ~ Sr.iso + d13C + d18O.VPDB, AFB_data, tol = 7.0e-06)
  tab <- table(AFB_data$family, CV.lda$class, dnn = c('Actual Group','Predicted Group'))
  prediction <- predict(lda, lda.data) # used in plotting
  # tab 2 isn't necessary. It doesn't include the cross-validation
  tab2 <- table(AFB_data$family, prediction$class, dnn = c('Actual Group','Predicted Group'))
  }

# Figure 3
{
  # I got this function from someone on the internet, I can't remember who. Thank you anonymous person. Note that I made them all the same length post-R. 
  lda.arrows <- function(x, myscale = 1, tex = 0.75, choices = c(1,2), ...){
    ## adds `biplot` arrows to an lda using the discriminant function values
    heads <- coef(x)
    arrows(x0 = 0, y0 = 0, 
           x1 = myscale * heads[,choices[1]], 
           y1 = myscale * heads[,choices[2]], ...)
    text(myscale * heads[,choices], labels = row.names(heads), 
         cex = tex)}
  
  # cairo_pdf("figure 3.pdf", width = 5, height = 4)
  # the dimensions matching the publication are width = 5, height = 4
  par(mfrow = c(1,1), mar = c(5,4,1,1))
  plot(prediction$x[,1], prediction$x[,2], col = AFB_data$color, pch = AFB_data$pch, xlab = "LD1 (93.4%)", ylab = "LD2 (0.06%)", xlim = c(-5,6), ylim = c(-4,4))
  lda.arrows(lda, col = "black", myscale = 1)
  Ellipses(prediction$x[,1:2], AFB_data$family, level = 0.95, centers = T, c.pch = 16, c.cex = 0.75, match.color = F)
  # dev.off()
  # legend is the same figure 2, added post-R
}

# Clementz et al. (2008) Semi-aquatic test
{
  # calculations
  {
    # convert MacFadden data from VPDB to VSMOW 
    VPDBtoVSMOW <- function(VPDB.values){VPDB.values * 1.03091 + 30.91}
    
    MacFadden_data$d18O.VSMOW <- round(VPDBtoVSMOW(MacFadden_data$d18O.VPDB), digits = 1)
    Oxygen.fig.data <- rbind(Oxygen_data, MacFadden_data[,c(1:3,5)]) 
    
    O.loc <- unique(Oxygen.fig.data$Locality)
    # locality name abbreviations:
      # HR = Hottel Ranch
      # MF = Myers Farm
      # NS = North Shore
      # PS = Pratt Slide
      # Cam = Cambridge
      # Clem = Clementz et al. (2008)'s 11.8 Ma Nebraska
      # Ash = this study (Ashfall)
      # Opt = Optima
      # LBB = Love Bone Bed
      # McG = McGehee Farm
      # Wth = Withlacoochie
      # BV = Bone Valley
    
    O.families <- unique(Oxygen.fig.data$family)
    O.families <- O.families[O.families != "T.rhinocerotid"]
    
    se <- function(x){sd(x)/length(x)^.5} # why is there not a standard error function?
    
    Oxygen.localities <- data.frame()
    # these are for Figure 5
    for(i in 1:length(O.loc)){
      Oxygen.localities[i, "Locality"] <- O.loc[i]
      Oxygen.localities[i, "Teleo.n"] <- length(Oxygen.fig.data[Oxygen.fig.data$Locality == O.loc[i] & Oxygen.fig.data$Taxon == "Teleoceras", "d18O.VSMOW"])
      Oxygen.localities[i, "Teleo.mean"] <- round(mean(Oxygen.fig.data[Oxygen.fig.data$Locality == O.loc[i] & Oxygen.fig.data$Taxon == "Teleoceras", "d18O.VSMOW"]), digits = 1)
      Oxygen.localities[i, "Teleo.se"] <- round(se(Oxygen.fig.data[Oxygen.fig.data$Locality == O.loc[i] & Oxygen.fig.data$Taxon == "Teleoceras", "d18O.VSMOW"]), digits = 1)
      Oxygen.localities[i, "Fauna.n"] <- length(Oxygen.fig.data[Oxygen.fig.data$Locality == O.loc[i] & Oxygen.fig.data$Taxon != "Teleoceras", "d18O.VSMOW"])
      Oxygen.localities[i, "Fauna.mean"] <- round(mean(Oxygen.fig.data[Oxygen.fig.data$Locality == O.loc[i] & Oxygen.fig.data$Taxon != "Teleoceras", "d18O.VSMOW"]), digits = 1)
      Oxygen.localities[i, "Fauna.se"] <- round(se(Oxygen.fig.data[Oxygen.fig.data$Locality == O.loc[i] & Oxygen.fig.data$Taxon != "Teleoceras", "d18O.VSMOW"]), digits = 1)
      # these are for Figure S1
      for(j in 1:length(O.families)){
        Oxygen.localities[i, paste(O.families[j],".n", sep = "")] <- length(Oxygen.fig.data[Oxygen.fig.data$Locality == O.loc[i] & Oxygen.fig.data$family == O.families[j], "d18O.VSMOW"])
        Oxygen.localities[i, paste(O.families[j],".percent", sep = "")] <- signif(( length(Oxygen.fig.data[Oxygen.fig.data$Locality == O.loc[i] & Oxygen.fig.data$family == O.families[j], "d18O.VSMOW"]) / Oxygen.localities[Oxygen.localities$Locality == O.loc[i], "Fauna.n"] ) * 100, digits = 3)
      }
    }
  }
  
  # Figure 5
  {
    # cairo_pdf("figure 5.pdf", width = 5, height = 5)
    # the dimensions matching the publication are width = 5, height = 5
    par(mfrow = c(1,1), mar = c(5,5,3,2))
    plot(Oxygen.localities$Fauna.mean, Oxygen.localities$Teleo.mean, xlim = c(24.5, 31.5), ylim = c(23.5, 32.5), pch = c(0:10), xlab = expression(bold(delta^18*"O"["Fauna, VSMOW"]* " (\u2030)")), ylab = expression(bold(delta^18*"O"["Teleoceras, VSMOW"]* " (\u2030)")))
    abline(a = 0, b = 1, lty = 3, col = "black")
    abline(a = -1.67, b = 0.96, col = "black")
    arrows(x0 = Oxygen.localities$Fauna.mean-Oxygen.localities$Fauna.se, x1 = Oxygen.localities$Fauna.mean+Oxygen.localities$Fauna.se, y0 = Oxygen.localities$Teleo.mean, length = 0.05, angle = 90, code = 3, col = "black")
    arrows(x0 = Oxygen.localities$Fauna.mean, y0 = Oxygen.localities$Teleo.mean-Oxygen.localities$Teleo.se, y1 = Oxygen.localities$Teleo.mean+Oxygen.localities$Teleo.se, length = 0.05, angle = 90, code = 3, col = "black")
    points(Oxygen.localities$Fauna.mean, Oxygen.localities$Teleo.mean, pch = c(0:10), col = c("blue", "blue", "green", "green", "green", "skyblue", "red","purple", "orange", "orange", "orange", "orange"))
    legend("topleft", legend = Oxygen.localities$Locality, pch = c(0:14), cex = .6, pt.cex = 1, inset = c(0,0), xpd = T, ncol = 4, col = c("blue", "blue", "green", "green", "green", "skyblue", "red","purple", "orange", "orange", "orange", "orange"))
    # dev.off()
    # Color designates source of data: 
      # blue = Nguy and Secord, 2022
      # green = Kita et al., 2014
      # skyblue = Clementz et al., 2008 (11.8 Ma Nebraska)
      # red = this study
      # purple = Frederickson et al., 2022
      # orange = MacFadden and Cerling, 1996; MacFadden, 1998
  }
  
  # Supplemental Figure S1
  {
    fam.abbr <- c("A", "B", "Ca", "Ch", "D", "E", "L", "Me", "P", "R", "Tp", "Ty")
    # family name abbreviations:
      # A = Antilocapridae
      # B = Blastomerycinae
      # Ca = Camelidae
      # Ch = Chalicotheriidae
      # D = Dromomerycidae
      # E = Equidae
      # L = Leptomerycidae
      # Me = Merycoidodontidae
      # P = Palaeomerycidae
      # R = non-Teleoceras Rhinocerotidae
      # Tp = Tapiridae
      # Ty = Tayassuidae
    
    family.percents <- data.frame()
    for(i in 1:length(O.loc)){
      family.percents[i, "aaLocality"] <- O.loc[i]     
      for(j in 1:length(O.families)){
        family.percents[i, paste(O.families[j],".percent", sep = "")] <- signif(( length(Oxygen.fig.data[Oxygen.fig.data$Locality == O.loc[i] & Oxygen.fig.data$family == O.families[j], "d18O.VSMOW"]) / Oxygen.localities[Oxygen.localities$Locality == O.loc[i], "Fauna.n"] ) * 100, digits = 3)
      }
    }
    
    # alphabetize
    family.percents <- family.percents[ , order(names(family.percents))]
    
    # figure S1 with x-axis labels
    {
      # cairo_pdf("figure S1_1.pdf", width = 4, height = 12)
      # the dimensions matching the publication are width = 4, height = 12
      par(mfrow = c(12, 1))
      for(i in c(12, 8, 5, 11, 10, 3, 4, 9, 6, 7, 2, 1)){
        par(mar = c(2.2, 2, 0.5, 0))
        barplot(as.matrix(family.percents[i, c(2:13)]), ylim = c(0, 65), col = "white", names.arg = fam.abbr)
      }
      # dev.off()
    }
    
    # figure S1 without x-axis labels
    {
      # cairo_pdf("figure S1_2.pdf", width = 4, height = 12)
      # the dimensions matching the publication are width = 4, height = 12
      par(mfrow = c(12, 1))
      for(i in c(12, 8, 5, 11, 10, 3, 4, 9, 6, 7, 2, 1)){
        par(mar = c(0.75, 2, 0.5, 0))
        barplot(as.matrix(family.percents[i, c(2:13)]), ylim = c(0, 65), col = "white", axisnames = F)
      }
      # dev.off()
    }
  }
}






