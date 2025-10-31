# Packages ----
needed_packages <- c("tidyverse", # package version 2.0.0
                     "ape", 
                     "geiger", 
                     "nlme", 
                     "cowplot",
                     "phytools",
                     "caper",
                     "RRPP",
                     "picante",
                     "rr2",
                     "patchwork",
                     "GGally",
                     "modelsummary",
                     "DHARMa",
                     "car",
                     "phylobase",
                     "phylosignal",
                     'foreach', # for looping construct (package version 1.5.2)
                     'doParallel' # for parallel computing (v. 1.0.17)
)

new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)

# Set your directory ----
setwd("D:/repos/hearing_allometry")
# Load ----
load("raw_data.RData")

# individual data
anura_individual_data
# species mean data
anura_species_data
# phylogeny data
anurantree <- read.nexus("phylogeny/tree_amphibia.tre") # phylogeny
# remove species duplicated
anurantree <- drop.tip(anurantree, tip = c(
  "Boana_albopunctata 2",
  "Boana_semilineata 2",
  "Crossodactylus_caramaschii 2",              
  "Hylodes_asper 2",            
   "Hylodes_cardosoi_Iporanga_SP",             
   "Hylodes_heyeri 2",
   "Hylodes_phyllodes 2",
  "Phantasmarana_boticariana 2",
  "Phantasmarana_tamuia",        
   "Pipa_pipa",                   
   "Rhinella_ornata 2",          
   "Thoropa_taophora 2"))

# Descriptive statistics ----
# Define custom labels
custom_labels <- c("SVL" = "Snout–vent length",
                   "TA" = "Tympanum area",
                   #"DF" = "Dominant frequency",
                   #"HearingPeak" = "Hearing Peak",
                   "DeltaF2" = "HFR") # O que eh o DeltaF2?

# Create the ggpairs plot with custom labels
p <- ggpairs(
  anura_individual_data, 
  columns = c(3:5), 
  upper = list(continuous = wrap("cor", method = "spearman")),
  lower = list(continuous = wrap("points", alpha = 0.5)),
  diag = list(continuous = wrap("barDiag", bins = 20, fill = "skyblue", color = "black")),
  labeller = as_labeller(custom_labels) # apply custom labels
) + 
  theme(
    axis.text = element_text(size = 8),        
    strip.text = element_text(size = 7, face = 'bold')
  ); p

# Correlacoes moderadas e altas entre as variaveis
cor(anura_individual_data[ , c("SVL", "DeltaF2", "TA")], 
    method = "spearman", use = "complete.obs")
# Altamente e significativamente correlacionadas
cor.test(anura_individual_data$SVL, anura_individual_data$TA, method = "spearman") # 0.91
usdm::vif(anura_individual_data[ , c("SVL", "DeltaF2", "TA")])

# Save the image
dir.create('figures') # create folder to store images
ggsave(paste0(getwd(), "/figures/FigS1.ResponseCorr.pdf"), plot=p,
       width=7, height=5, units="in", dpi = "print")
ggsave(paste0(getwd(), "/figures/FigS1.ResponseCorr.png"), plot=p,
       width=7, height=5, units="in", dpi = "print", bg = 'white')

# PGLS ----
anura_data <- anura_species_data %>% 
  mutate(
    TA_raz = scale(TA/SVL),
    DeltaF2 = scale(DeltaF2),
    HearingPeak = scale(HearingPeak),
    DF = scale(DF)) 

# H1
# The peak frequency of hearing  corresponds spectrally to the 
# call's dominant frequency, as predicted by the Matched Filter Hypothesis
anura_data_edit <- anura_data %>% remove_missing() # remove P. boticariana
anurantree_edit <- drop.tip(anurantree,
                            tip = "Phantasmarana_boticariana")

# PGLS nao funciona com Pagel, ta dificil entender
# o que ta rolando, mas eh algo na estrutura de correlacao de Pagel, mudei para
# brownian e funciona...
pgls1 <- gls(HearingPeak ~ DF,
             correlation = corBrownian(phy = anurantree_edit, form = ~species),
             data = anura_data_edit,
             method = "ML")

summary(pgls1)

# Residuos padronizados
par(mfrow = c(1, 2))
# Residuos vs preditor
res <- residuals(pgls1, type = "normalized")
plot(anura_data_edit$DF, res, xlab = "Fitted values", 
     ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)
# QQ-plot
qqnorm(res); qqline(res)

# H2 
# The frequency range of hearing is related positively to the 
# relative tympanum area in macroevolutionary (interspecific level) scenarios.
pgls2 <- gls(DeltaF2 ~ TA_raz,
             correlation = corBrownian(phy = anurantree, form = ~species),
             data = anura_data, method = "ML")
summary(pgls2)

# Residuos padronizados
par(mfrow = c(1, 2))
# Residuos vs preditor
res <- residuals(pgls2, type = "normalized")
plot(anura_data$TA_raz, res, xlab = "Fitted values", 
     ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)
# QQ-plot
qqnorm(res); qqline(res)

# H3: 
# The frequency range of hearing is related positively to the 
# relative tympanum area in micro (interspecific level) 
# e-PGLS ----
# transform variables if needed
anurandata <- anura_individual_data %>%
  mutate(
    TA_raz = scale(TA/SVL),
    TA_resid = residuals(lm(TA ~ SVL)),
    TA = scale(TA),
    SVL = scale(SVL),
    DeltaF2 = scale(DeltaF2))
glimpse(anurandata)
length(unique(anurandata$species))
cor(anurandata$TA_resid, anurandata$TA_raz)

# Model
#hist(anurandata$SVL)
#hist(anurandata$TA)
#hist(anurandata$DeltaF2)
rdf <- rrpp.data.frame(DeltaF2 = anurandata$DeltaF2,
                       species = as.factor(anurandata$species),
                       TA = anurandata$TA,
                       TA_resid = as.numeric(anurandata$TA_resid),
                       TA_raz = as.numeric(anurandata$TA_raz),
                       SVL = anurandata$SVL)
# TA: Efeito do tamanho ajustado do tímpano (resíduo)
# species: Efeito das diferenças entre espécies (efeito fixo ou 
# estrutura filogenética)
# TA:species Interação entre tímpano ajustado e espécie 
# (ou seja, se a relação muda por espécie
# Residuals: variação não explicada
# fit model with phylogeny
fit.model <- lm.rrpp.ws(DeltaF2 ~ TA_raz * species, data = rdf,
                        subjects = "species",
                        Cov = vcv(anurantree), # filogenia entra aqui para
                        # controlar a covariancia dado o parentesco entre as especies
                        # o que quebraria a premissa de um modelo linear. 
                        delta = 0.5, # covariancia entre individuos 
                        gamma = "sample", # dados desbalanceados
                        #turbo = FALSE, verbose = TRUE,
                        print.progress = FALSE
)
anova(fit.model)
summary(fit.model)
coef(fit.model)

par(mfrow = c(1, 2))
plot(fit.model)

## Interspecific variation ----
anurandata %>%
  filter(!species %in% c("Hylodes_heyeri", "Phantasmarana_boticariana",
                         "Boana_albopunctata","Rhinella_ornata")) %>%
  ggplot(aes(x = TA_raz, y = DeltaF2, color = species)) +  
  geom_point() +
  geom_smooth(aes(fill = species), method = "lm", se = TRUE, alpha = 0.3) +
  theme_bw() +
  scale_color_manual(values = c(
    "Crossodactylus_caramaschii" = "#8fb08d",
    "Hylodes_asper"     = "#58AD53",
    "Hylodes_cardosoi"  = "#40ED37",
    "Hylodes_phyllodes" = "#526E50",
    "Thoropa_taophora"  = "#b8336a", 
    "Boana_semilineata" = "#8799e0"
  )) +
  scale_fill_manual(values = c(   # mesmo vetor usado no color
    "Crossodactylus_caramaschii" = "#8fb08d",
    "Hylodes_asper"     = "#58AD53",
    "Hylodes_cardosoi"  = "#40ED37",
    "Hylodes_phyllodes" = "#526E50",
    "Thoropa_taophora"  = "#b8336a", 
    "Boana_semilineata" = "#8799e0"
  )) +
  xlab("Tympanum-body ratio ") + 
  ylab("Frequency Range of Hearing (HFR)") +
  theme(
    legend.position = c(0.72, 0.10),
    legend.justification = "center",
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

# Ancestral reconstructing ----
# We calculated the tympanum allometric residual as the deviation from the 
# expected tympanum size given body size, based on a log-log regression. 
# This variable represents the allometric scaling of the auditory system, 
# positive values indicate relatively larger tympana than expected for a given body size
# and negative values indicate relatively smaller ones.
anura_data <- anura_species_data %>% 
  mutate(
    TA_raz = log(TA/SVL), 
    DeltaF2 = log(DeltaF2)) 

# Tympanum allometry
species_ta <- setNames(anura_data$TA_raz, anura_data$species)

# residuals TA
map_ta <- contMap(anurantree, species_ta,
                  lims = c(-2.3437,-1.467), # para ajustar o plot
                  plot = FALSE)
map_ta <- setMap(map_ta, hcl.colors(n=100))

# Delta F2
species_deltaf2 <- setNames(anura_data$DeltaF2, anura_data$species)
map_deltaf2 <- contMap(anurantree, species_deltaf2,
                       lims = c(7.1426,8.6788),
                                plot = FALSE)
map_deltaf2 <- setMap(map_deltaf2, hcl.colors(n=100))

# — Bloco de plotagem —
dev.off()
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot(map_deltaf2, direction="rightwards", lwd=6, max.width=0.8,
     ftype="off", border="black", legend=60,
     leg.txt="Hearing Frequence \nRange (log)")
errorbar.contMap(map_deltaf2, lwd=8)

plot.new()
plot.window(xlim = c(-0.1, 0.1),
            ylim = get("last_plot.phylo", envir = .PlotPhyloEnv)$y.lim)
lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_positions <- lp$yy[1:Ntip(map_deltaf2$tree)]

text(
  x = rep(0, Ntip(map_deltaf2$tree)),
  y = tip_positions,
  labels = gsub("_", " ", map_deltaf2$tree$tip.label),
  font = 3, cex = 2
)

plot(map_ta, direction="leftwards", lwd=6, max.width=0.8,
     ftype="off", border="black", legend=100,
     leg.txt="Relative \ntympanum area (log)")
errorbar.contMap(map_ta, lwd=8)
dev.off() # "desliga" o mfrow linha 421

## Phylomorphospace ----
#par(mar=c(5.1,5.1,1.1,1.1),
#    cex.axis=0.7,cex.lab=0.9)
#anura_data <- data.frame(anura_data)
#rownames(anura_data) <- anura_data$species
#
### graph phylomorphospace projection
#phylomorphospace(anurantree,
#                 anura_data[,c("TA_raz", "DeltaF2")],
#                 colors=cols,bty="l",node.by.map=TRUE,
#                 node.size=c(0,1.2),
#                 xlab = "Relative tympanum area (TA)",
#                 ylab="Hearing Frequence \nRange (HFR)")
#?phylomorphospace
### overlay points onto the phylomorphospace plot
#points(anura_data$TA_raz,anura_data$DeltaF2,pch=21,bg="gray",cex=1.2)
### add gridlines
#grid()
### clip plot
#clip(min(anura_data$TA),max(anura_data$TA),
#     min(anura_data$DeltaF2),max(anura_data$DeltaF2))

