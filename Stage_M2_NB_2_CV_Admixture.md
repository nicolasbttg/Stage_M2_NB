Script Rmarkdown associé :
- [2_CV_Admixture_Rmd](scripts/Stage_M2_NB_2_CV_Admixture.Rmd)
- [2_CV_Admixture_html](scripts/Stage_M2_NB_2_CV_Admixture.html)

## Chargement des packages R

```{r, include=FALSE}
library(ggplot2)
library(stringr)
library(vcfR)
library(vioplot)
library(qqman)
library(dplyr)
library(tidyr)
library(reshape2)
library(gplots)
library(pheatmap)
library(MASS)
```

Lien vers les scripts bash et commandes utilisées issus de la publication de Wragg et al., 2022 : <https://github.com/avignal5/SeqApiPop/tree/v1.5>

### Création des fichiers des échantillons SeqApiPop - BeeMuSe

#### Garder seulement les populations de références SeqApiPop - 301 Samples

```{r}
# Charger les fichiers CSV
df_labels <- read.csv('~/Documents/Stage_NB/data/SeqApiPop_labels.csv')
df_samples <- read.table('~/Documents/Stage_NB/data/Diversity_Study_629_Samples.txt', header = F)
colnames(df_samples)[colnames(df_samples) == "V1"] <- "name"

# Fusionner les DataFrames sur la colonne 'name'
merged_df <- merge(df_samples, df_labels,  by = 'name')

# Sélectionner les lignes où 'GeneticOrigin' n'est pas égal à 'unknown' et autres pour avoir reference populations
filtered_df <- merged_df[merged_df$GeneticOrigin != 'Unknown' &
                           merged_df$Label != 'Ariege Conservatory' &
                           merged_df$Label != 'Brittany Conservatory' &
                           merged_df$UniqueInHive != 'Unknown' &
                           merged_df$UniqueInHive != 'Buckfast' &
                           merged_df$GeneticOrigin != 'Buckfast', ]

# Obtenir la liste finale des noms de samples
selected_samples <- filtered_df$name

# Exporter la liste vers un fichier texte
write.table(selected_samples, file = 'Diversity_Study_RefPop_Samples.txt', col.names = FALSE, row.names = FALSE)
```

```         
sed -i 's/"//g' Diversity_Study_RefPop_Samples.txt

less Diversity_Study_RefPop_Samples.txt \| wc -l
```

-   301 samples - reference populations

#### Retirer Mellifera Ouessant / Colonsay de SeqApiPop - 561 Samples

```{r}
# Charger les fichiers CSV
df_labels <- read.csv('~/Documents/Stage_NB/data/SeqApiPop_labels.csv')
df_samples <- read.table('~/Documents/Stage_NB/data/Diversity_Study_629_Samples.txt', header = F)
colnames(df_samples)[colnames(df_samples) == "V1"] <- "name"

# Fusionner les DataFrames sur la colonne 'name'
merged_df <- merge(df_samples, df_labels,  by = 'name')

# Sélectionner les lignes où 'GeneticOrigin' n'est pas égal à 'unknown' et autres pour avoir reference populations
filtered_df <- merged_df[merged_df$Label != 'Ouessant Conservatory' &
                           merged_df$Label != 'Colonsay Conservatory', ]

# Obtenir la liste finale des noms de samples
selected_samples <- filtered_df$name

# Exporter la liste vers un fichier texte
write.table(selected_samples, file = 'Diversity_Study_561_Samples.txt', col.names = FALSE, row.names = FALSE)
```

```         
sed -i 's/"//g' Diversity_Study_561_Samples.txt

less Diversity_Study_561_Samples.txt \| wc -l
```

=\> 561 samples - sans Colonsay / Ouessant

#### Création du fichier ind2pop_merged_data 629 Samples maf 001 LD03 - 3848 SNPs

```{r}
setwd("~/Documents/Stage_NB/data/merged_BeeMuSe_SeqApiPop_629_filtered_maf001_LD03")

# Chargement des données
ind2pop_merged_data <- read.table("ind2pop_merged_data.txt", header = FALSE)

# Renommer les colonnes
colnames(ind2pop_merged_data) <- c("Sample", "name")

# Convertir les colonnes en caractères
ind2pop_merged_data$Sample <- as.character(ind2pop_merged_data$Sample)
ind2pop_merged_data$name <- as.character(ind2pop_merged_data$name)

# Appliquer la modification sur chaque ligne
for (i in 1:nrow(ind2pop_merged_data)) {
  if (!grepl("Beemuse", ind2pop_merged_data[i, "Sample"]) && !grepl("Beemuse", ind2pop_merged_data[i, "name"])) {
    next  # Si "Beemuse" n'est pas présent dans la ligne, passez à la suivante
  } else {
    ind2pop_merged_data[i, "Sample"] <- ifelse(grepl("Beemuse", ind2pop_merged_data[i, "Sample"]), "Beemuse", ind2pop_merged_data[i, "Sample"])
    ind2pop_merged_data[i, "name"] <- ifelse(grepl("Beemuse", ind2pop_merged_data[i, "name"]), "Beemuse", ind2pop_merged_data[i, "name"])
  }
}

# Supprimer la deuxième colonne
ind2pop_merged_data <- ind2pop_merged_data[, -2]

# Écrire le résultat dans un nouveau fichier
write.table(ind2pop_merged_data, "ind2pop_merged_data_2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

ind2pop_merged_data_2 <- read.table("ind2pop_merged_data_2.txt", header = FALSE)
```

# Analyse en Composantes Principales (ACP)

## BeeMuSe

### 748 échantillons - SNPsBeeMuSe filtered - 10256 SNPs

```{r}
# Chargement des données
eigenvec <- read.table("BeeMuse_filtered.eigenvec", header = F) 
#plink2.eigenvec = non filtrées / BeeMuse_filtered.eigenvec : filtre liste marqueurs
eigenval <- read.table("BeeMuse_filtered.eigenval", header = F) 
#plink2.eigenval = non filtrées / BeeMuse_filtered.eigenval : filtre liste marqueurs
colnames(eigenvec)[colnames(eigenvec) == "V1"] <- "Sample"
colnames(eigenvec)[colnames(eigenvec) == "V2"] <- "Filename"

Corres_ID_E756 <- read.table("/home/nbettembourg/Documents/Stage_NB/data/BeeMuSe_Wageningen_E756_17_03_2023_Alain/BeeMuSe_Wageningen_E756_17_03_2023_Alain/Corres_ID_E756.txt",header=F)
input_pedigree_BeeMuSe <- read.table("/home/nbettembourg/Documents/Stage_NB/data/BeeMuSe_Wageningen_E756_17_03_2023_Alain/BeeMuSe_Wageningen_E756_17_03_2023_Alain/input-pedigree_BeeMuSe.csv",header=T,sep=";")
```

```{r}
#extraire 'Pool' et '-100' obtenir 'Pool-100'
eigenvec$name <- paste(sub("Beemuse_", "", eigenvec$Sample), 
                 sub("_(.*?)\\..*", "\\1", eigenvec$Filename), 
                 sep = "-")
eigenvec$name <- str_extract(eigenvec$name, "[A-Za-z0-9]+-[0-9]+")

colnames(Corres_ID_E756)[colnames(Corres_ID_E756) == "V1"] <- "name"
colnames(Corres_ID_E756)[colnames(Corres_ID_E756) == "V2"] <- "ID_1a"

Corres_ID_E756$ID_1a <- gsub("o", "_", Corres_ID_E756$ID_1a)

# Remplacement de tous les "o" par "_" sauf quand présent dans un mot
Corres_ID_E756$ID_1a <- gsub("Pers_", "Perso", Corres_ID_E756$ID_1a)
Corres_ID_E756$ID_1a <- gsub("L_c", "Loc", Corres_ID_E756$ID_1a)

Corres_ID_E756_eigenvec <- merge(eigenvec, Corres_ID_E756,  by = 'name')
merged_3 <- merge(Corres_ID_E756_eigenvec, input_pedigree_BeeMuSe,  by = 'ID_1a')
colnames(merged_3)[colnames(merged_3) == "V3.x"] <- "V3"

ind2pop_ID_2a = subset(merged_3, select = c(name, ID_2a))

id_counts <- as.data.frame(table(merged_3$ID_2a))
names(id_counts) <- c("ID_2a", "Occurrence")
```

```{r}
# Importation des données de matrice d'apparentement dans R et visualisation par clustering hiérarchique
matrice_app <- read.table("plink2.rel", header = FALSE)
dist_matrice <- dist(matrice_app)
hc <- hclust(dist_matrice, method = "ward.D2")
plot(hc)

#heatplot(as.matrix(dist(matrice_app,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')
```

```{r}
# Proportion de la variance expliquée
eigen_percent <- round((eigenval / (sum(eigenval) )*100),2)

lambda <- eigenval$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)
```

#### PC1/PC2

```{r}
  # ACP Sample : plaque de provenance génotypage échantillon
ggplot(eigenvec, aes(x = V3, y = V4, label = Sample, color = Sample)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(legend.position = "right")  
```

```{r}
  # ACP ID_breeder : identification de l'apiculteur de chez qui proviennent les abeilles
ggplot(merged_3, aes(x = V3, y = V4, label = ID_breeder, color = ID_breeder)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(legend.position = "right") 
```

```{r}
  #plot ID_2a : reines mères des reines génotypées
ggplot(merged_3, aes(x = V3, y = V4, label = ID_2a, color = ID_2a)) +
  geom_point(size = 1) +  
  scale_color_manual(values = c("#FF0000FF", "#FF3300FF", "#FF6600FF", "#FF9900FF", "#FFCC00FF", 
                                 "#FFFF00FF", "#CCFF00FF", "#99FF00FF", "#66FF00FF", "#33FF00FF", 
                                 "#00FF00FF", "#00FF33FF", "#00FF66FF", "#00FF99FF", "#00FFCCFF", 
                                 "#00FFFFFF", "#00CCFFFF", "#0099FFFF", "#0066FFFF", "#0033FFFF", 
                                 "#0000FFFF", "#3300FFFF", "#6600FFFF", "#9900FFFF", "#CC00FFFF", 
                                 "#FF00FFFF", "#FF00CCFF", "#FF0099FF", "#FF0066FF")) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(legend.position = "right") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = -0.1, xmax = -0.16, ymin = 0.05, ymax = -0.01)


ggplot(merged_3, aes(x = V3, y = V4, label = ID_2a)) +
  geom_point(size = 1, color = "red") + 
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(legend.position = "right") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = -0.1, xmax = -0.16, ymin = 0.05, ymax = -0.01)
```

```{r}
# Visualisaiton ACP PC1 PC2 avec ellipses de seuil de confiance de 0.97
ggplot(merged_3, aes(x = V3, y = V4, label = ID_2a, color = ID_2a)) +
  geom_point(size = 0.8) +  
  stat_ellipse(aes(group = ID_2a), geom = "polygon", level = 0.97, alpha = 0, size = 0.1, color = "black") +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(legend.position = "right") 

#15 familles bien regroupés entre elles
ggplot(subset(merged_3, ID_2a %in% c("BAH_20-19", "BBS_6-19","BER_11-19", "BLS_53-19","BHA_2-20",  "KLSU_14-19", "MM_31-20", "MP_10-20", "PersoLD_2021", "PersoLD_2022", "S_GZ_2-19", "SJ_16-20", "SJ_24-20", "TL_13-20", "TL_19-20")), aes(x = V3, y = V4, color = ID_2a)) +
  geom_point(size = 1) +
  stat_ellipse(aes(group = ID_2a, fill = ID_2a), geom = "polygon", level = 0.97, alpha = 0.15, size = 0.15, color = "black") +
  #stat_ellipse(aes(group = ID_2a), geom = "polygon", level = 0.97, alpha = 0, size = 0.1, color = "black") +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(legend.position = "right")

#14 familles moins bien regroupées entre elles
ggplot(subset(merged_3, ID_2a %in% c( "BH_44","BH_7-19", "ER_13-19","MM_37-20", "KBJ_1-19", "KLoc_37-19","KBru_6-20", "PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoUB_2021", "PersoUB_2022", "SJ_30-20","SBJ_3-19")), aes(x = V3, y = V4, color = ID_2a)) +
  geom_point(size = 1) +
  stat_ellipse(aes(group = ID_2a, fill = ID_2a), geom = "polygon", level = 0.97, alpha = 0.15, size = 0.15, color = "black") +
  #stat_ellipse(aes(group = ID_2a), geom = "polygon", level = 0.97, alpha = 0, size = 0.1, color = "black") +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(legend.position = "right")

```

```{r}
  # ACP selon variable Sample (plaques)
ggplot(merged_3, aes(x = V3, y = V4, label = Sample, color = Sample)) +
  geom_point() +
 # stat_ellipse(aes(group = Sample), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  # stat_ellipse(aes(group = Sample, fill = Sample), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(legend.position = "right") 

# 748 -> 681 -> 612
values_in_A_not_in_B <- eigenvec$name[!(eigenvec$name %in% merged_3$name)]
length(values_in_A_not_in_B) 
```

#### PC3/PC4

```{r}
 # ACP PC3 PC4 - ID_2a : reines mères des reines génotypées
ggplot(merged_3, aes(x = V5, y = V6, label = ID_2a, color = ID_2a)) +
  geom_point(size = 1) + 
  scale_color_manual(values = c("#FF0000FF", "#FF3300FF", "#FF6600FF", "#FF9900FF", "#FFCC00FF", 
                                 "#FFFF00FF", "#CCFF00FF", "#99FF00FF", "#66FF00FF", "#33FF00FF", 
                                 "#00FF00FF", "#00FF33FF", "#00FF66FF", "#00FF99FF", "#00FFCCFF", 
                                 "#00FFFFFF", "#00CCFFFF", "#0099FFFF", "#0066FFFF", "#0033FFFF", 
                                 "#0000FFFF", "#3300FFFF", "#6600FFFF", "#9900FFFF", "#CC00FFFF", 
                                 "#FF00FFFF", "#FF00CCFF", "#FF0099FF", "#FF0066FF")) +
  labs(title = "PCA Plot - BeeMuSe", x = "PC3", y = "PC4") +
  theme(legend.position = "right") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = -0.055, xmax = -0.09, ymin = -0.04, ymax = -0.1)

#15 bien regroupés
ggplot(subset(merged_3, ID_2a %in% c("BAH_20-19", "BBS_6-19","BER_11-19", "BLS_53-19","BHA_2-20",  "KLSU_14-19", "MM_31-20", "MP_10-20", "PersoLD_2021", "PersoLD_2022", "S_GZ_2-19", "SJ_16-20", "SJ_24-20", "TL_13-20", "TL_19-20")), aes(x = V5, y = V6, color = ID_2a)) +
  geom_point(size = 1) +
  stat_ellipse(aes(group = ID_2a, fill = ID_2a), geom = "polygon", level = 0.97, alpha = 0.15, size = 0.15, color = "black") +
  #stat_ellipse(aes(group = ID_2a), geom = "polygon", level = 0.97, alpha = 0, size = 0.1, color = "black") +
  labs(title = "PCA Plot - BeeMuSe", x = "PC3", y = "PC4") +
  theme(legend.position = "right")

#14  autres
ggplot(subset(merged_3, ID_2a %in% c( "BH_44","BH_7-19", "ER_13-19","MM_37-20", "KBJ_1-19", "KLoc_37-19","KBru_6-20", "PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoUB_2021", "PersoUB_2022", "SJ_30-20","SBJ_3-19")), aes(x = V5, y = V6, color = ID_2a)) +
  geom_point(size = 1) +
 stat_ellipse(aes(group = ID_2a, fill = ID_2a), geom = "polygon", level = 0.97, alpha = 0.15, size = 0.15, color = "black") +
  #stat_ellipse(aes(group = ID_2a), geom = "polygon", level = 0.97, alpha = 0, size = 0.1, color = "black") +
  labs(title = "PCA Plot - BeeMuSe", x = "PC3", y = "PC4") +
  theme(legend.position = "right")

ggplot(merged_3, aes(x = V5, y = V6, label = ID_2a)) +
  geom_point(size = 1, color = "red") +  
  labs(title = "PCA Plot - BeeMuSe", x = "PC3", y = "PC4") +
  theme(legend.position = "right") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = -0.055, xmax = -0.09, ymin = -0.04, ymax = -0.1)
```

#### Correspondance ID_2a - fichier .fam

```{r}
  #  BEEMUSE - ID_2a
beemuse_samples_ID_2a <- merged_3[, c("Sample", "Filename", "ID_2a")]
nom_fichier <- "beemuse_samples_ID_2a.txt"
write.table(beemuse_samples_ID_2a, file = nom_fichier, sep = "\t", row.names = FALSE,quote=FALSE)
```

```{r}
  # FAM MERGED SEQAPIPOP BEEMUSE 561
merged_fam <- read.table("merged_BeeMuSe_SeqApiPop_561_filtered_maf001_LD_default.fam", sep = " ", stringsAsFactors = FALSE)

merged_fam <- merged_fam[, 1:2]
colnames(merged_fam)[colnames(merged_fam) == "V1"] <- "Sample"
colnames(merged_fam)[colnames(merged_fam) == "V2"] <- "Filename"
```

```{r}
# Créer une nouvelle colonne dans le fichier 2 pour stocker les valeurs de la troisième colonne
merged_fam$ID_2a <- NA

# Parcourir les lignes du fichier 2
for (i in 1:nrow(merged_fam)) {
    # Vérifier si les valeurs des colonnes Sample et Filename du fichier 2 sont présentes dans le fichier 1
    correspondance <- beemuse_samples_ID_2a$ID_2a[beemuse_samples_ID_2a$Sample == merged_fam$Sample[i] & beemuse_samples_ID_2a$Filename == merged_fam$Filename[i]]
    
    # Si une correspondance exacte est trouvée, ajouter la valeur correspondante de la colonne "ID_2a" du fichier 1 au fichier 2
    if (length(correspondance) > 0) {
        merged_fam$ID_2a[i] <- as.character(correspondance[1]) # Utilisez la correspondance trouvée
    } else {
        # Sinon, ajouter "Beemuse" dans la troisième colonne du fichier 2
        merged_fam$ID_2a[i] <- "Beemuse"
    }
}

# Créer fichier de correspondance entre les échantillons de BeeMuSe ID_2a avec le fichier .fam
correspondance <- merged_fam[, c("Sample","Filename","ID_2a")]
nom_fichier <- "correspondance_samples.txt"
write.table(correspondance, file = nom_fichier, sep = "\t", row.names = FALSE,quote=FALSE)
```

## SeqApiPop

### 629 échantillons - MAF \> 0.01

#### LD pruning = 0.3 (fenêtre de 1749 SNPS et pas de 175 bp)

```{r}
setwd("~/Documents/Stage_NB/data/maf001_LD03")

# fichiers pour SeqApiPop 629 échantillons et filtre maf001
#eigenvec_refpop <- read.table("SeqApiPop_629_maf001_acp.eigenvec", header = F)
#eigenval_refpop <- read.table("SeqApiPop_629_maf001_acp.eigenval", header = F)

# fichiers pour SeqApiPop 629 échantillons et filtre maf001 LD pruning = 0.3 (fenêtre de 1749 SNPs et pas de 175 bp)
eigenvec_refpop <- read.table("SeqApiPop_629_maf001_LD03_acp.eigenvec", header = F)
eigenval_refpop <- read.table("SeqApiPop_629_maf001_LD03_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_refpop)[colnames(eigenvec_refpop) == "V2"] <- "name"
eigenvec_refpop_seq_api_labels <- merge(eigenvec_refpop, seq_api_labels, by = "name")

eigen_percent_refpop <- round((eigenval_refpop / (sum(eigenval_refpop) )*100),2)
```

```{r}
#Clustering hiérarchique
#Tree
setwd("~/Documents/Stage_NB/data/maf001_LD03")
matrice_app_refpop <- read.table("SeqApiPop_629_maf001_LD03_acp.rel", header = FALSE)

dist_matrice_refpop <- dist(matrice_app_refpop)
hc_refpop <- hclust(dist_matrice_refpop, method = "ward.D2")
plot(hc_refpop)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')
```

```{r}
#ACP
#filter 629 -> 301 RefPop
eigenvec_refpop_seq_api_labels <- eigenvec_refpop_seq_api_labels[eigenvec_refpop_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                   eigenvec_refpop_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                   eigenvec_refpop_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                   eigenvec_refpop_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                   eigenvec_refpop_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                   eigenvec_refpop_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label2 <- c("black", "mediumpurple4", "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_refpop$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)
```

##### PC1/PC2

```{r}
#PC1/PC2
# ACP avec variance expliquée
ggplot(data = eigenvec_refpop_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.44, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.083, xmax = 0.05, ymin = -0.11, ymax = -0.04)

```

```{r}
ggplot(data = eigenvec_refpop_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.44, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

```{r}
#ellipses autour des points selon Label
ggplot(data = eigenvec_refpop_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  geom_text(aes(x = 0.045, y = 0.04, label = "C-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.02, y = 0.05, label = "M-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.02, y = -0.16, label = "O-lineage"), size = 4, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.44, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 
```

```{r}
# ellipses avec couleur
ggplot(data = eigenvec_refpop_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  geom_text(aes(x = 0.045, y = 0.04, label = "M-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.02, y = 0.05, label = "C-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.02, y = -0.16, label = "O-lineage"), size = 4, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  theme(legend.position = c(0.44, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

##### PC3/PC4

```{r}
#ACP - plot variance - PC3 - PC4
ggplot(data = eigenvec_refpop_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.56, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = -0.04, xmax = -0.02, ymin = -0.16, ymax = -0.06)
```

#### LD pruning = 0.2 (fenêtre de 1749 SNPS et pas de 175 bp)

```{r}
#LD02
setwd("~/Documents/Stage_NB/data/maf001_LD02")

eigenvec_LD02 <- read.table("SeqApiPop_629_maf001_LD02_acp.eigenvec", header = F)
eigenval_LD02 <- read.table("SeqApiPop_629_maf001_LD02_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_LD02)[colnames(eigenvec_LD02) == "V2"] <- "name"
eigenvec_LD02_seq_api_labels <- merge(eigenvec_LD02, seq_api_labels, by = "name")

eigen_percent_LD02 <- round((eigenval_LD02 / (sum(eigenval_LD02) )*100),2)

#Clustering hiérarchique
#Tree
matrice_app_refpop <- read.table("SeqApiPop_629_maf001_LD02_acp.rel", header = FALSE)

dist_matrice_refpop <- dist(matrice_app_refpop)
hc_refpop <- hclust(dist_matrice_refpop, method = "ward.D2")
plot(hc_refpop)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')
```

```{r}
# filter 629 -> 301 RefPop
eigenvec_LD02_seq_api_labels <- eigenvec_LD02_seq_api_labels[eigenvec_LD02_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                   eigenvec_LD02_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                   eigenvec_LD02_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                   eigenvec_LD02_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                   eigenvec_LD02_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                   eigenvec_LD02_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label2 <- c("black", "mediumpurple4", "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_LD02$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)
```

##### PC1/PC2

```{r}
ggplot(data = eigenvec_LD02_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() + 
  theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.105, xmax = 0.06, ymin = -0.01, ymax = 0.07)

ggplot(data = eigenvec_LD02_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() + 
  theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))

# ellipses autour des points selon Label
ggplot(data = eigenvec_LD02_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() + 
  theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 

# ellipses avec couleur
ggplot(data = eigenvec_LD02_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  scale_y_reverse() + 
  theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))

```

##### PC3/PC4

```{r}
ggplot(data = eigenvec_LD02_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = -0.13, xmax = -0.18, ymin = 0.03, ymax = 0.08)

```

#### LD pruning = 0.1 (fenêtre de 1749 SNPS et pas de 175 bp)

##### PC1/PC2

```{r}
#LD01
setwd("~/Documents/Stage_NB/data/maf001_LD01")

eigenvec_LD01 <- read.table("SeqApiPop_629_maf001_LD01_acp.eigenvec", header = F)
eigenval_LD01 <- read.table("SeqApiPop_629_maf001_LD01_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_LD01)[colnames(eigenvec_LD01) == "V2"] <- "name"
eigenvec_LD01_seq_api_labels <- merge(eigenvec_LD01, seq_api_labels, by = "name")

eigen_percent_LD01 <- round((eigenval_LD01 / (sum(eigenval_LD01) )*100),2)

#Clustering hiérarchique
#Tree
matrice_app_refpop <- read.table("SeqApiPop_629_maf001_LD01_acp.rel", header = FALSE)

dist_matrice_refpop <- dist(matrice_app_refpop)
hc_refpop <- hclust(dist_matrice_refpop, method = "ward.D2")
plot(hc_refpop)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')


#ACP
#filter 629 -> 301 RefPop
eigenvec_LD01_seq_api_labels <- eigenvec_LD01_seq_api_labels[eigenvec_LD01_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                               eigenvec_LD01_seq_api_labels$Label != 'Ariege Conservatory' &
                                                               eigenvec_LD01_seq_api_labels$Label != 'Brittany Conservatory' &
                                                               eigenvec_LD01_seq_api_labels$UniqueInHive != 'Unknown' &
                                                               eigenvec_LD01_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                               eigenvec_LD01_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label2 <- c("black", "mediumpurple4", "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_LD01$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)

ggplot(data = eigenvec_LD01_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.05, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.13, xmax = 0.08, ymin = -0.02, ymax = 0.06)


ggplot(data = eigenvec_LD01_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.05, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))

#ellipses autour des points selon Label
ggplot(data = eigenvec_LD01_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.05, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 

# ellipses avec couleur
ggplot(data = eigenvec_LD01_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.05, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))

```

##### PC3/PC4

```{r}
#ACP - plot variance - PC3 - PC4
ggplot(data = eigenvec_LD01_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.5, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.09, xmax = 0.135, ymin = -0.05, ymax = 0.01)
```

#### LD pruning = 0.05 (fenêtre de 1749 SNPS et pas de 175 bp)

##### PC1/PC2

```{r}
#LD005
setwd("~/Documents/Stage_NB/data/maf001_LD005")

eigenvec_LD005 <- read.table("SeqApiPop_629_maf001_LD005_acp.eigenvec", header = F)
eigenval_LD005 <- read.table("SeqApiPop_629_maf001_LD005_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_LD005)[colnames(eigenvec_LD005) == "V2"] <- "name"
eigenvec_LD005_seq_api_labels <- merge(eigenvec_LD005, seq_api_labels, by = "name")

eigen_percent_LD005 <- round((eigenval_LD005 / (sum(eigenval_LD005) )*100),2)

#Clustering hiérarchique
#Tree
matrice_app_refpop <- read.table("SeqApiPop_629_maf001_LD005_acp.rel", header = FALSE)

dist_matrice_refpop <- dist(matrice_app_refpop)
hc_refpop <- hclust(dist_matrice_refpop, method = "ward.D2")
plot(hc_refpop)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')

#ACP
#filter 629 -> 301 RefPop
eigenvec_LD005_seq_api_labels <- eigenvec_LD005_seq_api_labels[eigenvec_LD005_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                               eigenvec_LD005_seq_api_labels$Label != 'Ariege Conservatory' &
                                                               eigenvec_LD005_seq_api_labels$Label != 'Brittany Conservatory' &
                                                               eigenvec_LD005_seq_api_labels$UniqueInHive != 'Unknown' &
                                                               eigenvec_LD005_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                               eigenvec_LD005_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label2 <- c("black", "mediumpurple4", "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_LD005$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)

ggplot(data = eigenvec_LD005_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.26, 0.64), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.077, xmax = 0.042, ymin = 0.015, ymax = 0.08)


ggplot(data = eigenvec_LD005_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.26, 0.64), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))

#ellipses autour des points selon Label
ggplot(data = eigenvec_LD005_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.25, 0.64), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 

# ellipses avec couleur
ggplot(data = eigenvec_LD005_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.25, 0.64), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

##### PC3/PC4

```{r}
#ACP - plot variance - PC3 - PC4
ggplot(data = eigenvec_LD005_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.06, xmax = 0.11, ymin = -0.06, ymax = 0)
```

#### LD pruning = 0.04 (fenêtre de 1749 SNPS et pas de 175 bp)

##### PC1/PC2

```{r}
#LD004
setwd("~/Documents/Stage_NB/data/maf001_LD004")

eigenvec_LD004 <- read.table("SeqApiPop_629_maf001_LD004_acp.eigenvec", header = F)
eigenval_LD004 <- read.table("SeqApiPop_629_maf001_LD004_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_LD004)[colnames(eigenvec_LD004) == "V2"] <- "name"
eigenvec_LD004_seq_api_labels <- merge(eigenvec_LD004, seq_api_labels, by = "name")

eigen_percent_LD004 <- round((eigenval_LD004 / (sum(eigenval_LD004) )*100),2)

#Clustering hiérarchique
#Tree
matrice_app_refpop <- read.table("SeqApiPop_629_maf001_LD004_acp.rel", header = FALSE)

dist_matrice_refpop <- dist(matrice_app_refpop)
hc_refpop <- hclust(dist_matrice_refpop, method = "ward.D2")
plot(hc_refpop)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')

#ACP
#filter 629 -> 301 RefPop
eigenvec_LD004_seq_api_labels <- eigenvec_LD004_seq_api_labels[eigenvec_LD004_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                 eigenvec_LD004_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                 eigenvec_LD004_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                 eigenvec_LD004_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                 eigenvec_LD004_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                 eigenvec_LD004_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label2 <- c("black", "mediumpurple4", "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_LD004$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)

ggplot(data = eigenvec_LD004_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.32, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.068, xmax = 0.042, ymin = -0.03, ymax = -0.11)

ggplot(data = eigenvec_LD004_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.32, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))

#ellipses autour des points selon Label
ggplot(data = eigenvec_LD004_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.42, 0.06), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 

# ellipses avec couleur
ggplot(data = eigenvec_LD004_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  theme(legend.position = c(0.42, 0.06), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

##### PC3/PC4

```{r}
#ACP - plot variance - PC3 - PC4
ggplot(data = eigenvec_LD004_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.03, xmax = 0.09, ymin = 0.06, ymax = 0.12)
```

#### LD pruning = 0.03 (fenêtre de 1749 SNPS et pas de 175 bp)

##### PC1/PC2

```{r}
#LD003
setwd("~/Documents/Stage_NB/data/maf001_LD003")

eigenvec_LD003 <- read.table("SeqApiPop_629_maf001_LD003_acp.eigenvec", header = F)
eigenval_LD003 <- read.table("SeqApiPop_629_maf001_LD003_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_LD003)[colnames(eigenvec_LD003) == "V2"] <- "name"
eigenvec_LD003_seq_api_labels <- merge(eigenvec_LD003, seq_api_labels, by = "name")

eigen_percent_LD003 <- round((eigenval_LD003 / (sum(eigenval_LD003) )*100),2)

#Clustering hiérarchique
#Tree
matrice_app_refpop <- read.table("SeqApiPop_629_maf001_LD003_acp.rel", header = FALSE)

dist_matrice_refpop <- dist(matrice_app_refpop)
hc_refpop <- hclust(dist_matrice_refpop, method = "ward.D2")
plot(hc_refpop)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')

#ACP
#filter 629 -> 301 RefPop
eigenvec_LD003_seq_api_labels <- eigenvec_LD003_seq_api_labels[eigenvec_LD003_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                 eigenvec_LD003_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                 eigenvec_LD003_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                 eigenvec_LD003_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                 eigenvec_LD003_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                 eigenvec_LD003_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label2 <- c("black", "mediumpurple4", "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_LD003$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)

ggplot(data = eigenvec_LD003_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.035, xmax = 0.065, ymin = -0.04, ymax = -0.11)


ggplot(data = eigenvec_LD003_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))

#ellipses autour des points selon Label
ggplot(data = eigenvec_LD003_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 

# ellipses avec couleur
ggplot(data = eigenvec_LD003_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

##### PC3/PC4

```{r}
#ACP - plot variance - PC3 - PC4
ggplot(data = eigenvec_LD003_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.05, xmax = 0.1, ymin = 0.08, ymax = 0.16)
```

#### LD pruning = 0.01 (fenêtre de 1749 SNPS et pas de 175 bp)

##### PC1/PC2

```{r}
#LD001
setwd("~/Documents/Stage_NB/data/maf001_LD001")

eigenvec_LD001 <- read.table("SeqApiPop_629_maf001_LD001_acp.eigenvec", header = F)
eigenval_LD001 <- read.table("SeqApiPop_629_maf001_LD001_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_LD001)[colnames(eigenvec_LD001) == "V2"] <- "name"
eigenvec_LD001_seq_api_labels <- merge(eigenvec_LD001, seq_api_labels, by = "name")

eigen_percent_LD001 <- round((eigenval_LD001 / (sum(eigenval_LD001) )*100),2)

#Clustering hiérarchique
#Tree
matrice_app_refpop <- read.table("SeqApiPop_629_maf001_LD001_acp.rel", header = FALSE)

dist_matrice_refpop <- dist(matrice_app_refpop)
hc_refpop <- hclust(dist_matrice_refpop, method = "ward.D2")
plot(hc_refpop)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')

#ACP
#filter 629 -> 301 RefPop
eigenvec_LD001_seq_api_labels <- eigenvec_LD001_seq_api_labels[eigenvec_LD001_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                 eigenvec_LD001_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                 eigenvec_LD001_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                 eigenvec_LD001_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                 eigenvec_LD001_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                 eigenvec_LD001_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label2 <- c("black", "mediumpurple4", "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_LD001$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)

ggplot(data = eigenvec_LD001_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.035, xmax = 0.065, ymin = -0.04, ymax = -0.11)


ggplot(data = eigenvec_LD001_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))

#ellipses autour des points selon Label
ggplot(data = eigenvec_LD001_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 

# ellipses avec couleur
ggplot(data = eigenvec_LD001_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))

```

##### PC3/PC4

```{r}
#ACP - plot variance - PC3 - PC4
ggplot(data = eigenvec_LD001_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = -0.04, xmax = 0, ymin = 0.08, ymax = 0.15)
```

### 561 échantillons - MAF \> 0.01

#### LD pruning = 0.3 (fenêtre de 1749 SNPS et pas de 175 bp)

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_561_maf001_LD03")

eigenvec_refpop <- read.table("SeqApiPop_561_maf001_LD03_acp.eigenvec", header = F)
eigenval_refpop <- read.table("SeqApiPop_561_maf001_LD03_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_refpop)[colnames(eigenvec_refpop) == "V2"] <- "name"
eigenvec_refpop_seq_api_labels <- merge(eigenvec_refpop, seq_api_labels, by = "name")

eigen_percent_refpop <- round((eigenval_refpop / (sum(eigenval_refpop) )*100),2)
```

```{r}
#Clustering hiérarchique
#Tree
setwd("~/Documents/Stage_NB/data/SeqApiPop_561_maf001_LD03")
matrice_app_refpop <- read.table("SeqApiPop_561_maf001_LD03_acp.rel", header = FALSE)

dist_matrice_refpop <- dist(matrice_app_refpop)
hc_refpop <- hclust(dist_matrice_refpop, method = "ward.D2")
plot(hc_refpop)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')
```

```{r}
#ACP
#filter 629 -> 301 RefPop
eigenvec_refpop_seq_api_labels <- eigenvec_refpop_seq_api_labels[eigenvec_refpop_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                   eigenvec_refpop_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                   eigenvec_refpop_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                   eigenvec_refpop_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                   eigenvec_refpop_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                   eigenvec_refpop_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label2 <- c("mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_refpop$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)
```

##### PC1/PC2

```{r}
ggplot(data = eigenvec_refpop_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.49, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.105, xmax = 0.07, ymin = -0.11, ymax = -0.04)

```

```{r}
#ellipses autour des points selon Label
ggplot(data = eigenvec_refpop_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  geom_text(aes(x = 0.05, y = 0.04, label = "M-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.01, y = 0.05, label = "C-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.02, y = -0.16, label = "O-lineage"), size = 4, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c( "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.49, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 
```

```{r}
# ellipses avec couleur
ggplot(data = eigenvec_refpop_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  geom_text(aes(x = 0.045, y = 0.04, label = "M-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.01, y = 0.05, label = "C-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.02, y = -0.16, label = "O-lineage"), size = 4, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c( "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c( "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  theme(legend.position = c(0.49, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

##### PC3/PC4

```{r}
#ACP - plot variance
ggplot(data = eigenvec_refpop_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = -0.11, xmax = -0.2, ymin = 0.07, ymax = 0.02)
```

```{r}
#ellipses autour des points selon Label
ggplot(data = eigenvec_refpop_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c( "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 

```

```{r}
# ellipses avec couleur
ggplot(data = eigenvec_refpop_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c( "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c( "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

### 629 échantillons - SNPsBeeMuSe filtered

#### No LD pruning - 10030 SNPS

##### PC1/PC2

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_629_SNPsBeeMuSe")

# SNPsBeeMuSe filtered
eigenvec_SNPsBeeMuSe <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_acp.eigenvec", header = F)
eigenval_SNPsBeeMuSe <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_SNPsBeeMuSe)[colnames(eigenvec_SNPsBeeMuSe) == "V2"] <- "name"
eigenvec_SNPsBeeMuSe_seq_api_labels <- merge(eigenvec_SNPsBeeMuSe, seq_api_labels, by = "name")

eigen_percent_SNPsBeeMuSe <- round((eigenval_SNPsBeeMuSe / (sum(eigenval_SNPsBeeMuSe) )*100),2)

#Clustering hiérarchique
#Tree
matrice_app_refpop <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_acp.rel", header = FALSE)

dist_matrice_refpop <- dist(matrice_app_refpop)
hc_refpop <- hclust(dist_matrice_refpop, method = "ward.D2")
plot(hc_refpop)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')

#ACP
#filter 629 -> 301 RefPop
eigenvec_SNPsBeeMuSe_seq_api_labels <- eigenvec_SNPsBeeMuSe_seq_api_labels[eigenvec_SNPsBeeMuSe_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label2 <- c("black", "mediumpurple4", "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_SNPsBeeMuSe$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)

ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.035, xmax = 0.065, ymin = -0.04, ymax = -0.11)

#ellipses autour des points selon Label
ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 

# ellipses avec couleur
ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

##### PC3/PC4

```{r}
#PC3/PC4
ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_y_reverse() +  # Reverse the y-axis
  theme(legend.position = c(0.53, 0.7), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = -0.046, xmax = -0.024, ymin = 0.08, ymax = 0.02)
```

#### MAF \> 0.01 - LD pruning = 0.3 (fenêtre de 1749 SNPs et pas de 175 bp) - 3848 SNPs

```{r}
  # SNPsBeeMuSe filtered maf001 LD03
setwd("~/Documents/Stage_NB/data/SeqApiPop_629_SNPsBeeMuSe")

eigenvec_SNPsBeeMuSe <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_acp.eigenvec", header = F)
eigenval_SNPsBeeMuSe <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_SNPsBeeMuSe)[colnames(eigenvec_SNPsBeeMuSe) == "V2"] <- "name"
eigenvec_SNPsBeeMuSe_seq_api_labels <- merge(eigenvec_SNPsBeeMuSe, seq_api_labels, by = "name")

eigen_percent_SNPsBeeMuSe <- round((eigenval_SNPsBeeMuSe / (sum(eigenval_SNPsBeeMuSe) )*100),2)

#Clustering hiérarchique
#Tree
matrice_app_refpop <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_acp.rel", header = FALSE)

dist_matrice_refpop <- dist(matrice_app_refpop)
hc_refpop <- hclust(dist_matrice_refpop, method = "ward.D2")
plot(hc_refpop)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')

#ACP
#filter 629 -> 301 RefPop
eigenvec_SNPsBeeMuSe_seq_api_labels <- eigenvec_SNPsBeeMuSe_seq_api_labels[eigenvec_SNPsBeeMuSe_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label2 <- c("black", "mediumpurple4", "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_SNPsBeeMuSe$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)
```

##### PC1/PC2

```{r}
ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.035, xmax = 0.065, ymin = -0.04, ymax = -0.11)

ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))

#ellipses autour des points selon Label
ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 

# ellipses avec couleur
ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

##### PC3/PC4

```{r}
#PC3/PC4
ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.5, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.025, xmax = 0.05, ymin = 0.02, ymax = 0.11)

```

#### MAF \> 0.01 - LD pruning = 0.1 (fenêtre de 50 SNPs et pas de 10 bp) - 1055 SNPs

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_629_SNPsBeeMuSe")

eigenvec_SNPsBeeMuSe <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_default_acp.eigenvec", header = F)
eigenval_SNPsBeeMuSe <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_default_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_SNPsBeeMuSe)[colnames(eigenvec_SNPsBeeMuSe) == "V2"] <- "name"
eigenvec_SNPsBeeMuSe_seq_api_labels <- merge(eigenvec_SNPsBeeMuSe, seq_api_labels, by = "name")

eigen_percent_SNPsBeeMuSe <- round((eigenval_SNPsBeeMuSe / (sum(eigenval_SNPsBeeMuSe) )*100),2)

#Clustering hiérarchique
#Tree
matrice_app_refpop <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_default_acp.rel", header = FALSE)

dist_matrice_refpop <- dist(matrice_app_refpop)
hc_refpop <- hclust(dist_matrice_refpop, method = "ward.D2")
plot(hc_refpop)

#ACP
#filter 629 -> 301 RefPop
eigenvec_SNPsBeeMuSe_seq_api_labels <- eigenvec_SNPsBeeMuSe_seq_api_labels[eigenvec_SNPsBeeMuSe_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                 eigenvec_SNPsBeeMuSe_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label2 <- c("black", "mediumpurple4", "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_SNPsBeeMuSe$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)
```

##### PC1/PC2

```{r}
ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.035, xmax = 0.065, ymin = -0.04, ymax = -0.11)

ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.34, 0.02), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))

#ellipses autour des points selon Label
ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.5, 0.04), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 

# ellipses avec couleur
ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label2, 
                    breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  theme(legend.position = c(0.5, 0.04), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1, "lines"),  
        legend.text = element_text(size = 10))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

##### PC3/PC4

```{r}
#PC3/PC4
ggplot(data = eigenvec_SNPsBeeMuSe_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label2, 
                     breaks = c("Ouessant Conservatory", "Colonsay Conservatory", "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Mellifera Ouessant", "Mellifera Colonsay", "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.5, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.025, xmax = 0.05, ymin = 0.02, ymax = 0.11)
```

### 561 échantillons - SNPsBeeMuSe filtered

#### MAF \> 0.01 - LD pruning = 0.3 (fenêtre de 1749 SNPs et pas de 175 bp) - 3848 SNPs

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_561_SNPsBeeMuse_LD03")

eigenvec_561_SNPsBeeMuSe <- read.table("SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD03_acp.eigenvec", header = F)
eigenval_561_SNPsBeeMuSe <- read.table("SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD03_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_561_SNPsBeeMuSe)[colnames(eigenvec_561_SNPsBeeMuSe) == "V2"] <- "name"
eigenvec_561_SNPsBeeMuSe_seq_api_labels <- merge(eigenvec_561_SNPsBeeMuSe, seq_api_labels, by = "name")

eigen_percent_561_SNPsBeeMuSe <- round((eigenval_561_SNPsBeeMuSe / (sum(eigenval_561_SNPsBeeMuSe) )*100),2)
```

```{r}
#Clustering hiérarchique
#Tree
  setwd("~/Documents/Stage_NB/data/SeqApiPop_561_SNPsBeeMuse_LD03")
matrice_app_561_default <- read.table("SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD03_acp.rel", header = FALSE)

dist_matrice_561_default <- dist(matrice_app_561_default)
hc_561_d <- hclust(dist_matrice_561_default, method = "ward.D2")
plot(hc_561_d)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')
```

```{r}
#ACP
eigenvec_561_SNPsBeeMuSe_seq_api_labels <- eigenvec_561_SNPsBeeMuSe_seq_api_labels[eigenvec_561_SNPsBeeMuSe_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                   eigenvec_561_SNPsBeeMuSe_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                   eigenvec_561_SNPsBeeMuSe_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                   eigenvec_561_SNPsBeeMuSe_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                   eigenvec_561_SNPsBeeMuSe_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                   eigenvec_561_SNPsBeeMuSe_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label <- c("mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_561_SNPsBeeMuSe$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)
```

##### PC1/PC2

```{r}
#ACP - plot variance
ggplot(data = eigenvec_561_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.45, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.085, xmax = 0.045, ymin = -0.12, ymax = -0.05)

```

##### PC3/PC4

```{r}
#ACP - plot variance
ggplot(data = eigenvec_561_SNPsBeeMuSe_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.48, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.025, xmax = 0.05, ymin = -0.04, ymax = -0.12)
```

#### MAF \> 0.01 - LD pruning = 0.1 (fenêtre de 50 SNPs et pas de 10 bp) - 1055 SNPs

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_561_SNPsBeeMuSe_LD_default")

eigenvec_561_SNPsBeeMuSe <- read.table("SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD03_default_pruned_acp.eigenvec", header = F)
eigenval_561_SNPsBeeMuSe <- read.table("SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD03_default_pruned_acp.eigenval", header = F)

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

colnames(eigenvec_561_SNPsBeeMuSe)[colnames(eigenvec_561_SNPsBeeMuSe) == "V2"] <- "name"
eigenvec_561_SNPsBeeMuSe_seq_api_labels <- merge(eigenvec_561_SNPsBeeMuSe, seq_api_labels, by = "name")

eigen_percent_561_SNPsBeeMuSe <- round((eigenval_561_SNPsBeeMuSe / (sum(eigenval_561_SNPsBeeMuSe) )*100),2)
```

```{r}
#Clustering hiérarchique
#Tree
  setwd("~/Documents/Stage_NB/data/SeqApiPop_561_SNPsBeeMuSe_LD_default")
matrice_app_561_default <- read.table("SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD03_default_pruned_acp.rel", header = FALSE)

dist_matrice_561_default <- dist(matrice_app_561_default)
hc_561_d <- hclust(dist_matrice_561_default, method = "ward.D2")
plot(hc_561_d)

#heatmap
#heatplot(as.matrix(dist(matrice_app_refpop,diag=T)), cols.default = F, lowcol = 'blue', highcol='yellow', dualScale = F, scale='none', method='ward.D2')
```

```{r}
#ACP
eigenvec_561_SNPsBeeMuSe_seq_api_labels <- eigenvec_561_SNPsBeeMuSe_seq_api_labels[eigenvec_561_SNPsBeeMuSe_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                   eigenvec_561_SNPsBeeMuSe_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                   eigenvec_561_SNPsBeeMuSe_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                   eigenvec_561_SNPsBeeMuSe_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                   eigenvec_561_SNPsBeeMuSe_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                   eigenvec_561_SNPsBeeMuSe_seq_api_labels$GeneticOrigin != 'Buckfast', ]

custom_colors_label <- c("mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen")

lambda <- eigenval_561_SNPsBeeMuSe$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)
```

##### PC1/PC2

```{r}
# ACP avec variance expliquée
ggplot(data = eigenvec_561_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.45, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.085, xmax = 0.045, ymin = -0.12, ymax = -0.05)

```

```{r}
#ellipses autour des points selon Label
ggplot(data = eigenvec_561_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  geom_text(aes(x = 0.05, y = 0.04, label = "M-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.01, y = 0.05, label = "C-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.01, y = -0.15, label = "O-lineage"), size = 4, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c( "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.45, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 
```

```{r}
# ellipses avec couleur
ggplot(data = eigenvec_561_SNPsBeeMuSe_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  geom_text(aes(x = 0.045, y = 0.04, label = "M-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.01, y = 0.05, label = "C-lineage"), size = 4, color = "black") +
  geom_text(aes(x = -0.01, y = -0.15, label = "O-lineage"), size = 4, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC1", y = "PC2") +
  scale_color_manual(values = custom_colors_label, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label, 
                    breaks = c( "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c( "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  theme(legend.position = c(0.45, 0.03), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

##### PC3/PC4

```{r}
ggplot(data = eigenvec_561_SNPsBeeMuSe_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.45, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.03, xmax = 0.055, ymin = -0.05, ymax = -0.15)
```

```{r}
#ellipses autour des points selon Label
ggplot(data = eigenvec_561_SNPsBeeMuSe_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label), geom = "polygon", level = 0.97, alpha = 0, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c( "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  theme(legend.position = c(0.45, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2)) 
```

```{r}
# ellipses avec couleur
ggplot(data = eigenvec_561_SNPsBeeMuSe_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point() +
  stat_ellipse(aes(group = Label, fill = Label), geom = "polygon", level = 0.97, alpha = 0.2, size = 0.2, color = "black") +
  labs(title = "PCA Plot - reference populations", x = "PC3", y = "PC4") +
  scale_color_manual(values = custom_colors_label, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +
  scale_fill_manual(values = custom_colors_label, 
                    breaks = c( "Iberiensis Spain", 
                               "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                    labels = c( "Iberiensis Spain", 
                               "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                               "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                               "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) +  
  theme(legend.position = c(0.45, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(1.2, "lines"),  
        legend.text = element_text(size = 11))  + 
  guides(color = guide_legend(override.aes = list(size = 3.5), ncol = 2))
```

## Merged Data

### 748 échantillons BeeMuSe - 561 échantillons SeqApiPop

#### MAF \> 0.01 - LD pruning = 0.1 (fenêtre de 50 SNPs et pas de 10 bp) - 1055 SNPs

```{r}
setwd("~/Documents/Stage_NB/data/merged_data_1055_561_not_supervised")

eigenvec_merged_maf001_LD03 <- read.table("merged_BeeMuSe_SeqApiPop_561_filtered_maf001_LD_default_acp.eigenvec", header = F) 
eigenval_merged_maf001_LD03 <- read.table("merged_BeeMuSe_SeqApiPop_561_filtered_maf001_LD_default_acp.eigenval", header = F) 

colnames(eigenvec_merged_maf001_LD03)[colnames(eigenvec_merged_maf001_LD03) == "V1"] <- "Sample"
colnames(eigenvec_merged_maf001_LD03)[colnames(eigenvec_merged_maf001_LD03) == "V2"] <- "name"

Corres_ID_E756 <- read.table("/home/nbettembourg/Documents/Stage_NB/data/BeeMuSe_Wageningen_E756_17_03_2023_Alain/BeeMuSe_Wageningen_E756_17_03_2023_Alain/Corres_ID_E756.txt",header=F)
input_pedigree_BeeMuSe <- read.table("/home/nbettembourg/Documents/Stage_NB/data/BeeMuSe_Wageningen_E756_17_03_2023_Alain/BeeMuSe_Wageningen_E756_17_03_2023_Alain/input-pedigree_BeeMuSe.csv",header=T,sep=";")

seq_api_labels <- read.csv("~/Documents/Stage_NB/data/SeqApiPop_labels.csv")

eigenvec_merged_maf001_LD03_seq_api_labels <- merge(eigenvec_merged_maf001_LD03, seq_api_labels, by = "name")

eigen_percent <- round((eigenval_merged_maf001_LD03 / (sum(eigenval_merged_maf001_LD03) )*100),2)
```

```{r}
lambda <- eigenval_merged_maf001_LD03$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)
```

```{r}
# filter 629 -> 301 RefPop Samples
eigenvec_merged_maf001_LD03_seq_api_labels <- eigenvec_merged_maf001_LD03_seq_api_labels[eigenvec_merged_maf001_LD03_seq_api_labels$GeneticOrigin != 'Unknown' &
                                                                   eigenvec_merged_maf001_LD03_seq_api_labels$Label != 'Ariege Conservatory' &
                                                                   eigenvec_merged_maf001_LD03_seq_api_labels$Label != 'Brittany Conservatory' &
                                                                   eigenvec_merged_maf001_LD03_seq_api_labels$UniqueInHive != 'Unknown' &
                                                                   eigenvec_merged_maf001_LD03_seq_api_labels$UniqueInHive != 'Buckfast' &
                                                                   eigenvec_merged_maf001_LD03_seq_api_labels$GeneticOrigin != 'Buckfast', ]

# Filtrer les lignes où les facteurs des colonnes Sample et name sont différents
eigenvec_BeeMuSe_748_Samples <- eigenvec_merged_maf001_LD03[as.character(eigenvec_merged_maf001_LD03$Sample) != as.character(eigenvec_merged_maf001_LD03$name), ]
```

```{r}
#extraire 'Pool' et '-100' obtenir 'Pool-100'
eigenvec_merged_maf001_LD03$name <- paste(sub("Beemuse_", "", eigenvec_merged_maf001_LD03$Sample), 
                 sub("_(.*?)\\..*", "\\1", eigenvec_merged_maf001_LD03$name), 
                 sep = "-")
eigenvec_merged_maf001_LD03$name <- str_extract(eigenvec_merged_maf001_LD03$name, "[A-Za-z0-9]+-[0-9]+")

colnames(Corres_ID_E756)[colnames(Corres_ID_E756) == "V1"] <- "name"
colnames(Corres_ID_E756)[colnames(Corres_ID_E756) == "V2"] <- "ID_1a"

Corres_ID_E756$ID_1a <- gsub("o", "_", Corres_ID_E756$ID_1a)

# Remplacement de "o" par "_" sauf quand présent dans un mot
Corres_ID_E756$ID_1a <- gsub("Pers_", "Perso", Corres_ID_E756$ID_1a)
Corres_ID_E756$ID_1a <- gsub("L_c", "Loc", Corres_ID_E756$ID_1a)

Corres_ID_E756_eigenvec <- merge(eigenvec_merged_maf001_LD03, Corres_ID_E756,  by = 'name')
eigenvec_merged_Corres_ID_E756_pedigree <- merge(Corres_ID_E756_eigenvec, input_pedigree_BeeMuSe,  by = 'ID_1a')
colnames(eigenvec_merged_Corres_ID_E756_pedigree)[colnames(eigenvec_merged_Corres_ID_E756_pedigree) == "V3.x"] <- "V3"
```

##### PC1/PC2

```{r}
# ACP - merged data
lambda <- eigenval_merged_maf001_LD03$V1
variance_proportion <- lambda / sum(lambda)
variance_df <- data.frame(PC = seq_along(variance_proportion), Variance = variance_proportion)

ggplot(eigenvec_merged_maf001_LD03, aes(x = V3, y = V4)) +
  geom_point() +
  scale_x_reverse() + 
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC1", y = "PC2") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.065, xmax = 0.04, ymin = -0.1, ymax = -0.05)
```

```{r}
# ACP plot - BeeMuSe 748 & SeqApiPop 301 RefPop
# ACP SeqApiPop
plot1 <- ggplot(eigenvec_merged_maf001_LD03_seq_api_labels, aes(x = V3, y = V4, color = "SeqApiPop")) +
  geom_point() +
  scale_x_reverse() 

# Ajout ACP pour BeeMuSe
plot2 <- plot1 + geom_point(data = eigenvec_BeeMuSe_748_Samples, aes(x = V3, y = V4, color = "BeeMuSe")) +
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop",x = "PC1", y = "PC2")

# Ajout de la légende
plot2 + scale_color_manual(name = "Dataset", 
                            values = c(SeqApiPop = "black", BeeMuSe = "red"),
                            labels = c(SeqApiPop = "SeqApiPop 301 RefPop", BeeMuSe = "BeeMuSe 748 Samples")) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.075, xmax = 0.04, ymin = -0.1, ymax = -0.05)

# ACP plot - BeeMuSe 612 Pedigree & SeqApiPop 301 RefPop
# ACP SeqApiPop
plot1 <- ggplot(eigenvec_merged_maf001_LD03_seq_api_labels, aes(x = V3, y = V4, color = "SeqApiPop")) +
  geom_point() +
  scale_x_reverse()

# Ajout ACP pour BeeMuSe
plot2 <- plot1 + geom_point(data = eigenvec_merged_Corres_ID_E756_pedigree, aes(x = V3, y = V4, color = "BeeMuSe")) +
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop",x = "PC1", y = "PC2")

# Ajout de la légende
plot2 + scale_color_manual(name = "Dataset", 
                            values = c(SeqApiPop = "black", BeeMuSe = "red"),
                            labels = c(SeqApiPop = "SeqApiPop 301 RefPop", BeeMuSe = "BeeMuSe 612 Pedigree")) +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.075, xmax = 0.04, ymin = -0.1, ymax = -0.05)

# Couleurs 3 lignées - Label
custom_colors_2 <- c( "black", "black", "black", "black", "orange", "orange", "orange", "orange", "orange", "orange", "mediumseagreen", "red")

# ACP SeqApiPop
plot1 <- ggplot(eigenvec_merged_maf001_LD03_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point(size = 1) +
  scale_x_reverse() +
  #scale_y_reverse() +
  scale_color_manual(values = custom_colors_2, 
                     breaks = c("Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) 

# Ajout ACP pour BeeMuSe
plot1 + geom_point(data = eigenvec_BeeMuSe_748_Samples, aes(x = V3, y = V4), color = "red", size =1) +
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC1", y = "PC2") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.075, xmax = 0.04, ymin = -0.1, ymax = -0.05)

# Couleurs des 11 populations de référence - Label
custom_colors <- c( "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen", "red")

# ACP SeqApiPop
plot1 <- ggplot(eigenvec_merged_maf001_LD03_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point(size=1) +
  scale_x_reverse() +
  scale_color_manual(values = custom_colors, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c( "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) 

# Ajout ACP pour BeeMuSe
plot1 + geom_point(data = eigenvec_BeeMuSe_748_Samples, aes(x = V3, y = V4), color = "red", size=1) +
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC1", y = "PC2") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.075, xmax = 0.04, ymin = -0.1, ymax = -0.05)
```

```{r}
# SeqApiPop - 301 Samples - BeeMuSe 612 Samples - ID_2a - colors
custom_colors <- c( "black", "black", "black", 
                   "black", "black", "black", "black", "black", "black", "black", 
                   "black", 
                   "#FF0000FF", "#FF3300FF", "#FF6600FF", "#FF9900FF", "#FFCC00FF", 
                   "#FFFF00FF", "#CCFF00FF", "#99FF00FF", "#66FF00FF", "#33FF00FF", 
                   "#00FF00FF", "#00FF33FF", "#00FF66FF", "#00FF99FF", "#00FFCCFF", 
                   "#00FFFFFF", "#00CCFFFF", "#0099FFFF", "#0066FFFF", "#0033FFFF", 
                   "#0000FFFF", "#3300FFFF", "#6600FFFF", "#9900FFFF", "#CC00FFFF", 
                   "#FF00FFFF", "#FF00CCFF", "#FF0099FF", "#FF0066FF")

# ACP SeqApiPop
plot1 <- ggplot(eigenvec_merged_maf001_LD03_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point(size=1) +
  scale_x_reverse() +
  scale_color_manual(values = custom_colors, 
                     breaks = c( "Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France", "BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20", "BLS_53-19", "ER_13-19", "KBJ_1-19", "KBru_6-20", "KLoc_37-19", "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", 
"PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", "PersoUB_2021", 
"PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", "SJ_24-20", "SJ_30-20", "TL_13-20", 
"TL_19-20"),
                     labels = c( "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia","BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20", "BLS_53-19", "ER_13-19", 
 "KBJ_1-19","KBru_6-20", "KLoc_37-19", "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", 
"PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", "PersoUB_2021", 
"PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", "SJ_24-20", "SJ_30-20", "TL_13-20", 
"TL_19-20")) +  
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC1", y = "PC2") +
  theme(legend.position = "right")

# Ajout ACP BeeMuSe
plot1 + geom_point(data = eigenvec_merged_Corres_ID_E756_pedigree, aes(x = V3, y = V4, color = ID_2a),size=1) +
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC1", y = "PC2") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.075, xmax = 0.04, ymin = -0.1, ymax = -0.05)

# SeqApiPop - 301 Samples - Label - colors - 3 lignées
custom_colors <- c( "black", "black", "black", 
                   "black", "orange", "orange", "orange", "orange", "orange", "orange", 
                   "mediumseagreen", 
                   "#FF0000FF", "#FF3300FF", "#FF6600FF", "#FF9900FF", "#FFCC00FF", 
                   "#FFFF00FF", "#CCFF00FF", "#99FF00FF", "#66FF00FF", "#33FF00FF", 
                   "#00FF00FF", "#00FF33FF", "#00FF66FF", "#00FF99FF", "#00FFCCFF", 
                   "#00FFFFFF", "#00CCFFFF", "#0099FFFF", "#0066FFFF", "#0033FFFF", 
                   "#0000FFFF", "#3300FFFF", "#6600FFFF", "#9900FFFF", "#CC00FFFF", 
                   "#FF00FFFF", "#FF00CCFF", "#FF0099FF", "#FF0066FF")

# ACP SeqApiPop
plot1 <- ggplot(eigenvec_merged_maf001_LD03_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point(size=1) +
  scale_x_reverse() +
  scale_color_manual(values = custom_colors, 
                     breaks = c("Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France", "BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20", "BLS_53-19", "ER_13-19", "KBJ_1-19", "KBru_6-20", "KLoc_37-19", "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", 
"PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", "PersoUB_2021", 
"PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", "SJ_24-20", "SJ_30-20", "TL_13-20", 
"TL_19-20"),
                     labels = c( "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia","BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20", "BLS_53-19", "ER_13-19", 
 "KBJ_1-19","KBru_6-20", "KLoc_37-19", "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", 
"PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", "PersoUB_2021", 
"PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", "SJ_24-20", "SJ_30-20", "TL_13-20", 
"TL_19-20")) +  
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC1", y = "PC2") +
  theme(legend.position = "right")

# Ajout ACP BeeMuSe
plot1 + geom_point(data = eigenvec_merged_Corres_ID_E756_pedigree, aes(x = V3, y = V4, color = ID_2a),size=1) +
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC1", y = "PC2") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.075, xmax = 0.04, ymin = -0.1, ymax = -0.05)

# SeqApiPop - 301 - Label - couleurs - 11 populations de référence
custom_colors <- c("mediumpurple3", "lightskyblue", "mediumvioletred", 
                   "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", 
                   "mediumseagreen", 
                   "#FF0000FF", "#FF3300FF", "#FF6600FF", "#FF9900FF", "#FFCC00FF", 
                   "#FFFF00FF", "#CCFF00FF", "#99FF00FF", "#66FF00FF", "#33FF00FF", 
                   "#00FF00FF", "#00FF33FF", "#00FF66FF", "#00FF99FF", "#00FFCCFF", 
                   "#00FFFFFF", "#00CCFFFF", "#0099FFFF", "#0066FFFF", "#0033FFFF", 
                   "#0000FFFF", "#3300FFFF", "#6600FFFF", "#9900FFFF", "#CC00FFFF", 
                   "#FF00FFFF", "#FF00CCFF", "#FF0099FF", "#FF0066FF")

# ACP SeqApiPop
plot1 <- ggplot(eigenvec_merged_maf001_LD03_seq_api_labels, aes(x = V3, y = V4, color = Label)) +
  geom_point(size = 1) +  
  scale_x_reverse() +
  scale_color_manual(values = custom_colors, 
                     breaks = c("Iberiensis Spain", "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France", "BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20", "BLS_53-19", "ER_13-19", "KBJ_1-19", "KBru_6-20", "KLoc_37-19", "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", 
"PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", "PersoUB_2021", 
"PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", "SJ_24-20", "SJ_30-20", "TL_13-20", 
"TL_19-20"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia","BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20", "BLS_53-19", "ER_13-19", 
 "KBJ_1-19","KBru_6-20", "KLoc_37-19", "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", 
"PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", "PersoUB_2021", 
"PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", "SJ_24-20", "SJ_30-20", "TL_13-20", 
"TL_19-20")) +  
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC1", y = "PC2") +
  theme(legend.position = "right")

# Ajout ACP BeeMuSe
plot1 + geom_point(data = eigenvec_merged_Corres_ID_E756_pedigree, aes(x = V3, y = V4, color = ID_2a), size = 1) + 
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC1", y = "PC2") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.075, xmax = 0.04, ymin = -0.1, ymax = -0.05)
```

##### PC3/PC4

```{r}
# ACP - PC3 / PC4

# ACP SeqApiPop
plot1 <- ggplot(eigenvec_merged_maf001_LD03_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point(size = 1) + 
  scale_x_reverse() +
  scale_color_manual(values = custom_colors_2, 
                     breaks = c("Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c( "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) 

# Ajout ACP BeeMuSe
plot1 + geom_point(data = eigenvec_BeeMuSe_748_Samples, aes(x = V5, y = V6), color = "red", size = 1) + 
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC3", y = "PC4") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.097, xmax = 0.034, ymin = -0.1, ymax = -0.05)

# Couleurs 11 populations de référence
custom_colors <- c( "mediumpurple3", "lightskyblue", "mediumvioletred", "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", "mediumseagreen", "red")

# ACP SeqApiPop
plot1 <- ggplot(eigenvec_merged_maf001_LD03_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point(size = 1) +
  scale_x_reverse() +
  scale_color_manual(values = custom_colors, 
                     breaks = c("Iberiensis Spain", 
                                "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France"),
                     labels = c( "Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia")) 

# Ajout ACP BeeMuSe
plot1 + geom_point(data = eigenvec_BeeMuSe_748_Samples, aes(x = V5, y = V6), color = "red", size = 1) +
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC3", y = "PC4") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.097, xmax = 0.034, ymin = -0.1, ymax = -0.05)

# Couleurs des 11 populations de référence SeqApiPop + 29 familles ID_2a 
custom_colors <- c("mediumpurple3", "lightskyblue", "mediumvioletred", 
                   "hotpink1", "yellow", "gold", "orange", "chocolate", "brown", "olivedrab3", 
                   "mediumseagreen", 
                   "#FF0000FF", "#FF3300FF", "#FF6600FF", "#FF9900FF", "#FFCC00FF", 
                   "#FFFF00FF", "#CCFF00FF", "#99FF00FF", "#66FF00FF", "#33FF00FF", 
                   "#00FF00FF", "#00FF33FF", "#00FF66FF", "#00FF99FF", "#00FFCCFF", 
                   "#00FFFFFF", "#00CCFFFF", "#0099FFFF", "#0066FFFF", "#0033FFFF", 
                   "#0000FFFF", "#3300FFFF", "#6600FFFF", "#9900FFFF", "#CC00FFFF", 
                   "#FF00FFFF", "#FF00CCFF", "#FF0099FF", "#FF0066FF")

# ACP SeqApiPop
plot1 <- ggplot(eigenvec_merged_maf001_LD03_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point(size = 1) +  
  scale_x_reverse() +
  scale_color_manual(values = custom_colors, 
                     breaks = c("Iberiensis Spain", "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France", "BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20", "BLS_53-19", "ER_13-19", "KBJ_1-19", "KBru_6-20", "KLoc_37-19", "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", 
"PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", "PersoUB_2021", 
"PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", "SJ_24-20", "SJ_30-20", "TL_13-20", 
"TL_19-20"),
                     labels = c("Iberiensis Spain", 
                                "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia","BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20", "BLS_53-19", "ER_13-19", 
 "KBJ_1-19","KBru_6-20", "KLoc_37-19", "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", 
"PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", "PersoUB_2021", 
"PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", "SJ_24-20", "SJ_30-20", "TL_13-20", 
"TL_19-20")) +  
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC3", y = "PC4") +
  theme(legend.position = "right")

# Ajouter ACP BeeMuSe
plot1 + geom_point(data = eigenvec_merged_Corres_ID_E756_pedigree, aes(x = V5, y = V6, color = ID_2a), size = 1) + 
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC3", y = "PC4") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.097, xmax = 0.034, ymin = -0.1, ymax = -0.05)

# Couleurs 3 lignées SeqApiPop + 29 familles ID_2a 
custom_colors <- c( "black", "black", "black", 
                   "black", "orange", "orange", "orange", "orange", "orange", "orange", 
                   "mediumseagreen", 
                   "#FF0000FF", "#FF3300FF", "#FF6600FF", "#FF9900FF", "#FFCC00FF", 
                   "#FFFF00FF", "#CCFF00FF", "#99FF00FF", "#66FF00FF", "#33FF00FF", 
                   "#00FF00FF", "#00FF33FF", "#00FF66FF", "#00FF99FF", "#00FFCCFF", 
                   "#00FFFFFF", "#00CCFFFF", "#0099FFFF", "#0066FFFF", "#0033FFFF", 
                   "#0000FFFF", "#3300FFFF", "#6600FFFF", "#9900FFFF", "#CC00FFFF", 
                   "#FF00FFFF", "#FF00CCFF", "#FF0099FF", "#FF0066FF")

# ACP SeqApiPop
plot1 <- ggplot(eigenvec_merged_maf001_LD03_seq_api_labels, aes(x = V5, y = V6, color = Label)) +
  geom_point(size = 1) +
  scale_x_reverse() +
  scale_color_manual(values = custom_colors, 
                     breaks = c( "Iberiensis Spain", "Savoy Conservatory", "Porquerolles Conservatory", "Sollies Conservatory", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia France", "BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20", "BLS_53-19", "ER_13-19", "KBJ_1-19", "KBru_6-20", "KLoc_37-19", "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", 
"PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", "PersoUB_2021", 
"PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", "SJ_24-20", "SJ_30-20", "TL_13-20", 
"TL_19-20"),
                     labels = c( "Iberiensis Spain", "Mellifera Savoy", "Mellifera Porquerolles", "Mellifera Sollies", 
                                "Ligustica Italy", "Carnica Slovenia", "Carnica Germany", 
                                "Carnica Switzerland", "Carnica France", "Carnica Poland", "Caucasia","BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20", "BLS_53-19", "ER_13-19", 
 "KBJ_1-19","KBru_6-20", "KLoc_37-19", "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", 
"PersoBC_2021", "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", "PersoUB_2021", 
"PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", "SJ_24-20", "SJ_30-20", "TL_13-20", 
"TL_19-20"))

# Ajout ACP BeeMuSe
plot1 + geom_point(data = eigenvec_merged_Corres_ID_E756_pedigree, aes(x = V5, y = V6, color = ID_2a), size = 1 ) +
  labs(title = "PCA Plot - BeeMuSe & SeqApiPop", x = "PC3", y = "PC4") +
  annotation_custom(
    ggplotGrob(
      ggplot(variance_df, aes(x = factor(PC), y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ), xmin = 0.097, xmax = 0.034, ymin = -0.1, ymax = -0.05)
```

## Analyse d'ADMIXTURE - CV error

### Admixture non supervisée

### CV error plot - SeqApiPop - 629 échantillons - MAF \> 0.01

#### CV error plot - LD pruning = 0.3 (fenêtre de 1749 SNPs et pas de 175 bp)

```{r}
 #LD03 - CV Error plot
setwd("~/Documents/Stage_NB/data/maf001_LD03")

cv_error <- read.table("SeqApiPop_629_maf001_LD03_1.cv.error", header = F)

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_629_maf001_LD03_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

  # 26/1
merge_cv_error <- read.table("merge_cv_error", header = F)

point_min <- merge_cv_error[which.min(merge_cv_error[, 2]), ]

#box plot
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.475, 0.480, 0.485, 0.490, 0.495, 0.500, 0.505),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.475, 0.505))

#jitter plot
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.475, 0.480, 0.485, 0.490, 0.495, 0.500, 0.505),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.475, 0.505))
```

#### CV error plot - LD pruning = 0.2

```{r}
#####LD02
setwd("~/Documents/Stage_NB/data/maf001_LD02")

cv_error <- read.table("SeqApiPop_629_maf001_LD02_1.cv.error", header = F)

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:10) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_629_maf001_LD02_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

point_min <- merge_cv_error[which.min(merge_cv_error[, 2]), ]

#box plot LD02
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.430, 0.435, 0.440, 0.445, 0.450),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.430, 0.450))

#jitter plot LD02
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.430, 0.435, 0.440, 0.445, 0.450),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.430, 0.450))
```

#### CV error plot - LD pruning = 0.1

```{r}
#####LD01
setwd("~/Documents/Stage_NB/data/maf001_LD01")

cv_error <- read.table("SeqApiPop_629_maf001_LD01_1.cv.error", header = F)

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:10) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_629_maf001_LD01_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

point_min <- merge_cv_error[which.min(merge_cv_error[, 2]), ]

#box plot LD01
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.370, 0.375, 0.380, 0.385, 0.390),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.370, 0.390))

#jitter plot LD01
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.370, 0.375, 0.380, 0.385, 0.390),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.370, 0.390))
```

#### CV error plot - LD pruning = 0.05

```{r}
#####LD005
setwd("~/Documents/Stage_NB/data/maf001_LD005")

cv_error <- read.table("SeqApiPop_629_maf001_LD005_1.cv.error", header = F)

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:10) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_629_maf001_LD005_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD01
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.460, 0.465, 0.470, 0.475, 0.480),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.460, 0.480))

#jitter plot LD01
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.460, 0.465, 0.470, 0.475, 0.480),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.460, 0.480))
```

#### CV error plot - LD pruning = 0.04

```{r}
#####LD004
setwd("~/Documents/Stage_NB/data/maf001_LD004")

cv_error <- read.table("SeqApiPop_629_maf001_LD004_1.cv.error", header = F)

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:10) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_629_maf001_LD004_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD004
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.575, 0.58, 0.585, 0.59, 0.595, 0.6),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.575, 0.6))

#jitter plot LD004
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.575, 0.58, 0.585, 0.59, 0.595, 0.6),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.575, 0.6))
```

#### CV error plot - LD pruning = 0.03

```{r}
#####LD003
setwd("~/Documents/Stage_NB/data/maf001_LD003")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:10) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_629_maf001_LD003_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD01
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.78, 0.785, 0.79, 0.795, 0.8, 0.805, 0.81, 0.815),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.78, 0.815))

#jitter plot LD01
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.78, 0.785, 0.79, 0.795, 0.8, 0.805, 0.81, 0.815),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.78, 0.815))
```

#### CV error plot - LD pruning = 0.01

```{r}
#####LD001
setwd("~/Documents/Stage_NB/data/maf001_LD001")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:10) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_629_maf001_LD001_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD01
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.895, 0.9, 0.905, 0.91, 0.915, 0.92, 0.925, 0.93, 0.935, 0.94, 0.945),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.895, 0.945))

#jitter plot LD01
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.895, 0.9, 0.905, 0.91, 0.915, 0.92, 0.925, 0.93, 0.935, 0.94, 0.945),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.895, 0.945))
```

### CV error plot - SeqApiPop - 561 échantillons - MAF \> 0.01

#### CV error plot - LD pruning = 0.3 (fenêtre de 1749 SNPs et pas de 175 bp)

```{r}
#LD03 - CV Error plot
setwd("~/Documents/Stage_NB/data/SeqApiPop_561_maf001_LD03")

cv_error <- read.table("SeqApiPop_561_maf001_LD03_1.cv.error", header = F)

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_561_maf001_LD03_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

  # 26/1
merge_cv_error <- read.table("merge_cv_error", header = F)

point_min <- merge_cv_error[which.min(merge_cv_error[, 2]), ]

#box plot
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.46, 0.465, 0.47,0.475, 0.480, 0.485, 0.490),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.46, 0.49))

#jitter plot
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.46, 0.465, 0.47,0.475, 0.480, 0.485, 0.490),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.46, 0.49))
```

### 

### CV error plot - SeqApiPop - 629 échantillons - SNPsBeeMuSe filtered

#### CV error plot - filtered

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_629_SNPsBeeMuSe")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_629_SNPsBeeMuSe_filtered_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.72,0.725,0.73,0.735,0.74,0.745,0.75,0.755,0.76,0.765,0.77,0.775,0.78,0.785,0.79,0.795,0.8),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.72, 0.8))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.72,0.725,0.73,0.735,0.74,0.745,0.75,0.755,0.76,0.765,0.77,0.775,0.78,0.785,0.79,0.795,0.8),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.72, 0.8))
```

#### CV error plot - MAF \> 0.01 - LD pruning = 0.3 (fenêtre de 1749 SNPs et pas de 175 bp)

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_629_SNPsBeeMuSe")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.765, 0.77, 0.775, 0.78, 0.785, 0.79, 0.795, 0.8, 0.805, 0.81, 0.815,0.82,0.825,0.83,0.835,0.84),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.765, 0.84))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.765, 0.77, 0.775, 0.78, 0.785, 0.79, 0.795, 0.8, 0.805, 0.81, 0.815,0.82,0.825,0.83,0.835,0.84),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.765, 0.84))
```

#### CV error plot - MAF \> 0.01 - LD pruning = 0.1 (fenêtre de 50 SNPs et pas de 10 bp)

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_629_SNPsBeeMuSe")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_default_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.755, 0.76, 0.765, 0.77, 0.775, 0.78, 0.785, 0.79, 0.795, 0.8, 0.805, 0.81, 0.815,0.82),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.755, 0.82))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.755, 0.76, 0.765, 0.77, 0.775, 0.78, 0.785, 0.79, 0.795, 0.8, 0.805, 0.81, 0.815,0.82),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.755, 0.82))
```

### CV error plot - SeqApiPop - 561 échantillons - SNPsBeeMuSe filtered

#### CV error plot - MAF \> 0.01 - LD pruning = 0.3 (fenêtre de 1749 SNPS et pas de 175 bp)

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_561_SNPsBeeMuse_LD03")

liste_de_donnees <- list()

for (i in 1:30) {
  merge_cv_error <- paste0('SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD03_pruned_', i, '.cv.error')
  donnees <- read.table(merge_cv_error, header = FALSE)
  liste_de_donnees[[i]] <- donnees
}

donnees_combinees <- do.call(rbind, liste_de_donnees)
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)
merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.775, 0.78, 0.785, 0.79,0.795,0.8,0.805,0.81,0.815,0.82,0.825,0.83,0.835,0.84,0.845,0.85),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.775, 0.85))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.775, 0.78, 0.785, 0.79,0.795,0.8,0.805,0.81,0.815,0.82,0.825,0.83,0.835,0.84,0.845,0.85),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.775, 0.85))
```

#### CV error plot - MAF \> 0.01 - LD pruning = 0.1 (fenêtre de 50 SNPS et pas de 10 bp)

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_561_SNPsBeeMuSe_LD_default")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD03_default_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.765,0.77, 0.775, 0.78, 0.785, 0.79,0.795,0.8,0.805,0.81,0.815,0.82,0.825,0.83),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.765, 0.83))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.765,0.77, 0.775, 0.78, 0.785, 0.79,0.795,0.8,0.805,0.81,0.815,0.82,0.825,0.83),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.765, 0.83))
```

### CV error plot - BeeMuSe 12000 SNPs

```{r}
K_new <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70)
CV_new <- c(0.48943, 0.46307, 0.44338, 0.42781, 0.41962, 0.40411, 0.40055, 0.38682, 0.38044, 0.37802, 0.37102, 0.36619, 0.36130, 0.35863, 0.35431, 0.35162, 0.34960, 0.34742, 0.34640, 0.34555, 0.34281, 0.34187, 0.34055, 0.34081, 0.33861, 0.33930, 0.33562, 0.33588, 0.33543, 0.33518, 0.33262, 0.33176, 0.33216, 0.33019, 0.33013, 0.32908, 0.33166, 0.32875, 0.32961, 0.32977, 0.32796, 0.32878, 0.32828, 0.32772, 0.32533, 0.32954, 0.32862, 0.32522, 0.32771, 0.32770, 0.32558, 0.33150, 0.33390, 0.32462, 0.32938, 0.32893, 0.33254, 0.32727, 0.32883, 0.33015, 0.33254, 0.33260, 0.33633, 0.33744, 0.33756, 0.33277, 0.33579)

# Trouver l'indice de la valeur la plus basse de CV
indice_min_new <- which.min(CV_new)

# Créer le graphique
plot(K_new, CV_new, type="l", col="black", xlab="K", ylab="CV", main="CV error plot - BeeMuSe  3848 SNPs")

# Ajouter la ligne avec le trait hachuré bleu pour la valeur la plus basse
abline(h=CV_new[indice_min_new], col="blue", lty=2)

# Ajouter la droite verticale rouge pour la valeur la plus basse de CV
abline(v=K_new[indice_min_new], col="red", lty=2)

# Ajouter les lignes de grille horizontales à intervalles de 0.01
abline(h=seq(0, 1, by=0.01), col="lightgray")

# Ajouter la légende avec la valeur de K correspondant à la plus basse erreur CV
legend("topright", legend=sprintf("lowest CV error: %.5f (K = %d)", CV_new[indice_min_new], K_new[indice_min_new]), col="blue", lty=2, cex=0.8)
```

### 

### CV error plot - Merged Data - BeeMuSe SeqApiPop

#### CV error plot - MAF \> 0.01 - LD pruning = 0.3 (fenêtre de 1749 SNPS et pas de 175 bp)

```{r}
# Valeurs CV - Admixture non supervisée
K <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48 ,49, 50)
CV <- c(0.66344, 0.62032, 0.60790, 0.60015, 0.58943, 0.58458, 0.58509, 0.57812, 0.57116, 0.56933, 0.56023, 0.56145, 0.55645, 0.55416, 0.55222, 0.54907, 0.54576, 0.54370, 0.54270, 0.54220, 0.54162, 0.54375, 0.54147, 0.53958, 0.54064, 0.53829, 0.54184, 0.53940, 0.53965, 0.54062, 0.53549, 0.54011, 0.53806, 0.53622, 0.53630, 0.53986, 0.54062, 0.53707, 0.54252, 0.53593, 0.53901, 0.54513, 0.54334, 0.54278, 0.54128, 0.54482, 0.54184, 0.55066, 0.55151
)

# Trouver l'indice de la valeur la plus basse de CV
indice_min <- which.min(CV)

# Créer le graphique
plot(K, CV, type="l", col="black", xlab="K", ylab="CV", main="CV error plot - Merged BeeMuSe SeqApiPop 3848 SNPs")

# Ajouter la ligne avec le trait hachuré bleu pour la valeur la plus basse
abline(h=CV[indice_min], col="blue", lty=2)

# Ajouter la droite verticale rouge pour la valeur la plus basse de CV
abline(v=K[indice_min], col="red", lty=2)

# Ajouter les lignes de grille horizontales à intervalles de 0.01
abline(h=seq(0, 1, by=0.01), col="lightgray")

# Ajouter la légende
legend("topright", legend=sprintf("lowest CV error: %.5f (K = %d)", CV[indice_min], K[indice_min]), col="blue", lty=2, cex=0.8)

# Nouvelles données
K <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)
CV <- c(0.66350, 0.62020, 0.60778, 0.59670, 0.58968, 0.58467, 0.58258, 0.57804, 0.57563, 0.56626, 0.56028, 0.55808, 0.55537, 0.55140, 0.55279, 0.54843, 0.54649, 0.54583, 0.54412, 0.54363, 0.54187, 0.54172, 0.54028, 0.53822, 0.54265, 0.53707, 0.53598, 0.53746, 0.54359, 0.53759, 0.53651, 0.53563, 0.53614, 0.53701, 0.53993, 0.53787, 0.53929, 0.53656, 0.53982, 0.54181, 0.54237, 0.53793, 0.54851, 0.54505, 0.54506, 0.54203, 0.54622, 0.55210, 0.55177)

# Trouver l'indice de la valeur la plus basse de CV
indice_min <- which.min(CV)

# Créer le graphique
plot(K, CV, type="l", col="black", xlab="K", ylab="CV", main="CV error plot - Merged BeeMuSe SeqApiPop 3848 SNPs")

# Ajouter la ligne avec le trait hachuré bleu pour la valeur la plus basse
abline(h=CV[indice_min], col="blue", lty=2)

# Ajouter la droite verticale rouge pour la valeur la plus basse de CV
abline(v=K[indice_min], col="red", lty=2)

# Ajouter les lignes de grille horizontales à intervalles de 0.01
abline(h=seq(0, 1, by=0.01), col="lightgray")

# Ajouter la légende
legend("topright", legend=sprintf("lowest CV error: %.5f (K = %d)", CV[indice_min], K[indice_min]), col="blue", lty=2, cex=0.8)

```

##### **CV error plot - 629 échantillons - K2 à K9 - 30 exécutions**

```{r}
setwd("~/Documents/Stage_NB/data/merged_BeeMuSe_SeqApiPop_629_filtered_maf001_LD03")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('merged_BeeMuSe_SeqApiPop_629_filtered_maf001_LD03_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.55, 0.67))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.55, 0.67))
```

##### **CV error plot - 561 échantillons - K2 à K9 - 30 exécutions**

```{r}
setwd("~/Documents/Stage_NB/data/merged_data_3848_561_not_supervised")

liste_de_donnees <- list()
for (i in 1:30) {
  merge_cv_error <- paste0('merged_BeeMuSe_SeqApiPop_561_filtered_maf001_LD03_', i, '.cv.error')
  donnees <- read.table(merge_cv_error, header = FALSE)
  liste_de_donnees[[i]] <- donnees
}

donnees_combinees <- do.call(rbind, liste_de_donnees)
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)
merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.54,0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65,0.66,0.67),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.55, 0.67))

#jitter plot filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.54,0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65,0.66,0.67),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.55, 0.67))
```

#### CV error plot - MAF \> 0.01 - LD pruning = 0.1 (fenêtre de 50 SNPS et pas de 10 bp)

##### CV error plot - 561 échantillons - K2 à K9 - 30 exécutions

```{r}
setwd("~/Documents/Stage_NB/data/merged_data_1055_561_not_supervised")

liste_de_donnees <- list()
for (i in 1:30) {
  merge_cv_error <- paste0('merged_BeeMuSe_SeqApiPop_561_filtered_MAF001_LD_default_', i, '.cv.error')
  donnees <- read.table(merge_cv_error, header = FALSE)
  liste_de_donnees[[i]] <- donnees
}

donnees_combinees <- do.call(rbind, liste_de_donnees)
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)
merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.54,0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.54, 0.65))

#jitter plot filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.54,0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.54, 0.65))
```

### Admixture supervisée

### CV error plot - Merge Data - BeeMuSe SeqApiPop

#### CV error plot - 629 échantillons - LD pruning = 0.1 (fenêtre de 50 et pas de 10 bp)

##### **K = 5**

```{r}
setwd("~/Documents/Stage_NB/data/merged_data_3848_629_supervised")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('merged_BeeMuSe_SeqApiPop_629_filtered_MAF001_LD_default_K5_95_supervised_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.616,0.617,0.618,0.619,0.62,0.621,0.622),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.616, 0.622))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.616,0.617,0.618,0.619,0.62,0.621,0.622),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.616, 0.622))
```

##### **K = 6**

```{r}
setwd("~/Documents/Stage_NB/data/merged_data_3848_629_supervised")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('merged_BeeMuSe_SeqApiPop_629_filtered_MAF001_LD_default_K6_95_supervised_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.611,0.612,0.613,0.614,0.615,0.616,0.617),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.611, 0.617))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.611,0.612,0.613,0.614,0.615,0.616,0.617),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.611, 0.617))
```

#### CV error plot - 561 échantillons - LD pruning = 0.3 (fenêtre de 1749 et pas de 175 bp) **K = 5**

```{r}
setwd("~/Documents/Stage_NB/data/merged_data_3848_561_supervised")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('merged_BeeMuSe_SeqApiPop_561_filtered_MAF001_LD_default_K5_95_supervised_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.616,0.617,0.618,0.619,0.62),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.616, 0.62))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.616,0.617,0.618,0.619,0.62),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.616, 0.62))
```

##### K = 6

```{r}
setwd("~/Documents/Stage_NB/data/merged_data_3848_561_supervised")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('merged_BeeMuSe_SeqApiPop_561_filtered_MAF001_LD_default_K6_95_supervised_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.61,0.611,0.612,0.613,0.614,0.615),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.61, 0.615))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.61,0.611,0.612,0.613,0.614,0.615),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.61, 0.615))
```

#### CV error plot - 561 échantillons - LD pruning = 0.1 (fenêtre de 50 et pas de 10 bp)

##### **K = 5**

```{r}
setwd("~/Documents/Stage_NB/data/merged_data_1055_561_supervised")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('merged_BeeMuSe_SeqApiPop_561_filtered_MAF001_LD_default_K5_90_supervised_', i, '.cv.error')
  #merge_cv_error <- paste0('merged_BeeMuSe_SeqApiPop_561_filtered_MAF001_LD_default_K5_95_supervised_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.59, 0.591, 0.592, 0.593,0.594,0.595,0.596,0.597,0.598,0.599),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.59, 0.599))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.59, 0.591, 0.592, 0.593,0.594,0.595,0.596,0.597,0.598,0.599),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.59, 0.599))
```

##### **K = 6**

```{r}
setwd("~/Documents/Stage_NB/data/merged_data_1055_561_supervised")

# Étape 1: Créer une liste vide pour stocker les données
liste_de_donnees <- list()

# Étape 2: Parcourir les fichiers
for (i in 1:30) {
  # Générer le nom du fichier
  merge_cv_error <- paste0('merged_BeeMuSe_SeqApiPop_561_filtered_MAF001_LD_default_K6_90_supervised_', i, '.cv.error')
  #merge_cv_error <- paste0('merged_BeeMuSe_SeqApiPop_561_filtered_MAF001_LD_default_K5_95_supervised_', i, '.cv.error')
  # Lire les données du fichier
  donnees <- read.table(merge_cv_error, header = FALSE)
  # Ajouter les données à la liste
  liste_de_donnees[[i]] <- donnees
}

# Étape 3: Combinez les données en une seule structure (par exemple, data frame)
donnees_combinees <- do.call(rbind, liste_de_donnees)

# Étape 4: Enregistrez ou affichez le résultat final sans numéro de lignes
write.table(donnees_combinees, "merge_cv_error", sep = "\t", col.names = FALSE, row.names = FALSE)

merge_cv_error <- read.table("merge_cv_error", header = F)

#box plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.586,0.587,0.588,0.589,0.59, 0.591, 0.592, 0.593,0.594,0.595),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.587, 0.595))

#jitter plot LD03 filtered
ggplot(merge_cv_error, aes(x = factor(merge_cv_error[,1]), y = merge_cv_error[,2])) +
  geom_hline(
    yintercept = c(0.587,0.588,0.589,0.59, 0.591, 0.592, 0.593,0.594,0.595),
    color = "black",
    linetype = "solid",
    size = 0.5
  ) +
  geom_boxplot(width = 0.5, fill = "yellow", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "red") +
  labs(title = "Cross-validation Error Plot",
       x = "K",
       y = "CV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0.587, 0.595))
```

### CV plot pour légende de figure d'Admixture

#### SeqApiPop 629 échantillons - MAF \> 0.01 - LD pruning = 0.3 (fenêtre de 1749 SNPs et pas de 175 bp)

```{r}
# Définir les valeurs et les noms des groupes
valeurs <- c(0.475215, 0.4747, 0.47582, 0.4794, 0.4888, 0.5012)

# Créer le barplot avec des barres plus espacées et moins larges
bp <- barplot(valeurs, horiz = TRUE, xlim = c(0.47, 0.51), width = 0.1,
        xlab = "CV", las = 1, col = "darkblue", border = NA)

# Ajouter un cadrillage
abline(v = seq(0.47, 0.51, by = 0.005), col = "lightgray", lty = 2)
```

#### SeqApiPop 629 échantillons - SNPsBeeMuSe filtered

```{r}
# Définir les valeurs et les noms des groupes
valeurs <- c(0.7218, 0.7219, 0.7242, 0.7282, 0.7405, 0.7955)

# Créer le barplot avec des barres plus espacées et moins larges
bp <- barplot(valeurs, horiz = TRUE, xlim = c(0.72, 0.8), width = 0.1,
        xlab = "CV", las = 1, col = "darkblue", border = NA)

# Ajouter un cadrillage
abline(v = seq(0.72, 0.8, by = 0.01), col = "lightgray", lty = 2)
```

#### SeqApiPop 629 échantillons - SNPsBeeMuSe filtered - MAF \> 0.01 - LD pruning = 0.1 (fenêtre de 50 SNPs et pas de 10 bp)

```{r}
# Définir les valeurs et les noms des groupes
valeurs <- c(0.7593, 0.7595, 0.7669, 0.7815, 0.8179)

# Créer le barplot avec des barres plus espacées et moins larges
bp <- barplot(valeurs, horiz = TRUE, xlim = c(0.75, 0.82), width = 0.1,
        xlab = "CV", las = 1, col = "darkblue", border = NA)

# Ajouter un cadrillage
abline(v = seq(0.75, 0.82, by = 0.01), col = "lightgray", lty = 2)
```

#### SeqApiPop 561 échantillons - SNPsBeeMuSe filtered - MAF \> 0.01 - LD pruning = 0.1 (fenêtre de 50 SNPs et pas de 10 bp)

```{r}
# Définir les valeurs et les noms des groupes
valeurs <- c(0.7734 ,0.7709, 0.7680, 0.7715, 0.7871, 0.8271)

# Créer le barplot avec des barres plus espacées et moins larges
bp <- barplot(valeurs, horiz = TRUE, xlim = c(0.76, 0.83), width = 0.1,
        xlab = "CV", las = 1, col = "darkblue", border = NA)

# Ajouter un cadrillage
abline(v = seq(0.76, 0.83, by = 0.01), col = "lightgray", lty = 2)
```

## Admixture supervisé : Création fichier liste individu - population

### 561 échantillons - MAF \> 0.01 - LD pruning = 0.3 (fenêtre de 1749 SNPs et pas de 175 bp) - 3848 SNPs

##### **K = 3**

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles/SeqApiPop_561_maf001_LD03")

Q_3_561 <- read.table("SeqApiPop_561_maf001_LD03_pruned.3.r10.Q", header = FALSE)

colnames(Q_3_561)[colnames(Q_3_561) == "V1"] <- "Vert"
colnames(Q_3_561)[colnames(Q_3_561) == "V2"] <- "Noir"
colnames(Q_3_561)[colnames(Q_3_561) == "V3"] <- "Orange"

# Create an empty vector to store the category for each row
categories <- character(nrow(Q_3_561))

# Initialisation du vecteur de catégories
categories <- rep("-", nrow(Q_3_561))

# Itérer à travers chaque ligne
for (i in 1:nrow(Q_3_561)) {
  # Vérifier si aucune valeur dans la ligne ne dépasse 0.9
  if (all(Q_3_561[i,] <= 0.9)) {
    categories[i] <- "-"
  } else {
    # Vérifier quelle colonne a la valeur supérieure à 0.9
    if (Q_3_561[i,1] > 0.9) {
      categories[i] <- "Vert"
    } else if (Q_3_561[i,2] > 0.9) {
      categories[i] <- "Noir"
    } else if (Q_3_561[i,3] > 0.9) {
      categories[i] <- "Orange"
    }
  }
}

# Write the categories to a single list file
write(categories, file = "output_list_K3_561_90.txt")
```

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles")

output_list_K3_561_95_merged <- readLines("output_list_K3_561_95_merged.txt")
output_list_K3_561_95_merged <- gsub("Beemuse", "-", output_list_K3_561_95_merged)
writeLines(output_list_K3_561_95_merged, "merged_data_K3_561_95.pop")

output_list_K3_561_90_merged <- readLines("output_list_K3_561_90_merged.txt")
output_list_K3_561_90_merged <- gsub("Beemuse", "-", output_list_K3_561_90_merged)
writeLines(output_list_K3_561_90_merged, "merged_data_K3_561_90.pop")
```

```{r}
texte_complet <- paste(output_list_K3_561_95_merged, collapse = " ")
K3_95 <- unlist(strsplit(texte_complet, "\\s+"))
nombre_apparitions <- table(K3_95)
print(nombre_apparitions)

texte_complet <- paste(output_list_K3_561_90_merged, collapse = " ")
K3_90 <- unlist(strsplit(texte_complet, "\\s+"))
nombre_apparitions <- table(K3_90)
print(nombre_apparitions)
```

##### K = 5

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_561_maf001_LD03")

labels <- read.csv('~/Documents/Stage_NB/data/SeqApiPop_labels.csv')
samples_561 <- read.table("SeqApiPop_561_maf001_LD03_pruned.fam", header = FALSE)

samples_561 <- samples_561[, 1:2] # Keep only the first two columns
colnames(samples_561)[colnames(samples_561) == "V1"] <- "name"

merged_labels_samples_561 <- merge(labels, samples_561,  by = 'name')

Label_samples_561 <- subset(merged_labels_samples_561, select = "Label") 

writeLines(as.character(Label_samples_561$Label), "Label_samples_561.txt")
```

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles/SeqApiPop_561_maf001_LD03")

Q_5_561 <- read.table("SeqApiPop_561_maf001_LD03_pruned.5.r10.Q", header = FALSE)

colnames(Q_5_561)[colnames(Q_5_561) == "V1"] <- "Bleu"
colnames(Q_5_561)[colnames(Q_5_561) == "V2"] <- "Jaune"
colnames(Q_5_561)[colnames(Q_5_561) == "V3"] <- "Vert"
colnames(Q_5_561)[colnames(Q_5_561) == "V4"] <- "Orange"
colnames(Q_5_561)[colnames(Q_5_561) == "V5"] <- "Noir"

# Create an empty vector to store the category for each row
categories <- character(nrow(Q_5_561))

# Initialisation du vecteur de catégories
categories <- rep("-", nrow(Q_5_561))

# Itérer à travers chaque ligne
for (i in 1:nrow(Q_5_561)) {
  # Vérifier si aucune valeur dans la ligne ne dépasse 0.95
  if (all(Q_5_561[i,] <= 0.95)) {
    categories[i] <- "-"
  } else {
    # Vérifier quelle colonne a la valeur supérieure à 0.95
    if (Q_5_561[i,1] > 0.95) {
      categories[i] <- "Bleu"
    } else if (Q_5_561[i,2] > 0.95) {
      categories[i] <- "Jaune"
    } else if (Q_5_561[i,3] > 0.95) {
      categories[i] <- "Vert"
    } else if (Q_5_561[i,4] > 0.95) {
      categories[i] <- "Orange"
    } else if (Q_5_561[i,5] > 0.95) {
      categories[i] <- "Noir"
    }
  }
}

# Write the categories to a single list file
write(categories, file = "output_list_K5_561.txt")
```

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles")

# Read the file content
output_list_K5_561_merged <- readLines("output_list_K5_561_merged.txt")

# Replace "Beemuse" with "-"
output_list_K5_561_merged <- gsub("Beemuse", "-", output_list_K5_561_merged)

# Write the modified content back to the file
writeLines(output_list_K5_561_merged, "merged_data_K5_561.pop")
```

##### K = 6

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles/SeqApiPop_561_maf001_LD03")

Q_6_561 <- read.table("SeqApiPop_561_maf001_LD03_pruned.6.r13.Q", header = FALSE)

colnames(Q_6_561)[colnames(Q_6_561) == "V1"] <- "Rouge"
colnames(Q_6_561)[colnames(Q_6_561) == "V2"] <- "Bleu"
colnames(Q_6_561)[colnames(Q_6_561) == "V3"] <- "Orange"
colnames(Q_6_561)[colnames(Q_6_561) == "V4"] <- "Vert"
colnames(Q_6_561)[colnames(Q_6_561) == "V5"] <- "Jaune"
colnames(Q_6_561)[colnames(Q_6_561) == "V6"] <- "Noir"

# Create an empty vector to store the category for each row
categories <- character(nrow(Q_6_561))

# Initialisation du vecteur de catégories
categories <- rep("-", nrow(Q_6_561))

# Itérer à travers chaque ligne
for (i in 1:nrow(Q_6_561)) {
  # Vérifier si aucune valeur dans la ligne ne dépasse 0.95
  if (all(Q_6_561[i,] <= 0.95)) {
    categories[i] <- "-"
  } else {
    # Vérifier quelle colonne a la valeur supérieure à 0.95
    if (Q_6_561[i,1] > 0.95) {
      categories[i] <- "Rouge"
    } else if (Q_6_561[i,2] > 0.95) {
      categories[i] <- "Bleu"
    } else if (Q_6_561[i,3] > 0.95) {
      categories[i] <- "Orange"
    } else if (Q_6_561[i,4] > 0.95) {
      categories[i] <- "Vert"
    } else if (Q_6_561[i,5] > 0.95) {
      categories[i] <- "Jaune"
    } else if (Q_6_561[i,6] > 0.95) {
      categories[i] <- "Noir"
    }
  }
}

# Write the categories to a single list file
write(categories, file = "output_list_K6_561.txt")
```

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles")

# Read the file content
output_list_K6_561_merged <- readLines("output_list_K6_561_merged.txt")

# Replace "Beemuse" with "-"
output_list_K6_561_merged <- gsub("Beemuse", "-", output_list_K6_561_merged)

# Write the modified content back to the file
writeLines(output_list_K6_561_merged, "merged_data_K6_561.pop")
```

### 561 échantillons - MAF \> 0.01 - LD pruning = 0.1 (fenêtre de 50 SNPs et pas de 10 bp) - 1055 SNPs

##### **K = 3**

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles/SeqApiPop_561_LD03_default_1055")

Q_3_561_default_default <- read.table("SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD03_default_pruned.3.r0.Q", header = FALSE)

colnames(Q_3_561_default)[colnames(Q_3_561_default) == "V1"] <- "Noir"
colnames(Q_3_561_default)[colnames(Q_3_561_default) == "V2"] <- "Vert"
colnames(Q_3_561_default)[colnames(Q_3_561_default) == "V3"] <- "Orange"

# Create an empty vector to store the category for each row
categories <- character(nrow(Q_3_561_default))

# Initialisation du vecteur de catégories
categories <- rep("-", nrow(Q_3_561_default))

# Itérer à travers chaque ligne
for (i in 1:nrow(Q_3_561_default)) {
  # Vérifier si aucune valeur dans la ligne ne dépasse 0.9
  if (all(Q_3_561_default[i,] <= 0.9)) {
    categories[i] <- "-"
  } else {
    # Vérifier quelle colonne a la valeur supérieure à 0.9
    if (Q_3_561_default[i,1] > 0.9) {
      categories[i] <- "Noir"
    } else if (Q_3_561_default[i,2] > 0.9) {
      categories[i] <- "Vert"
    } else if (Q_3_561_default[i,3] > 0.9) {
      categories[i] <- "Orange"
    } 
  }
}

# Write the categories to a single list file
write(categories, file = "output_list_K3_561_LD_default_90.txt")
```

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles")

output_list_K3_561_default_95_merged <- readLines("output_list_K3_561_LD_default_95_merged.txt")
output_list_K3_561_default_95_merged <- gsub("Beemuse", "-", output_list_K3_561_default_95_merged)
writeLines(output_list_K3_561_default_95_merged, "merged_data_K3_561_LD_default_95.pop")

output_list_K3_561_default_90_merged <- readLines("output_list_K3_561_LD_default_90_merged.txt")
output_list_K3_561_default_90_merged <- gsub("Beemuse", "-", output_list_K3_561_default_90_merged)
writeLines(output_list_K3_561_default_90_merged, "merged_data_K3_561_LD_default_90.pop")
```

```{r}
texte_complet <- paste(output_list_K3_561_default_95_merged, collapse = " ")
K3_95 <- unlist(strsplit(texte_complet, "\\s+"))
nombre_apparitions <- table(K3_95)
print(nombre_apparitions)

texte_complet <- paste(output_list_K3_561_default_90_merged, collapse = " ")
K3_90 <- unlist(strsplit(texte_complet, "\\s+"))
nombre_apparitions <- table(K3_90)
print(nombre_apparitions)
```

##### K = 5

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles/SeqApiPop_561_LD03_default_1055")

Q_5_561_default <- read.table("SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD03_default_pruned.5.r22.Q", header = FALSE)

colnames(Q_5_561_default)[colnames(Q_5_561_default) == "V1"] <- "Rouge"
colnames(Q_5_561_default)[colnames(Q_5_561_default) == "V2"] <- "Vert"
colnames(Q_5_561_default)[colnames(Q_5_561_default) == "V3"] <- "Noir"
colnames(Q_5_561_default)[colnames(Q_5_561_default) == "V4"] <- "Orange"
colnames(Q_5_561_default)[colnames(Q_5_561_default) == "V5"] <- "Jaune"

categories <- character(nrow(Q_5_561_default))
categories <- rep("-", nrow(Q_5_561_default))

for (i in 1:nrow(Q_5_561_default)) {
  # Vérifier si aucune valeur dans la ligne ne dépasse 0.95
  if (all(Q_5_561_default[i,] <= 0.95)) {
    categories[i] <- "-"
  } else {
    # Vérifier quelle colonne a la valeur supérieure à 0.95
    if (Q_5_561_default[i,1] > 0.95) {
      categories[i] <- "Rouge"
    } else if (Q_5_561_default[i,2] > 0.95) {
      categories[i] <- "Vert"
    } else if (Q_5_561_default[i,3] > 0.95) {
      categories[i] <- "Noir"
    } else if (Q_5_561_default[i,4] > 0.95) {
      categories[i] <- "Orange"
    } else if (Q_5_561_default[i,5] > 0.95) {
      categories[i] <- "Jaune"
    } 
  }
}

write(categories, file = "output_list_K5_561_LD_default_95.txt")
```

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles")

output_list_K5_561_default_95_merged <- readLines("output_list_K5_561_LD_default_95_merged.txt")
output_list_K5_561_default_95_merged <- gsub("Beemuse", "-", output_list_K5_561_default_95_merged)
writeLines(output_list_K5_561_default_95_merged, "merged_data_K5_561_LD_default_95.pop")

output_list_K5_561_default_90_merged <- readLines("output_list_K5_561_LD_default_90_merged.txt")
output_list_K5_561_default_90_merged <- gsub("Beemuse", "-", output_list_K5_561_default_90_merged)
writeLines(output_list_K5_561_default_90_merged, "merged_data_K5_561_LD_default_90.pop")
```

```{r}
texte_complet <- paste(output_list_K5_561_default_95_merged, collapse = " ")
K5_95 <- unlist(strsplit(texte_complet, "\\s+"))
nombre_apparitions <- table(K5_95)
print(nombre_apparitions)

texte_complet <- paste(output_list_K5_561_default_90_merged, collapse = " ")
K5_90 <- unlist(strsplit(texte_complet, "\\s+"))
nombre_apparitions <- table(K5_90)
print(nombre_apparitions)
```

##### **K = 6**

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles/SeqApiPop_561_LD03_default_1055")

Q_6_561_default <- read.table("SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD03_default_pruned.6.r23.Q", header = FALSE)

colnames(Q_6_561_default)[colnames(Q_6_561_default) == "V1"] <- "Bleu"
colnames(Q_6_561_default)[colnames(Q_6_561_default) == "V2"] <- "Rouge"
colnames(Q_6_561_default)[colnames(Q_6_561_default) == "V3"] <- "Noir"
colnames(Q_6_561_default)[colnames(Q_6_561_default) == "V4"] <- "Jaune"
colnames(Q_6_561_default)[colnames(Q_6_561_default) == "V5"] <- "Orange"
colnames(Q_6_561_default)[colnames(Q_6_561_default) == "V6"] <- "Vert"

categories <- character(nrow(Q_6_561_default))
categories <- rep("-", nrow(Q_6_561_default))

for (i in 1:nrow(Q_6_561_default)) {
  # Vérifier si aucune valeur dans la ligne ne dépasse 0.9
  if (all(Q_6_561_default[i,] <= 0.9)) {
    categories[i] <- "-"
  } else {
    # Vérifier quelle colonne a la valeur supérieure à 0.9
    if (Q_6_561_default[i,1] > 0.9) {
      categories[i] <- "Bleu"
    } else if (Q_6_561_default[i,2] > 0.9) {
      categories[i] <- "Rouge"
    } else if (Q_6_561_default[i,3] > 0.9) {
      categories[i] <- "Noir"
    } else if (Q_6_561_default[i,4] > 0.9) {
      categories[i] <- "Jaune"
    } else if (Q_6_561_default[i,5] > 0.9) {
      categories[i] <- "Orange"
    } else if (Q_6_561_default[i,6] > 0.9) {
      categories[i] <- "Vert"
    }
  }
}

write(categories, file = "output_list_K6_561_LD_default_90.txt")
```

```{r}
setwd("~/Documents/Stage_NB/data/Qfiles")

output_list_K6_561_default_95_merged <- readLines("output_list_K6_561_LD_default_95_merged.txt")
output_list_K6_561_default_95_merged <- gsub("Beemuse", "-", output_list_K6_561_default_95_merged)
writeLines(output_list_K6_561_default_95_merged, "merged_data_K6_561_LD_default_95.pop")

output_list_K6_561_default_90_merged <- readLines("output_list_K6_561_LD_default_90_merged.txt")
output_list_K6_561_default_90_merged <- gsub("Beemuse", "-", output_list_K6_561_default_90_merged)
writeLines(output_list_K6_561_default_90_merged, "merged_data_K6_561_LD_default_90.pop")
```

```{r}
texte_complet <- paste(output_list_K6_561_default_95_merged, collapse = " ")
K6_95 <- unlist(strsplit(texte_complet, "\\s+"))
nombre_apparitions <- table(K6_95)
print(nombre_apparitions)

texte_complet <- paste(output_list_K6_561_default_90_merged, collapse = " ")
K6_90 <- unlist(strsplit(texte_complet, "\\s+"))
nombre_apparitions <- table(K6_90)
print(nombre_apparitions)
```

## FST plot

### Création du fichier FST ind2pop - Merged Data - SeqApiPop BeeMuSe

```{r}
#output_list_1_merged.txt
setwd("~/Documents/Stage_NB/data/merged_BeeMuSe_SeqApiPop_629_filtered_maf001_LD03") #2 premières colonnes .fam

ind2pop_merged_data <- read.table("ind2pop_merged_data.txt", header = FALSE)

# Remplacement de la deuxième colonne par les valeurs spécifiées
ind2pop_merged_data$V2 <- ind_pop

write.table(ind2pop_merged_data, "ind2pop_merged_data_3.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
```

```{r}
setwd("~/Documents/Stage_NB/data/merged_BeeMuSe_SeqApiPop_629_filtered_maf001_LD03") 
fst_file <- read.table("ind2pop_seqapipop_629_fst.txt", header = FALSE)
print(unique(fst_file$V2))

fst_file <- cbind(fst_file[, "V1"], fst_file, fst_file[, "V1"])

# Remove the original V1 column
fst_file <- fst_file[, -4]

# Rename the duplicated V1 column to V2
names(fst_file)[1] <- "V1"
names(fst_file)[2] <- "V2"
names(fst_file)[3] <- "V3"

write.table(fst_file, "ind2pop_seqapipop_629_fst_2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
```

## 

#### **SeqApiPop 629 échantillons**

```{r}
setwd("~/Documents/Stage_NB/data/fst_files") 
# Charger les données
fst_data_SeqApiPop_629 <- read.table("SeqApiPop_629_fst.fst", header=TRUE)

fstsubset <- fst_data_SeqApiPop_629[complete.cases(fst_data_SeqApiPop_629),]
SNP <- c(1:(nrow(fstsubset)))
mydf <- data.frame(SNP, fstsubset)
 
manhattan(mydf,chr="CHR",bp="POS",p="FST",snp="SNP",logp=FALSE,ylab="Fst")

# Nuage de points
plot(fst_data_SeqApiPop_629$POS, fst_data_SeqApiPop_629$FST, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
# Histogramme
barplot(fst_data_SeqApiPop_629$FST, names.arg=fst_data_SeqApiPop_629$POS, xlab="Position (POS)", ylab=" Fst", main="Graphique FST")

# Créer un violin plot des valeurs de FST
vioplot(fst_data_SeqApiPop_629$FST, names="FST", col="#FF0000",ylim=c(0, 1), horizontal=FALSE, main="Violin Plot FST - SeqApiPop_629")
```

#### SeqApiPop 629 échantillons - MAF \> 0.01

```{r}
setwd("~/Documents/Stage_NB/data/fst_files") 
# Charger les données
fst_data_SeqApiPop_629_maf001 <- read.table("SeqApiPop_629_maf001_fst.fst", header=TRUE)

fstsubset2 <- fst_data_SeqApiPop_629_maf001[complete.cases(fst_data_SeqApiPop_629_maf001),]
SNP2 <- c(1:(nrow(fstsubset2)))
mydf2 <- data.frame(SNP2, fstsubset2)
 
manhattan(mydf2,chr="CHR",bp="POS",p="FST",snp="SNP",logp=FALSE,ylab="Fst")

# Nuage de points
plot(fst_data_SeqApiPop_629_maf001$POS, fst_data_SeqApiPop_629_maf001$FST, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
# Histogramme
barplot(fst_data_SeqApiPop_629_maf001$FST, names.arg=fst_data_SeqApiPop_629_maf001$POS, xlab="Position (POS)", ylab=" Fst", main="Graphique FST")

# Créer un violin plot des valeurs de FST
vioplot(fst_data_SeqApiPop_629_maf001$FST, names="FST", col="#FFA500",ylim=c(0, 1), horizontal=FALSE, main="Violin Plot FST - SeqApiPop_629_maf001")
```

#### SeqApiPop 629 échantillons - MAF \> 0.01 - LD pruning = 0.3 (fenêtre de 1749 SNPs et pas de 175 bp)

```{r}
setwd("~/Documents/Stage_NB/data/fst_files") 
# Charger les données
fst_data_SeqApiPop_629_maf001_LD03 <- read.table("SeqApiPop_629_maf001_LD03_pruned_fst.fst", header=TRUE)

fstsubset3 <- fst_data_SeqApiPop_629_maf001_LD03[complete.cases(fst_data_SeqApiPop_629_maf001_LD03),]
SNP3 <- c(1:(nrow(fstsubset3)))
mydf3 <- data.frame(SNP3, fstsubset3)
 
manhattan(mydf3,chr="CHR",bp="POS",p="FST",snp="SNP",logp=FALSE,ylab="Fst")

# Nuage de points
plot(fst_data_SeqApiPop_629_maf001_LD03$POS, fst_data_SeqApiPop_629_maf001_LD03$FST, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
# Histogramme
barplot(fst_data_SeqApiPop_629_maf001_LD03$FST, names.arg=fst_data_SeqApiPop_629_maf001_LD03$POS, xlab="Position (POS)", ylab=" Fst", main="Graphique FST")

# Créer un violin plot des valeurs de FST
vioplot(fst_data_SeqApiPop_629_maf001_LD03$FST, names="FST", col="#FFFF00",ylim=c(0, 1), horizontal=FALSE, main="Violin Plot FST - SeqApiPop_629_maf001_LD03")
```

```{r}
# Création d'une variable pour chaque jeu de données
fst_data_SeqApiPop_629$Group <- "seqapipop_629"
fst_data_SeqApiPop_629_maf001$Group <- "seqapipop_629_maf001"
fst_data_SeqApiPop_629_maf001_LD03$Group <- "seqapipop_629_maf001_LD03"

# Combiner les jeux de données en un seul dataframe
combined_data <- rbind(fst_data_SeqApiPop_629, fst_data_SeqApiPop_629_maf001, fst_data_SeqApiPop_629_maf001_LD03)

# Changer le type de la variable 'Group' en factor avec l'ordre spécifié
combined_data$Group <- factor(combined_data$Group, levels = group_order)

# Créer le violin plot avec ggplot
ggplot(combined_data, aes(x = Group, y = FST, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = c("seqapipop_629" = "#FF0000", 
                               "seqapipop_629_maf001" = "#FFA500", 
                               "seqapipop_629_maf001_LD03" = "#FFFF00")) +
  labs(title = "Violin Plot of FST values",
       x = "Group",
       y = "FST") +
  theme_minimal()
```

#### SeqApiPop 629 échantillons - SNPsBeeMuSe filtered

##### 10030 SNPs

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_629_SNPsBeeMuSe") 

# Charger les données
fst_data_10030 <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_fst.fst", header=TRUE)

fstsubset <- fst_data_10030[complete.cases(fst_data_10030),]
SNP <- c(1:(nrow(fstsubset)))
mydf <- data.frame(SNP, fstsubset)
 
manhattan(mydf,chr="CHR",bp="POS",p="FST"
,snp="SNP",logp=FALSE,ylab="Fst")

# Nuage de points
#plot(fst_data_10030$POS, fst_data_10030$FST, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
# Histogramme
#barplot(fst_data_10030$FST, names.arg=fst_data_10030$POS, xlab="Position (POS)", ylab=" Fst", main="Graphique FST")

# Créer un violin plot des valeurs de FST
vioplot(fst_data_10030$FST, names="FST", col="#1f77b4",ylim=c(0, 1), horizontal=FALSE, main="Violin Plot FST - 10030 SNPs")

# Filtrer les données pour inclure uniquement les positions du chromosome 11
fst_data_10030_chr11 <- subset(fst_data_10030, CHR == "11")

SNP <- seq_len(nrow(fst_data_10030_chr11))
mydf_11 <- data.frame(SNP = SNP, fst_data_10030_chr11)
manhattan(mydf_11, chr = "CHR", bp = "POS", p = "FST", snp = "SNP", logp = FALSE, ylab = "Fst")

#plot(fst_data_10030_chr11$POS, fst_data_10030_chr11$FST, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
#barplot(fst_data_10030_chr11$FST, names.arg=fst_data_10030_chr11$POS, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
vioplot(fst_data_10030_chr11$FST, names="FST", col="#aec7e8", ylim=c(0, 1),horizontal=FALSE, main="Violin Plot FST - Chr 11 - 10030 SNPs")
```

##### 3848 SNPs

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_629_SNPsBeeMuSe") 

# Charger les données
fst_data_3848 <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_pruned_fst.fst", header=TRUE)


fstsubset <- fst_data_3848[complete.cases(fst_data_3848),]
SNP <- c(1:(nrow(fstsubset)))
mydf <- data.frame(SNP, fstsubset)
 
manhattan(mydf,chr="CHR",bp="POS",p="FST"
,snp="SNP",logp=FALSE,ylab="Fst")

# Nuage de points
#plot(fst_data_3848$POS, fst_data_3848$FST, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
# Histogramme
#barplot(fst_data_3848$FST, names.arg=fst_data_3848$POS, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")

vioplot(fst_data_3848$FST, names="FST", col="#d62728", ylim=c(0, 1),horizontal=FALSE, main="Violin Plot FST - 3848 SNPs")

# Filtrer les données pour inclure uniquement les positions du chromosome 11
fst_data_3848_chr11 <- subset(fst_data_3848, CHR == "11")

SNP <- seq_len(nrow(fst_data_3848_chr11))
mydf_11 <- data.frame(SNP = SNP, fst_data_3848_chr11)
manhattan(mydf_11, chr = "CHR", bp = "POS", p = "FST", snp = "SNP", logp = FALSE, ylab = "Fst")

#plot(fst_data_3848_chr11$POS, fst_data_3848_chr11$FST, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
#barplot(fst_data_3848_chr11$FST, names.arg=fst_data_3848_chr11$POS, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
vioplot(fst_data_3848_chr11$FST, names="FST", col="#ff9896",ylim=c(0, 1), horizontal=FALSE, main="Violin Plot FST - Chr 11 - 3848 SNPs")
```

##### 1055 SNPs

```{r}
setwd("~/Documents/Stage_NB/data/SeqApiPop_629_SNPsBeeMuSe") 

# Charger les données
fst_data_1055 <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_default_pruned_fst.fst", header=TRUE)


fstsubset <- fst_data_1055[complete.cases(fst_data_1055),]
SNP <- c(1:(nrow(fstsubset)))
mydf <- data.frame(SNP, fstsubset)
 
manhattan(mydf,chr="CHR",bp="POS",p="FST"
,snp="SNP",logp=FALSE,ylab="Fst")

# Nuage de points
#plot(fst_data_1055$POS, fst_data_1055$FST, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
# Histogramme
#barplot(fst_data_1055$FST, names.arg=fst_data_1055$POS, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")

vioplot(fst_data_1055$FST, names="FST", col="#FFA500", ylim=c(0, 1),horizontal=FALSE, main="Violin Plot FST - 3848 SNPs")

# Filtrer les données pour inclure uniquement les positions du chromosome 11
fst_data_1055_chr11 <- subset(fst_data_1055, CHR == "11")

SNP <- seq_len(nrow(fst_data_1055_chr11))
mydf_11 <- data.frame(SNP = SNP, fst_data_1055_chr11)
manhattan(mydf_11, chr = "CHR", bp = "POS", p = "FST", snp = "SNP", logp = FALSE, ylab = "Fst")

#plot(fst_data_1055_chr11$POS, fst_data_1055_chr11$FST, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
#barplot(fst_data_1055_chr11$FST, names.arg=fst_data_1055_chr11$POS, xlab="Position (POS)", ylab="Valeurs FST", main="Graphique FST")
vioplot(fst_data_1055_chr11$FST, names="FST", col="#ffcc7a",ylim=c(0, 1), horizontal=FALSE, main="Violin Plot FST - Chr 11 - 3848 SNPs")
```

```{r}
# Création d'un facteur pour distinguer les différentes données
fst_data_10030$Dataset <- "fst_data_10030"
fst_data_10030_chr11$Dataset <- "fst_data_10030_chr11"
fst_data_3848$Dataset <- "fst_data_3848"
fst_data_3848_chr11$Dataset <- "fst_data_3848_chr11"
fst_data_1055$Dataset <- "fst_data_1055"
fst_data_1055_chr11$Dataset <- "fst_data_1055_chr11"

fst_dataset <- rbind(fst_data_10030, fst_data_10030_chr11, fst_data_3848, fst_data_3848_chr11, fst_data_1055, fst_data_1055_chr11)

# Define the order of groups
group_order <- c("fst_data_10030", "fst_data_10030_chr11", 
                 "fst_data_3848", "fst_data_3848_chr11", 
                 "fst_data_1055", "fst_data_1055_chr11")
# Convert Group column to factor with custom levels
fst_dataset$Dataset <- factor(fst_dataset$Dataset, levels = group_order)

# Create the violin plot with y-axis limits set from 0 to 1
ggplot(fst_dataset, aes(x = Dataset, y = FST, fill = Dataset)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = c("fst_data_10030" = "#1f77b4", 
                               "fst_data_10030_chr11" = "#aec7e8", 
                               "fst_data_3848" = "#d62728",
                               "fst_data_3848_chr11" = "#ff9896",
                               "fst_data_1055" = "#FFA500",
                               "fst_data_1055_chr11" = "#ffcc7a")) +
  labs(title = "Violin Plot of FST values",
       x = "Group",
       y = "FST") +
  theme_minimal() +
  ylim(0, 1) 


# Calculate mean values for each dataset
mean_values <- aggregate(FST ~ Dataset, data = fst_dataset, FUN = mean)

# Print mean values
print(mean_values)
```

## Degré d'apparentement

### IBS - SeqApiPop 629 - SNPsBeeMuse - maf001 LD03

```{r}
setwd("~/Documents/Stage_NB/data/IBD") 

# Charger les données des coefficients de parenté
relatedness <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_pruned_ibd.genome", header=TRUE) 

hist(relatedness$DST, breaks=20, col="lightgrey", main="Histogramme de la distance génétique entre les 629 individus SeqApiPop", xlab="Distance génétique", cex.main=0.8)

# Créer un violin plot pour la distance génétique entre les individus SeqApiPop
ggplot(relatedness, aes(x = "", y = DST, fill = "SeqApiPop_629")) +
  geom_violin(trim = FALSE, color = "black") +
  scale_fill_manual(values = c("SeqApiPop_629" = "lightgrey")) +
  labs(title = "Violin Plot de la distance génétique entre les 629 individus SeqApiPop",
       x = "", y = "Distance génétique") +
  theme_minimal()

relatedness_car_ger_pol<- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_pruned_ibd_carnica_germany_poland.genome", header=TRUE) 

hist(relatedness_car_ger_pol$DST, breaks=20, col="orange", main="Histogramme de la distance génétique entre individus Carnica Germany - Carnica Poland", xlab="Distance génétique", cex.main=0.8)

# Créer un violin plot pour la distance génétique entre individus Carnica Germany - Carnica Poland
ggplot(relatedness_car_ger_pol, aes(x = "", y = DST, fill = "Carnica Germany - Carnica Poland")) +
  geom_violin(trim = FALSE, color = "black") +
  scale_fill_manual(values = c("Carnica Germany - Carnica Poland" = "orange")) +
  labs(title = "Violin Plot de la distance génétique entre individus Carnica Germany - Poland",
       x = "", y = "Distance génétique") +
  theme_minimal()

relatedness_mel_oue_col <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_pruned_ibd_mellifera_ouessant_colonsay.genome", header=TRUE) 

hist(relatedness_mel_oue_col$DST, breaks=20, col="black", main="Histogramme de la distance génétique entre individus Mellifera Ouessant - Mellifera Colonsay", xlab="Distance génétique", cex.main=0.8)

# Créer un violin plot pour la distance génétique entre individus Mellifera Ouessant - Mellifera Colonsay
ggplot(relatedness_mel_oue_col, aes(x = "", y = DST, fill = "Mellifera Ouessant - Mellifera Colonsay")) +
  geom_violin(trim = FALSE, color = "black") +
  scale_fill_manual(values = c("Mellifera Ouessant - Mellifera Colonsay" = "black")) +
  labs(title = "Violin Plot de la distance génétique entre individus Mellifera Ouessant - Colonsay",
       x = "", y = "Distance génétique") +
  theme_minimal()

# Charger les données des coefficients de parenté
relatedness_rep <- read.table("SeqApiPop_629_SNPsBeeMuSe_filtered_maf001_LD03_pruned_ibd_representants.genome", header=TRUE) 

hist(relatedness_rep$DST, breaks=20, col="grey", main="Histogramme de la distance génétique entre les individus représentants SeqApiPop", xlab="Distance génétique", cex.main=0.8)

# Créer un violin plot pour la distance génétique entre les individus SeqApiPop
ggplot(relatedness_rep, aes(x = "", y = DST, fill = "SeqApiPop_629")) +
  geom_violin(trim = FALSE, color = "black") +
  scale_fill_manual(values = c("SeqApiPop_629" = "grey")) +
  labs(title = "Violin Plot de la distance génétique entre les individus représentants SeqApiPop",
       x = "", y = "Distance génétique") +
  theme_minimal()
```

```{r}
matrice_dist <- acast(relatedness, IID1 ~ IID2, value.var = "DST")
matrice_dist2 <- acast(relatedness_car_ger_pol, IID1 ~ IID2, value.var = "DST")
matrice_dist3 <- acast(relatedness_mel_oue_col, IID1 ~ IID2, value.var = "DST")
matrice_dist4 <- acast(relatedness_rep, IID1 ~ IID2, value.var = "DST")

heatmap(matrice_dist, symm = TRUE, scale = "none", col = rev(heat.colors(256)),
        margins = c(10, 10), main = "Heatmap de la matrice de distance génétique",
        cexRow = 0.3, cexCol = 0.3)

heatmap(matrice_dist2, symm = TRUE, scale = "none", col = rev(heat.colors(256)),
        margins = c(10, 10), main = "Heatmap de la matrice de distance génétique",
        cexRow = 0.5, cexCol = 0.5)

heatmap(matrice_dist3, symm = TRUE, scale = "none", col = rev(heat.colors(256)),
        margins = c(10, 10), main = "Heatmap de la matrice de distance génétique",
        cexRow = 0.5, cexCol = 0.5)

heatmap(matrice_dist4, symm = TRUE, scale = "none", col = rev(heat.colors(256)),
        margins = c(10, 10), main = "Heatmap de la matrice de distance génétique",
        cexRow = 0.4, cexCol = 0.4)
```

```{r}
ggplot(data = melt(matrice_dist), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap de la matrice de distance génétique") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 2)) +
  theme(axis.text.y = element_text(hjust = 1, size = 2))

ggplot(data = melt(matrice_dist2), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap de la matrice de distance génétique") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
  theme(axis.text.y = element_text(hjust = 1, size = 5))

ggplot(data = melt(matrice_dist3), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap de la matrice de distance génétique") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
  theme(axis.text.y = element_text(hjust = 1, size = 5))

ggplot(data = melt(matrice_dist4), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap de la matrice de distance génétique") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
  theme(axis.text.y = element_text(hjust = 1, size = 5))
```

```{r}
# Créer un heatmap avec pheatmap
pheatmap(matrice_dist, 
         color = colorRampPalette(c("blue", "red"))(256), 
         fontsize_row = 3, fontsize_col = 3,
         main = "Heatmap de la matrice de distance génétique")

pheatmap(matrice_dist2, 
         color = colorRampPalette(c("blue", "red"))(256), 
         fontsize_row = 5, fontsize_col = 5,
         main = "Heatmap de la matrice de distance génétique")

pheatmap(matrice_dist3, 
         color = colorRampPalette(c("blue", "red"))(256), 
         fontsize_row = 5, fontsize_col = 5,
         main = "Heatmap de la matrice de distance génétique")

pheatmap(matrice_dist4, 
         color = colorRampPalette(c("blue", "red"))(256), 
         fontsize_row = 5, fontsize_col = 5,
         main = "Heatmap de la matrice de distance génétique")
```

### **IBD - BeeMuSe**

```{r}
setwd("~/Documents/Stage_NB/data/IBD") 

beemuse_genome <- read.table("BeeMuse_genome_full.genome", header=TRUE) 
beemuse_BAH_20_19_genome <- read.table("BeeMuse_BAH_20-19_genome.genome", header=TRUE) 
beemuse_BBS_6_19_genome <- read.table("BeeMuse_BBS_6-19_genome.genome", header=TRUE) 
beemuse_BER_11_19_genome <- read.table("BeeMuse_BER_11-19_genome.genome", header=TRUE) 
beemuse_BH_44_genome <- read.table("BeeMuse_BH_44_genome.genome", header=TRUE) 
beemuse_BH_7_19_genome <- read.table("BeeMuse_BH_7-19_genome.genome", header=TRUE) 
beemuse_BHA_2_20_genome <- read.table("BeeMuse_BHA_2-20_genome.genome", header=TRUE) 
beemuse_BLS_53_19_genome <- read.table("BeeMuse_BLS_53-19_genome.genome", header=TRUE) 
beemuse_ER_13_19_genome <- read.table("BeeMuse_ER_13-19_genome.genome", header=TRUE) 
beemuse_KBJ_1_19_genome <- read.table("BeeMuse_KBJ_1-19_genome.genome", header=TRUE) 
beemuse_KBru_6_20_genome <- read.table("BeeMuse_KBru_6-20_genome.genome", header=TRUE) 
beemuse_KLoc_37_19_genome <- read.table("BeeMuse_KLoc_37-19_genome.genome", header=TRUE) 
beemuse_KLSU_14_19_genome <- read.table("BeeMuse_KLSU_14-19_genome.genome", header=TRUE) 
beemuse_MM_31_20_genome <- read.table("BeeMuse_MM_31-20_genome.genome", header=TRUE) 
beemuse_MM_37_20_genome <- read.table("BeeMuse_MM_37-20_genome.genome", header=TRUE) 
beemuse_MP_10_20_genome <- read.table("BeeMuse_MP_10-20_genome.genome", header=TRUE) 
beemuse_PersoBC_2021_genome <- read.table("BeeMuse_PersoBC_2021_genome.genome", header=TRUE) 
beemuse_PersoJLL_2021_genome <- read.table("BeeMuse_PersoJLL_2021_genome.genome", header=TRUE) 
beemuse_PersoJLL_2022_genome <- read.table("BeeMuse_PersoJLL_2022_genome.genome", header=TRUE) 
beemuse_PersoLD_2021_genome <- read.table("BeeMuse_PersoLD_2021_genome.genome", header=TRUE) 
beemuse_PersoLD_2022_genome <- read.table("BeeMuse_PersoLD_2022_genome.genome", header=TRUE) 
beemuse_PersoUB_2021_genome <- read.table("BeeMuse_PersoUB_2021_genome.genome", header=TRUE) 
beemuse_PersoUB_2022_genome <- read.table("BeeMuse_PersoUB_2022_genome.genome", header=TRUE) 
beemuse_S_GZ_2_19_genome <- read.table("BeeMuse_S_GZ_2-19_genome.genome", header=TRUE) 
beemuse_SBJ_3_19_genome <- read.table("BeeMuse_SBJ_3-19_genome.genome", header=TRUE) 
beemuse_SJ_16_20_genome <- read.table("BeeMuse_SJ_16-20_genome.genome", header=TRUE) 
beemuse_SJ_24_20_genome <- read.table("BeeMuse_SJ_24-20_genome.genome", header=TRUE) 
beemuse_SJ_30_20_genome <- read.table("BeeMuse_SJ_30-20_genome.genome", header=TRUE) 
beemuse_TL_13_20_genome <- read.table("BeeMuse_TL_13-20_genome.genome", header=TRUE) 
beemuse_TL_19_20_genome <- read.table("BeeMuse_TL_19-20_genome.genome", header=TRUE) 
beemuse_unknown_genome <- read.table("BeeMuse_Unknown_genome.genome", header=TRUE) 

# 1
x_min <- min(c(beemuse_genome$DST, beemuse_BBS_6_19_genome$DST, beemuse_KBJ_1_19_genome$DST, beemuse_BAH_20_19_genome$DST))
x_max <- max(c(beemuse_genome$DST, beemuse_BBS_6_19_genome$DST, beemuse_KBJ_1_19_genome$DST, beemuse_BAH_20_19_genome$DST))

par(mfrow=c(2,2))

hist(beemuse_genome$DST, breaks=20, col="lightblue", main="Histogramme de la distance génétique (IBS) entre les 748 échantilons BeeMuSe", xlab="Distance génétique", cex.main=0.8, xlim=c(x_min, x_max))

hist(beemuse_BBS_6_19_genome$DST, breaks=20, col="lightgreen", main="Histogramme de la distance génétique (IBS)  - BBS_6-19", xlab="Distance génétique", cex.main=0.8, xlim=c(x_min, x_max))

hist(beemuse_KBJ_1_19_genome$DST, breaks=20, col="lightpink", main="Histogramme de la distance génétique (IBS) - KBJ_1-19", xlab="Distance génétique", cex.main=0.8, xlim=c(x_min, x_max))

hist(beemuse_BAH_20_19_genome$DST, breaks=20, col="lightyellow", main="Histogramme de la distance génétique (IBS) - BAH_20-19 ", xlab="Distance génétique", cex.main=0.8, xlim=c(x_min, x_max))

par(mfrow=c(1,1))

# 2
x_min <- min(c(beemuse_genome$DST, beemuse_BBS_6_19_genome$DST, beemuse_KBJ_1_19_genome$DST, beemuse_BAH_20_19_genome$DST))
x_max <- max(c(beemuse_genome$DST, beemuse_BBS_6_19_genome$DST, beemuse_KBJ_1_19_genome$DST, beemuse_BAH_20_19_genome$DST))

hist(beemuse_genome$DST, breaks=20, col=rgb(0,0,1,0.2), main="Histogramme de la distance génétique (IBS)", xlab="Distance génétique", cex.main=0.8, xlim=c(x_min, x_max), ylim=c(0, 25), probability=TRUE)

hist(beemuse_BBS_6_19_genome$DST, breaks=20, col=rgb(0,1,0,0.2), add=TRUE, probability=TRUE)
hist(beemuse_KBJ_1_19_genome$DST, breaks=20, col=rgb(1,0,0,0.2), add=TRUE, probability=TRUE)
hist(beemuse_BAH_20_19_genome$DST, breaks=20, col=rgb(1,1,0,0.2), add=TRUE, probability=TRUE)

legend("topleft", legend=c("BeeMuSe", "BBS_6-19", "KBJ_1-19", "BAH_20-19"), fill=c(rgb(0,0,1,0.2), rgb(0,1,0,0.2), rgb(1,0,0,0.2), rgb(1,1,0,0.2)))
```

```{r}
# 1
x_min <- min(c(beemuse_genome$PI_HAT, beemuse_BBS_6_19_genome$PI_HAT, beemuse_KBJ_1_19_genome$PI_HAT, beemuse_BAH_20_19_genome$PI_HAT))
x_max <- max(c(beemuse_genome$PI_HAT, beemuse_BBS_6_19_genome$PI_HAT, beemuse_KBJ_1_19_genome$PI_HAT, beemuse_BAH_20_19_genome$PI_HAT))

par(mfrow=c(2,2))

hist(beemuse_genome$PI_HAT, breaks=20, col="lightblue", main="Histogramme du degré d'apparentement (IBD) entre les 748 échantilons BeeMuSe", xlab="IBD", cex.main=0.8, xlim=c(x_min, x_max))

hist(beemuse_BBS_6_19_genome$PI_HAT, breaks=20, col="lightgreen", main="Histogramme du degré d'apparentement (IBD)- BBS_6-19", xlab="IBD", cex.main=0.8, xlim=c(x_min, x_max))

hist(beemuse_KBJ_1_19_genome$PI_HAT, breaks=20, col="lightpink", main="Histogramme du degré d'apparentement (IBD) - KBJ_1-19", xlab="IBD", cex.main=0.8, xlim=c(x_min, x_max))

hist(beemuse_BAH_20_19_genome$PI_HAT, breaks=20, col="lightyellow", main="Histogramme du degré d'apparentement (IBD) - BAH_20-19 ", xlab="IBD", cex.main=0.8, xlim=c(x_min, x_max))

par(mfrow=c(1,1))

# 2
x_min <- min(c(beemuse_genome$PI_HAT, beemuse_BBS_6_19_genome$PI_HAT, beemuse_KBJ_1_19_genome$PI_HAT, beemuse_BAH_20_19_genome$PI_HAT))
x_max <- max(c(beemuse_genome$PI_HAT, beemuse_BBS_6_19_genome$PI_HAT, beemuse_KBJ_1_19_genome$PI_HAT, beemuse_BAH_20_19_genome$PI_HAT))

hist(beemuse_genome$PI_HAT, breaks=20, col=rgb(0,0,1,0.2), main="Histogramme du degré d'apparentement (IBD)", xlab="IBD", cex.main=0.8, xlim=c(x_min, x_max), ylim=c(0, 25), probability=TRUE)

hist(beemuse_BBS_6_19_genome$PI_HAT, breaks=20, col=rgb(0,1,0,0.2), add=TRUE, probability=TRUE)
hist(beemuse_KBJ_1_19_genome$PI_HAT, breaks=20, col=rgb(1,0,0,0.2), add=TRUE, probability=TRUE)
hist(beemuse_BAH_20_19_genome$PI_HAT, breaks=20, col=rgb(1,1,0,0.2), add=TRUE, probability=TRUE)

legend("topleft", legend=c("BeeMuSe", "BBS_6-19", "KBJ_1-19", "BAH_20-19"), fill=c(rgb(0,0,1,0.2), rgb(0,1,0,0.2), rgb(1,0,0,0.2), rgb(1,1,0,0.2)))
```

```{r}
# Création d'un facteur pour distinguer les différentes données
beemuse_genome$Dataset <- "BeeMuSe_748"
beemuse_BBS_6_19_genome$Dataset <- "BBS_6-19"
beemuse_KBJ_1_19_genome$Dataset <- "KBJ_1-19"
beemuse_BAH_20_19_genome$Dataset <- "BAH_20-19"

#all_data <- rbind(beemuse_genome, beemuse_BBS_6_19_genome, beemuse_KBJ_1_19_genome, beemuse_BAH_20_19_genome)

#group_order <- c("BeeMuSe_748", "BBS_6-19", "BAH_20-19", "KBJ_1-19")
#group_colors <- c("BeeMuSe_748" = "lightblue", "BBS_6-19" = "lightgreen", "KBJ_1-19" = "lightpink", "BAH_20-19" = "lightyellow")

all_data <- rbind(beemuse_BBS_6_19_genome, beemuse_KBJ_1_19_genome, beemuse_BAH_20_19_genome)

group_order <- c ("KBJ_1-19" ,"BBS_6-19","BAH_20-19")
group_colors <- c("BBS_6-19" = "lightgreen", "KBJ_1-19" = "lightpink", "BAH_20-19" = "lightyellow")

all_data$Dataset <- factor(all_data$Dataset, levels = group_order)

# vertical
ggplot(all_data, aes(x = Dataset, y = DST, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = group_colors) + 
  labs(title = "Violin Plot - IBS",
       x = "", y = "IBS") +
  theme_minimal() 

ggplot(all_data, aes(x = Dataset, y = PI_HAT, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.2) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = group_colors) + 
  labs(title = "Violin Plot - IBD",
       x = "", y = "IBD") +
  theme_minimal()

#horizontal
ggplot(all_data, aes(x = Dataset, y = DST, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = group_colors) + 
  labs(title = "Violin Plot - IBS",
       x = "", y = "IBS") +
  theme_minimal() +
  coord_flip()

ggplot(all_data, aes(x = Dataset, y = PI_HAT, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.2) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = group_colors) + 
  labs(title = "Violin Plot - IBD",
       x = "", y = "IBD") +
  theme_minimal() +
  coord_flip()
```

##### 29 groups ID_2a - plink --genome --keep

```{r}
# Définir la palette de couleurs pastel
pastel_colors <- c(
    "#FF0000", "#FF3300", "#FF6600", "#FF9900", "#FFCC00",
    "#FFFF00", "#CCFF00", "#99FF00", "#66FF00", "#33FF00",
    "#00FF00", "#00FF33", "#00FF66", "#00FF99", "#00FFCC",
    "#00FFFF", "#00CCFF", "#0099FF", "#0066FF", "#0033FF",
    "#0000FF", "#3300FF", "#6600FF", "#9900FF", "#CC00FF",
    "#FF00FF", "#FF00CC", "#FF0099", "#FF0066", "#FF0033"
)

# Création d'un facteur pour distinguer les différentes données
beemuse_BAH_20_19_genome$Dataset <- "BAH_20-19"
beemuse_BBS_6_19_genome$Dataset <- "BBS_6-19"
beemuse_BER_11_19_genome$Dataset <- "BER_11-19"
beemuse_BH_44_genome$Dataset <- "BH_44"
beemuse_BH_7_19_genome$Dataset <- "BH_7-19"
beemuse_BHA_2_20_genome$Dataset <- "BHA_2-20"
beemuse_BLS_53_19_genome$Dataset <- "BLS_53-19"
beemuse_ER_13_19_genome$Dataset <- "ER_13-19"
beemuse_KBJ_1_19_genome$Dataset <- "KBJ_1-19"
beemuse_KBru_6_20_genome$Dataset <- "KBru_6-20"
beemuse_KLoc_37_19_genome$Dataset <- "KLoc_37-19"
beemuse_KLSU_14_19_genome$Dataset <- "KLSU_14-19"
beemuse_MM_31_20_genome$Dataset <- "MM_31-20"
beemuse_MM_37_20_genome$Dataset <- "MM_37-20"
beemuse_MP_10_20_genome$Dataset <- "MP_10-20"
beemuse_PersoBC_2021_genome$Dataset <- "PersoBC_2021"
beemuse_PersoJLL_2021_genome$Dataset <- "PersoJLL_2021"
beemuse_PersoJLL_2022_genome$Dataset <- "PersoJLL_2022"
beemuse_PersoLD_2021_genome$Dataset <- "PersoLD_2021"
beemuse_PersoLD_2022_genome$Dataset <- "PersoLD_2022" 
beemuse_PersoUB_2021_genome$Dataset <- "PersoUB_2021"
beemuse_PersoUB_2022_genome$Dataset <- "PersoUB_2022"
beemuse_S_GZ_2_19_genome$Dataset <- "S_GZ_2-19"
beemuse_SBJ_3_19_genome$Dataset <- "SBJ_3-19"
beemuse_SJ_16_20_genome$Dataset <- "SJ_16-20"
beemuse_SJ_24_20_genome$Dataset <- "SJ_24-20"
beemuse_SJ_30_20_genome$Dataset <- "SJ_30-20"
beemuse_TL_13_20_genome$Dataset <- "TL_13-20"
beemuse_TL_19_20_genome$Dataset <- "TL_19-20"
beemuse_unknown_genome$Dataset <- "Unknown"

# Combiner les ensembles de données
all_data <- rbind(beemuse_BAH_20_19_genome, beemuse_BBS_6_19_genome, beemuse_BER_11_19_genome, beemuse_BH_44_genome, beemuse_BH_7_19_genome, beemuse_BHA_2_20_genome, beemuse_BLS_53_19_genome, beemuse_ER_13_19_genome, beemuse_KBJ_1_19_genome, beemuse_KBru_6_20_genome, beemuse_KLoc_37_19_genome, beemuse_KLSU_14_19_genome, beemuse_MM_31_20_genome, beemuse_MM_37_20_genome, beemuse_MP_10_20_genome, beemuse_PersoBC_2021_genome, beemuse_PersoJLL_2021_genome, beemuse_PersoJLL_2022_genome, beemuse_PersoLD_2021_genome, beemuse_PersoLD_2022_genome, beemuse_PersoUB_2021_genome, beemuse_PersoUB_2022_genome, beemuse_S_GZ_2_19_genome, beemuse_SBJ_3_19_genome, beemuse_SJ_16_20_genome, beemuse_SJ_24_20_genome, beemuse_SJ_30_20_genome, beemuse_TL_13_20_genome, beemuse_TL_19_20_genome, beemuse_unknown_genome)

# Define the order of the groups and reverse it
group_order <- c(
    "BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20",
    "BLS_53-19", "ER_13-19", "KBJ_1-19", "KBru_6-20", "KLoc_37-19", 
    "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", "PersoBC_2021", 
    "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", 
    "PersoUB_2021", "PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", 
    "SJ_24-20", "SJ_30-20", "TL_13-20", "TL_19-20", "Unknown"
)

group_order <- rev(group_order)

# Convert the Dataset variable to a factor with the reversed order
all_data$Dataset <- factor(all_data$Dataset, levels = group_order)

# Créer le graphique de violon pour IBD
ggplot(all_data, aes(x = PI_HAT, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.3) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) + 
  labs(title = "Violin Plot - IBD",
       x = "IBD", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))

# Create the violin plot for IBS
ggplot(all_data, aes(x = DST, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.2) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) +  
  labs(title = "Violin Plot - IBS",
       x = "IBS", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))
```

```{r}
# x15 - 1
# Définir la palette de couleurs pastel
pastel_colors <- c(
    "#FF0000", "#FF3300", "#FF6600", "#FF9900", "#FFCC00",
    "#FFFF00", "#CCFF00", "#99FF00", "#66FF00", "#33FF00",
    "#00FF00", "#00FF33", "#00FF66", "#00FF99", "#00FFCC")

# Création d'un facteur pour distinguer les différentes données
beemuse_BAH_20_19_genome$Dataset <- "BAH_20-19"
beemuse_BBS_6_19_genome$Dataset <- "BBS_6-19"
beemuse_BER_11_19_genome$Dataset <- "BER_11-19"
beemuse_BH_44_genome$Dataset <- "BH_44"
beemuse_BH_7_19_genome$Dataset <- "BH_7-19"
beemuse_BHA_2_20_genome$Dataset <- "BHA_2-20"
beemuse_BLS_53_19_genome$Dataset <- "BLS_53-19"
beemuse_ER_13_19_genome$Dataset <- "ER_13-19"
beemuse_KBJ_1_19_genome$Dataset <- "KBJ_1-19"
beemuse_KBru_6_20_genome$Dataset <- "KBru_6-20"
beemuse_KLoc_37_19_genome$Dataset <- "KLoc_37-19"
beemuse_KLSU_14_19_genome$Dataset <- "KLSU_14-19"
beemuse_MM_31_20_genome$Dataset <- "MM_31-20"
beemuse_MM_37_20_genome$Dataset <- "MM_37-20"
beemuse_MP_10_20_genome$Dataset <- "MP_10-20"

# Combiner les ensembles de données
all_data <- rbind(beemuse_BAH_20_19_genome, beemuse_BBS_6_19_genome, beemuse_BER_11_19_genome, beemuse_BH_44_genome, beemuse_BH_7_19_genome, beemuse_BHA_2_20_genome, beemuse_BLS_53_19_genome, beemuse_ER_13_19_genome, beemuse_KBJ_1_19_genome, beemuse_KBru_6_20_genome, beemuse_KLoc_37_19_genome, beemuse_KLSU_14_19_genome, beemuse_MM_31_20_genome, beemuse_MM_37_20_genome, beemuse_MP_10_20_genome)

# Définir l'ordre des groupes
group_order <- c(
    "BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20",
    "BLS_53-19", "ER_13-19", "KBJ_1-19", "KBru_6-20", "KLoc_37-19", 
    "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20"
)

group_order <- rev(group_order)

# Convertir la variable Dataset en un facteur avec l'ordre spécifié
all_data$Dataset <- factor(all_data$Dataset, levels = group_order)

# Créer le graphique de violon pour IBS
ggplot(all_data, aes(x = DST, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.3) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) + 
  labs(title = "Violin Plot - IBS",
       x = "IBS", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))

# Créer le graphique de violon pour IBD
ggplot(all_data, aes(x = PI_HAT, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.1) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) + 
  labs(title = "Violin Plot - IBD",
       x = "IBD", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))
```

```{r}
# x15 - 2
# Définir la palette de couleurs pastel
pastel_colors <- c(
    "#00FFFF", "#00CCFF", "#0099FF", "#0066FF", "#0033FF",
    "#0000FF", "#3300FF", "#6600FF", "#9900FF", "#CC00FF",
    "#FF00FF", "#FF00CC", "#FF0099", "#FF0066", "#FF0033")

# Création d'un facteur pour distinguer les différentes données
beemuse_PersoBC_2021_genome$Dataset <- "PersoBC_2021"
beemuse_PersoJLL_2021_genome$Dataset <- "PersoJLL_2021"
beemuse_PersoJLL_2022_genome$Dataset <- "PersoJLL_2022"
beemuse_PersoLD_2021_genome$Dataset <- "PersoLD_2021"
beemuse_PersoLD_2022_genome$Dataset <- "PersoLD_2022" 
beemuse_PersoUB_2021_genome$Dataset <- "PersoUB_2021"
beemuse_PersoUB_2022_genome$Dataset <- "PersoUB_2022"
beemuse_S_GZ_2_19_genome$Dataset <- "S_GZ_2-19"
beemuse_SBJ_3_19_genome$Dataset <- "SBJ_3-19"
beemuse_SJ_16_20_genome$Dataset <- "SJ_16-20"
beemuse_SJ_24_20_genome$Dataset <- "SJ_24-20"
beemuse_SJ_30_20_genome$Dataset <- "SJ_30-20"
beemuse_TL_13_20_genome$Dataset <- "TL_13-20"
beemuse_TL_19_20_genome$Dataset <- "TL_19-20"
beemuse_unknown_genome$Dataset <- "Unknown"

# Combiner les ensembles de données
all_data <- rbind(beemuse_PersoBC_2021_genome, beemuse_PersoJLL_2021_genome, beemuse_PersoJLL_2022_genome, beemuse_PersoLD_2021_genome, beemuse_PersoLD_2022_genome, beemuse_PersoUB_2021_genome, beemuse_PersoUB_2022_genome, beemuse_S_GZ_2_19_genome, beemuse_SBJ_3_19_genome, beemuse_SJ_16_20_genome, beemuse_SJ_24_20_genome, beemuse_SJ_30_20_genome, beemuse_TL_13_20_genome, beemuse_TL_19_20_genome, beemuse_unknown_genome)

# Définir l'ordre des groupes
group_order <- c("PersoBC_2021", 
    "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", 
    "PersoUB_2021", "PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", 
    "SJ_24-20", "SJ_30-20", "TL_13-20", "TL_19-20", "Unknown"
)

group_order <- rev(group_order)

# Convertir la variable Dataset en un facteur avec l'ordre spécifié
all_data$Dataset <- factor(all_data$Dataset, levels = group_order)

# Créer le graphique de violon pour IBS
ggplot(all_data, aes(x = DST, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.4) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) + 
  labs(title = "Violin Plot - IBS",
       x = "IBS", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))

# Créer le graphique de violon pour IBD
ggplot(all_data, aes(x = PI_HAT, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.6) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) + 
  labs(title = "Violin Plot - IBD",
       x = "IBD", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))
```

##### 748 Samples

```{r}
#library(ggplot2)
#library(reshape2)
similarity_matrix <- acast(beemuse_genome, IID1 ~ IID2, value.var = "PI_HAT")

ggplot(beemuse_genome, aes(x = "", y = DST, fill = "IBS")) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = "lightblue") + 
  labs(title = "Violin Plot - IBS",
       x = "", y = "IBS") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1)) +
  coord_flip()

ggplot(beemuse_genome, aes(x = "", y = IBS0/(IBS0+IBS1+IBS2), fill = "IBS0")) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = "lightblue") + 
  labs(title = "Violin Plot - IBS0",
       x = "", y = "IBS0") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1)) +
  coord_flip()

ggplot(beemuse_genome, aes(x = "", y = PI_HAT, fill = "IBD")) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = "lightblue") + 
  labs(title = "Violin Plot - IBD",
       x = "", y = "IBD") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1)) +
  coord_flip()

# Tracer la heatmap
ggplot(melt(similarity_matrix), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap - IBD - BeeMuSe",
       x = "Individu 1",
       y = "Individu 2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calcul de la matrice de similarité
similarity_matrix <- as.matrix(similarity_matrix)

# Calcul de la distance euclidienne
distance <- dist(1 - similarity_matrix)

# Création du dendrogramme
dendrogram <- hclust(distance)
plot(dendrogram, main = "Cluster Dendrogram - IBD Similarity")

#library(MASS)
# MDS
mds <- cmdscale(distance)

# Création du plot MDS avec noms d'individus
plot(mds, xlab = "Dimension 1", ylab = "Dimension 2", main = "MDS Plot - IBD Similarity")
text(mds, labels = rownames(similarity_matrix), pos = c(3, 4), col = "black", cex = 0.6)

# Distribution de la similarité
ggplot(beemuse_genome, aes(x = PI_HAT)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  labs(title = "Distribution Plot - IBD Similarity", x = "Similarity", y = "Density")
```

##### BAH_20-19

```{r}
#BAH_20-19
#library(ggplot2)
#library(reshape2)
similarity_matrix4 <- acast(beemuse_BAH_20_19_genome, IID1 ~ IID2, value.var = "PI_HAT")

# Tracer la heatmap
ggplot(melt(similarity_matrix4), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap - IBD - BeeMuSe",
       x = "Individu 1",
       y = "Individu 2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calcul de la matrice de similarité
similarity_matrix <- as.matrix(similarity_matrix4)

# Calcul de la distance euclidienne
distance <- dist(1 - similarity_matrix)

# Création du dendrogramme
dendrogram <- hclust(distance)
plot(dendrogram, main = "Cluster Dendrogram - IBD Similarity")

#library(MASS)
# MDS
mds <- cmdscale(distance)

# Création du plot MDS avec noms d'individus
plot(mds, xlab = "Dimension 1", ylab = "Dimension 2", main = "MDS Plot - IBD Similarity")
text(mds, labels = rownames(similarity_matrix), pos = c(1, 3), col = "black", cex = 0.6)


# Distribution de la similarité
ggplot(beemuse_BAH_20_19_genome, aes(x = PI_HAT)) +
  geom_density(fill = "lightyellow", alpha = 0.7) +
  labs(title = "Distribution Plot - IBD Similarity", x = "Similarity", y = "Density")
```

##### BBS_6-19

```{r}
#library(ggplot2)
#library(reshape2)
similarity_matrix2 <- acast(beemuse_BBS_6_19_genome, IID1 ~ IID2, value.var = "PI_HAT")

# Tracer la heatmap
ggplot(melt(similarity_matrix2), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap - IBD - BeeMuSe",
       x = "Individu 1",
       y = "Individu 2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calcul de la matrice de similarité
similarity_matrix <- as.matrix(similarity_matrix2)

# Calcul de la distance euclidienne
distance <- dist(1 - similarity_matrix)

# Création du dendrogramme
dendrogram <- hclust(distance)
plot(dendrogram, main = "Cluster Dendrogram - IBD Similarity")

#library(MASS)
# MDS
mds <- cmdscale(distance)

# Création du plot MDS avec noms d'individus
plot(mds, xlab = "Dimension 1", ylab = "Dimension 2", main = "MDS Plot - IBD Similarity")
text(mds, labels = rownames(similarity_matrix), pos = c(3, 4), col = "black", cex = 0.6)


# Distribution de la similarité
ggplot(beemuse_BBS_6_19_genome, aes(x = PI_HAT)) +
  geom_density(fill = "lightgreen", alpha = 0.7) +
  labs(title = "Distribution Plot - IBD Similarity", x = "Similarity", y = "Density")
```

##### KBJ_1-19

```{r}
#KBJ_1-19
#library(ggplot2)
#library(reshape2)
similarity_matrix3 <- acast(beemuse_KBJ_1_19_genome, IID1 ~ IID2, value.var = "PI_HAT")

# Tracer la heatmap
ggplot(melt(similarity_matrix3), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap - IBD - BeeMuSe",
       x = "Individu 1",
       y = "Individu 2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calcul de la matrice de similarité
similarity_matrix <- as.matrix(similarity_matrix3)

# Calcul de la distance euclidienne
distance <- dist(1 - similarity_matrix)

# Création du dendrogramme
dendrogram <- hclust(distance)
plot(dendrogram, main = "Cluster Dendrogram - IBD Similarity")

#library(MASS)
# MDS
mds <- cmdscale(distance)

# Création du plot MDS avec noms d'individus
plot(mds, xlab = "Dimension 1", ylab = "Dimension 2", main = "MDS Plot - IBD Similarity")
text(mds, labels = rownames(similarity_matrix), pos = c(1, 3), col = "black", cex = 0.6)


# Distribution de la similarité
ggplot(beemuse_KBJ_1_19_genome, aes(x = PI_HAT)) +
  geom_density(fill = "lightpink", alpha = 0.7) +
  labs(title = "Distribution Plot - IBD Similarity", x = "Similarity", y = "Density")
```

##### 29 groups ID_2a - subset plink --genome 748 samples

```{r}
setwd("~/Documents/Stage_NB/data/IBD") 

beemuse_BAH_20_19_genome <- read.table("BAH_20-19.genome", header=TRUE) 
beemuse_BBS_6_19_genome <- read.table("BBS_6-19.genome", header=TRUE) 
beemuse_BER_11_19_genome <- read.table("BER_11-19.genome", header=TRUE) 
beemuse_BH_44_genome <- read.table("BH_44.genome", header=TRUE) 
beemuse_BH_7_19_genome <- read.table("BH_7-19.genome", header=TRUE) 
beemuse_BHA_2_20_genome <- read.table("BHA_2-20.genome", header=TRUE) 
beemuse_BLS_53_19_genome <- read.table("BLS_53-19.genome", header=TRUE) 
beemuse_ER_13_19_genome <- read.table("ER_13-19.genome", header=TRUE) 
beemuse_KBJ_1_19_genome <- read.table("KBJ_1-19.genome", header=TRUE) 
beemuse_KBru_6_20_genome <- read.table("KBru_6-20.genome", header=TRUE) 
beemuse_KLoc_37_19_genome <- read.table("KLoc_37-19.genome", header=TRUE) 
beemuse_KLSU_14_19_genome <- read.table("KLSU_14-19.genome", header=TRUE) 
beemuse_MM_31_20_genome <- read.table("MM_31-20.genome", header=TRUE) 
beemuse_MM_37_20_genome <- read.table("MM_37-20.genome", header=TRUE) 
beemuse_MP_10_20_genome <- read.table("MP_10-20.genome", header=TRUE) 
beemuse_PersoBC_2021_genome <- read.table("PersoBC_2021.genome", header=TRUE) 
beemuse_PersoJLL_2021_genome <- read.table("PersoJLL_2021.genome", header=TRUE) 
beemuse_PersoJLL_2022_genome <- read.table("PersoJLL_2022.genome", header=TRUE) 
beemuse_PersoLD_2021_genome <- read.table("PersoLD_2021.genome", header=TRUE) 
beemuse_PersoLD_2022_genome <- read.table("PersoLD_2022.genome", header=TRUE) 
beemuse_PersoUB_2021_genome <- read.table("PersoUB_2021.genome", header=TRUE) 
beemuse_PersoUB_2022_genome <- read.table("PersoUB_2022.genome", header=TRUE) 
beemuse_S_GZ_2_19_genome <- read.table("S_GZ_2-19.genome", header=TRUE) 
beemuse_SBJ_3_19_genome <- read.table("SBJ_3-19.genome", header=TRUE) 
beemuse_SJ_16_20_genome <- read.table("SJ_16-20.genome", header=TRUE) 
beemuse_SJ_24_20_genome <- read.table("SJ_24-20.genome", header=TRUE) 
beemuse_SJ_30_20_genome <- read.table("SJ_30-20.genome", header=TRUE) 
beemuse_TL_13_20_genome <- read.table("TL_13-20.genome", header=TRUE) 
beemuse_TL_19_20_genome <- read.table("TL_19-20.genome", header=TRUE) 
beemuse_unknown_genome <- read.table("Unknown.genome", header=TRUE) 
```

```{r}
# Définir la palette de couleurs pastel
pastel_colors <- c(
    "#FF0000", "#FF3300", "#FF6600", "#FF9900", "#FFCC00",
    "#FFFF00", "#CCFF00", "#99FF00", "#66FF00", "#33FF00",
    "#00FF00", "#00FF33", "#00FF66", "#00FF99", "#00FFCC",
    "#00FFFF", "#00CCFF", "#0099FF", "#0066FF", "#0033FF",
    "#0000FF", "#3300FF", "#6600FF", "#9900FF", "#CC00FF",
    "#FF00FF", "#FF00CC", "#FF0099", "#FF0066"
)

# Création d'un facteur pour distinguer les différentes données
beemuse_BAH_20_19_genome$Dataset <- "BAH_20-19"
beemuse_BBS_6_19_genome$Dataset <- "BBS_6-19"
beemuse_BER_11_19_genome$Dataset <- "BER_11-19"
beemuse_BH_44_genome$Dataset <- "BH_44"
beemuse_BH_7_19_genome$Dataset <- "BH_7-19"
beemuse_BHA_2_20_genome$Dataset <- "BHA_2-20"
beemuse_BLS_53_19_genome$Dataset <- "BLS_53-19"
beemuse_ER_13_19_genome$Dataset <- "ER_13-19"
beemuse_KBJ_1_19_genome$Dataset <- "KBJ_1-19"
beemuse_KBru_6_20_genome$Dataset <- "KBru_6-20"
beemuse_KLoc_37_19_genome$Dataset <- "KLoc_37-19"
beemuse_KLSU_14_19_genome$Dataset <- "KLSU_14-19"
beemuse_MM_31_20_genome$Dataset <- "MM_31-20"
beemuse_MM_37_20_genome$Dataset <- "MM_37-20"
beemuse_MP_10_20_genome$Dataset <- "MP_10-20"
beemuse_PersoBC_2021_genome$Dataset <- "PersoBC_2021"
beemuse_PersoJLL_2021_genome$Dataset <- "PersoJLL_2021"
beemuse_PersoJLL_2022_genome$Dataset <- "PersoJLL_2022"
beemuse_PersoLD_2021_genome$Dataset <- "PersoLD_2021"
beemuse_PersoLD_2022_genome$Dataset <- "PersoLD_2022" 
beemuse_PersoUB_2021_genome$Dataset <- "PersoUB_2021"
beemuse_PersoUB_2022_genome$Dataset <- "PersoUB_2022"
beemuse_S_GZ_2_19_genome$Dataset <- "S_GZ_2-19"
beemuse_SBJ_3_19_genome$Dataset <- "SBJ_3-19"
beemuse_SJ_16_20_genome$Dataset <- "SJ_16-20"
beemuse_SJ_24_20_genome$Dataset <- "SJ_24-20"
beemuse_SJ_30_20_genome$Dataset <- "SJ_30-20"
beemuse_TL_13_20_genome$Dataset <- "TL_13-20"
beemuse_TL_19_20_genome$Dataset <- "TL_19-20"
#beemuse_unknown_genome$Dataset <- "Unknown"

# Combiner les ensembles de données
all_data <- rbind(beemuse_BAH_20_19_genome, beemuse_BBS_6_19_genome, beemuse_BER_11_19_genome, beemuse_BH_44_genome, beemuse_BH_7_19_genome, beemuse_BHA_2_20_genome, beemuse_BLS_53_19_genome, beemuse_ER_13_19_genome, beemuse_KBJ_1_19_genome, beemuse_KBru_6_20_genome, beemuse_KLoc_37_19_genome, beemuse_KLSU_14_19_genome, beemuse_MM_31_20_genome, beemuse_MM_37_20_genome, beemuse_MP_10_20_genome, beemuse_PersoBC_2021_genome, beemuse_PersoJLL_2021_genome, beemuse_PersoJLL_2022_genome, beemuse_PersoLD_2021_genome, beemuse_PersoLD_2022_genome, beemuse_PersoUB_2021_genome, beemuse_PersoUB_2022_genome, beemuse_S_GZ_2_19_genome, beemuse_SBJ_3_19_genome, beemuse_SJ_16_20_genome, beemuse_SJ_24_20_genome, beemuse_SJ_30_20_genome, beemuse_TL_13_20_genome, beemuse_TL_19_20_genome)

# Define the order of the groups and reverse it
group_order <- c(
    "BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20",
    "BLS_53-19", "ER_13-19", "KBJ_1-19", "KBru_6-20", "KLoc_37-19", 
    "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", "PersoBC_2021", 
    "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", 
    "PersoUB_2021", "PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", 
    "SJ_24-20", "SJ_30-20", "TL_13-20", "TL_19-20"
)

group_order <- rev(group_order)

# Convert the Dataset variable to a factor with the reversed order
all_data$Dataset <- factor(all_data$Dataset, levels = group_order)

# Créer le graphique de violon pour IBD
ggplot(all_data, aes(x = PI_HAT, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.1) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) + 
  labs(title = "Violin Plot - IBD",
       x = "IBD", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))

# Create the violin plot for IBS
ggplot(all_data, aes(x = DST, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.2) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) +  
  labs(title = "Violin Plot - IBS",
       x = "IBS", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))
```

```{r}
# x15 - 1
# Définir la palette de couleurs pastel
pastel_colors <- c(
    "#FF0000", "#FF3300", "#FF6600", "#FF9900", "#FFCC00",
    "#FFFF00", "#CCFF00", "#99FF00", "#66FF00", "#33FF00",
    "#00FF00", "#00FF33", "#00FF66", "#00FF99", "#00FFCC")

# Création d'un facteur pour distinguer les différentes données
beemuse_BAH_20_19_genome$Dataset <- "BAH_20-19"
beemuse_BBS_6_19_genome$Dataset <- "BBS_6-19"
beemuse_BER_11_19_genome$Dataset <- "BER_11-19"
beemuse_BH_44_genome$Dataset <- "BH_44"
beemuse_BH_7_19_genome$Dataset <- "BH_7-19"
beemuse_BHA_2_20_genome$Dataset <- "BHA_2-20"
beemuse_BLS_53_19_genome$Dataset <- "BLS_53-19"
beemuse_ER_13_19_genome$Dataset <- "ER_13-19"
beemuse_KBJ_1_19_genome$Dataset <- "KBJ_1-19"
beemuse_KBru_6_20_genome$Dataset <- "KBru_6-20"
beemuse_KLoc_37_19_genome$Dataset <- "KLoc_37-19"
beemuse_KLSU_14_19_genome$Dataset <- "KLSU_14-19"
beemuse_MM_31_20_genome$Dataset <- "MM_31-20"
beemuse_MM_37_20_genome$Dataset <- "MM_37-20"
beemuse_MP_10_20_genome$Dataset <- "MP_10-20"

# Combiner les ensembles de données
all_data <- rbind(beemuse_BAH_20_19_genome, beemuse_BBS_6_19_genome, beemuse_BER_11_19_genome, beemuse_BH_44_genome, beemuse_BH_7_19_genome, beemuse_BHA_2_20_genome, beemuse_BLS_53_19_genome, beemuse_ER_13_19_genome, beemuse_KBJ_1_19_genome, beemuse_KBru_6_20_genome, beemuse_KLoc_37_19_genome, beemuse_KLSU_14_19_genome, beemuse_MM_31_20_genome, beemuse_MM_37_20_genome, beemuse_MP_10_20_genome)

# Définir l'ordre des groupes
group_order <- c(
    "BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20",
    "BLS_53-19", "ER_13-19", "KBJ_1-19", "KBru_6-20", "KLoc_37-19", 
    "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20"
)

group_order <- rev(group_order)

# Convertir la variable Dataset en un facteur avec l'ordre spécifié
all_data$Dataset <- factor(all_data$Dataset, levels = group_order)

# Créer le graphique de violon pour IBS
ggplot(all_data, aes(x = DST, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.3) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) + 
  labs(title = "Violin Plot - IBS",
       x = "IBS", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))

# Créer le graphique de violon pour IBD
ggplot(all_data, aes(x = PI_HAT, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.3) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) + 
  labs(title = "Violin Plot - IBD",
       x = "IBD", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))

# x15 - 2
# Définir la palette de couleurs pastel
pastel_colors <- c(
    "#00FFFF", "#00CCFF", "#0099FF", "#0066FF", "#0033FF",
    "#0000FF", "#3300FF", "#6600FF", "#9900FF", "#CC00FF",
    "#FF00FF", "#FF00CC", "#FF0099", "#FF0066", "#FF0033")

# Création d'un facteur pour distinguer les différentes données
beemuse_PersoBC_2021_genome$Dataset <- "PersoBC_2021"
beemuse_PersoJLL_2021_genome$Dataset <- "PersoJLL_2021"
beemuse_PersoJLL_2022_genome$Dataset <- "PersoJLL_2022"
beemuse_PersoLD_2021_genome$Dataset <- "PersoLD_2021"
beemuse_PersoLD_2022_genome$Dataset <- "PersoLD_2022" 
beemuse_PersoUB_2021_genome$Dataset <- "PersoUB_2021"
beemuse_PersoUB_2022_genome$Dataset <- "PersoUB_2022"
beemuse_S_GZ_2_19_genome$Dataset <- "S_GZ_2-19"
beemuse_SBJ_3_19_genome$Dataset <- "SBJ_3-19"
beemuse_SJ_16_20_genome$Dataset <- "SJ_16-20"
beemuse_SJ_24_20_genome$Dataset <- "SJ_24-20"
beemuse_SJ_30_20_genome$Dataset <- "SJ_30-20"
beemuse_TL_13_20_genome$Dataset <- "TL_13-20"
beemuse_TL_19_20_genome$Dataset <- "TL_19-20"
beemuse_unknown_genome$Dataset <- "Unknown"

# Combiner les ensembles de données
all_data <- rbind(beemuse_PersoBC_2021_genome, beemuse_PersoJLL_2021_genome, beemuse_PersoJLL_2022_genome, beemuse_PersoLD_2021_genome, beemuse_PersoLD_2022_genome, beemuse_PersoUB_2021_genome, beemuse_PersoUB_2022_genome, beemuse_S_GZ_2_19_genome, beemuse_SBJ_3_19_genome, beemuse_SJ_16_20_genome, beemuse_SJ_24_20_genome, beemuse_SJ_30_20_genome, beemuse_TL_13_20_genome, beemuse_TL_19_20_genome, beemuse_unknown_genome)

# Définir l'ordre des groupes
group_order <- c("PersoBC_2021", 
    "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", 
    "PersoUB_2021", "PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", 
    "SJ_24-20", "SJ_30-20", "TL_13-20", "TL_19-20","Unknown"
)

group_order <- rev(group_order)

# Convertir la variable Dataset en un facteur avec l'ordre spécifié
all_data$Dataset <- factor(all_data$Dataset, levels = group_order)

# Créer le graphique de violon pour IBS
ggplot(all_data, aes(x = DST, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.4) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) + 
  labs(title = "Violin Plot - IBS",
       x = "IBS", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))

# Créer le graphique de violon pour IBD
ggplot(all_data, aes(x = PI_HAT, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.2) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) + 
  labs(title = "Violin Plot - IBD",
       x = "IBD", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))
```

```{r}
#library(ggplot2)
#library(reshape2)
similarity_matrix2 <- acast(beemuse_BBS_6_19_genome, IID1 ~ IID2, value.var = "PI_HAT")

# Tracer la heatmap
ggplot(melt(similarity_matrix2), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap - IBD - BeeMuSe",
       x = "Individu 1",
       y = "Individu 2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calcul de la matrice de similarité
similarity_matrix <- as.matrix(similarity_matrix2)

# Calcul de la distance euclidienne
distance <- dist(1 - similarity_matrix)

# Création du dendrogramme
dendrogram <- hclust(distance)
plot(dendrogram, main = "Cluster Dendrogram - IBD Similarity")

#library(MASS)
# MDS
mds <- cmdscale(distance)

# Création du plot MDS avec noms d'individus
plot(mds, xlab = "Dimension 1", ylab = "Dimension 2", main = "MDS Plot - IBD Similarity")
text(mds, labels = rownames(similarity_matrix), pos = c(3, 4), col = "black", cex = 0.6)


# Distribution de la similarité
ggplot(beemuse_BBS_6_19_genome, aes(x = PI_HAT)) +
  geom_density(fill = "lightgreen", alpha = 0.7) +
  labs(title = "Distribution Plot - IBD Similarity", x = "Similarity", y = "Density")
```

```{r}
#KBJ_1-19
#library(ggplot2)
#library(reshape2)
similarity_matrix3 <- acast(beemuse_KBJ_1_19_genome, IID1 ~ IID2, value.var = "PI_HAT")

# Tracer la heatmap
ggplot(melt(similarity_matrix3), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap - IBD - BeeMuSe",
       x = "Individu 1",
       y = "Individu 2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calcul de la matrice de similarité
similarity_matrix <- as.matrix(similarity_matrix3)

# Calcul de la distance euclidienne
distance <- dist(1 - similarity_matrix)

# Création du dendrogramme
dendrogram <- hclust(distance)
plot(dendrogram, main = "Cluster Dendrogram - IBD Similarity")

#library(MASS)
# MDS
mds <- cmdscale(distance)

# Création du plot MDS avec noms d'individus
plot(mds, xlab = "Dimension 1", ylab = "Dimension 2", main = "MDS Plot - IBD Similarity")
text(mds, labels = rownames(similarity_matrix), pos = c(1, 3), col = "black", cex = 0.6)


# Distribution de la similarité
ggplot(beemuse_KBJ_1_19_genome, aes(x = PI_HAT)) +
  geom_density(fill = "lightpink", alpha = 0.7) +
  labs(title = "Distribution Plot - IBD Similarity", x = "Similarity", y = "Density")
```

```{r}
#library(ggplot2)
#library(reshape2)
similarity_matrix2 <- acast(beemuse_MM_31_20_genome, IID1 ~ IID2, value.var = "PI_HAT")

# Tracer la heatmap
ggplot(melt(similarity_matrix2), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap - IBD - BeeMuSe",
       x = "Individu 1",
       y = "Individu 2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calcul de la matrice de similarité
similarity_matrix <- as.matrix(similarity_matrix2)

# Calcul de la distance euclidienne
distance <- dist(1 - similarity_matrix)

# Création du dendrogramme
dendrogram <- hclust(distance)
plot(dendrogram, main = "Cluster Dendrogram - IBD Similarity")

#library(MASS)
# MDS
mds <- cmdscale(distance)

# Création du plot MDS avec noms d'individus
plot(mds, xlab = "Dimension 1", ylab = "Dimension 2", main = "MDS Plot - IBD Similarity")
text(mds, labels = rownames(similarity_matrix), pos = c(3, 4), col = "black", cex = 0.6)


# Distribution de la similarité
ggplot(beemuse_BBS_6_19_genome, aes(x = PI_HAT)) +
  geom_density(fill = "lightyellow", alpha = 0.7) +
  labs(title = "Distribution Plot - IBD Similarity", x = "Similarity", y = "Density")
```

#### PLINK2 --make-king-table

```{r}
setwd("~/Documents/Stage_NB/data/IBD") 

beemuse_p2_genome <- read.table("BeeMuse_plink2_genome.kin0", header=FALSE, sep="\t", 
                   col.names=c("FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP"))

ggplot(beemuse_p2_genome, aes(x = "", y = IBS0, fill = "IBS0")) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = "lightblue") + 
  labs(title = "Violin Plot - IBS0",
       x = "", y = "IBS0") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1)) +
  coord_flip()

ggplot(beemuse_p2_genome, aes(x = "", y = KINSHIP, fill = "KINSHIP")) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = "lightblue") + 
  labs(title = "Violin Plot - KINSHIP",
       x = "", y = "KINSHIP") +
  theme_minimal() +
  scale_y_continuous(limits = c(-2, 0.5)) +
  coord_flip()

#library(ggplot2)
#library(reshape2)
similarity_matrix <- acast(beemuse_p2_genome, IID2 ~ IID1, value.var = "KINSHIP")

# Tracer la heatmap
ggplot(melt(similarity_matrix), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap - KINSHIP - BeeMuSe",
       x = "Individu 2",
       y = "Individu 1") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calcul de la matrice de similarité
similarity_matrix <- as.matrix(similarity_matrix)

# Calcul de la distance euclidienne
distance <- dist(1 - similarity_matrix)

# Création du dendrogramme
dendrogram <- hclust(distance)
plot(dendrogram, main = "Cluster Dendrogram - KINSHIP Similarity")

#library(MASS)
# MDS
mds <- cmdscale(distance)

# Création du plot MDS avec noms d'individus
plot(mds, xlab = "Dimension 1", ylab = "Dimension 2", main = "MDS Plot - KINSHIP Similarity")
text(mds, labels = rownames(similarity_matrix), pos = c(3, 4), col = "black", cex = 0.6)


# Distribution de la similarité
ggplot(beemuse_p2_genome, aes(x = KINSHIP)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  labs(title = "Distribution Plot - KINSHIP Similarity", x = "Similarity", y = "Density")
```

##### 29 groups ID_2a - plink2 --make-king-table

```{r}
setwd("~/Documents/Stage_NB/data/IBD") 

beemuse_BAH_20_19_genome <- read.table("BAH_20-19_plink2.genome", header=TRUE) 
beemuse_BBS_6_19_genome <- read.table("BBS_6-19_plink2.genome", header=TRUE) 
beemuse_BER_11_19_genome <- read.table("BER_11-19_plink2.genome", header=TRUE) 
beemuse_BH_44_genome <- read.table("BH_44_plink2.genome", header=TRUE) 
beemuse_BH_7_19_genome <- read.table("BH_7-19_plink2.genome", header=TRUE) 
beemuse_BHA_2_20_genome <- read.table("BHA_2-20_plink2.genome", header=TRUE) 
beemuse_BLS_53_19_genome <- read.table("BLS_53-19_plink2.genome", header=TRUE) 
beemuse_ER_13_19_genome <- read.table("ER_13-19_plink2.genome", header=TRUE) 
beemuse_KBJ_1_19_genome <- read.table("KBJ_1-19_plink2.genome", header=TRUE) 
beemuse_KBru_6_20_genome <- read.table("KBru_6-20_plink2.genome", header=TRUE) 
beemuse_KLoc_37_19_genome <- read.table("KLoc_37-19_plink2.genome", header=TRUE) 
beemuse_KLSU_14_19_genome <- read.table("KLSU_14-19_plink2.genome", header=TRUE) 
beemuse_MM_31_20_genome <- read.table("MM_31-20_plink2.genome", header=TRUE) 
beemuse_MM_37_20_genome <- read.table("MM_37-20_plink2.genome", header=TRUE) 
beemuse_MP_10_20_genome <- read.table("MP_10-20_plink2.genome", header=TRUE) 
beemuse_PersoBC_2021_genome <- read.table("PersoBC_2021_plink2.genome", header=TRUE) 
beemuse_PersoJLL_2021_genome <- read.table("PersoJLL_2021_plink2.genome", header=TRUE) 
beemuse_PersoJLL_2022_genome <- read.table("PersoJLL_2022_plink2.genome", header=TRUE) 
beemuse_PersoLD_2021_genome <- read.table("PersoLD_2021_plink2.genome", header=TRUE) 
beemuse_PersoLD_2022_genome <- read.table("PersoLD_2022_plink2.genome", header=TRUE) 
beemuse_PersoUB_2021_genome <- read.table("PersoUB_2021_plink2.genome", header=TRUE) 
beemuse_PersoUB_2022_genome <- read.table("PersoUB_2022_plink2.genome", header=TRUE) 
beemuse_S_GZ_2_19_genome <- read.table("S_GZ_2-19_plink2.genome", header=TRUE) 
beemuse_SBJ_3_19_genome <- read.table("SBJ_3-19_plink2.genome", header=TRUE) 
beemuse_SJ_16_20_genome <- read.table("SJ_16-20_plink2.genome", header=TRUE) 
beemuse_SJ_24_20_genome <- read.table("SJ_24-20_plink2.genome", header=TRUE) 
beemuse_SJ_30_20_genome <- read.table("SJ_30-20_plink2.genome", header=TRUE) 
beemuse_TL_13_20_genome <- read.table("TL_13-20_plink2.genome", header=TRUE) 
beemuse_TL_19_20_genome <- read.table("TL_19-20_plink2.genome", header=TRUE) 
beemuse_unknown_genome <- read.table("Unknown_plink2.genome", header=TRUE) 
```

```{r}
# Définir la palette de couleurs pastel
pastel_colors <- c(
    "#FF0000", "#FF3300", "#FF6600", "#FF9900", "#FFCC00",
    "#FFFF00", "#CCFF00", "#99FF00", "#66FF00", "#33FF00",
    "#00FF00", "#00FF33", "#00FF66", "#00FF99", "#00FFCC",
    "#00FFFF", "#00CCFF", "#0099FF", "#0066FF", "#0033FF",
    "#0000FF", "#3300FF", "#6600FF", "#9900FF", "#CC00FF",
    "#FF00FF", "#FF00CC", "#FF0099", "#FF0066", "#FF0033"
)

# Création d'un facteur pour distinguer les différentes données
beemuse_BAH_20_19_genome$Dataset <- "BAH_20-19"
beemuse_BBS_6_19_genome$Dataset <- "BBS_6-19"
beemuse_BER_11_19_genome$Dataset <- "BER_11-19"
beemuse_BH_44_genome$Dataset <- "BH_44"
beemuse_BH_7_19_genome$Dataset <- "BH_7-19"
beemuse_BHA_2_20_genome$Dataset <- "BHA_2-20"
beemuse_BLS_53_19_genome$Dataset <- "BLS_53-19"
beemuse_ER_13_19_genome$Dataset <- "ER_13-19"
beemuse_KBJ_1_19_genome$Dataset <- "KBJ_1-19"
beemuse_KBru_6_20_genome$Dataset <- "KBru_6-20"
beemuse_KLoc_37_19_genome$Dataset <- "KLoc_37-19"
beemuse_KLSU_14_19_genome$Dataset <- "KLSU_14-19"
beemuse_MM_31_20_genome$Dataset <- "MM_31-20"
beemuse_MM_37_20_genome$Dataset <- "MM_37-20"
beemuse_MP_10_20_genome$Dataset <- "MP_10-20"
beemuse_PersoBC_2021_genome$Dataset <- "PersoBC_2021"
beemuse_PersoJLL_2021_genome$Dataset <- "PersoJLL_2021"
beemuse_PersoJLL_2022_genome$Dataset <- "PersoJLL_2022"
beemuse_PersoLD_2021_genome$Dataset <- "PersoLD_2021"
beemuse_PersoLD_2022_genome$Dataset <- "PersoLD_2022" 
beemuse_PersoUB_2021_genome$Dataset <- "PersoUB_2021"
beemuse_PersoUB_2022_genome$Dataset <- "PersoUB_2022"
beemuse_S_GZ_2_19_genome$Dataset <- "S_GZ_2-19"
beemuse_SBJ_3_19_genome$Dataset <- "SBJ_3-19"
beemuse_SJ_16_20_genome$Dataset <- "SJ_16-20"
beemuse_SJ_24_20_genome$Dataset <- "SJ_24-20"
beemuse_SJ_30_20_genome$Dataset <- "SJ_30-20"
beemuse_TL_13_20_genome$Dataset <- "TL_13-20"
beemuse_TL_19_20_genome$Dataset <- "TL_19-20"
beemuse_unknown_genome$Dataset <- "Unknown"

# Combiner les ensembles de données
all_data <- rbind(beemuse_BAH_20_19_genome, beemuse_BBS_6_19_genome, beemuse_BER_11_19_genome, beemuse_BH_44_genome, beemuse_BH_7_19_genome, beemuse_BHA_2_20_genome, beemuse_BLS_53_19_genome, beemuse_ER_13_19_genome, beemuse_KBJ_1_19_genome, beemuse_KBru_6_20_genome, beemuse_KLoc_37_19_genome, beemuse_KLSU_14_19_genome, beemuse_MM_31_20_genome, beemuse_MM_37_20_genome, beemuse_MP_10_20_genome, beemuse_PersoBC_2021_genome, beemuse_PersoJLL_2021_genome, beemuse_PersoJLL_2022_genome, beemuse_PersoLD_2021_genome, beemuse_PersoLD_2022_genome, beemuse_PersoUB_2021_genome, beemuse_PersoUB_2022_genome, beemuse_S_GZ_2_19_genome, beemuse_SBJ_3_19_genome, beemuse_SJ_16_20_genome, beemuse_SJ_24_20_genome, beemuse_SJ_30_20_genome, beemuse_TL_13_20_genome, beemuse_TL_19_20_genome, beemuse_unknown_genome)

# Define the order of the groups and reverse it
group_order <- c(
    "BAH_20-19", "BBS_6-19", "BER_11-19", "BH_44", "BH_7-19", "BHA_2-20",
    "BLS_53-19", "ER_13-19", "KBJ_1-19", "KBru_6-20", "KLoc_37-19", 
    "KLSU_14-19", "MM_31-20", "MM_37-20", "MP_10-20", "PersoBC_2021", 
    "PersoJLL_2021", "PersoJLL_2022", "PersoLD_2021", "PersoLD_2022", 
    "PersoUB_2021", "PersoUB_2022", "S_GZ_2-19", "SBJ_3-19", "SJ_16-20", 
    "SJ_24-20", "SJ_30-20", "TL_13-20", "TL_19-20", "Unknown"
)

group_order <- rev(group_order)

# Convert the Dataset variable to a factor with the reversed order
all_data$Dataset <- factor(all_data$Dataset, levels = group_order)

ggplot(all_data, aes(x = KINSHIP, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.3) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) + 
  labs(title = "Violin Plot - KINSHIP",
       x = "KINSHIP", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(-0.5, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(-0.5, 1))

ggplot(all_data, aes(x = IBS0, y = Dataset, fill = Dataset)) +
  geom_violin(trim = FALSE, color = "black", width = 1.2) +
  geom_boxplot(width = 0.1, fill = alpha("white", 0), color = "black", position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = pastel_colors, limits = rev(group_order)) +  
  labs(title = "Violin Plot - IBS0",
       x = "IBS0", y = "") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  coord_cartesian(xlim = c(0, 1))
```

#### KING -b \_ --related

```{r}
setwd("~/Documents/Stage_NB/data/IBD") 

# Charger les données du fichier .kin
data <- read.table("votre_fichier.kin", header=TRUE)

# Filtrer les paires d'individus appartenant à des familles différentes
pairs_diff_families <- subset(data, FID1 != FID2)

# Comparer les coefficients de parenté entre les familles différentes
boxplot(Kinship ~ interaction(FID1, FID2), data=pairs_diff_families, xlab="Familles", ylab="Coefficient de parenté KING", main="Comparaison des coefficients de parenté entre familles")
```
