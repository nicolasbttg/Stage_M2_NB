Certains scripts et commandes sont issus depuis le GitHub pour [SeqApiPop](https://github.com/avignal5/SeqApiPop)

(Alain Vignal, David Wragg, Sonia E. Eynard, Benjamin Basso, Kamila Canale-Tabet, Emmanuelle Labarthe, Olivier Bouchez, Kaspar Bienefeld, Małgorzata Bieńkowska, Cecilia Costa, Aleš Gregorc, Per Kryger, Melanie Parejo, M. Alice Pinto, Jean-Pierre Bidanel, Bertrand Servin, & Yves Le Conte. (2021). Data from: Complex population structure and haplotype patterns in Western Europe honey bee from sequencing a large panel of haploid drones [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5592452)

Les scripts R pour réaliser certaines analyses et obtenir les figures sont accessibles au format Rmarkdown :
- [1_ACP](Stage_M2_NB_1_ACP.md)
- [2_CV_Admixture](Stage_M2_NB_2_CV_Admixture.md)
- [3_FST](Stage_M2_NB_3_FST.md)
- [4_IBS_IBD_KINSHIP](Stage_M2_NB_4_IBS_IBD_KINSHIP.md)


# BeeMuSe

```
wc -l E756_BeeMuSe.map
```
Le jeu de données BeeMuSe correspond à 12000 SNPs pour 748 échantillons

## Conversion au format plink
- convertnumchrmap.sh

```
#!/bin/bash
# Fichier de référence pour faire correspondre les ID de chromosome aux numéros
chromosome_list="HAv3_1_Chromosomes.list"

# Fichier d'entrée et de sortie
input_map="E756_BeeMuSe.map"
output_map="E756_BeeMuSe_num_chr.map"

# Créer un tableau associatif pour stocker la correspondance entre les ID de chromosome et les numéros
declare -A chromosome_mapping

# Lire le fichier de la liste des chromosomes et remplir le tableau
while read -r line; do
    chromosome_id=$(echo "$line" | awk '{print $1}')
    chromosome_number=$(echo "$line" | awk '{print $2}')
    chromosome_mapping["$chromosome_id"]=$chromosome_number
done < "$chromosome_list"

# Lire le fichier d'entrée et remplacer les ID de chromosome
while read -r line; do
    # Vérifier si la ligne est une ligne de commentaires (commençant par ##)
    if [[ $line == "##"* ]]; then
        echo "$line"  # Laisser les lignes de commentaires inchangées
    else
        # Extraire l'ID de chromosome de la ligne
        chromosome_id=$(echo "$line" | awk '{print $1}')
        rest_of_line=$(echo "$line" | cut -f2-)

        # Vérifier si l'ID de chromosome existe dans la correspondance
        if [ -n "${chromosome_mapping[$chromosome_id]}" ]; then
            chromosome_number="${chromosome_mapping[$chromosome_id]}"
            echo -e "$chromosome_number\t$rest_of_line"
        else
            # Si l'ID de chromosome n'a pas de correspondance, laisser la ligne inchangée
            echo "$line"
        fi
    fi
done < "$input_map" > "$output_map"
```

## Conversion en fichiers .bim, .bed et .fam

-makebedbeemuse.sh
```
#!/bin/bash
module load bioinfo/PLINK/1.90b7

plink --file E756_BeeMuSe_num_chr --make-bed --no-parents --no-sex --no-pheno

less -S +376 E756_BeeMuSe_num_chr.ped
```
2 colonnes puis 1 à cause de '-' au lieu de ' '
```
sed 's/-/ /g' E756_BeeMuSe_num_chr.ped > E756_BeeMuSe_num_chr_2col.ped

module load bioinfo/PLINK/1.90b7
plink --file E756_BeeMuSe_num_chr_2col --make-bed --no-parents --no-sex --no-pheno -out BeeMuse
```
## ACP
- acpbeemuse.sh
```
#!/bin/bash
module load bioinfo/PLINK/2.00a4

plink2 --bfile BeeMuse --make-rel square --nonfounders
plink2 --bfile BeeMuse --pca  --nonfounders
```
# SeqApiPop
## Obtention des données
- Récupérer les fichiers de référence SeqApiPop (Wragg et al., 2021), avec le fichier VCF de 7 023 689 SNPs des 870 échantillons d'abeilles domestiques et la liste des 629 échantillons de référence pour les analyses de structure des populations à extraire depuis [zenodo](https://zenodo.org/records/5592452)
> Diversity_Study_629_Samples.txt

> MetaGenotypesCalled870_raw_snps_allfilter.vcf.gz

(Wragg D, Eynard SE, Basso B, Canale-Tabet K, Labarthe E, Bouchez O, et al. Complex population structure and haplotype patterns in Western Europe honey bee from sequencing a large panel of haploid drones [Internet]. Genetics; 2021 Sep. Available from: http://biorxiv.org/lookup/doi/10.1101/2021.09.20.460798)
```
wget https://zenodo.org/record/5592452/files/Diversity_Study_629_Samples.txt
wget https://zenodo.org/record/5592452/files/MetaGenotypesCalled870_raw_snps_allfilter.vcf.gz
```
```
zcat MetaGenotypesCalled870_raw_snps_allfilter.vcf.gz | head
```

## Conversion du fichier VCF au format plink (.bed, .bim, .fam)

- vcftoplink.sh
```
#!/bin/bash
#SBATCH --mem=8G

module load bioinfo/Bcftools/1.9
module load bioinfo/PLINK/2.00a4

bcftools norm -m-any --output-type z -o subset_501_samples_ref_split.vcf.gz subset_501_samples_ref.vcf
plink2 --vcf subset_501_samples_ref_split.vcf.gz --make-bed --out subset_501_samples_ref --allow-extra-chr

echo "done"
```

## ACP

```
#!/bin/bash
module load bioinfo/PLINK/2.00a4

plink2 --bfile subset_501_samples_ref --make-rel square --nonfounders --out subset_501_samples_ref
plink2 --bfile subset_501_samples_ref --pca --nonfounders --out subset_501_samples_ref

echo "done"
```
sbatch acp501samples.bash
squeue -u nbettembour

Refaire la même chose seulement reference populations
```
filtered_data <- subset(eigenvec_501_seq_api_labels, Label != 'Ariege Conservatory' & Label != 'Brittany Conservatory' & UniqueInHive != 'Unknown' & UniqueInHive != 'Buckfast' & GeneticOrigin != 'Buckfast')
```
lien E756_BeeMuSe_num_chr_2col.map et list_markers_ID_to_keep.csv

.map : Marker ID ($2)
.csv : probeset_id ($1)

puis refaire ACP - beemuse 
```
awk -F';' 'NR>1 {print $1}' list_markers_ID_to_keep.csv > list_markers_ID_to_keep.txt
```
```
plink --file E756_BeeMuSe_num_chr_2col --make-bed --no-parents --no-sex --no-pheno -out BeeMuse
```
- extractlistmarkersID.bash
```
module load bioinfo/PLINK/1.90b7
# Supprimer les marqueurs indésirables du fichier BED
plink --bfile BeeMuse --extract list_markers_ID_to_keep.txt --make-bed --out BeeMuse_filtered
```

Effectuer l'ACP sur les données filtrées

- acpBeeMusesamples.sh


# Merged Data - Fusion des deux jeux de données
```
 head BeeMuse_filtered.bim
```
1	AX-643871392	0	14449	A	G
1	AX-643872730	0	31950	C	T
1	AX-643872782	0	36727	T	C

```
 head subset_RefPop_samples.bim
```
NC_037638.1	.	0	5671	C	T
NC_037638.1	.	0	5698	C	T
NC_037638.1	.	0	6621	G	A

```
head BeeMuse_filtered.fam
```
Beemuse_Pool 100_C2.CEL 0 0 0 -9
Beemuse_Pool 102_G2.CEL 0 0 0 -9
Beemuse_Pool 103_I2.CEL 0 0 0 -9

```
head subset_RefPop_samples.fam
```
0	Ab-PacBio	0	0	0	-9
0	BER10	0	0	0	-9
0	BER11	0	0	0	-9

- subsetforPlink.sh
```
#!/bin/bash

sed -i 's/NC_037638.1/1/g' E756_BeeMuSe.vcf
sed -i 's/NC_037639.1/2/g' E756_BeeMuSe.vcf
sed -i 's/NC_037640.1/3/g' E756_BeeMuSe.vcf
sed -i 's/NC_037641.1/4/g' E756_BeeMuSe.vcf
sed -i 's/NC_037642.1/5/g' E756_BeeMuSe.vcf
sed -i 's/NC_037643.1/6/g' E756_BeeMuSe.vcf
sed -i 's/NC_037644.1/7/g' E756_BeeMuSe.vcf
sed -i 's/NC_037645.1/8/g' E756_BeeMuSe.vcf
sed -i 's/NC_037646.1/9/g' E756_BeeMuSe.vcf
sed -i 's/NC_037647.1/10/g' E756_BeeMuSe.vcf
sed -i 's/NC_037648.1/11/g' E756_BeeMuSe.vcf
sed -i 's/NC_037649.1/12/g' E756_BeeMuSe.vcf
sed -i 's/NC_037650.1/13/g' E756_BeeMuSe.vcf
sed -i 's/NC_037651.1/14/g' E756_BeeMuSe.vcf
sed -i 's/NC_037652.1/15/g' E756_BeeMuSe.vcf
sed -i 's/NC_037653.1/16/g' E756_BeeMuSe.vcf
sed -i 's/NC_001566.1/17/g' E756_BeeMuSe.vcf

echo "done"
```
```
bgzip E756_BeeMuSe.vcf
tabix -p vcf E756_BeeMuSe.vcf.gz
```

```
grep -c '^NC_' subset_RefPop_samples_ref.vcf
```
7 023 976 SNPS
```
grep -c '^NC_' E756_BeeMuSe.vcf
```
12000 SNPs
```
head E756_BeeMuSe.vcf
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Beemuse_Pool 100_C2.CEL Beemuse_Pool 102_G2.CEL Beemuse_Pool 103_I2.CEL Beemuse_Pool 104_K2.CEL Beemuse_Pool 105_M2.CEL Beemuse_Po>
17      752     AX-643870442    N       .       .       FAIL    CR=98;ConversionType=Other      GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>
17      1454    AX-643891109    N       .       .       FAIL    CR=100;ConversionType=Other     GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>

```
head subset_RefPop_samples_ref.vcf
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Ab-PacBio       BER10   BER11   BER12   BER13   BER14   BER15   BER16   BER17   BER18   BER19   BER2    BER21   BER4    BER5    BE>
1       5671    .       T       C       332.46  .       AC=0;AF=0.001196;AN=576;BaseQRankSum=-0.969;DP=13376;ExcessHet=0.013;FS=4.15;InbreedingCoeff=0.0188;MLEAC=1;MLEAF=0.0005981;MQ=52.1;MQRankSum=1.51>
1       5698    .       T       C       3633.8  .       AC=0;AF=0.004779;AN=576;BaseQRankSum=0.544;DP=12793;ExcessHet=0.0011;FS=0;InbreedingCoeff=0.0802;MLEAC=9;MLEAF=0.005376;MQ=50.23;MQRankSum=3.5;QD=>

On a ici deux formats différents pour les fichiers VCF de SeqApiPop et BeeMuSe, dont un issu de génotypage par séquençage et l'autre de génotypage par puce de 12000 SNPs.
On va créer les fichiers fusionnés .bim, .bed, .fam en ayant attribué les mêmes identifiants de marqueurs de BeeMuSe à ceux en commun avec SeqApiPop.

## Extraction des SNPs de bonne qualité à garder - list_markers_ID_to_keep.csv (SeqApiPop - BeeMuse)

- script Python dico_CHR_POS_ID.py
``` 
module load devel/python/Python-3.11.1
python dico_CHR_POS_ID.py
```
```
import re

def dico_CHR_POS_ID(nom_fichier_1):
    # Initialiser le dictionnaire
    dico = {}

    # Ouvrir le fichier en mode lecture
    with open(nom_fichier_1, 'r') as fichier:
        # Lire chaque ligne du fichier
        for ligne in fichier:
            # Vérifier si la ligne commence par "NC_"
            if ligne.startswith('NC_'):
                # Utiliser une expression régulière pour extraire CHR, POS et ID
                match = re.match(r'(\S+)\s+(\d+)\s+(\S+)\s+', ligne)

                # Vérifier si la correspondance a réussi
                if match:
                    # Construire la clé CHR_POS
                    chr_pos = f"{match.group(1)}_{match.group(2)}"

                    # Ajouter la paire clé-valeur au dictionnaire
                    dico[chr_pos] = match.group(3)
    return dico

def replace_ID_values(dico, nom_fichier_2):
    # Ouvrir le deuxième fichier en mode lecture
    with open(nom_fichier_2, 'r') as fichier:
        # Lire chaque ligne du fichier
        for ligne in fichier:
            # Vérifier si la ligne commence par "NC_"
            if ligne.startswith('NC_'):
                # Utiliser une expression régulière pour extraire CHR, POS et ID
                match = re.match(r'(\S+)\s+(\d+)\s+(\S+)\s+', ligne)

                # Vérifier si la correspondance a réussi
                if match:
                    # Construire la clé CHR_POS
                    chr_pos = f"{match.group(1)}_{match.group(2)}"

                    # Vérifier si la clé existe dans le premier dictionnaire
                    if chr_pos in dico:
                        # Remplacer la valeur dans le deuxième dictionnaire
                        dico[chr_pos] = match.group(3)

# Premier fichier .vcf importé
vcf_file_1 = 'E756_BeeMuSe.vcf'
dictionnaire_resultats = dico_CHR_POS_ID(vcf_file_1)

# Deuxième fichier .vcf
vcf_file_2 = 'subset_RefPop_samples_ref.vcf'
replace_ID_values(dictionnaire_resultats, vcf_file_2)

# Imprimer les 10 premiers éléments du dictionnaire
for i, (cle, valeur) in enumerate(dictionnaire_resultats.items()):
    print(f"{cle}: {valeur}")

    # Arrêter après les 10 premiers résultats
    if i == 9:
        break
```

On obtient : 

NC_001566.1_752: .

NC_001566.1_1454: .

NC_001566.1_3721: .

NC_001566.1_3816: .

NC_001566.1_5219: AX-643870441

NC_001566.1_10111: AX-643891057

NC_001566.1_10902: .

NC_037638.1_14449: .

NC_037638.1_31950: .

NC_037638.1_36727: .

- Combien y-a-t-il de marqueurs en commun ? 
```
grep 'AX-' list_markers_ID_to_keep.txt | wc -l
```
=> 10256
```
grep 'AX-' E756_BeeMuSe.vcf | wc -l
```
=> 12000
```
grep 'AX-' subset_RefPop_samples_ref_2.vcf | wc -l
```
=> 11709

On a 11709 marqueurs en commun entre E756_BeeMuSe.vcf et subset_RefPop_samples_ref_2.vcf

- vcftoplinkRefPop.sh
```
#!/bin/bash
#SBATCH --mem=8G

module load bioinfo/Bcftools/1.9
module load bioinfo/PLINK/2.00a4

bcftools norm -m-any --output-type z -o subset_501_samples_ref_split.vcf.gz subset_501_samples_ref.vcf
plink2 --vcf subset_501_samples_ref_split.vcf.gz --make-bed --out subset_501_samples_ref --allow-extra-chr

echo "done"
```

- extractlistmarkersID.sh
```
#!/bin/bash
module load bioinfo/PLINK/1.90b7

# Supprimer les marqueurs indésirables du fichier BED
plink --bfile subset_RefPop_samples_ref_2 --extract list_markers_ID_to_keep.txt --make-bed --out subset_RefPop_samples_filtered --allow-extra-chr
```
```
grep 'AX-' list_markers_ID_to_keep.txt | wc -l
```
**10256** SNPs
```
grep 'AX-' BeeMuse_filtered.bim | wc -l
```
12000 => **10256** SNPs
```
grep 'AX-' subset_RefPop_samples_filtered.bim | wc -l
```
7023976 => 11709 => **10030** SNPs

Soit **10030** SNPs filtrés en commun entre les 7023976 SNPs du jeu de données de référence des 870 échantillons SeqApiPop et les 10256 SNPs de bonne qualité à garder de la puce BeeMuSe.

- maf001.sh
```
#! /bin/bash
module load bioinfo/PLINK/1.90b7

NAME=subset_RefPop_samples_filtered
plink --bfile ${NAME} \
  --maf 0.01 \
  --out ${NAME}_maf001 \
  --make-bed
```

 Filtre LD par défault
 
- filterLD.sh
```
#! /bin/bash

module load bioinfo/PLINK/1.90b7
module load bioinfo/PLINK/2.00a4

NAME1=subset_RefPop_samples_filtered_maf001
NAME2=subset_RefPop_samples_filtered_maf001_LD_default

plink --bfile ${NAME1} \
  --out ${NAME2} \
  --indep-pairwise 50 10 0.1

plink --bfile ${NAME1} \
  --out ${NAME2}_pruned \
  --extract ${NAME2}.prune.in \
  --make-bed

plink2 --bfile ${NAME2}_pruned --make-rel square --nonfounders --out ${NAME2}_acp --allow-extra-chr
plink2 --bfile ${NAME2}_pruned --pca --nonfounders --out ${NAME2}_acp --allow-extra-chr
```
