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
Le jeu de données BeeMuSe correspond à 12000 SNPs pour 748 échantillons d'abeilles domestiques

## Conversion au format plink

- convertplinkbeemuse.sh

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

- makebedbeemuse.bash
```
#!/bin/bash
module load bioinfo/PLINK/1.90b7

sed 's/-/ /g' E756_BeeMuSe_num_chr.ped > E756_BeeMuSe_num_chr_2.ped
plink --file E756_BeeMuSe_num_chr_2 --make-bed --no-parents --no-sex --no-pheno -out BeeMuse
```

## ACP

- acpBeeMuSe.bash
```
#!/bin/bash
module load bioinfo/PLINK/2.00a4

plink2 --bfile BeeMuse --make-rel square  --out BeeMuse_acp --allow-extra-chr
plink2 --bfile BeeMuse --pca --out BeeMuse_acp --allow-extra-chr

echo "done"
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

- substForPlinkWrite.bash
```
#!/bin/bash

#substForPlinkWrite.bash

printf "#!/bin/bash\n\n"

printf "cp MetaGenotypesCalled870_raw_snps_allfilter.vcf MetaGenotypesCalled870_raw_snps_allfilter_plink.vcf\n\n"

j=0
for i in `cut -f1 HAv3_1_Chromosomes.list`
do
j=$((j+1))
printf "sed -i \'s/${i}/${j}/g\' MetaGenotypesCalled870_raw_snps_allfilter_plink.vcf\n"
done
```
```
bash substForPlinkWrite.bash > substForPlink.bash
```
```
grep -v '^#' MetaGenotypesCalled870_raw_snps_allfilter_plink.vcf | wc -l
```
=> **7 023 976** SNPs

- select629.bash
```
#! /bin/bash
#select629.bash

module load bioinfo/PLINK/1.90b7
module load bioinfo/Bcftools/1.9

bcftools view -S Diversity_Study_629_Samples.txt MetaGenotypesCalled870_raw_snps_allfilter_plink.vcf > MetaGenotypesCalled870_raw_snps_allfilter_plink_629.vcf 

NAME=SeqApiPop_629

VCFin=MetaGenotypesCalled870_raw_snps_allfilter_plink_629.vcf
VCFout=${NAME}

plink --vcf ${VCFin} \
  --keep-allele-order \
  --a2-allele ${VCFin} 4 3 '#' \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --set-missing-var-ids @:#[HAV3.1]\$1\$2 \
  --chr 1-16 \
  --mind 0.1 \
  --geno 0.1 \
  --out ${VCFout} \
  --make-bed \
  --missing

plink --bfile ${VCFout} \
  --maf 0.01 \
  --out ${VCFout}_maf001 \
  --make-bed
```

Une liste de 629 échantillons est sélectionnée pour les analyses de structure des populations en enlevant :
- Les échantillons dupliqués d'une même ruche
- Les échantillons des sous-populations expérimentales
- Les échantillons d'une autre étude
- Les 15 échantillons avec > 0,1 de données manquantes.

- select629_LD.bash
```
#! /bin/bash
#select629_LD.bash

module load bioinfo/PLINK/1.90b7

NAME1=SeqApiPop_629_maf001
NAME2=SeqApiPop_629_maf001_LD03

plink --bfile ${NAME1} \
  --out ${NAME2} \
  --indep-pairwise 1749 175 0.3

plink --bfile ${NAME1} \
  --out ${NAME2}_pruned \
  --extract ${NAME2}.prune.in \
  --make-bed

plink --bfile ${NAME2}_pruned \
  --out ${NAME2}_acp \
  --pca

plink --bfile ${NAME2}_pruned \
  --out ${NAME2}_acp \
  --make-rel square
```

- launchAdmixtureRunsWriteScriptsMAF001.bash
```
#!/bin/bash
#launchAdmixtureRunsWriteScriptsMAF001.bash

LD=LD03

for i in $(seq 00 29)
do
mkdir SeqApiPop_629_MAF001_${LD}rep${i}

echo \#!/bin/bash > SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo
echo \#admixtureAnalysis_multiThread.sh >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo module load bioinfo/ADMIXTURE/1.3.0
 >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo IN=SeqApiPop_629_maf001_${LD}_pruned.bed >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo "for K in 2 3 4 5 6 7 8 9 10 11 12;" >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo do >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo -e sbatch --cpus-per-task=4 --mem-per-cpu=4G \\ >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo -e "	-J \${K}admixt -o \${IN}.\${K}.o -e \${IN}.\${K}.e" \\ >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo -en '	--wrap="admixture --cv'  >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo -en " -s ${RANDOM}"  >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo  -e ' ../${IN} ${K} -j4 | tee ${IN}.log${K}"' >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo done >> SeqApiPop_629_MAF001_${LD}rep${i}/admixtureAnalysis.sh

done
```

- obtainCVerrorAllSeqApiPop.bash
```
#!/bin/bash

#obtainCVerrorAllSeqApiPop.bash

for i in `ls | grep ^SeqApiPop`
do
grep CV ${i}/*log* | \
	awk -v var="$i" 'BEGIN{OFS="\t"}{print $3,$4, var}' | \
	sed 's/(K=//' | \
	sed 's/)://' | \
	awk 'BEGIN{FS="_";OFS="\t"}{print $1,$4,$3}' | \
	awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$5}'
done
```

- remaneQmatrixes.bash
```
#!/bin/bash
#renameQmatrixes.bash

#copies and renames the Q matrix outputs in the directory Qfiles

#Edit

MAF=MAF001

for h in $(seq 3 3)
do
    	LD=LD0${h}
    	for i in $(seq 0 49)
    	do
            	for j in `ls SeqApiPop_629_${MAF}_${LD}rep${i}/ | grep Q$`
            	do
                    	cp SeqApiPop_629_${MAF}_${LD}rep${i}/${j}  Qfiles/${j%.*}.r${i}.Q
            	done
            	#echo ${i}
    	done
done
```

Obtenir les erreurs de CV
```
grep "CV" *log* | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > SeqApiPop_629_maf001_LD03_0.cv.error

```
Créer le fichier pong avec les matrices Q Admixture
```
ls  ../Qfiles/* | \
	grep LD03 | \
	awk 'BEGIN{FS="/";OFS="\t"}{print $3, $0}' | \
	awk 'BEGIN{FS=".";OFS="\t"}{print $1"_"$2"_"$3"_"$4, $2, ".."$6"."$7"."$8"."$9}' | \
	awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' > pong_filemap_629_maf001_LD03_K2K12
```

### PONG

```
awk '{print $1}' SeqApiPop_629_maf001_LD03_pruned.fam > ind2pop_629.list

awk -F',' 'NR==FNR{a[$1]=$5; next} {print a[$1]}' SeqApiPop_labels.csv ind2pop_629.list > ind2pop_629_Label.list
awk -F',' 'NR==FNR{a[$1]=$6; next} {print a[$1]}' SeqApiPop_labels.csv ind2pop_629.list > ind2pop_629_GeographicOrigin.list

sed -i 's/ /_/g' ind2pop_629_Label.list
sed -i 's/ /_/g' ind2pop_629_GeographicOrigin.list
```

- popOrder_629_Label.list

```
Iberiensis_Spain IberiensisSpain
Ouessant_Conservatory MelliferaOuessant
Colonsay_Conservatory MelliferaColonsay
Porquerolles_Conservatory MelliferaPorquerolles
Sollies_Conservatory MelliferaSollies
Savoy_Conservatory SavoyConservatory
Ligustica_Italy LigusticaItaly
Carnica_Slovenia CarnicaSlovenia
Carnica_Germany CarnicaGermany
Carnica_France CarnicaFrance
Carnica_Switzerland CarnicaSwitzerland
Carnica_Poland CarnicaPoland
Caucasia_France CaucasiaFrance
China China
Royal_Jelly_France RoyalJellyFrance
Corsica_Breeder Corsica
Buckfast_France BuckfastFrance
Buckfast_Switzerland BuckfastSwitzerland
Tarn_1_Breeder Tarn1
Tarn_2_Breeder Tarn2
Hautes_Pyrenees_Breeder HautesPyrenees
Ariege_Breeder AriegeBreeder
Ariege_Conservatory AriegeConservatory
Brittany_Conservatory BrittanyConservatory
Brittany_Breeder BrittanyBreeder
Isere_1_Breeder Isere1Breeder
Isere_2_Breeder Isere2Breeder
Herault_Breeder HeraultBreeder
Sarthe_Breeder Sarthe
Vaucluse_Breeder VaucluseBreeder
Unknown Unknown
```

Création du fichier des couleurs pour la visualisation Admixture

- colors
```
black
#FE9001
#018a16
#EDFE01
#0603a6
#FE0101
#9B9B9B
#A15D19
#DEDAD7
#FEB7F7
#D09CFE
#9CF4FE
```
On lance PONG pour K2 à K12 :
```
pong -m pong_filemap_629_maf001_LD03_K2K12 -n popOrder_629_Label.list -i ind2pop_629_Label.list -l colors -s 0.98
```
Pour K2 à K9 :
```
pong -m pong_filemap_629_maf001_LD03_K2K9 -n popOrder_629_Label.list -i ind2pop_629_Label.list -l colors -s 0.98
```

## ACP

- acpseqapipop.bash
```
#!/bin/bash
#acpseqapipop.bash
module load bioinfo/PLINK/2.00a4

plink2 --bfile subset_501_samples_ref --make-rel square --nonfounders --out subset_501_samples_ref
plink2 --bfile subset_501_samples_ref --pca --nonfounders --out subset_501_samples_ref

echo "done"
```

lien E756_BeeMuSe_num_chr_2col.map et list_markers_ID_to_keep.csv

.map : Marker ID ($2)
.csv : probeset_id ($1)

```
awk -F';' 'NR>1 {print $1}' list_markers_ID_to_keep.csv > list_markers_ID_to_keep.txt
```
```
plink --file E756_BeeMuSe_num_chr_2col --make-bed --no-parents --no-sex --no-pheno -out BeeMuse
```

- extractlistmarkersID.bash
```
#!/bin/bash
#extractlistmarkersID.bash
module load bioinfo/PLINK/1.90b7

# Supprimer les marqueurs indésirables du fichier BED
plink --bfile BeeMuse --extract list_markers_ID_to_keep.txt --make-bed --out BeeMuse_filtered
```

Effectuer l'ACP sur les données filtrées

- acpRefPop.bash
```
#!/bin/bash
#acpRefPop.bash
module load bioinfo/PLINK/2.00a4

plink2 --bfile subset_RefPop_samples_filtered --make-rel square --nonfounders --out subset_RefPop_samples_filtered_acp --allow-extra-chr
plink2 --bfile subset_RefPop_samples_filtered --pca  --nonfounders --out subset_RefPop_samples_filtered_acp --allow-extra-chr

echo "done"
```


## SeqApiPop - MAF > 0.01 + LD pruning = 0.1 (fenêtre de 50 et pas de 10 bp)

- maf001.bash
```
#! /bin/bash
#maf001.bash
module load bioinfo/PLINK/1.90b7

NAME=subset_RefPop_samples_filtered
plink --bfile ${NAME} \
  --maf 0.01 \
  --out ${NAME}_maf001 \
  --make-bed
```

 Filtre LD par défault et  ACP
 
- filterLD.bash
```
#! /bin/bash
#filterLD.bash

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

## Jeu de données de référence sans Mellifera Colonsay / Ouessant - 561 échantillons

On enlève les échantillons SeqApiPop des populaitons Mellifera Colonsay / Ouessant

On extrait la liste des noms des 561 échantillons à garder dans : Diversity_Study_561_Samples.txt, puis on relance les scripts permettant les analyses d'ACP et d'Admixture.

- select561.bash
```
#!/bin/bash
#select561.bash

module load bioinfo/PLINK/1.90b7
module load bioinfo/Bcftools/1.9

bcftools view -S Diversity_Study_561_Samples.txt MetaGenotypesCalled870_raw_snps_allfilter_plink.vcf.gz > MetaGenotypesCalled870_raw_snps_allfilter_plink_561.vcf.gz

NAME=SeqApiPop_561

VCFin=MetaGenotypesCalled870_raw_snps_allfilter_plink_561.vcf.gz
VCFout=${NAME}

plink --vcf ${VCFin} \
  --keep-allele-order \
  --a2-allele ${VCFin} 4 3 '#' \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --set-missing-var-ids @:#[HAV3.1]\$1\$2 \
  --chr 1-16 \
  --mind 0.1 \
  --geno 0.1 \
  --out ${VCFout} \
  --make-bed \
  --missing

plink --bfile ${VCFout} \
  --maf 0.01 \
  --out ${VCFout}_maf001 \
  --make-bed
```

- select561_LD.bash
```
#!/bin/bash
#select561_LD.bash

module load bioinfo/PLINK/1.90b7

NAME1=SeqApiPop_561_maf001
NAME2=SeqApiPop_561_maf001_LD03

plink --bfile ${NAME1} \
  --out ${NAME2} \
  --indep-pairwise 1749 175 0.3

plink --bfile ${NAME1} \
  --out ${NAME2}_pruned \
  --extract ${NAME2}.prune.in \
  --make-bed

plink --bfile ${NAME2}_pruned \
  --out ${NAME2}_acp \
  --pca

plink --bfile ${NAME2}_pruned \
  --out ${NAME2}_acp \
  --make-rel square
```

- launchAdmixtureRunsWriteScripts561MAF001.bash
```
#!/bin/bash
#launchAdmixtureRunsWriteScripts561MAF001.bash

LD=LD03

for i in $(seq 00 29)
do
mkdir SeqApiPop_561_MAF001_${LD}rep${i}

echo \#!/bin/bash > SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo
echo \#admixtureAnalysis_multiThread.sh >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo module load bioinfo/ADMIXTURE/1.3.0
 >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo IN=SeqApiPop_561_maf001_${LD}_pruned.bed >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo "for K in 2 3 4 5 6 7 8 9 10 11 12;" >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo do >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo -e sbatch --cpus-per-task=4 --mem-per-cpu=4G \\ >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo -e "	-J \${K}admixt -o \${IN}.\${K}.o -e \${IN}.\${K}.e" \\ >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo -en '	--wrap="admixture --cv'  >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo -en " -s ${RANDOM}"  >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo  -e ' ../${IN} ${K} -j4 | tee ${IN}.log${K}"' >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh
echo done >> SeqApiPop_561_MAF001_${LD}rep${i}/admixtureAnalysis.sh

done
```
- Obtenir les erreurs de CV d'Admixture : 
```
grep "CV" *log* | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > SeqApiPop_561_maf001_LD03_1.cv.error
```
- Créer le fichier PONG des matrices Q d'Admixture : 
```
ls ../SeqApiPop_561_maf001_LD03/* | grep LD03 | grep -vE '\.10\.|\.11\.|\.12\.' | awk 'BEGIN{FS="/";OFS="\t"}{print $3, $0}' | awk 'BEGIN{FS=".";OFS="\t"}{print $1"_"$2"_"$3"_"$4, $2, ".."$6"."$7"."$8"."$9}' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' > pong_filemap_561_maf001_LD03_K2K9
```

- popOrder_561_Label.list
```
Iberiensis_Spain IberiensisSpain
Porquerolles_Conservatory MelliferaPorquerolles
Sollies_Conservatory MelliferaSollies
Savoy_Conservatory SavoyConservatory
Ligustica_Italy LigusticaItaly
Carnica_Slovenia CarnicaSlovenia
Carnica_Germany CarnicaGermany
Carnica_France CarnicaFrance
Carnica_Switzerland CarnicaSwitzerland
Carnica_Poland CarnicaPoland
Caucasia_France CaucasiaFrance
China China
Royal_Jelly_France RoyalJellyFrance
Corsica_Breeder Corsica
Buckfast_France BuckfastFrance
Buckfast_Switzerland BuckfastSwitzerland
Tarn_1_Breeder Tarn1
Tarn_2_Breeder Tarn2
Hautes_Pyrenees_Breeder HautesPyrenees
Ariege_Breeder AriegeBreeder
Ariege_Conservatory AriegeConservatory
Brittany_Conservatory BrittanyConservatory
Brittany_Breeder BrittanyBreeder
Isere_1_Breeder Isere1Breeder
Isere_2_Breeder Isere2Breeder
Herault_Breeder HeraultBreeder
Sarthe_Breeder Sarthe
Vaucluse_Breeder VaucluseBreeder
Unknown Unknown
```
- Visualisation d'Admixture avec PONG :
```
pong -m pong_filemap_561_maf001_LD03_K2K9 -n popOrder_561_Label.list -i ind2pop_561_Label.list -l colors_K9 -s 0.98
```

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
#dico_CHR_POS_ID.py
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
=> **10256** SNPs
```
grep 'AX-' E756_BeeMuSe.vcf | wc -l
```
=> **12000** SNPs
```
grep 'AX-' subset_RefPop_samples_ref_2.vcf | wc -l
```
=> **11709** SNPs

On a 11709 marqueurs en commun entre E756_BeeMuSe.vcf et subset_RefPop_samples_ref_2.vcf

- extractlistmarkersID.sh
```
#!/bin/bash
#extractlistmarkersID.sh
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

On extrait les 561 échantillons de référence d'intérêt et on applique les filtres MAF > 0.01 et LD pruning = 0.1 avec une fenêtre de 50 SNPs et un pas de 10 bp. 

Ensuite on récupère les identifiants des marqueurs obtenu après le filtre LD 
```
awk '{print $2}' SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD_default_pruned.bim > marker_ids_filtered_maf001_LD_default.txt
```
```
less -S marker_ids_filtered_maf001_LD_default.txt | wc -l
```
=> **1055** SNPs après filtre MAF > 0.01 et LD pruning = 0.1 (fenêtre de 50 SNPs et pas de 10 bp)

Extraction des  1055 marqueurs du fichier VCF BeeMuSe

- extractlistmarkersIDbeemuse.bash
```
#!/bin/bash
#extractlistmarkersIDbeemuse.bash
module load bioinfo/PLINK/1.90b7

plink --bfile BeeMuse_filtered --extract marker_ids_filtered_maf001_LD_default.txt --make-bed --out BeeMuSe_filtered_maf001_LD_default --allow-extra-chr
```
```
less BeeMuSe_filtered_maf001_LD_default.bim | wc -
```
=> **1055** SNPs

On a bien les mêmes marqueurs dans les deux jeux de données, on peut alors les fusionner :

- bmerge.bash
```
#!/bin/bash
#bmerge.bash
module load bioinfo/PLINK/1.90b7

plink --bfile BeeMuSe_filtered_maf001_LD_default \
   	--bmerge SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD_default_pruned.bed \
            	SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD_default_pruned.bim \
            	SeqApiPop_561_SNPsBeeMuSe_filtered_maf001_LD_default_pruned.fam \
   	--make-bed \
   	--out merged_BeeMuSe_SeqApiPop_561_filtered_maf001_LD_default
```

## ACP

- acpmergedBeeMuSeSeqApiPop.bash
```
#!/bin/bash
#acpmergedBeeMuSeSeqApiPop.bash
module load bioinfo/PLINK/2.00a4

plink2 --bfile merged_BeeMuSe_SeqApiPop_561_filtered_maf001_LD_default --make-rel square  --out merged_BeeMuSe_SeqApiPop_561_filtered_maf001_LD_default_acp --allow-extra-chr
plink2 --bfile merged_BeeMuSe_SeqApiPop_561_filtered_maf001_LD_default --pca --out merged_BeeMuSe_SeqApiPop_561_filtered_maf001_LD_default_acp --allow-extra-chr
```
