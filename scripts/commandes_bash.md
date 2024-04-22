Certains scripts et commandes sont issus depuis le GitHub pour [SeqApiPop](https://github.com/avignal5/SeqApiPop)

(Alain Vignal, David Wragg, Sonia E. Eynard, Benjamin Basso, Kamila Canale-Tabet, Emmanuelle Labarthe, Olivier Bouchez, Kaspar Bienefeld, Małgorzata Bieńkowska, Cecilia Costa, Aleš Gregorc, Per Kryger, Melanie Parejo, M. Alice Pinto, Jean-Pierre Bidanel, Bertrand Servin, & Yves Le Conte. (2021). Data from: Complex population structure and haplotype patterns in Western Europe honey bee from sequencing a large panel of haploid drones [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5592452)

- Récupérer les fichiers de référence SeqApiPop depuis [zenodo](https://zenodo.org/records/5592452)
- Diversity_Study_629_Samples.txt
- MetaGenotypesCalled870_raw_snps_allfilter.vcf.gz

(Wragg D, Eynard SE, Basso B, Canale-Tabet K, Labarthe E, Bouchez O, et al. Complex population structure and haplotype patterns in Western Europe honey bee from sequencing a large panel of haploid drones [Internet]. Genetics; 2021 Sep. Available from: http://biorxiv.org/lookup/doi/10.1101/2021.09.20.460798)

# BeeMuse
```
touch convertnumchrmap.bash
nano convertnumchrmap
```
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
```
sbatch convertnumchrmap.bash
```
#Chromosome     Marker ID       Genetic distance        Physical position
17      AX-643870442    0.000   752
17      AX-643891109    0.000   1454
```
module load bioinfo/PLINK/1.90b7
plink --file E756_BeeMuSe_num_chr --make-bed --no-fid --no-parents --no-sex --no-pheno
```
wc -l E756_BeeMuSe.map 
#12000
head -n 2 E756_BeeMuSe.ped | awk '{print NF}'
#24000

plink --file E756_BeeMuSe_num_chr --make-bed --no-parents --no-sex --no-pheno

less -S +376 E756_BeeMuSe_num_chr.ped

# pb 2 colonnes puis 1 à cause de '-' au lieu de ' '

sed 's/-/ /g' E756_BeeMuSe_num_chr.ped > E756_BeeMuSe_num_chr_2col.ped

plink --file E756_BeeMuSe_num_chr_2col --make-bed --no-parents --no-sex --no-pheno -out BeeMuse

# ACP
module load bioinfo/PLINK/2.00a4
plink2 --bfile BeeMuse --make-rel square --nonfounders
plink2 --bfile BeeMuse --pca  --nonfounders

	#15/1
#Zenodo vcf + liste 629 samples
wget https://zenodo.org/record/5592452/files/Diversity_Study_629_Samples.txt
wget https://zenodo.org/record/5592452/files/MetaGenotypesCalled870_raw_snps_allfilter.vcf.gz

#R -> liste 501 samples connus : Diversity_Study_501_Samples.txt

zcat MetaGenotypesCalled870_raw_snps_allfilter.vcf.gz | head

rm E756_BeeMuSe_num_chr_2col
touch subset501ref.bash

# script subset501ref.bash

#!/bin/bash
module load bioinfo/Bcftools/1.9

# Extraire les individus de référence du VCF d'origine
bcftools view -S Diversity_Study_501_Samples.txt MetaGenotypesCalled870_raw_snps_allfilter.vcf.gz > subset_501_samples_ref.vcf

sbatch subset501ref.bash
squeue -u nbettembour

# VCF -> plink .bed .bim .fam

touch vcftoplink.bash

	#script vcftoplink
#!/bin/bash
#SBATCH --mem=8G

module load bioinfo/Bcftools/1.9
module load bioinfo/PLINK/2.00a4

bcftools norm -m-any --output-type z -o subset_501_samples_ref_split.vcf.gz subset_501_samples_ref.vcf
plink2 --vcf subset_501_samples_ref_split.vcf.gz --make-bed --out subset_501_samples_ref --allow-extra-chr

echo "done"

sbatch vcftoplink.bash
squeue -u nbettembour

sbatch vcftoplinkRefPop.bash
squeue -u nbettembour

#ACP

#!/bin/bash
module load bioinfo/PLINK/2.00a4

plink2 --bfile subset_501_samples_ref --make-rel square --nonfounders --out subset_501_samples_ref
plink2 --bfile subset_501_samples_ref --pca --nonfounders --out subset_501_samples_ref

echo "done"

sbatch acp501samples.bash
squeue -u nbettembour

	# 16/1

#Refaire même chose seulement reference populations
filtered_data <- subset(eigenvec_501_seq_api_labels, Label != 'Ariege Conservatory' & Label != 'Brittany Conservatory' & UniqueInHive != 'Unknown' & UniqueInHive != 'Buckfast' & GeneticOrigin != 'Buckfast')

#lien E756_BeeMuSe_num_chr_2col.map et list_markers_ID_to_keep.csv

.map : Marker ID ($2)
.csv : probeset_id ($1)

# puis refaire ACP - beemuse ?

awk -F';' 'NR>1 {print $1}' list_markers_ID_to_keep.csv > list_markers_ID_to_keep.txt

#plink --file E756_BeeMuSe_num_chr_2col --make-bed --no-parents --no-sex --no-pheno -out BeeMuse

touch extractlistmarkersID.bash

module load bioinfo/PLINK/1.90b7
# Supprimer les marqueurs indésirables du fichier BED
plink --bfile BeeMuse --extract list_markers_ID_to_keep.txt --make-bed --out BeeMuse_filtered

sbatch extractlistmarkersID.bash

# Effectuer l'ACP sur les données filtrées
# script acpBeeMusesamples.bash
touch acpBeeMusesamples.bash

sbatch acpBeeMusesamples.bash
squeue -u nbettembour


	# 17/1
#ACP merged data : RefPop + BeeMuse puce 12000 SNPs
BeeMuse_filtered
subset_RefPop_samples

 head BeeMuse_filtered.bim 
1	AX-643871392	0	14449	A	G
1	AX-643872730	0	31950	C	T
1	AX-643872782	0	36727	T	C

 head subset_RefPop_samples.bim
NC_037638.1	.	0	5671	C	T
NC_037638.1	.	0	5698	C	T
NC_037638.1	.	0	6621	G	A

[nbettembour@genobioinfo1 data]$ head BeeMuse_filtered.fam
Beemuse_Pool 100_C2.CEL 0 0 0 -9
Beemuse_Pool 102_G2.CEL 0 0 0 -9
Beemuse_Pool 103_I2.CEL 0 0 0 -9

[nbettembour@genobioinfo1 data]$ head subset_RefPop_samples.fam
0	Ab-PacBio	0	0	0	-9
0	BER10	0	0	0	-9
0	BER11	0	0	0	-9

	# script mergebimbedfam.bash
	
#plink --bfile BeeMuse_filtered --bmerge subset_RefPop_samples.bed subset_RefPop_samples.bim subset_RefPop_samples.fam --make-bed --out merged_BeeMuse_RefPop --allow-extra-chr
#plink --bfile BeeMuse_filtered --bmerge subset_RefPop_samples.bed subset_RefPop_samples.bim subset_RefPop_samples.fam --merge-mode 6 --make-bed --out merged_BeeMuse_RefPop --allow-extra-chr

sbatch mergebimbedfam.bash
squeue -u nbettembour

#plink --bfile BeeMuse_filtered --merge-list merged_BeeMuse_RefPop-merge.missnp --make-bed --out merged_BeeMuse_RefPop_updated

module load bioinfo/Bcftools/1.9

bcftools merge E756_BeeMuSe.vcf.gz subset_RefPop_samples_split.vcf.gz -o merged_BeeMuse_RefPop.vcf.gz

bcftools sort merged.vcf.gz -o merged_sorted.vcf.gz
bcftools index merged_sorted.vcf.gz

#head BeeMuse_filtered.bim
1	AX-643871392	0	14449	A	G
1	AX-643872730	0	31950	C	T

#head subset_RefPop_samples.bim
1	.	0	5671	C	T
1	.	0	5698	C	T

	#script mergebimbedfam
plink --bfile BeeMuse_filtered --bmerge subset_RefPop_samples.bed subset_RefPop_samples.bim subset_RefPop_samples.fam --make-bed --out merged_BeeMuse_RefPop
#plink --bfile BeeMuse_filtered --merge-list empty_merge_list.txt --make-bed --out merged_BeeMuse_RefPop_updated

	#subsetforPlink.bash
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

bgzip E756_BeeMuSe.vcf
tabix -p vcf E756_BeeMuSe.vcf.gz

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Beemuse_Pool 100_C2.CEL Beemuse_Pool 102_G2.CEL Beemuse_Pool 103_I2.CEL Beemuse_Pool 104_K2.CEL Beemuse_Pool 105_M2.CEL Beemuse_Po>
17      752     AX-643870442    N       .       .       FAIL    CR=98;ConversionType=Other      GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>
17      1454    AX-643891109    N       .       .       FAIL    CR=100;ConversionType=Other     GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Ab-PacBio       BER10   BER11   BER12   BER13   BER14   BER15   BER16   BER17   BER18   BER19   BER2    BER21   BER4    BER5    BE>
1       5671    .       T       C       332.46  .       AC=0;AF=0.001196;AN=576;BaseQRankSum=-0.969;DP=13376;ExcessHet=0.013;FS=4.15;InbreedingCoeff=0.0188;MLEAC=1;MLEAF=0.0005981;MQ=52.1;MQRankSum=1.51>
1       5698    .       T       C       3633.8  .       AC=0;AF=0.004779;AN=576;BaseQRankSum=0.544;DP=12793;ExcessHet=0.0011;FS=0;InbreedingCoeff=0.0802;MLEAC=9;MLEAF=0.005376;MQ=50.23;MQRankSum=3.5;QD=>

	#script bash
bcftools merge E756_BeeMuSe.vcf subset_RefPop_samples_ref.vcf.gz -o merged_BeeMuse_RefPop.vcf.gz

#refaire bim bed fam en ayant changé ID chromosomes même que BeeMuse ?
sbatch vcftoplink.bash

	#concaténation
cat BeeMuse_filtered_cat.bim subset_RefPop_samples_ref_num_chr.bim > merged_BeeMuse_RefPop_cat.bim
cat BeeMuse_filtered.fam subset_RefPop_samples_ref_cat.fam > merged_BeeMuse_RefPop_cat.fam
#cat BeeMuse_filtered.bed subset_RefPop_samples_ref.bed > merged_BeeMuse_RefPop_cat.bed

awk 'BEGIN{OFS=FS="\t"} {$2="."; print}' BeeMuse_filtered.bim > BeeMuse_filtered_cat.bim

awk '{printf "%s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6}' subset_RefPop_samples_ref.fam > subset_RefPop_samples_ref_cat.fam

	# ACP
#!/bin/bash
module load bioinfo/PLINK/2.00a4

plink2 --bfile merged_BeeMuse_RefPop_cat --make-rel square --nonfounders --out merged_BeeMuse_RefPop_cat_pca --allow-extra-chr
plink2 --bfile merged_BeeMuse_RefPop_cat --pca --nonfounders --out merged_BeeMuse_RefPop_cat_pca --allow-extra-chr

echo "done"

#!/bin/bash
module load bioinfo/PLINK/2.00a4

# Sort the variants and create a new set of files
plink2 --bfile merged_BeeMuse_RefPop_cat --make-pgen --sort-vars --out merged_BeeMuse_RefPop_sorted

# Perform PCA using the sorted files
plink2 --pfile merged_BeeMuse_RefPop_sorted --pca --nonfounders --out merged_BeeMuse_RefPop_cat_pca

echo "done"

	# 18/1
grep -c '^NC_' subset_RefPop_samples_ref.vcf
7023976
grep -c '^NC_' E756_BeeMuSe.vcf
12000

# refaire cat
# comment faire lien list_markers_ID_to_keep et RefPop csv ? - BeeMuse.csv ?

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Beemuse_Pool 100_C2.CEL Beemuse_Pool 102_G2.CEL Beemuse_Pool 103_I2.CEL Beemuse_Pool 104_K2.CEL Beemuse_Pool 105_M2.CEL Beemuse_Po>
NC_001566.1     752     AX-643870442    N       .       .       FAIL    CR=98;ConversionType=Other      GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>
NC_001566.1     1454    AX-643891109    N       .       .       FAIL    CR=100;ConversionType=Other     GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>
NC_001566.1     3721    AX-643891061    N       .       .       PASS    CR=100;ConversionType=MonoHighResolution        GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>
NC_001566.1     3816    AX-643870915    N       .       .       PASS    CR=100;ConversionType=NoMinorHom        GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>
NC_001566.1     5219    AX-643870441    N       .       .       PASS    CR=100;ConversionType=NoMinorHom        GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>
NC_001566.1     10111   AX-643891057    N       .       .       FAIL    CR=100;ConversionType=OTV       GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>
NC_001566.1     10902   AX-643891108    N       .       .       FAIL    CR=99;ConversionType=OTV        GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>
NC_037638.1     14449   AX-643871392    N       .       .       PASS    CR=100;ConversionType=PolyHighResolution        GT      ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./>

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Ab-PacBio       BER10   BER11   BER12   BER13   BER14   BER15   BER16   BER17   BER18   BER19   BER2    BER21   BER4    BER5    BE>
NC_037638.1     5671    .       T       C       332.46  .       AC=0;AF=0.001196;AN=576;BaseQRankSum=-0.969;DP=13376;ExcessHet=0.013;FS=4.15;InbreedingCoeff=0.0188;MLEAC=1;MLEAF=0.0005981;MQ=52.1;MQRank>
NC_037638.1     5698    .       T       C       3633.8  .       AC=0;AF=0.004779;AN=576;BaseQRankSum=0.544;DP=12793;ExcessHet=0.0011;FS=0;InbreedingCoeff=0.0802;MLEAC=9;MLEAF=0.005376;MQ=50.23;MQRankSum>

NC_037638.1     14449   .       G       A       219798  .       AC=314;AF=0.446;AN=602;BaseQRankSum=1.25;DP=15313;ExcessHet=0;FS=0;InbreedingCoeff=0.966;MLEAC=778;MLEAF=0.453;MQ=60;MQRankSum=0;QD=32.2;R>
=> AX-643871392

	# script Python dico_CHR_POS_ID
module load devel/python/Python-3.11.1

#dico 2ème vcf
#regarder mêmes clés CHR_POS
#remplacer valeurs clés 2 par valeurs clés 1
#print vcf

python dico_CHR_POS_ID.py

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

# Exemple d'utilisation
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

=> NC_001566.1_752: .
NC_001566.1_1454: .
NC_001566.1_3721: .
NC_001566.1_3816: .
NC_001566.1_5219: AX-643870441
NC_001566.1_10111: AX-643891057
NC_001566.1_10902: .
NC_037638.1_14449: .
NC_037638.1_31950: .
NC_037638.1_36727: .

import re

def dico_CHR_POS_ID(nom_fichier):
    dico = {}
    with open(nom_fichier, 'r') as fichier:
        for ligne in fichier:
            if ligne.startswith('NC_'):
                match = re.match(r'(\S+)\s+(\d+)\s+(\S+)\s+', ligne)
                if match:
                    chr_pos = f"{match.group(1)}_{match.group(2)}"
                    dico[chr_pos] = match.group(3)
    return dico

def replace_ID_values(dico, nom_fichier_entree, nom_fichier_sortie):
    with open(nom_fichier_entree, 'r') as fichier_entree, open(nom_fichier_sortie, 'w') as fichier_sortie:
        for ligne in fichier_entree:
            if ligne.startswith('NC_'):
                match = re.match(r'(\S+)\s+(\d+)\s+(\S+)\s+', ligne)
                if match:
                    chr_pos = f"{match.group(1)}_{match.group(2)}"
                    if chr_pos in dico:
                        nouvelle_ligne = re.sub(r'(\S+\s+\d+\s+)(\S+)(\s+)', f"\\1{dico[chr_pos]}\\3", ligne)
                        fichier_sortie.write(nouvelle_ligne)
                    else:
                        fichier_sortie.write(ligne)
                else:
                    fichier_sortie.write(ligne)
            else:
                fichier_sortie.write(ligne)

# Exemple d'utilisation
vcf_file_entree = 'E756_BeeMuSe.vcf'
dictionnaire_resultats = dico_CHR_POS_ID(vcf_file_entree)

# Deuxième fichier .vcf
vcf_file_sortie = 'subset_RefPop_samples_ref_2.vcf'
replace_ID_values(dictionnaire_resultats, 'subset_RefPop_samples_ref.vcf', vcf_file_sortie)

NC_037638.1     54359   .       C       T       2836.52 .       AC=4;AF=0.004037;AN=602;BaseQRankSum=1.98;DP=14732;ExcessHet=0;FS=1.688;InbreedingCoeff=0.3221;MLEAC=7;MLEAF=0.004037;MQ=60;MQRankSum=0;QD>
NC_037638.1     54377   AX-643872963    G       A       142698  AX-643872963    AC=206;AF=0.314;AN=600;BaseQRankSum=-0.319;DP=14743;ExcessHet=-0;FS=0;InbreedingCoeff=0.9591;MLEAC=547;MLEAF=0.315;MQ=60;M>
NC_037638.1     54447   .       G       A       283.27  .       AC=0;AF=0.001152;AN=602;DP=14853;ExcessHet=0.0449;FS=0;InbreedingCoeff=0.1199;MLEAC=1;MLEAF=0.000576;MQ=60;QD=28.33;SOR=0.693   GT:AD:DP:G>

	# 19/1

# compter cb en communs de marqueurs ?
grep 'AX-' list_markers_ID_to_keep.txt | wc -l
=> 10256
grep 'AX-' E756_BeeMuSe.vcf | wc -l
=> 12000
grep 'AX-' subset_RefPop_samples_ref_2.vcf | wc -l
=> 11709

# refaire bim bed fam 
sbatch vcftoplinkRefPop

	#script extractlistmarkersID.bash
#!/bin/bash
module load bioinfo/PLINK/1.90b7

# Supprimer les marqueurs indésirables du fichier BED
plink --bfile subset_RefPop_samples_ref_2 --extract list_markers_ID_to_keep.txt --make-bed --out subset_RefPop_samples_filtered --allow-extra-chr

grep 'AX-' list_markers_ID_to_keep.txt | wc -l
=> 10256  
grep 'AX-' BeeMuse_filtered.bim | wc -l
12000 => 10256
grep 'AX-' subset_RefPop_samples_filtered.bim | wc -l
7023976 => 11709 => 10030

	#scrip extractlistmarkersIDvcf.bash
#!/bin/bash
module load bioinfo/VCFtools/0.1.16

vcftools --vcf subset_RefPop_samples_ref_2.vcf --snps list_markers_ID_to_keep.txt --recode --recode-INFO-all --out subset_RefPop_samples_filtered

grep 'AX-' subset_RefPop_samples_filtered.recode.vcf | wc -l
=> 10030

bgzip subset_RefPop_samples_filtered.recode.vcf
tabix subset_RefPop_samples_filtered.recode.vcf.gz

	#script mergevcf.bash
module load bioinfo/Bcftools/1.9

bcftools merge E756_BeeMuSe.vcf subset_RefPop_samples_filtered.recode.vcf -o merged_BeeMuse_RefPop.vcf.gz

#x
[W::hts_idx_load2] The index file is older than the data file: subset_RefPop_samples_filtered.recode.vcf.gz.tbi
[W::hts_idx_load2] The index file is older than the data file: subset_RefPop_samples_filtered.recode.vcf.gz.tbi
[W::vcf_parse] FILTER 'AX-643871392' is not defined in the header
[W::vcf_parse] FILTER 'AX-643872730' is not defined in the header
The REF prefixes differ: G vs N (1,1)
Failed to merge alleles at NC_037638.1:14449 in subset_RefPop_samples_filtered.recode.vcf.gz

	#essayer avec marqueurs envoyés ~18000
awk -F',' 'NR>1 {print $3}' chip_list_hav31_corrected.csv > chip_list_hav31_corrected.txt
=> 18123

	#format CHR:POS
sed -i 's/_/:/2' chip_list_hav31_corrected.txt
sed -i 's/:/\t/g' chip_list_hav31_corrected.txt

	#script extractchiplisthav31corrected.bash
#!/bin/bash

bcftools filter -R chip_list_hav31_corrected.txt subset_RefPop_samples_ref.vcf -O v -o subset_RefPop_samples_filtered_18000.vcf

grep 'NC_' subset_RefPop_samples_filtered_18000.vcf | wc -l
=> 122

awk '!/^#/{print $1"\t"$2}' subset_RefPop_samples_ref.vcf > test.txt

comm -12 <(sort chip_list_hav31_corrected.txt) <(sort test.txt) | wc -l

#!/bin/bash
module load bioinfo/Bcftools/1.9

# Input VCF file
vcf_file="subset_RefPop_samples_ref.vcf.gz"

# Output filtered VCF file
output_vcf="subset_RefPop_samples_filtered_18000.vcf"

# File containing CHR and POS
positions_file="chip_list_hav31_corrected.txt"

bcftools view -R "$positions_file" "$vcf_file" >> "$output_vcf"

echo "Filtered VCF file created: $output_vcf"


module load devel/python/Python-3.11.1

nano filter_vcf.py

def filter_vcf(input_vcf, output_vcf, chr_pos_list):
    # Open the input VCF file for reading
    with open(input_vcf, 'r') as input_file:
        # Open the output VCF file for writing
        with open(output_vcf, 'w') as output_file:
            # Write VCF header lines to the output file
            for line in input_file:
                if line.startswith('#'):
                    output_file.write(line)
                else:
                    break

            # Iterate over each record in the input VCF file
            for line in input_file:
                fields = line.strip().split('\t')
                chrom, pos = fields[0], int(fields[1])

                # Construct the CHR_POS string for the current record
                chr_pos = f"{chrom}_{pos}"

                # Check if the CHR_POS is in the list
                if chr_pos in chr_pos_list:
                    # Write the record to the output VCF file
                    output_file.write(line)

if __name__ == "__main__":
    # Define input VCF file, output VCF file, and CHR_POS list file
    input_vcf_file = "subset_RefPop_samples_ref.vcf"
    output_vcf_file = "subset_RefPop_samples_filtered_18000.vcf"
    chr_pos_list_file = "chip_list_hav31_corrected.txt"

    # Read CHR_POS values from the list file
    with open(chr_pos_list_file, 'r') as f:
        chr_pos_list = [line.strip() for line in f]

    # Call the filter_vcf function
    filter_vcf(input_vcf_file, output_vcf_file, chr_pos_list)
	
		#scrip extractlistmarkersIDvcf.bash

#!/bin/bash
module load bioinfo/Bcftools/1.9

input_vcf="subset_RefPop_samples_ref.vcf.gz"
output_vcf="subset_RefPop_samples_filtered_18000.vcf"
position_file="chip_list_hav31_corrected.txt"
temp_dir="temp_dir"

# Créer le répertoire temporaire
mkdir -p "$temp_dir"

# Créer un tableau pour stocker les PID des processus en arrière-plan
pids=()

# Boucle de lecture du fichier position.txt
while IFS= read -r line; do
	# Utilisation de bcftools view pour extraire la région spécifiée en arrière-plan
	bcftools view "$input_vcf" -t "$line" -O v -o "$temp_dir/$(basename "$line").vcf" &

	# Stocker le PID du processus en arrière-plan
	pids+=($!)

done < "$position_file"

# Attendre que tous les processus en arrière-plan soient terminés
wait "${pids[@]}"

# Concaténer tous les fichiers temporaires dans un seul fichier VCF de sortie
bcftools concat -o "$output_vcf" "$temp_dir"/*.vcf

# Supprimer le répertoire temporaire
rm -r "$temp_dir"


	#maf001
#! /bin/bash
module load bioinfo/PLINK/1.90b7

plink --bfile ${VCFout} \
  --maf 0.01 \
  --out ${VCFout}_maf001 \
  --make-bed


	# Filtre LD03 
	# script filterLD03.bash
#! /bin/bash

module load bioinfo/PLINK/1.90b7
module load bioinfo/PLINK/2.00a4

NAME1=subset_RefPop_samples_maf001
NAME2=subset_RefPop_samples_maf001_LD03

plink --bfile ${NAME1} \
  --out ${NAME2} \
  --indep-pairwise 1749 175 0.3

plink --bfile ${NAME1} \
  --out ${NAME2}_pruned \
  --extract ${NAME2}.prune.in \
  --make-bed

plink2 --bfile subset_RefPop_samples_maf001_LD03_pruned --make-rel square --nonfounders --out subset_RefPop_samples_maf001_LD03_acp --allow-extra-chr
plink2 --bfile subset_RefPop_samples_maf001_LD03_pruned --pca --nonfounders --out subset_RefPop_samples_ref_maf001_acp --allow-extra-chr

	#23/1
	
	#script merge vcf python ?
def concatenate_vcf(vcf_file1, vcf_file2, output_file):
    # Open the first VCF file for reading
    with open(vcf_file1, 'r') as file1:
        # Read the content of the first VCF file
        content1 = file1.read()

    # Open the second VCF file for reading
    with open(vcf_file2, 'r') as file2:
        # Read the content of the second VCF file
        content2 = file2.read()

    # Concatenate the contents of the two VCF files
    concatenated_content = content1 + content2

    # Open the output file for writing
    with open(output_file, 'w') as output_file:
        # Write the concatenated content to the output file
        output_file.write(concatenated_content)

# Example usage
vcf_file1 = 'path/to/your/first.vcf'
vcf_file2 = 'path/to/your/second.vcf'
output_file = 'path/to/your/output.vcf'

concatenate_vcf(vcf_file1, vcf_file2, output_file)


###### other method

# merging requires that the files be indexed
bcftools index E756_BeeMuSe.vcf.gz
bcftools index subset_RefPop_samples_ref.vcf.gz

# merge those into a file with 6 samples
bcftools merge -Oz E756_BeeMuSe.vcf.gz subset_RefPop_samples_ref.vcf.gz > merged_data.vcf.gz


###bcftools norm --multiallelics '-any' a.vcf | bcftools norm -f '/path/to/genome.fa' > a.normed.vcf
