cd /home/proj/MDW_genomics/steepale/genome_landscape

# Assign the sample variable
Tbird=$1

# Combine the SNV and INDEL files
(grep "^#" '/home/proj/MDW_genomics/steepale/mutated_gene_annotation/data/somaticseq_vcf/'$Tbird'_somaticseq_snv_vep.vcf'; \
grep -v "^#" '/home/proj/MDW_genomics/steepale/mutated_gene_annotation/data/somaticseq_vcf/'$Tbird'_somaticseq_snv_vep.vcf'; \
grep -v "^#" '/home/proj/MDW_genomics/steepale/mutated_gene_annotation/data/somaticseq_vcf/'$Tbird'_somaticseq_indel_vep.vcf') > \
'/home/proj/MDW_genomics/steepale/mutated_gene_annotation/data/somaticseq_vcf/'$Tbird'_somaticseq_snv_indel_vep.vcf'

# Create an empty file or python script will error
touch './data/all_somatic_snv_indel_'$Tbird'.txt'

# Run python filtering script
python ./scripts/custom_filter_somatic_snvs_indels.py \
'/home/proj/MDW_genomics/steepale/mutated_gene_annotation/data/somaticseq_vcf/'$Tbird'_somaticseq_snv_indel_vep.vcf' \
'./data/all_somatic_snv_indel_'$Tbird'.txt' \
$Tbird

# Sort the sites that passed filtering to reduce redundency
(grep "^#" './data/all_somatic_snv_indel_'$Tbird'.txt'; \
grep -v "^#" './data/all_somatic_snv_indel_'$Tbird'.txt' | sort | uniq) > \
'./data/all_somatic_snv_indel_'$Tbird'_final.txt'
