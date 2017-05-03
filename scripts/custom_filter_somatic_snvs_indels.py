import sys
import os
import re
import subprocess
from subprocess import check_output

# General Strategy of script:
# Script will take a single vep file from one tumor sample and search for each called variant
# across it's own tumor as well as all germline bam files with samtools mpileup. Then filters
# will be applied.

# Filters include:
# Variant must be in tumor bam file from which it was originally called by a variant caller
# Variants must have a variant allele frequency greater than or equal to 0.05
# Tumor file and germline bam pairs much each have a coverage of 4 at variant site
# Variant must not be found in more than 0 germline files (in any)

# Input SNV and INDEL file
vep_file = open(sys.argv[1], 'r')

# Output file
outfile = open(sys.argv[2], 'w')

# Tumor sample
tbird = sys.argv[3]

# Reference files
tumor_birds_file = "/home/users/a.steep/databases/samples/tumor_sample_dnaseq_list_NNN-N_SN.txt"
germline_birds_file = "/home/users/a.steep/databases/samples/germline_sample_dnaseq_list_NNN-N_SN.txt"

# Create a dictionary of germline bam files
germline_bam_file = {}
for gbird in open(germline_birds_file):
	gbird = gbird.rstrip()
	germline_bam_file[gbird] = '/home/proj/MDW_genomics/xu/final_bam/' + gbird + '_Bwa_RG_dedupped_realigned.bam'

# Create a dictionary of tumor bam files
tum_bam_file = {}
for tum_bird in open(tumor_birds_file):
	tum_bird = tum_bird.rstrip()
	tum_bam_file[tum_bird] = '/home/proj/MDW_genomics/xu/final_bam/' + tum_bird + '_Bwa_RG_dedupped_realigned.bam'

# Create another dictionary with the exact tumor bam
tumor_bam_file = {}
tumor_bam_file[tbird] = '/home/proj/MDW_genomics/xu/final_bam/' + tbird + '_Bwa_RG_dedupped_realigned.bam'

# Write header to outfile
outfile.write('#CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'MUT' + '\t' + 'IMPACT' + '\t' + 'SYMBOL' + '\t' + 'GENE_ID' + '\t' + 'TSN' + '\t' + 'SAMPLE' + '\t' + 'VAC' + '\t' + 'VAF' + '\t' + 'CALLER_CONSENSUS' + '\n')

# Iterate over lines in input file and perform custom filtering
for vep_line in vep_file:
	if vep_line[0] != '#':
		vep_line = vep_line.rstrip()
		vep_cols = vep_line.split('\t')
		vep_chr = vep_cols[0]
		vep_pos = vep_cols[1]
		vep_ref = vep_cols[3]
		vep_alt = vep_cols[4]
		vep_snv = vep_chr + '\t' + vep_pos + '\t' + vep_ref + '\t' + vep_alt
		vep_info = vep_cols[7]
		vep_format = vep_cols[8]
		vep_normal = vep_cols[9]
		vep_tumor = vep_cols[10]
		vep_sample = tbird.rstrip()
		if len(vep_alt) > len(vep_ref):
			var_type = 'INS'
		elif len(vep_alt) < len(vep_ref):
			var_type = 'DEL'
		elif len(vep_alt) == len(vep_ref):
			var_type = 'SNV'
		# Create samtools annotation for insertion
		if var_type == 'INS':
			vep_ins_length = len(vep_alt) - len(vep_ref) 
			vep_ins_samtools_str = '+' + str(vep_ins_length) + vep_alt[-vep_ins_length:]
		# Create samtools annotation for deletion
		elif var_type == 'DEL':
			vep_del_length = len(vep_ref) - len(vep_alt)
			vep_del_samtools_str = 'N' * vep_del_length
		# If the variants are somatic
		if re.search('SOMATIC', vep_info.split(';')[0]):
			vep_somat = vep_info.split(';')[0]
		else:
			vep_somat = 'NA'
		if re.search('MVJSDU', vep_info.split(';')[0]):
			vep_tools = vep_info.split(';')[0]
		elif re.search('MVJSDU', vep_info.split(';')[1]):
			vep_tools = vep_info.split(';')[1]
		if re.search('NUM_TOOLS', vep_info.split(';')[1]):
			vep_tool_num = vep_info.split(';')[1]
		elif re.search('NUM_TOOLS', vep_info.split(';')[2]):
			vep_tool_num = vep_info.split(';')[2]
		if re.search('CSQ=', vep_info.split(';')[2]):
			vep_info_ann = vep_info.split(';')[2].split('SQ=')[1]
			vep_info_ann_num = vep_info_ann.count(',') + 1
			info_ann = []
		elif re.search('CSQ=', vep_info.split(';')[3]):
			vep_info_ann = vep_info.split(';')[3].split('SQ=')[1]
			vep_info_ann_num = vep_info_ann.count(',') + 1
			info_ann = []
		for n in range(vep_info_ann_num):
			info_ann.append(n)
			info_ann[n] = vep_info_ann.split(',')[n]
			info_ann2read = info_ann[n]
			info_cols = info_ann2read.split('|')
			vep_allele = info_cols[0]
			vep_cons = info_cols[1]
			vep_impact = info_cols[2]
			vep_symbol = info_cols[3]
			vep_geneid = info_cols[4]
			vep_feat_type = info_cols[5]
			vep_feature = info_cols[6]
			vep_biotype = info_cols[7]
			vep_exon = info_cols[8]
			vep_intron = info_cols[9]
			vep_HGVSc = info_cols[10]
			vep_HGVSp = info_cols[11]
			vep_cDNA_pos = info_cols[12]
			vep_CDS_pos = info_cols[13]
			vep_protein_pos = info_cols[14]
			vep_aminos = info_cols[15]
			vep_codons = info_cols[16]
			vep_existing_var = info_cols[17]
			vep_distance = info_cols[18]
			vep_strand = info_cols[19]
			vep_flags = info_cols[20]
			vep_symbol_source = info_cols[21]
			vep_HGNC_ID = info_cols[22]
			vep_tsl = info_cols[23]
			vep_appris = info_cols[24]
			vep_ccds = info_cols[25]
			vep_ensp = info_cols[26]
			vep_swissprot = info_cols[27]
			vep_trembl = info_cols[28]
			vep_uniparc = info_cols[29]
			vep_sift = info_cols[30]
			vep_domains = info_cols[31]
			vep_hgvs_offset = info_cols[32]
			# Create counter for each germline and tumor bam file with variant at sufficient VAF and set to zero
			germ_sam_set = set()
			# Create coverage variables
			same_tumor_cov = "no"
			gleich_germline_cov = "no"
			# Reset all variables
			gleich_tumor_status = 'no'
			same_germline_status = 'no'
			tumor_in_germline_out = 'no'
			# Search each germline bam file for the variant
			for g_bird, germline_bam in germline_bam_file.items():
				# Set all counting variables to zero
				g_mpu_bases = ''
				g_mpu_depth = 0
				g_VAC = 0
				g_VAF = 0
				# Use samtools mpileup to show the actually mapped bases in the original BAM files to downcheck the accuracy of the callers
				# Each base needs to have a base quality of atleast 20 and a mapping quality of atleast 20
				g_samtools_cmd = 'samtools mpileup --min-MQ 20 --min-BQ 20 -r ' + vep_chr+':'+vep_pos+'-'+vep_pos+' '+germline_bam
				# Use subprocess.Popen to ellicit shell commands 
				g_samtools_proc = subprocess.Popen([g_samtools_cmd], stdout=subprocess.PIPE, shell=True)
				# Use communicate to capture the output in a 'bytes' object
				(g_out, g_err) = g_samtools_proc.communicate()
				# Decode the 'bytes' object to a string
				g_mpu_out = g_out.decode("utf-8")
				g_mpu = g_mpu_out.rstrip()
				# If the germline sample and tumor samples match and there is an output from samtools mpileup
				if g_bird[0:3] == vep_sample[0:3] and g_mpu != '' and g_mpu.split('\t')[3] != '0':
					# Collect variables on matching germline sample
					same_bird = vep_sample
					same_mpu_chr = g_mpu.split('\t')[0]
					same_mpu_pos = g_mpu.split('\t')[1]
					same_mpu_ref = g_mpu.split('\t')[2]
					same_mpu_depth = int(g_mpu.split('\t')[3])
					same_mpu_bases = g_mpu.split('\t')[4].upper()
					if var_type == 'INS':
						same_VAC = same_mpu_bases.count(vep_ins_samtools_str)
					elif var_type == 'DEL':
						same_VAC = same_mpu_bases.count(vep_del_samtools_str)
					elif var_type == 'SNV':
						same_VAC = same_mpu_bases.count(vep_alt)
					same_VAF = same_VAC/same_mpu_depth
					# If the Variant allele frequency and coverage of the germline is atleast 0.05 and 4, respectively, then consider site
					if same_VAF >= 0.05 and same_mpu_depth >= 4:
						germ_sam_set.add(g_bird)
						same_germline_status = 'yes'
					# Create a variable to express adequate coverage for paired germline sample
					if same_mpu_depth >= 4:
						same_tumor_cov = "yes"
					else:
						same_tumor_cov = "no"
				# Else if the germline sample has no coverage, then pass
				elif g_mpu == '' or int(g_mpu.split('\t')[3]) == 0:
					pass
				# For germline samples that displayed coverage but did not match to the tumor sample as a paired sample
				else:
					g_mpu_chr = g_mpu.split('\t')[0]
					g_mpu_pos = g_mpu.split('\t')[1]
					g_mpu_ref = g_mpu.split('\t')[2]
					g_mpu_depth = int(g_mpu.split('\t')[3])
					g_mpu_bases = g_mpu.split('\t')[4].upper()
					if var_type == 'INS':
						g_VAC = g_mpu_bases.count(vep_ins_samtools_str)
					elif var_type == 'DEL':
						g_VAC = g_mpu_bases.count(vep_del_samtools_str)
					elif var_type == 'SNV':
						g_VAC = g_mpu_bases.count(vep_alt)
					g_VAF = g_VAC/g_mpu_depth
					# Add to counter for each germline file with variant at sufficient VAF
					if g_VAF >= 0.10:
						germ_sam_set.add(g_bird)
			# Create empty dictionary for tumor values
			# Create a dictionary with samples as keys and a list with VAC and VAF as value for input tumor sample
			sam2VACVAF = {}
			# Create empty list variables
			tumor_sample_list = []
			VAC_list = []
			VAF_list = []
			# Search input tumor bam for each somatic called variant
			for t_bird, tumor_bam in tum_bam_file.items():
				# Use samtools mpileup to show the actually mapped bases in the original BAM files to downcheck the accuracy of the callers
				# Each base needs to have a base quality of atleast 20 and a mapping quality of atleast 20
				t_samtools_cmd = 'samtools mpileup --min-MQ 20 --min-BQ 20 -r ' + vep_chr+':'+vep_pos+'-'+vep_pos+' '+tumor_bam
				t_samtools_proc = subprocess.Popen([t_samtools_cmd], stdout=subprocess.PIPE, shell=True)
				(t_out, t_err) = t_samtools_proc.communicate()
				t_mpu_out = t_out.decode("utf-8")
				t_mpu = t_mpu_out.rstrip()
				# If the tumor is the input tumor and there is sufficient coverage
				if t_bird == vep_sample and t_mpu != '' and t_mpu.split('\t')[3] != '0':
					gleich_bird = vep_sample
					gleich_mpu_chr = t_mpu.split('\t')[0]
					gleich_mpu_pos = t_mpu.split('\t')[1]
					gleich_mpu_ref = t_mpu.split('\t')[2]
					gleich_mpu_depth = int(t_mpu.split('\t')[3])
					gleich_mpu_bases = t_mpu.split('\t')[4].upper()
					if var_type == 'INS':
						gleich_VAC = gleich_mpu_bases.count(vep_ins_samtools_str)
					elif var_type == 'DEL':
						gleich_VAC = gleich_mpu_bases.count(vep_del_samtools_str)
					elif var_type == 'SNV':
						gleich_VAC = gleich_mpu_bases.count(vep_alt)
					gleich_VAF = gleich_VAC/gleich_mpu_depth
					# If there variant allele frequency and coverage are sufficient
					if gleich_VAF >= 0.05 and gleich_mpu_depth >= 4:
						gleich_tumor_status = 'yes'
						sam2VACVAF[t_bird] = [str(gleich_VAC), str(gleich_VAF)[0:5]]
					if gleich_mpu_depth >= 4:
						gleich_germline_cov = "yes"
					else:
						gleich_germline_cov = "no"
				# Else if tumor is not the input tumor but there is sufficient coverage
				elif t_bird != vep_sample and t_mpu != '' and t_mpu.split('\t')[3] != '0':
					anders_bird = t_bird
					anders_mpu_chr = t_mpu.split('\t')[0]
					anders_mpu_pos = t_mpu.split('\t')[1]
					anders_mpu_ref = t_mpu.split('\t')[2]
					anders_mpu_depth = int(t_mpu.split('\t')[3])
					anders_mpu_bases = t_mpu.split('\t')[4].upper()
					if var_type == 'INS':
						anders_VAC = anders_mpu_bases.count(vep_ins_samtools_str)
					elif var_type == 'DEL':
						anders_VAC = anders_mpu_bases.count(vep_del_samtools_str)
					elif var_type == 'SNV':
						anders_VAC = anders_mpu_bases.count(vep_alt)
					anders_VAF = anders_VAC/anders_mpu_depth
					# If there variant allele frequency and coverage are sufficient
					if anders_VAF >= 0.05 and anders_mpu_depth >= 4:
						#anders_tumor_status = 'yes'
						sam2VACVAF[t_bird] = [str(anders_VAC), str(anders_VAF)[0:5]]
			# If stats pertain to the input tumor with site of adequate variant allele frequency and coverage,
			# AND if the compared germline samples are not the matching samples of the tumor
			if gleich_tumor_status == 'yes' and same_germline_status == 'no':
				tumor_in_germline_out = 'yes'
			else:
				tumor_in_germline_out = 'no'
			if vep_symbol == '':
				vep_symbol = 'NA'
			# print dictionary values is correct order
			for sample, VACVAF in sam2VACVAF.items():
				tumor_sample_list.append(sample)
				VAC_list.append(VACVAF[0])
				VAF_list.append(VACVAF[1])
			# Turn the lists into proper variables
			tumor_samples = ';'.join(map(str,tumor_sample_list))
			tumor_VAC = ';'.join(map(str,VAC_list))
			tumor_VAF = ';'.join(map(str,VAF_list))
			# Perform final filters if..., 
			# found at relevant variant allele freq in tumors
			# not found in germline samples
			# the input tumor sample shows adequate coverage
			# the paired germline sample shows adequate coverage
			if len(sam2VACVAF) > 0 and len(germ_sam_set) == 0 and gleich_germline_cov == 'yes' and same_tumor_cov == 'yes' and tumor_in_germline_out == 'yes':
				outfile.write(vep_chr + '\t' + vep_pos + '\t' + vep_ref + '\t' + vep_alt + '\t' + vep_cons + '\t' + vep_impact + '\t' + vep_symbol + '\t' + vep_geneid + '\t' + str(len(sam2VACVAF)) + '\t' + tumor_samples + '\t' + tumor_VAC + '\t' + tumor_VAF + '\t' + vep_tool_num + '\n')
outfile.close()
