### Post analysis in R


setwd('../PROSPECT')
library(tidyverse)
library(stringr)
library(VennDiagram)
library(CNTools)
library(MutationalPatterns)
library(BSgenome)
library(pheatmap)
library(RColorBrewer)

############### Functions ################## 

my_sum_table = function(data){
		data$patientID = factor(data$patientID)
		x1 = as.numeric(table(data$patientID))
		y1 = cbind.data.frame(names(table(data$patientID)), x1)
		colnames(y1) = c('tumour_name', 'mutation_counts')
		mb.dat = mutate(y1, mb = mutation_counts / 44100000 * 1000000)
}

 
Corner_text = function(text, location="topright"){
  legend(location,legend=text, bty ="n", pch=NA, cex = .9) 
}

add_key = function (x) {
	transform (x, key=paste0(read, "-", variant, "-", gene, "-", amino_acid))
}

############### Analysis starts here ################## 
dat = read.table("data", sep = "\t", header = T, quote = "")

patientID = str_sub(dat[,10], start=1, end=8)
dat[c("patientID")] = patientID

# Add unique key
dat = transform(dat, unique_key=paste0(patientID, ":", key), SAMPLE = paste0(sample_name, "-", ref_name))

# Additional columns
dat[c("t_total_count", "n_total_count", "normal_f")] = NA  
dat$t_total_count = dat$t_ref_count + dat$t_alt_count 
dat$n_total_count = dat$n_ref_count + dat$n_alt_count 
dat$normal_f = (dat$n_alt_count/dat$n_total_count)  


### New filtering criteria Set1  ### FINAL decision#####
filtered_dat = subset(dat, t_total_count>50 & n_total_count>25 & tumor_f>0.05 & normal_f <0.01)					
					

# Mutations per MB for exome sequencing
mb_dat = my_sum_table(filtered_dat)

# $ of shared mutations
venn = draw.pairwise.venn(area1 = length(filtered_dat$unique_key), 
                          area2 = length(filtered_dat2$unique_key), 
                          cross.area = length(intersect(filtered_dat$unique_key, filtered_dat2$unique_key)), 
                          category = c("NB", "TB"), 
                          cat.fontfamily = rep("sans", 2), 
                          fontfamily = rep("sans", 3), 
                          fill = c("#3288bd", "#d53e4f"), 
                          cat.pos = c(0, 0),
                          cat.dist = rep(0.025, 2),
                          cex = 1, 
                          cat.cex = 0.7,
                          lwd = rep(2, 2))

### Non-silent mutations
dat.nonsilent = subset(rescued_filtered_dat, exonicfunc=="nonsynonymous SNV" | exonicfunc=="stoploss" | exonicfunc=="stopgain")   


### Search for canonical mutations
# load the known genes
tsg_oncogenes = read.table('../Cancer_genes/tier1_tsg_oncogenes.tsv', sep ='\t', header = T, quote = "")
cancer.genes = read.table("../Cancer_genes/cancer_gene.tsv", sep = "\t", header =T)


# Sets of curated variants from DoCM 
curated_variants = read.table("../Cancer_genes/wustl_curated_variants.tsv", sep = "\t", header = TRUE) 
CIVic_Knowledgebase =read.table("../Cancer_genes/CIVic_Knowledgebase.tsv", sep = "\t", header=TRUE)
Drug_gene_knowledgebase_db = read.table("../Cancer_genes/Drug_gene_knowledgebase_db.tsv", sep = "\t", header=TRUE)
Kin_driver = read.table("../Cancer_genes/Kin-driver.tsv", sep = "\t", header=TRUE)
Literature = read.table("../Cancer_genes/Literature.tsv", sep = "\t", header=TRUE)
My_cancer_genome = read.table("../Cancer_genes/My_cancer_genome.tsv", sep = "\t", header=TRUE)
Oncomap_variants = read.table("../Cancer_genes/Oncomap_variants.tsv", sep = "\t", header=TRUE)
Pan_cancer_recurrent_hotspots = read.table("../Cancer_genes/Pan-cancer_recurrent_hotspots.tsv", sep = "\t", header=TRUE)
WashU_malignancy_list = read.table("../Cancer_genes/WashU_malignancy_list.tsv", sep = "\t", header=TRUE)

# Assign the type and key in the database
curated_variants$type = "curated"	
CIVic_Knowledgebase$type = "CIVic_Knowledgebase"
Drug_gene_knowledgebase_db$type = "Drug_gene_knowledgebase_db"
Kin_driver$type = "Kin-driver"
Literature$type = "Literature"
My_cancer_genome$type = "My_cancer_genome"
Oncomap_variants$type = "Oncomap_variants"
Pan_cancer_recurrent_hotspots$type = "Pan_cancer_recurrent_hotspots"
WashU_malignancy_list$type = "WashU_malignancy"

# Insert key
curated_variants = add_key(curated_variants)
CIVic_Knowledgebase = add_key(CIVic_Knowledgebase)
Drug_gene_knowledgebase_db = add_key(Drug_gene_knowledgebase_db)
Kin_driver = add_key(Kin_driver)
Literature = add_key(Literature)
My_cancer_genome = add_key(My_cancer_genome)
Oncomap_variants = add_key(Oncomap_variants)
Pan_cancer_recurrent_hotspots = add_key(Pan_cancer_recurrent_hotspots)
WashU_malignancy_list = add_key(WashU_malignancy_list)


# key for overlapping cancer canonical mutation
dat.nonsilent = transform(dat.nonsilent, cancerGenes_key=paste0(ref_allele, "-", alt_allele, "-", gene.knowngene, "-", aaannotation))

# Any canonical mutations - using the wustl table
# annotate whether if the genes are oncogenes or not?
dat.nonsilent$type = tsg_oncogenes$Type[match(dat.nonsilent$entrez_gene_id, tsg_oncogenes$Entrez_GeneId)]


# annotate whether if the genes are curated or not?
dat.nonsilent$curated = curated_variants$type[match(dat.nonsilent$cancerGenes_key, curated_variants$key)]  # example = C-G-LOC399753-p.Q27H
# annotate with the respective database
dat.nonsilent$CIVic_Knowledgebase = CIVic_Knowledgebase$type[match(dat.nonsilent$cancerGenes_key, CIVic_Knowledgebase$key)]
dat.nonsilent$Drug_gene_knowledgebase_db = Drug_gene_knowledgebase_db$type[match(dat.nonsilent$cancerGenes_key, Drug_gene_knowledgebase_db$key)]
dat.nonsilent$Kin_driver = Kin_driver$type[match(dat.nonsilent$cancerGenes_key, Kin_driver$key)]
dat.nonsilent$Literature = Literature$type[match(dat.nonsilent$cancerGenes_key, Literature$key)]
dat.nonsilent$My_cancer_genome = My_cancer_genome$type[match(dat.nonsilent$cancerGenes_key, My_cancer_genome$key)]
dat.nonsilent$Oncomap_variants = Oncomap_variants$type[match(dat.nonsilent$cancerGenes_key, Oncomap_variants$key)]
dat.nonsilent$Pan_cancer_recurrent_hotspots = Pan_cancer_recurrent_hotspots$type[match(dat.nonsilent$cancerGenes_key, Pan_cancer_recurrent_hotspots$key)]
dat.nonsilent$WashU_malignancy_list = WashU_malignancy_list$type[match(dat.nonsilent$cancerGenes_key, WashU_malignancy_list$key)]  


dat_canonical_mutations = subset(dat.nonsilent, curated == "curated" | (type == "tsg" & exonicfunc == "stopgain"))

# Measuring mutational signatures contribution
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

vcf_files = list.files("../TB", pattern = ".vcf", full.names = TRUE)
sample_names = read.table('../TB/sampleName.txt', sep = '\t', header = T)

# load a list of VCF files:
vcfs = read_vcfs_as_granges(vcf_files, sample_names$sampleName, genome = "hg19")
# mutation matrix
mut_mat = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

# Fit 96 mutation profiles to known signatures
sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"   # COSMIC 30 signatures
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)

# Reorder (to make the order of the trinucleotide changes the same)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]

# Only signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])

# select cancer signature of interest
cancer_signatures_of_interest = cbind(cancer_signatures[,1],cancer_signatures[,2], cancer_signatures[,3], cancer_signatures[,4], cancer_signatures[,6], cancer_signatures[,13])
colnames(cancer_signatures_of_interest) = c("Signature.1", "Signature.2", "Signature.3", "Signature.4", "Signature.6", "Signature.13")
fit_res = fit_to_signatures(mut_mat, cancer_signatures_of_interest)

# Convert to prop
df = as.data.frame(fit_res$contribution)
prop_df = prop.table(as.matrix(df),2)

# copy number analysis
cnv_dat = read.table('../sequenza/all.segments.NB.tsv', sep = '\t', header = T, quote ="")
cnv_dat = cnv_dat %>% select(ID = ID, chrom=chromosome, loc.start = start.pos, loc.end = end.pos, seg.mean = depth.ratio)

data("geneInfo")
datSeg = CNSeg(cnv_dat[which(is.element(cnv_dat[, "ID"], unique(cnv_dat[, "ID"]))),])
rdByGene.dat = getRS(datSeg , by = "gene", imput = FALSE, XY = FALSE, geneMap = geneInfo, what="mean")
rdByGene.dat.Out = rs(rdByGene.NB)

##### CNV unsupervise clustering
all_cnv = cbind(rdByGene.dat.Out[5:81], rrdByGene.dat2.Out[5:81])
all_cnv = all_cnv[apply(all_cnv[c(2:153)],1,function(z) any(z!=0)),]
mat_cnv = as.matrix(all_cnv)

colours = brewer.pal(8, "Set1")

pheatmap(mat_cnv, scale="none" , cluster_rows = FALSE, cluster_cols = TRUE,
         annotation_col = clab, annotation_legend = TRUE, legend = TRUE, 
         clustering_distance_cols = "euclidean", clustering_method = "ward.D",
         show_colnames = T, show_rownames = F, use_raster = T, fontsize = 8, filename = "all_CNV_no_scaled.pdf")
         

############### Plots ################## 
pdf('mutational_burden.pdf')
boxplot(mb.dat$mb, mb.df2$mb, ylab = "Mutations/MB", names= c("NB", "TB"), main="All mutations", staplelty=0)
tmp = wilcox.test(mb.dat$mb, mb.df2$mb)
Corner_text(text=paste("p-value:", round(tmp$p.value, 4)),location= "bottomright")
dev.off()
