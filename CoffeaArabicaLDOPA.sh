
#NCBI annot
makeblastdb -dbtype prot -in coffea_arabica.protein.faa -out coffea_arabica.protein

#sub canephora
makeblastdb -dbtype prot -in predicted.C.canephora.aa -out predicted.C.canephora
#sub eugenioides
makeblastdb -dbtype prot -in predicted.C.eugenioides.aa -out predicted.C.eugenioides

#V. NCBI annot
blastp -query c.canephora.1.10.3.1.fasta -db myAnnot/predicted.C.canephora -evalue 0.05  -outfmt 7 > subCanephora.putative.1.10.3.1.tab

#V. sub.canephora
blastp -query c.canephora.1.10.3.1.fasta -db myAnnot/predicted.C.eugenioides -evalue 0.05  -outfmt 7 > subEugenioides.putative.1.10.3.1.tab

#V. sub.eugenioides

#######################################
#DOPA descarboxylase

blastp -query PSO.DPAD.fasta -db ../myAnnot/predicted.C.canephora -evalue 0.05  -outfmt 7 > subCanephora.putative.4.1.1.28.tab

#V. sub.canephora
blastp -query PSO.DPAD.fasta -db ../myAnnot/predicted.C.eugenioides -evalue 0.05  -outfmt 7 > subEugenioides.putative.4.1.1.28.tab

#######################################
######## PPO alignment #################
######################################

#screen -S PPO
screen -r PPO

cd /home/lfmp/PaperDOPA/PPO

mafft --maxiterate 1000 --globalpair ALL.PPO.TO.Align.fasta > aligned.ALL.PPO.fasta

python3 ~/bin/ElConcatenero/ElConcatenero.py -of phylip -in aligned.ALL.PPO.fasta -o aligned.ALL.PPO

phylip seqboot

aligned.ALL.PPO.phy

mv outfile infile

phylip protdist

mv outfile infile

phylip neighbor

cp outtree intree

cp outfile infile

phylip consense

mv outtree PPO.nwk

cp PPO.nwk PPO.nwk.bk

#
blastp -task blastp-fast -remote -db nr -query ALL.PPO.TO.Align.fasta -outfmt 6 -max_target_seqs 1 -out used.PPO.blast.results.tab


#######################################
######## DDC alignment #################
######################################
#screen -S DDC
screen -r DDC

cd /home/lfmp/PaperDOPA/DDC

mafft --maxiterate 1000 --globalpair ALL.DDC.TO.align.fasta > aligned.ALL.DDC.fasta

python3 ~/bin/ElConcatenero/ElConcatenero.py -of phylip -in aligned.ALL.DDC.fasta -o aligned.ALL.DDC

phylip seqboot

aligned.ALL.DDC.phy

mv outfile infile

phylip protdist

mv outfile infile

phylip neighbor

cp outtree intree

cp outfile infile

phylip consense

mv outtree DDC.nwk

cp DDC.nwk DDC.nwk.bk
#########################################
###########Built pfam database##########
#######################################
#screen -S pfam
screen -r pfam

#wget http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.full.gz

wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz

#######################################
####### PPO hmmer domain search #######
######################################
gzip -d Pfam-A.seed.gz
ls -lh

#hmmbuild builds a profile from an alignment:
#hmmbuild Pfam-A.hmm Pfam-A.full
#hmmsearch searches a profile against a sequence database.

cd /home/lfmp/PaperDOPA/PPO

phmmer ALL.PPO.TO.Align.fasta ../../pfam/Pfam-A.seed  > putative.PPO.phmmsearch.out

#sreformat stockholm aligned.ALL.PPO.phy > aligned.ALL.PPO.sto

#hmmbuild putative.PPO.hmm aligned.ALL.PPO.sto

#hmmsearch  putative.PPO.hmm ../../pfam/Pfam-A.seed > putative.PPO.hmmsearch.out

#go to dashi and download the respective SRA libraries

############################
#screen -S PPOandDDC
screen -r PPOandDDC

mkdir PPOandDDC

cd PPOandDDC

mkdir libraries

cd libraries

#OpT_CA_5
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196520 &

#OpT_CA_4
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196523 &

#OpT_CA_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196521 &

#OpT_CA_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196522 &

#OpT_CA_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196524 &

#OpT_AC_5
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196525 &

#OpT_AC_4
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196530 &

#OpT_AC_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196526 &

#WaT_CA_5
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196529 &

#WaT_CA_4
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196527 &

#WaT_CA_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196528 &

#WaT_CA_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196531 &

#WaT_CA_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196532 &

#WaT_AC_5
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196537 &

#WaT_AC_4
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196533 &

#WaT_AC_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196536 &

#WaT_AC_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196534 &

#WaT_AC_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196535 &

#OpT_AC_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196538 &

#OpT_AC_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR11196539 &

###################
#####Field#########
###################

#PCO_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203959 &

#PCA_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203960 &

#PCA_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203961 &

 #PCA_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203962 &

#PAO_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203963 &

#PAO_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203964 &

#PAO_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203965 &

#VCO_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203966 &

#VCO_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203967 &

#VCO_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203968 &

#VCA_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203969 &

#PAA_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203970 &

#VCA_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203971 &

#VCA_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203972 &

#VAO_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203973 &

#VAO_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203974 &

#VAO_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203975 &

#VAA_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203976 &

#VAA_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203977 &

#VAA_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203978 &

#PCO_3
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203979 &

#PCO_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203980 &

#PAA_2
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203981 &

#PAA_1
~/bin/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR20203982 &

#################################
###########Rename################
#################################

#OpT_CA_5
mv SRR11196520_1.fastq OpT_CA_5_1.fq
mv SRR11196520_2.fastq OpT_CA_5_2.fq

#OpT_CA_4
mv SRR11196523_1.fastq OpT_CA_4_1.fq
mv SRR11196523_2.fastq OpT_CA_4_2.fq

#OpT_CA_3
mv SRR11196521_1.fastq OpT_CA_3_1.fq
mv SRR11196521_2.fastq OpT_CA_3_2.fq

#OpT_CA_2
mv SRR11196522_1.fastq OpT_CA_2_1.fq
mv SRR11196522_2.fastq OpT_CA_2_2.fq

#OpT_CA_1
mv SRR11196524_1.fastq OpT_CA_1_1.fq
mv SRR11196524_2.fastq OpT_CA_1_2.fq

#OpT_AC_5
mv SRR11196525_1.fastq OpT_AC_5_1.fq
mv SRR11196525_2.fastq OpT_AC_5_2.fq

#OpT_AC_4
mv SRR11196530_1.fastq OpT_AC_4_1.fq
mv SRR11196530_2.fastq OpT_AC_4_2.fq

#OpT_AC_3
mv SRR11196526_1.fastq OpT_AC_3_1.fq
mv SRR11196526_2.fastq OpT_AC_3_2.fq

#WaT_CA_5
mv SRR11196529_1.fastq WaT_CA_5_1.fq
mv SRR11196529_2.fastq WaT_CA_5_2.fq

#WaT_CA_4
mv SRR11196527_1.fastq WaT_CA_4_1.fq
mv SRR11196527_2.fastq WaT_CA_4_2.fq

#WaT_CA_3
mv SRR11196528_1.fastq WaT_CA_3_1.fq
mv SRR11196528_2.fastq WaT_CA_3_2.fq

#WaT_CA_2
mv SRR11196531_1.fastq WaT_CA_2_1.fq
mv SRR11196531_2.fastq WaT_CA_2_2.fq

#WaT_CA_1
mv SRR11196532_1.fastq WaT_CA_1_1.fq
mv SRR11196532_2.fastq WaT_CA_1_2.fq

#WaT_AC_5
mv SRR11196537_1.fastq WaT_AC_5_1.fq
mv SRR11196537_2.fastq WaT_AC_5_2.fq

#WaT_AC_4
mv SRR11196533_1.fastq WaT_AC_4_1.fq
mv SRR11196533_2.fastq WaT_AC_4_2.fq


#WaT_AC_3
mv SRR11196536_1.fastq WaT_AC_3_1.fq
mv SRR11196536_2.fastq WaT_AC_3_2.fq

#WaT_AC_2
mv SRR11196534_1.fastq WaT_AC_2_1.fq
mv SRR11196534_2.fastq WaT_AC_2_2.fq

#WaT_AC_1
mv SRR11196535_1.fastq WaT_AC_1_1.fq
mv SRR11196535_2.fastq WaT_AC_1_2.fq

#OpT_AC_2
mv SRR11196538_1.fastq OpT_AC_2_1.fq
mv SRR11196538_2.fastq OpT_AC_2_2.fq

#OpT_AC_1
mv SRR11196539_1.fastq OpT_AC_1_1.fq
mv SRR11196539_2.fastq OpT_AC_1_2.fq

###################
#####Field#########
###################

#PCO_1
mv SRR20203959_1.fastq PCO_1_1.fq
mv SRR20203959_2.fastq PCO_1_2.fq

#PCA_3
mv SRR20203960_1.fastq PCA_3_1.fq
mv SRR20203960_2.fastq PCA_3_2.fq

#PCA_2
mv SRR20203961_1.fastq PCA_2_1.fq
mv SRR20203961_2.fastq PCA_2_2.fq

 #PCA_1
mv SRR20203962_1.fastq PCA_1_1.fq
mv SRR20203962_2.fastq PCA_1_2.fq

#PAO_3
mv SRR20203963_1.fastq PAO_3_1.fq
mv SRR20203963_2.fastq PAO_3_2.fq

#PAO_2
mv SRR20203964_1.fastq PAO_2_1.fq
mv SRR20203964_2.fastq PAO_2_2.fq

#PAO_1
mv SRR20203965_1.fastq PAO_1_1.fq
mv SRR20203965_2.fastq PAO_1_2.fq


#VCO_3
mv SRR20203966_1.fastq VCO_3_1.fq
mv SRR20203966_2.fastq VCO_3_2.fq

#VCO_2
mv SRR20203967_1.fastq VCO_2_1.fq
mv SRR20203967_2.fastq VCO_2_2.fq


#VCO_1
mv SRR20203968_1.fastq VCO_1_1.fq
mv SRR20203968_2.fastq VCO_1_2.fq

#VCA_3
mv SRR20203969_1.fastq VCA_3_1.fq
mv SRR20203969_2.fastq VCA_3_2.fq

#PAA_3
mv SRR20203970_1.fastq PAA_3_1.fq
mv SRR20203970_2.fastq PAA_3_2.fq

#VCA_2
mv SRR20203971_1.fastq VCA_2_1.fq
mv SRR20203971_2.fastq VCA_2_2.fq

#VCA_1
mv SRR20203972_1.fastq VCA_1_1.fq
mv SRR20203972_2.fastq VCA_1_2.fq

#VAO_3
mv SRR20203973_1.fastq VAO_3_1.fq
mv SRR20203973_2.fastq VAO_3_2.fq

#VAO_2
mv SRR20203974_1.fastq VAO_2_1.fq
mv SRR20203974_2.fastq VAO_2_2.fq

#VAO_1
mv SRR20203975_1.fastq VAO_1_1.fq
mv SRR20203975_2.fastq VAO_1_2.fq

#VAA_3
mv SRR20203976_1.fastq VAA_3_1.fq
mv SRR20203976_2.fastq VAA_3_2.fq

#VAA_2
mv SRR20203977_1.fastq VAA_2_1.fq
mv SRR20203977_2.fastq VAA_2_2.fq

#VAA_1
mv SRR20203978_1.fastq VAA_1_1.fq
mv SRR20203978_2.fastq VAA_1_2.fq

#PCO_3
mv SRR20203979_1.fastq PCO_3_1.fq
mv SRR20203979_2.fastq PCO_3_2.fq

#PCO_2
mv SRR20203980_1.fastq PCO_2_1.fq
mv SRR20203980_2.fastq PCO_2_2.fq

#PAA_2
mv SRR20203981_1.fastq PAA_2_1.fq
mv SRR20203981_2.fastq PAA_2_2.fq

#PAA_1
mv SRR20203982_1.fastq PAA_1_1.fq
mv SRR20203982_2.fastq PAA_1_2.fq


##########################
#Use the quantification from previous experiments

mkdir expression
cd expression

cp -r ~/doc/WGCNA/counts ./ -r
cp ~/doc/WGCNA/Targets.txt ./

R
library("edgeR")
targets <- readTargets()
names <- targets$description
matrix_input <- readDGE(targets, comment.char = "!")

#remove meta Tags
MetaTags <- grep("^__", rownames(matrix_input))
matrix_input <- matrix_input[-MetaTags, ]


reads_before <- sum(matrix_input$counts)
#remove low expressed genes
rnaseqmatrix <- matrix_input$counts
rnaseqmatrix <- rnaseqmatrix[rowMeans(rnaseqmatrix) >=5,]

sum(rnaseqmatrix)/reads_before

conditions = matrix_input$samples[,2]
analysis_matrix <- DGEList(counts = rnaseqmatrix,group = conditions)
colnames(analysis_matrix$counts) <- names
design <- model.matrix(~0+group, data=analysis_matrix$samples)
colnames(design) <- levels(analysis_matrix$samples$group)


#NORMALIZATIONS
analysis_matrix <- calcNormFactors(analysis_matrix)

#To estimate common dispersion:
analysis_matrix <- estimateGLMCommonDisp(analysis_matrix, design)
#To estimate trended dispersions:
analysis_matrix <- estimateGLMTrendedDisp(analysis_matrix, design)
#To estimate tagwise dispersions:
analysis_matrix <- estimateGLMTagwiseDisp(analysis_matrix, design)

#Fit a negative binomial generalized log-linear model to the read counts for each gene.
fit <- glmFit(analysis_matrix,design)


rownames(analysis_matrix$counts)

DDCs <- c("SubC.c_24313","SubC.c_24311","SubC.e_31316","SubC.e_38178", "SubC.e_38125","SubC.e_3477")

names(DDCs) <- c("DDC.CAR1","DDC.CAR2","DDC.CAR3","DDC.CAR4","DDC.CAR5","DDC.CAR6")

PPOs <- c("SubC.c_14465","SubC.c_14404","SubC.c_14467", "SubC.c_14406","SubC.c_14369","SubC.c_14464","SubC.c_4204","SubC.e_5336")

names(PPOs) <- c("PPO.CAR1","PPO.CAR2","PPO.CAR3","PPO.CAR4", "PPO.CAR5","PPO.CAR6", "PPO.CAR7", "PPO.CAR8")

selectedSeqs <- c(PPOs,DDCs )

heatmap.cpm.matrix <- cpm(analysis_matrix)

selected.hm.cpm.matrix <- heatmap.cpm.matrix[rownames(heatmap.cpm.matrix) %in% selectedSeqs,]

selected.hm.cpm.matrix <- t(selected.hm.cpm.matrix)

colnames(selected.hm.cpm.matrix) <- c("PPO.CAR1-5","PPO.CAR6","DDC.CAR6","PPO.CAR7/8")

library(gplots)


colors<- colorRampPalette(c("white","#BF7F58"))(16)

pdf("heatmaptest.pdf",w=10/2.54,h=36/2.54)

heatmap.2(log(selected.hm.cpm.matrix+1),
Rowv=FALSE,
Colv=FALSE,
trace = "none",
margins = c(8, 8),
key=FALSE,
keysize = 0.3,
key.title="",
col =colors,
density.info="none",
dendrogram="none",
colsep = c(1,2,3),
rowsep = c(seq(1,24,1)),
sepwidth=c(0.001,0.001),
sepcolor="black",
cexRow=1.2,
cexCol=1.3
)

dev.off()

#FOR the key!!!
pdf("keymaptest.pdf",w=9,h=5)

heatmap.2(log(selected.hm.cpm.matrix+1),
Rowv=FALSE,
Colv=FALSE,
trace = "none",
margins = c(5, 5),
key=TRUE,
keysize = 1.5,
key.title="",
col =colors,
density.info="none",
dendrogram="none",
colsep = c(1,2,3),
rowsep = c(seq(1,24,1)),
sepwidth=c(0.001,0.001),
sepcolor="black",
cexRow=1.2,
cexCol=1.3
)

dev.off()
