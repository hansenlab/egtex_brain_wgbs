# Download data for use with LDSC (https://github.com/bulik/ldsc) with hg38
# co-ordinates.
# Peter Hickey
# 2020-04-20

# Here, we download the 'baseline-LD model v2.2' data, as recommended in
# https://data.broadinstitute.org/alkesgroup/LDSCORE/readme_baseline_versions,
# against the GRCh38 reference genome.
# Much of the data is in common with the data used in
# https://github.com/hansenlab/BrainEpigenomeNN and so does not need to be
# re-downloaded.
# We can re-use the already downloaded data; data may be re-used 'as is' or by
# lifting over from hg19 to hg38, depending on whether physical co-ordinates
# are needed or just rs-IDs.

# Shell variables --------------------------------------------------------------

LDSC_DIR="/dcl01/hansen/data/ldsc_GRCh38"
mkdir -p ${LDSC_DIR}

# NOTE: Files in these directories can be re-used (either 'as is' or by lifting
#       data over from hg19 to hg38).
PHASE1=${LDSC_DIR}/Phase1
ln -s /dcl01/hansen/data/ldsc/Phase1 ${PHASE1}
GWASSS=${LDSC_DIR}/GWAS_summary_stats
ln -s /dcl01/hansen/data/ldsc/GWAS_summary_stats ${GWASSS}
MUNGEDSS=${LDSC_DIR}/munge_sumstats
ln -s /dcl01/hansen/data/ldsc/munge_sumstats ${MUNDGEDSS}

# Download data for use with LDSC ----------------------------------------------

# Baseline model LD scores; this contains the `baseline.*` files
wget -O ${LDSC_DIR}/baselineLD_v2.2.tgz \
  https://data.broadinstitute.org/alkesgroup/LDSCORE/GRCh38/baselineLD_v2.2.tgz

# Allele frequencies; this contains the `1000G.EUR.hg38.*` files and PLINK files
wget -O ${LDSC_DIR}/plink_files.tgz \
  https://data.broadinstitute.org/alkesgroup/LDSCORE/GRCh38/plink_files.tgz

# Regression weights; this contains the `weights.*` files
wget -O ${LDSC_DIR}/weights.tgz \
  https://data.broadinstitute.org/alkesgroup/LDSCORE/GRCh38/weights.tgz

# Download baseline BED files --------------------------------------------------

# NOTE: These use hg19 co-ordinates
wget -O ${LDSC_DIR}/baselineLD_v2.2_bedfiles.tgz \
     https://data.broadinstitute.org/alkesgroup/LDSCORE/baselineLD_v2.2_bedfiles.tgz

# Inflate and extract downloaded data ------------------------------------------

tar xvfz ${LDSC_DIR}/baselineLD_v2.2.tgz -C ${LDSC_DIR}
tar xvfz ${LDSC_DIR}/plink_files.tgz -C ${LDSC_DIR}
tar xvfz ${LDSC_DIR}/weights.tgz -C ${LDSC_DIR}
# NOTE: This is really a bzip2 compressed tar archive (despite the extension).
tar xfvj ${LDSC_DIR}/baselineLD_v2.2_bedfiles.tgz -C ${LDSC_DIR}

# Process downloaded data ------------------------------------------------------

# NOTE: This advice came from Steven Gazal (2020-04-17). `list.txt` contains
#       the SNPs used to created the `baselineLD_v2.2` files.
awk '{if ($1!="SNP") {print $1} }' ${PHASE1}/w_hm3.snplist > ${LDSC_DIR}/list.txt
