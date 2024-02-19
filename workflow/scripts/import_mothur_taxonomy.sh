INPUTDIR="../imap-bioinformatics-mothur/mothur_process/otu_analysis/final.opti_mcc.0.03.cons.taxonomy"
OUTDIR="data/mothur"

echo PROGRESS: Importing MOTHUR taxonomy

mkdir -p "${OUTDIR}"

cp "${INPUTDIR}" "${OUTDIR}/mothur_taxonomy.tsv"