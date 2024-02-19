INPUTDIR="../imap-bioinformatics-mothur/mothur_process/otu_analysis/final.opti_mcc.shared"
OUTDIR="data/mothur"

echo PROGRESS: Importing MOTHUR shared file

mkdir -p "${OUTDIR}"

cp "${INPUTDIR}" "${OUTDIR}/mothur_otutable.tsv"
