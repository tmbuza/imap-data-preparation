INPUTDIR="../imap-qiime2-bioinformatics/qiime2_process/export/taxonomy.tsv"
OUTDIR="data/qiime2"

echo PROGRESS: Importing QIIME2 taxonomy

mkdir -p "${OUTDIR}"

cp "${INPUTDIR}" "${OUTDIR}/qiime2_taxonomy.tsv"