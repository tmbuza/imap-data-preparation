INPUTDIR="../imap-bioinformatics-qiime2/qiime2_process/export/feature-table.tsv"
OUTDIR="data/qiime2"

echo PROGRESS: Importing QiiME2 features

mkdir -p "${OUTDIR}"

cp "${INPUTDIR}" "${OUTDIR}/qiime2_features.tsv"
