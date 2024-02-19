METADATA="../imap-mothur-bioinformatics/resources/metadata/mothur_test_sample_metadata.tsv"
DESIGN="../imap-mothur-bioinformatics/resources/metadata/mothur_test_design_file.tsv"
OUTDIR="data/mothur"

echo PROGRESS: Importing MOTHUR metadata

mkdir -p "${OUTDIR}"

cp "${METADATA}" "${OUTDIR}/mothur_sample_metadata.tsv"
cp "${DESIGN}" "${OUTDIR}/mothur_design.tsv"