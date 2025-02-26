#!/bin/bash

SCRIPTPATH="$(dirname "$(readlink -f "$0")")"
AV_DB="/opt/annovar/humandb/hg38/"

if [ -f "clinvar.vcf.gz" ]; then
    rm clinvar.vcf.gz clinvar.vcf.gz.tbi
fi

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
echo "#Chr	Start	End	Ref	Alt	CLNALLELEID	CLNDN	CLNDISDB	CLNREVSTAT	CLNSIG" > $AV_DB/Hs38DH_clinvar_latest.txt
bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\t%ALLELEID\t%CLNDN\t%CLNDISDB\t%CLNREVSTAT\t%CLNSIG\n' clinvar.vcf.gz >> $AV_DB/Hs38DH_clinvar_latest.txt
perl $SCRIPTPATH/compileAnnovarIndex.pl $AV_DB/Hs38DH_clinvar_latest.txt 1000 > $AV_DB/Hs38DH_clinvar_latest.txt.idx
rm clinvar.vcf.gz clinvar.vcf.gz.tbi