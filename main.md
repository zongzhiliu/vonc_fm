# extract the variant info from FoundationMedicine data from dbGAP, calculate the allele frequencys.
* indication: #v/#total samples
* variantAF: #v/#total variant in the gene
* geneAF: #anyMutInTheGene/#total samples

## check out
```sh
cd ~/work/Foundation_Medicine/Test_Genotype_Files
find -name '*.gz' | xargs gunzip
:w
```

## liftover from hg38 to hg19
```
# install requirements
conda install -c bioconda ucsc-liftover ucsc-twobittofa picard
cd ~/data
# download chain file
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
    #liftOver oldFile map.chain newFile unMapped
# download genome file
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
twoBitToFa hg19.2bit hg19.fa

# lift over with picard
picard -Xmx10g CreateSequenceDictionary R=hg19.fa O=hg19.dict
    #! liftOver failed because of pos+1 and mismatching reference, should be a bug of picard?

# lift over with liftover: bed [0,1), vcf,igv [1,1], 
cat test1.vcf | sed '/^#/d' | awk '{ print $1, $2-1, $2, $8}' > test1.tmp
liftOver test1.tmp ~/data/hg38ToHg19.over.chain.gz test1.o test1.e
```
Turn out that the FM liftover to hg38 make a position+1 shift, so that a liftover from the downloaded vcf to hg19 give errors of mismatching reference. Decided to simply replace the Chrom and POS with the originalContig and INFO.originalStart.

## test run snpEff
```
conda install -c bioconda snpeff
snpEff databases | grep -i human
snpEff download -v GRCh38.p7.RefSeq
snpEff -Xmx4g GRCh38.p7.RefSeq test1.vcf > test1.hg38.eff.vcf
    # note that the eff annotation is inconsistent with FM
snpEff -Xmx4g -v GRCh37.p13.RefSeq test1.hg19.vcf > test1.hg19.ann.vcf
```
SnpEff on 18K file will be very time consuming, thinking of generate a combined vcf file.

## fix and merge the vcf files
- Fix: by replace CHOME with INFO.orginalChrom, and POS with INFO.originalPos, then replace the head with hg19
- Merge: merge the 18 vcf files by add a column sampleId with the Values defined by FORMAT

## run snpEff
```
#test run
vepVcf='/scratch/Foundation_Medicine/VCFS/combined_variants.vep.vcf' 
head $vepVcf -n 15 > test2.vep.vcf
snpEff -Xmx4g -v GRCh37.p13.RefSeq test2.vep.vcf > test2.vep.eff.vcf
#manual check INFO.ANN: it is consistence with vep in all cases for genename, c. p., effect, etc
#run

snpEff -Xmx10g -v GRCh37.p13.RefSeq $vepVcf > vep.eff.vcf
```
## check the snpEff results
* snpEff_summary.html
    * warnings and errors
    * number of processed: yes
    * variant rate details: NM?  >> check the snpEff ref
    * number of variants by type
    * number of effects by impact: most moderate, then High, modifier (as expected, as they are prefiltered by FM)
    * Number of effects by functional class: Misssense, nonsense, silent?  >> check for inconsistence
    * Number of effects by type and region: Exon, upstream, intron, downstream, splicing, UTR, ... >> check the up and downstream
    * ins/del lengths
    * base changes for SNPs: G2A/C2T is highly enriched?  >> research
    * Ts/Tv ratio = 1.4
    * Amino acid changes: E2K, R2H enriched
* snpEff_genes.txt
    todo: summarise by gene name
* vep.eff.vcf

## error and warning messages
```
ERRORS: Some errors were detected
Error type      Number of errors
ERROR_CHROMOSOME_NOT_FOUND      16


WARNINGS: Some warning were detected
Warning type    Number of warnings
INFO_REALIGN_3_PRIME    6238
WARNING_TRANSCRIPT_INCOMPLETE   2197
WARNING_TRANSCRIPT_INCOMPLETE&INFO_REALIGN_3_PRIME      2
WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS 8950
WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS&INFO_REALIGN_3_PRIME    30
WARNING_TRANSCRIPT_NO_START_CODON       3052
WARNING_TRANSCRIPT_NO_STOP_CODON        186
```

* report the vcf headers: anno of the fields and INFO subfields

## normalize the table
:wq

