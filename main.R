# read the clinical data
#conda install r-essentials
#MAF/clinical.tsv
- case_id: 17998, project_id: 18004  #there are a few patients with more than one samples
- classification_of_tumor: met: 9040, primary: 6500, Unknow:2464  #more mets than pris
- tissue_of_organ_of_origin: lung(3769), breast, colon, Unknown, ...
- primary_diagnosis: Adenocarcinoma, Carcinoma, IDC, SCC, ...
- site_of_resection_or_biopsy: Lung(2622), Liver, Lymph node, Brain, ...

tmp= read.delim('clinical.tsv')
summary(tmp)
clin = tmp
res = sort(table(clin$tissue_or_organ_of_origin))

#MAF/sample_sheet
tmp = read.delim('gdc_sample_sheet.2018-08-31.tsv', as.is=T)
tmp$nCase = sapply(tmp$Case.ID, function(x) length(strsplit(x, split=',')[[1]]))
tmp$nSample = sapply(tmp$Sample.ID, function(x) length(strsplit(x, split=',')[[1]]))
tmp$type = sapply(tmp$File.Name, function(x) strsplit(x, split='\\.')[[1]][2])
res = cbind(tmp$type, tmp$nCase, tmp$nSample)
res = as.data.frame(res)
sum(tmp$nCase)
sum(tmp$nSample)
dim(tmp)

sampleSheet = tmp
write.csv(sampleSheet, 'sampleSheet.csv', row.names=F)

# write a csv with diseaseType, CaseId
res = list()
for (i in 1:nrow(sampleSheet)) {
    for (x in strsplit(sampleSheet$Case.ID[[i]], split=', ')[[1]]) {
        res[[x]] = sampleSheet$type[[i]]
    }
}
mafInfo = data.frame(CaseID=names(res), TumorType=unlist(res))
write.csv(mafInfo, 'mafInfo.csv')

colnames(clin)[2] = 'CaseID'
res = sqldf('select * from mafInfo join clin using (CaseId)')
write.csv(res, 'sampleInfo.csv')

variant_call_tumortype = sqldf('select * from mafInfo m join variant_call v using (CaseID)')

# prepare a total patient for each type
maf_count = sqldf('select TumorType, count(CaseID) total_patient from mafInfo group by TumorType order by total_patient desc')
write.csv(maf_count, 'mafInfo.count.csv')

# read vcf from database

# check the Effect discrepancies
```R
effects = table(tmp$FMI_effect, tmp$Effect)
res = as.matrix(effects)
write.csv(t(res), 'snpEff_FM.csv')

tmp = read.csv('variant_call.csv')
table(tmp$CSQ_Consequence, tmp$FMI_effect)
res = table(tmp$CSQ_Consequence, tmp$FMI_effect)
write.csv(res, 'Vep_FMI.csv')
```
# calculate global AI
```sqldf
#source('http://bioconductor.org/biocLite.R')
#biocLite('sqldf')
require(sqldf)
require(dplyr)
colnames(tmp)[9] = 'CSQ_impact'
variant_call = tmp
total_patient = 18004 #17798 with any variant called

#gene level: the prob of seeing this gene mutated in the patients
#using variant_call, total_patient
calc_on_gene_level = function() {
    tmp = sqldf("select distinct FMI_gene, CaseID from variant_call where FMI_effect not in ('', 'nonsynonymous')")
    res = sqldf("select FMI_gene, count(distinct CaseID) count_patient from tmp group by FMI_gene order by count_patient desc")
    res$freq = round(res$count_patient/total_patient * 100, 2)
    res
}
gene_frequencies = calc_on_gene_level()
write.csv(gene_frequencies, 'gene_frequencies.csv')
    #todo: put filtered variant type into config

#indication level: the prob of see this variant in the patients
#using variant_call, total_patient
calc_on_indication_level = function() {
    tmp = sqldf("select distinct FMI_gene, HGVS_p, CaseID from variant_call where FMI_effect not in ('', 'nonsynonymous')")
    res = sqldf("select FMI_gene, HGVS_p, count(distinct CaseID) count_patient from tmp group by FMI_gene, HGVS_p order by count_patient desc")
    #res$freq = round(res$count_patient/total_patient * 100, 2)
    res %>% mutate(freq=round(count_patient/total_patient * 100, 2))
}
indication_frequencies = calc_on_indication_level()
write.csv(indication_frequencies, 'indication_frequencies.csv')
    #todo: against using HGVS_p here: 1) many to many with variant, 2) variant without influencing protein varies
    #todo: try a version using dplyr

#variant in gene: the prob of seeing this variant if the patient has the gene mutated.
#using indication_frequencies, gene_frequencies
calc_variant_in_gene = function() {
    tmp = sqldf("select FMI_gene, HGVS_p, i.count_patient i_count_patient, g.count_patient g_count_patient 
        from indication_frequencies i join gene_frequencies  g using (FMI_gene)
        order by g_count_patient desc")
    res = tmp %>% mutate(freq = round(i_count_patient / g_count_patient * 100, 2))
}

variant_with_gene = calc_variant_in_gene()
write.csv(variant_with_gene, 'variant_with_gene.csv')


## calc AI for each disease type
#using variant_call, mafInfo, maf_count
tmp = sqldf('select * from maf_count join mafInfo using (TumorType)')
masterView = sqldf("select distinct t.*, FMI_gene, HGVS_p 
    from tmp t join variant_call using(CaseID)
    where FMI_effect not in ('', 'nonsynonymous')")
dim(masterView)    #144074

#indication level: the prob of see this variant in the patients
res = sqldf("select TumorType, total_patient, FMI_gene, HGVS_p,
        count(distinct CaseID) count_patient
    from masterView
    group by TumorType, total_patient, FMI_gene, HGVS_p
    order by total_patient desc, count_patient desc")
dim(res)    #115338
i_res = sqldf('select *, cast(count_patient as float)/total_patient*100 freq
    from res')
head(i_res)
write.csv(i_res, 'indication_frequency.byTumorType.csv')

#gene level: the prob of seeing the gene mutated in the patients
res = sqldf("select TumorType, total_patient, FMI_gene,
        count(distinct CaseID) count_patient
    from masterView 
    group by TumorType, total_patient, FMI_gene
    order by total_patient desc, count_patient desc")
dim(res) #8002
g_res = sqldf('select *, cast(count_patient as float)/total_patient*100 freq
    from res')
head(g_res)
write.csv(g_res, 'gene_frequency.byTumorType.csv')

indication_frequencies = i_res
gene_frequencies = g_res
#variant in gene: prob of seeing the variant if the gene is mutated.
#using indication_frequencies, gene_frequencies
tmp = sqldf("select TumorType, FMI_gene, HGVS_p,
            cast(i.count_patient as float)/g.count_patient*100 freq
        from indication_frequencies i join gene_frequencies g using (TumorType, FMI_gene)
        order by g.total_patient desc, g.count_patient desc")
dim(tmp) #115338
head(tmp)
variant_with_gene = tmp
write.csv(variant_with_gene, 'variant_with_gene.byTumorType.csv')


## recalc AI for all disease type
#using variant_call, mafInfo, maf_count
#cohort_count = data.frame(cohort_name='FM_AD', total_patient=sum(maf_count$total_patient))
#cohort_count
masterView = sqldf("select distinct FMI_gene, HGVS_p, CaseID, 18004 as total_patient 
    from variant_call
    where FMI_effect not in ('', 'nonsynonymous')")
dim(masterView)    #144074
head(masterView

#indication level: the prob of see this variant in the patients
res = sqldf("select total_patient, FMI_gene, HGVS_p, 
        count(distinct CaseID) count_patient
    from masterView
    group by FMI_gene, HGVS_p
    order by count_patient desc")
dim(res)    #92889
head(res)
i_res = sqldf('select *, 100.*count_patient/total_patient as freq
    from res')
head(i_res)
write.csv(i_res, 'indication_frequency.cohort.csv')

#gene level: the prob of seeing the gene mutated in the patients
res = sqldf("select total_patient, FMI_gene,
        count(distinct CaseID) count_patient
    from masterView 
    group by total_patient, FMI_gene
    order by total_patient desc, count_patient desc")
dim(res) #8002
g_res = sqldf('select *, cast(count_patient as float)/total_patient*100 freq
    from res')
head(g_res)
write.csv(g_res, 'gene_frequency.cohort.csv')

indication_frequencies = i_res
gene_frequencies = g_res
#variant in gene: prob of seeing the variant if the gene is mutated.
#using indication_frequencies, gene_frequencies
tmp = sqldf("select FMI_gene, HGVS_p,
            cast(i.count_patient as float)/g.count_patient*100 freq
        from indication_frequencies i join gene_frequencies g using (FMI_gene)
        order by g.total_patient desc, g.count_patient desc")
dim(tmp) #92889
head(tmp)
variant_with_gene = tmp
write.csv(variant_with_gene, 'variant_with_gene.cohort.csv')




