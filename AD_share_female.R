################################################################################
##  Among 1702, 874 share & AD GWAS
################################################################################
.libPaths("/home/jiseon/R/x86_64-redhat-linux-gnu-library/4.3")
library(data.table)
library(dplyr)
library(coloc)



#outd="/home/jiseon/1702coloc/"
joblst<-read.csv("/home/jiseon/significant_sentinel_pqtl_anno_diff.csv") %>% data.frame() %>% filter(sig_group2 == 'share'& type == 'cis') 
analyteslist <- joblst$Analytes

#AD belleng gwas
belleng.data <- fread("/03-DryLab/01-References_Software/01-References/publication_refs/AD/AD_Risk/Bellenguez_etal_2022/GRCh38/GCST90027158_buildGRCh38.tsv.gz")

malefrq<-fread("/03-DryLab/04-Analyses/2022_Sex-specific_YS/2023_sex-pQTL/1_get-pQTL-results/input-files/csf_male.afreq") %>%
  mutate(MAF=ifelse(ALT_FREQS>0.5,1-ALT_FREQS,ALT_FREQS)) %>%
  dplyr::select(ID,MAF)

femalefrq<-fread("/03-DryLab/04-Analyses/2022_Sex-specific_YS/2023_sex-pQTL/1_get-pQTL-results/input-files/csf_female.afreq") %>%
  mutate(MAF=ifelse(ALT_FREQS>0.5,1-ALT_FREQS,ALT_FREQS)) %>%
  dplyr::select(ID,MAF)

joblst <- joblst %>%
  dplyr::mutate(index_pos=as.numeric(sapply(strsplit(ID,":"),function(x) x[2])),
                locus_start= index_pos - 200000,
                locus_end=index_pos + 200000,
                locus_startpos=locus_start,locus_endpos=locus_end)

analyteslist <- joblst$Analytes
myfinres<-c()

#unique(joblst$aptamer) use this. already list no need for loop
for (tr in 1:nrow(joblst)) {
  #for (tr in 91:91){
  file_path <- paste0("/03-DryLab/04-Analyses/2022_Sex-specific_YS/2023_sex-pQTL/1_get-pQTL-results/pQTL-female-results/pqtl.", analyteslist[tr], ".glm.linear.gz")
  if (file.exists(file_path)){
    #for (analyte in analyteslist) {
    infemale <- fread(file_path) %>%
      left_join(femalefrq, by = "ID") %>%
      mutate(
        chr = `#CHROM`,
        snp = paste0("chr", chr, "_", POS),
        position = POS,
        beta = BETA,
        se = SE,
        P = 10^(-LOG10_P),
        lenref = nchar(REF),
        lenalt = nchar(ALT),
        N = OBS_CT,
        alleles = paste0(REF, "_", ALT)
      ) %>%
      dplyr::filter(!is.na(chr) & lenref == 1 & lenalt == 1) %>%
      dplyr::filter(MAF > 0.01 & !alleles %in% c("A_T", "T_A", "C_G", "G_C")) %>%   # palindromic SNP remove
      dplyr::select(chr, snp, position, beta, se, MAF, P, REF, ALT, N)
    
    
    
    #male first
    locus_chr=joblst$chr[tr]
    locus_start = joblst$locus_start[tr]
    locus_end = joblst$locus_end[tr]
    myres.female<-infemale %>% dplyr::filter(chr==locus_chr & position>locus_start & position<locus_end)
    
    bellengtype <- "cc"
    
    belleng.data <- belleng.data %>%
      mutate(p_value = as.numeric(p_value))
    
    
    belleng.female <- belleng.data %>%
      mutate(chr = chromosome, 
             position = base_pair_location, 
             snp = paste0("chr", chr, "_", position), 
             A1 = effect_allele, 
             A2 = other_allele,
             beta = beta, 
             MAF=ifelse(effect_allele_frequency>0.5,1-effect_allele_frequency,effect_allele_frequency),
             P = p_value,  
             se = standard_error, 
             N = n_cases+n_controls) %>%
      select(chr, position, snp, MAF, P, A1, A2, beta, se, N) %>%
      distinct(snp, .keep_all = TRUE)
    
    overlap.snp <- myres.female %>%
      inner_join(belleng.female,by="snp") %>%
      mutate(compA1=case_when(
        A1=="A" ~ "T",
        A1=="T" ~ "A",
        A1=="C" ~ "G",
        A1=="G" ~ "C"
      ),
      compA2=case_when(
        A2=="A" ~ "T",
        A2=="T" ~ "A",
        A2=="C" ~ "G",
        A2=="G" ~ "C"
      ),A1_A2=paste0(A1,"_",A2),
      A2_A1=paste0(A2,"_",A1),
      compA1_compA2=paste0(compA1,"_",compA2),
      compA2_compA1=paste0(compA2,"_",compA1),
      ref_alt=paste0(REF,"_",ALT)
      ) %>%
      dplyr::filter(ref_alt==A1_A2 | ref_alt==A2_A1 | ref_alt==compA1_compA2 | ref_alt==compA2_compA1)
    
    nsnp=dim(overlap.snp)[1]
    p1=1e-4
    p2=1e-4
    p12=1e-5
    
    
    #==========coloc
    myres.female <- myres.female %>%
      dplyr::filter(snp %in% overlap.snp$snp)
    
    myres.female$varbeta <- (myres.female$se)^2
    myres.female<-as.list(myres.female)
    myres.female$N <- max(myres.female$N) # N is sample size information. thsis the sample size in GWAS
    myres.female$type = "quant"
    
    ingwas1<-belleng.female %>%
      dplyr::filter(snp %in% overlap.snp$snp)
    
    ingwas1$varbeta <- (ingwas1$se)^2
    ingwas1<-as.list(ingwas1)
    ingwas1$type=bellengtype
    
    check_dataset(myres.female)
    check_dataset(ingwas1)
    
    
    res1<-coloc.abf(myres.female, ingwas1, MAF = NULL, p1 = p1, p2 = p2, p12 = p12)
    
    
    res1.sum<-as.data.frame(t(as.data.frame(res1$summary)))
    res1.sum$pqtlgene<-analyteslist[tr]
    res1.sum$Chrom=locus_chr
    res1.sum$sex="female"
    res1.sum$startpos=locus_start
    res1.sum$endpos=locus_end
    res1.sum$gwastrait = joblst$belleng_f[tr]
    myfinres<-rbind(myfinres,res1.sum)
  } 
  write.csv(as.data.frame(myfinres), "/home/jiseon/share_AD_female_2nd.csv",quote=TRUE,row.names=FALSE)
  
  
}










