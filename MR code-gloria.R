#install.packages("usethis")
usethis::edit_r_environ()
Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJuaW1idXMyMDAwZUAxNjMuY29tIiwiaWF0IjoxNzM1MjAxNjY1LCJleHAiOjE3MzY0MTEyNjV9.YQgNR3VCRA1tbBN4KVL3N-_LvqQgtdzyJOS6spdzj4qtpjfDMQ36yHKHU36GnxmRClU9Gg-ZATBpGX5ejozP1CSkeYh5TVqffY_-OGFgfZZJ_M9b-ZxOsy4JTejLGFjJG2-jhVbOJHfYLi02-w4Fu_viSneswsV4HzU9-tiCfMI8M8OJ-7Vy7wO1bgVlEa8vfFsMKU8O-y5a6lVvPZJIguBjmsl4onNTR7o63N16uVE-P8vVUu6fGzISKB-4PF5M95NMguPE_CaIXSyd4swUaqdI9it7UpW1pC6s7P8ZWSsBpB1qRflnDHPoR_C5iiC8xOnUE9UQFsldtaJDEp1KWg")
#install.packages("ieugwasr")
library(ieugwasr)
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")
get_opengwas_jwt()
user()
#rm(list=ls())
wrdir <- ("/Users/liutiange/Documents/文章修改/孟德尔随机化/近视与眼病")
setwd(wrdir)
#ukb-b-6353
install.packages("openxlsx")
library(openxlsx)
#install.packages("TwoSampleMR")
library(TwoSampleMR)
#install.packages ("MRInstruments")
library(MRInstruments)
library(ggplot2)
###两个样本来自相同人群
##读取数据（获取工具变量及工具变量筛选）
##默认P<5E-8,r2=0.001，kb=10000
myopia<-extract_instruments(outcomes='ukb-b-6353', clump=T)

#save(myopia,file="myopia_exposure.Rdata")
load(file="myopia_exposure.Rdata")
#去掉不含EAF的
#library(tidyverse)
#data <- drop_na(data,eaf.exposure)

###计算MAF
eaf2maf <- function(eaf = NULL) {
  if (any(is.infinite(eaf))) {
    warning("The 'eaf' vector contains infinite values. 
            These will be turned into NAs.")
    is.na(eaf) <- is.infinite(eaf)}
  if (is.null(eaf)) {stop("No 'eaf' vector provided.")}
  if (is.character(eaf)) {stop("'eaf' must be a numeric vector.")}
  if (all(is.na(eaf))) {stop("All values in 'eaf' are NA.")}
  if (is.numeric(eaf)) {
    maf <- eaf
    ind <- which(eaf > 0.5)
    maf[ind] <- 1 - eaf[ind]
    return(maf)
  }
}
# 先加载上面的Function，在进行转化
exp_dat=TwoSampleMR::extract_instruments(outcomes = "ukb-b-6353")
exp_dat$maf <- eaf2maf(eaf=exp_dat$eaf.exposure)
exp_dat$maf

#write.xlsx(exp_dat,"exp_dat.xlsx",rowNames=F,colNames=T)

####计算R2和F统计量
get_f<-function(dat,F_value=10){
  log<-is.na(dat$eaf.exposure)
  log<-unique(log)
  if(length(log)==1)
  {if(log==TRUE){
    print("数据不包含eaf，无法计算F统计量")
    return(dat)}
  }
  if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("数据不包含beta，无法计算F统计量")
    return(dat)}
  if(is.null(dat$se.exposure[1])==T || is.na(dat$se.exposure[1])==T){print("数据不包含se，无法计算F统计量")
    return(dat)}
  if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("数据不包含samplesize(样本量)，无法计算F统计量")
    return(dat)}
  
  
  if("FALSE"%in%log && is.null(dat$beta.exposure[1])==F && is.na(dat$beta.exposure[1])==F && is.null(dat$se.exposure[1])==F && is.na(dat$se.exposure[1])==F && is.null(dat$samplesize.exposure[1])==F && is.na(dat$samplesize.exposure[1])==F){
    R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
    F<- (dat$samplesize.exposure-2)*R2/(1-R2)
    dat$R2<-R2
    dat$F<-F
    dat<-subset(dat,F>F_value)
    return(dat)
  }
}

exposue <-  get_f(exp_dat, F_value=0)
#save(exposue,file="exposue.Rdata")
load(file="exposue.Rdata")
##筛选
###F要大于10
#write.xlsx(exposue,"expose.xlsx",rowNames=F,colNames=T)

##去除混杂
library(MendelianRandomization)

confounder <- MendelianRandomization::phenoscanner(
  snpquery = unique(exposue$SNP),
  catalogue = "GWAS",
  pvalue = 5e-8,
  proxies = "None",
  r2 = 0.8,
  build = 37)

conf <- as.data.frame(cbind(confounder[["results"]][["snp"]],confounder[["results"]][["trait"]]))
conf$V2 <- tolower(conf$V2)##小写转化
conf_count <- dplyr::count(conf,V2,sort =TRUE)

c <-conf[conf$V2%in%c("basal metabolic rate","age-related macular degeneration","age of onset of myopia"),]
s <- unique(c$V1)

con2 <- conf[!conf$V1%in%s,]
conf_count2 <- dplyr::count(con2,V2,sort =TRUE)
exposue_conf <- myopia[myopia$SNP%in%con2$V1,]
#save(exposue_conf,file="exposue_conf.Rdata")
load(file="exposue_conf.Rdata")
#write.xlsx(exposue_conf,"expose2.xlsx",rowNames=F,colNames=T)

##进行筛选(F9)(89个)
outcome2 <- extract_outcome_data(
  snps = exposue_conf$SNP,
  ##结局----
  outcomes =c("finn-b-H7_CORNEALDYSTROPHY",
              "finn-b-H7_CORNEALNAS",
              "finn-b-H7_CORNEALOTH",
              "finn-b-H7_CORNEALSCAR",
              "finn-b-H7_CORNULCER",
              "finn-b-H7_CONJUHAEMOR",
              "finn-b-H7_CONJUNAS",
              "finn-b-H7_CONJUNCTDEGDEPOT",
              "finn-b-H7_CONJUNCTIVA",
              "finn-b-H7_GLAUCCLOSEPRIM",
              "finn-b-H7_GLAUCNAS",
              "finn-b-H7_GLAUCOMA",
              "finn-b-H7_GLAUCOMA_NTG",
              "finn-b-H7_GLAUCOMA_POAG",
              "finn-b-H7_GLAUCOMA_XFG",
              "finn-b-H7_GLAUCPRIMOPEN",
              "finn-b-H7_GLAUCSECDISEASE",
              "finn-b-H7_GLAUCSECINFLAM",
              "finn-b-H7_GLAUCSECTRAUMA",
              "finn-b-H7_GLAUCSUSP",
              "ebi-a-GCST009722",
              "finn-b-DRUGADVERS_CATAR",
              "finn-b-H7_CATARACTOTHER",
              "finn-b-H7_CATARACTSENILE",
              "finn-b-H7_DISORDERLENSNAS",
              "finn-b-H7_LENS",
              "finn-b-H7_LENSDISLOCATIO",
              "finn-b-H7_GLOBE",
              "finn-b-H7_GLOBEDEGENERATED",
              "finn-b-H7_VITRBODYGLOBE",
              "finn-b-H7_VITRCRYSTAL",
              "finn-b-H7_VITREOUS",
              "finn-b-H7_VITRHAEMORR",
              "finn-b-H7_VITROPACIT",
              "finn-b-H7_VITROTH",
              "finn-b-H7_CENTRRETARTOCC",
              "finn-b-H7_CHALAZION",
              "finn-b-H7_CHORIORETINFLAM",
              "finn-b-H7_CHOROIDOTH",
              "finn-b-H7_CHOROIDRETINA",
              "finn-b-H7_RETINAHAEMORR",
              "finn-b-H7_RETINALBREAK",
              "finn-b-H7_RETINALDETACH",
              "finn-b-H7_RETINALDETACHBREAK",
              "finn-b-H7_RETINALDETACHOTH",
              "finn-b-H7_RETINALDETACHTRACTION",
              "finn-b-H7_RETINALDISOTH",
              "finn-b-H7_RETINANAS",
              "finn-b-H7_RETINASEPAR",
              "finn-b-H7_RETINAVASC",
              "finn-b-H7_RETINOCHISISCYST",
              "finn-b-H7_RETINOPATHYDIAB",
              "finn-b-H7_RETINOPATHYDIAB_BKG",
              "finn-b-H7_RETINOPATHYDIAB_BKG_SEVERE",
              "finn-b-H7_RETINOPATHYDIAB_NAS",
              "finn-b-H7_RETINOPATHYDIAB_PROLIF",
              "finn-b-H7_RETIVASCOCCLUSION",
              "finn-b-H7_RETVASCNAS",
              "finn-b-H7_BCKRNDRETINOPAT",
              "finn-b-H7_HEREDRETINADYST",
              "finn-b-H7_OPTATROPHY",
              "finn-b-H7_OPTICDISCOTH",
              "finn-b-H7_OPTNERVE",
              "finn-b-H7_OPTNERVNAS",
              "finn-b-H7_OPTNEURITIS",
              "finn-b-H7_OTHRETARTOCC",
              "finn-b-H7_PAPILLOEDEMA",
              "finn-b-MACULA_NAS",
              "finn-b-MACULAR_CYST",
              "finn-b-MACULAR_HOLE",
              "finn-b-MACULAR_PUCKER",
              "finn-b-WET_AMD",
              "finn-b-H7_MACULADEGEN",
              "finn-b-H7_MACULOPATHYDIAB",
              "finn-b-DRY_AMD",
              "finn-b-H7_AMD",
              "finn-b-H7_BLINDMONOCULAR",
              "finn-b-DM_BCKGRND_RETINA",
              "finn-b-DM_BCKGRND_RETINA_NONPROLIF",
              "finn-b-DM_MACULOPATHY",
              "finn-b-DM_NEOVASCULAR_GLAUCOMA",
              "finn-b-DM_OPTHALMIC_COMORB",
              "finn-b-DM_RETINA_NOS",
              "finn-b-DM_RETINA_PROLIF",
              "finn-b-DM_RETINOPATHY",
              "finn-b-DM_VITREOUS_BLEEDING",
              "finn-b-H7_EXOPTHALMUS",
              "finn-b-H7_PTERYGIUM",
              "finn-b-H7_OCUPAIN"),
  proxies = FALSE,
  maf_threshold = 0.01,
  access_token = NULL)
####结局----


mydata2 <- harmonise_data(
  exposure_dat=exposue_conf,
  outcome_dat=outcome2,
  action= 2
)
res <- mr(mydata2)
##共89个结局
#save(res,file="res.Rdata")
load(file="res.Rdata")
###附件2
##write.xlsx(res,"89.xlsx",rowNames=F,colNames=T)

OR_89<-generate_odds_ratios(res)
res_scr_89 <- res[res$method=="Inverse variance weighted",]
#OR值计算
OR_scr<-generate_odds_ratios(res_scr_89)

##根据P值筛选
p_scr <- OR_scr[OR_scr$pval<0.05,]

####### FDR校正
res_scr_89$q <- p.adjust((res_scr_89$pval),  # P值列表
  method ="BH")     # FDR校正的方法

table <- merge(res_scr_89[,c(2,10)],res,by="id.outcome")
for(i in 1:nrow(table)){
  if(table$method[i]!="Inverse variance weighted")
  table$q[i]<- NA
}
##write.xlsx(table,"S3.xlsx",rowNames=F,colNames=T)

q_scr <- res_scr_89[res_scr_89$q<0.05,]


method <- data.frame()
for(i in 1:nrow(q_scr)){
  a<- OR_89[OR_89$id.outcome==q_scr$id.outcome[i],]
  method <- rbind(method,a)
  }
method <- method[-c(6:10),]

##write.xlsx(method,"筛选结果qvalue.xlsx",rowNames=F,colNames=T)
##write.xlsx(res,"孟德尔89结果.xlsx",rowNames=F,colNames=T)

########我们感兴趣的眼病###############
outcome_intere<-extract_outcome_data(
  snps = exposue_conf$SNP,
  outcomes =c("finn-b-H7_VITROTH","finn-b-H7_OPTICDISCOTH","ebi-a-GCST009722","finn-b-H7_CATARACTOTHER","finn-b-H7_GLAUCCLOSEPRIM","finn-b-H7_RETINALDETACHBREAK"),
  proxies = FALSE,
  maf_threshold = 0.01)

##"finn-b_H7 VITROTH","finn-b-H7_OPTICDISCOTH","ebi-a-GCST009722","finn-b-H7_CATARACTOTHER","finn-b-H7_GLAUCCLOSEPRIM","finn-b-H7_RETINALDETACHBREAK","finn-b-H7_VITRBODYGLOBE"
mydata_intere <- harmonise_data(
  exposure_dat=exposue_conf,
  outcome_dat=outcome_intere,
  action= 2
)

res <- mr(mydata_intere)
OR_intere<-generate_odds_ratios(res)
#write.xlsx(OR_intere,"6_table.1.xlsx",rowNames=F,colNames=T)
#read.xlsx("6_table.1.xlsx")
#######

##############
##library(RadialMR)
moren <- mr(mydata_intere)
ivw <- mr(mydata_intere,method_list=c('mr_ivw'))
suiji <- mr(mydata_intere,method_list=c('mr_ivw_mre'))
OR_suiji<-generate_odds_ratios(suiji)
gud <- mr(mydata_intere,method_list="mr_ivw_fe")
racial <- mr(mydata_intere,method_list="mr_ivw_radial")
plot_radial(racial)

library(MRPRESSO)
presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =mydata_intere, NbDistribution = 1000, SignifThreshold = 0.05) 

##异质性检验(主要是检验各个工具变量之间的差异,也可以用MR-presso)
#这些IV之间存在很强的异质性（Q_pval远小于0.05),这时候我们需要剔除某些outcome的SNP，或者直接使用随机效应模型
#如果不存在异质性，那就是固定效应模型
###使用随机效应模型:mr(mydata,method_list=c('mr_ivw_mre'))
het <- mr_heterogeneity(mydata_intere)
het
presso <- run_mr_presso(mydata_intere, NbDistribution = 5000)

#mr(mydata,method_list=c('mr_ivw_mre'))
###计算异质性的I2
library(MendelianRandomization)
MRInputObject <- MendelianRandomization::mr_input(
  bx = mydata_intere$beta.exposure,
  bxse = mydata_intere$se.exposure,
  by = mydata_intere$beta.outcome,
  byse = mydata_intere$se.outcome,
  snps = mydata_intere$SNP)
  #
  MendelianRandomization::mr_ivw(
  object = MRInputObject, model = "fixed" )


#多效性检验(主要检验多个工具变量是否存在水平多效性，可以理解为是否存在混杂因素)
#egger_intercept与0进行对比,p>0.05,说明不存在水平多效性。
pleio <- mr_pleiotropy_test(mydata_intere)
pleio
write.xlsx(het,"het_intere.xlsx",rowNames=F,colNames=T)
write.xlsx(pleio,"pleio_intere.xlsx",rowNames=F,colNames=T)
#逐个剔除检验-（leave-one-out）#主要是逐个剔除IV后计算剩下IV的MR结果
#如果无论去除哪个SNP都不会对结果产生根本影响（所有的线条均在0的一侧），那就是说这个MR结果实际是稳健的。

pdata <- mydata_intere
pdata$exposure <- "myopia"
pdata[1:14,]$originalname.outcome <- "disorders of the optic disc"
pdata$outcome <- pdata$originalname.outcome
single <- mr_leaveoneout(pdata)
 mr_leaveoneout_plot(single)
p[[1]]

#visualization：
#
p.5 <- mr_scatter_plot(res,pdata)
p.5
#
res_single <- mr_singlesnp(mydata_intere)
p <-mr_funnel_plot(res_single)
p1 <- p[[1]] + labs(title="Funnel plot for myopia on POAG")+
                          theme(plot.title=element_text(#
                                color="black", #
                                size=10,  #
                                hjust=0.5, #
                                vjust=0.5,
                                angle=360)) # 
p1 


#forest plot
res_single$exposure <- "myopia"
res_single$outcome[1:12] <- "Primary open angle glucoma"
res_single$outcome[13:27] <- "Other cataract"
res_single$outcome[28:42] <- "Primary angle−closure glaucoma"
res_single$outcome[43:57] <- "Disorders of optic disc"
res_single$outcome[58:72] <- "Retinal detachment with retinal break"
res_single$outcome[73:87] <- "Disorders of vitreous body"
a <- mr_forest_plot(res_single)
a
