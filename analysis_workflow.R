library(TwoSampleMR)
library(data.table)
library(MendelianRandomization)
library(ieugwasr)
library(plinkbinr)
library(gwasglue)
library(S4Vectors)
library(rio)
library(ggplot2)
library(ggsci)
library(MRPRESSO)
library("ggpubr")

setwd('/MR_head.and.neck_Gut.microbiome')

########################################################################################
mbg <- fread("MBG.allHits.p1e4.txt",head=TRUE,sep='\t')
mbg$samplesize <- "18340"
export(mbg,"MBG.allHits.p1e4.samplesize.txt",format = "\t")

exp_out <- read_exposure_data(filename = "MBG.allHits.p1e4.samplesize.txt",
                          sep='\t',
                          phenotype_col = "phenotype",
                          snp_col = "rsID",
                          beta_col = "beta",
                          se_col = "SE",
                          effect_allele_col = "eff.allele",
                          other_allele_col = "ref.allele",
                          pval_col = "P.weightedSumZ",
                          samplesize_col = "samplesize",
                          chr_col = "chr",
                          pos_col = "bp",
						            id_col="bac")

###############################################################################
exp_dat <- exp_out[exp_out$pval.exposure < 1e-05,]

################################################################################
r2 <- 2*exposure_dat_delLD$eaf.exposure*(1-exposure_dat_delLD$eaf)*(exposure_dat_delLD$beta.exposure/(sqrt(exposure_dat_delLD$samplesize.exposure)*exposure_dat_delLD$se.exposure))^2
#sqrt(exposure_dat_delLD$samplesize.exposure)*exposure_dat_delLD$se.exposure
F <- r2*(exposure_dat_delLD$samplesize.exposure-2)/(1-r2)
exposure_dat_delLD$F <-F

###############################################################################
exposure_dat <- clump_data(exp_dat, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, pop = "EUR")

##############################################################################
exposure_dat<- ld_clump(
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop ="EUR",
  dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, id=exp_dat$id.exposure),
  plink_bin ="/R/x86_64-pc-linux-gnu-library/4.3/plinkbinr/bin/plink_Linux",
  bfile ="/R/x86_64-pc-linux-gnu-library/4.3/plinkbinr/EUR")

export(exposure_dat,"exposure_delete_LD_snp_list.txt",format = "\t")

##########
exposure_dat_delLD<-subset(exp_dat,SNP %in% exposure_dat$rsid)

#########
vcf <- VariantAnnotation::readVcf("ieu-b-4912.vcf", "hg19")
ukb <- gwasvcf_to_TwoSampleMR(vcf,type = "outcome")
export(ukb,"ukb.txt",format = "\t") 

#############harmonise 数据
mydata <- harmonise_data(exposure_dat=exposure_dat,outcome_dat=ukb,action= 2)

export(mydata,"Harmonising_exposure_outcome.UKB.txt",format = "\t") 
mydata <- fread('Harmonising_exposure_outcome.UKB.txt',header = T)

##################################################################################
result <- mr(mydata, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median"))

export(result,"mr_result.txt",format = "\t") 

MR<-generate_odds_ratios(result)

export(MR,"mr_result_get_odds.txt",format = "\t")

################################################################################
mr_results <- read.table("mr_result_get_odds.txt",sep='\t',header=TRUE)

mr_results <- mr_results %>%
  mutate(category = sub("^(phylum|class|order|family|genus).*", "\\1", id.exposure))

unique(mr_results$category)

mr_results <- mr_results %>%
  group_by(category, method) %>%
  mutate(fdr = p.adjust(pval, method = "fdr")) %>%
  ungroup()

mr_results <- mr_results %>%
  relocate(fdr, .after = last_col())

export(mr_results,"mr_result_get_odds_addFDR.txt",format = "\t")

########################################################################################################
het <- mr_heterogeneity(mydata)
export(het,"mr_result_heterogeneity.txt",format = "\t")

res_single <- mr_singlesnp(mydata,parameters = default_parameters(),
                           single_method = "mr_wald_ratio",
                           all_method = c("mr_ivw", "mr_egger_regression"))

export(res_single,"mr_result_singleSNP.txt",format = "\t")

pleio <- mr_pleiotropy_test(mydata) 
export(pleio,"mr_result_pleiotropy.txt",format = "\t")

################################################################################################
single <- mr_leaveoneout(mydata)
export(single,"mr_result_Leave.one.out.txt",format = "\t")

############
set.seed(123)
presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = my.data, NbDistribution = 2000,  SignifThreshold = 0.05)
save(presso,file='mr_gut_UKB_Presso.Rdata')


presso <- run_mr_presso(data,NbDistribution = 2000)

############
mr_steiger_direction <-directionality_test(mydata)
export(mr_steiger_direction,"mr_Gut.vs.UKB_result_steiger_direction.txt",format = "\t")

##########

pdf('Result.pdf',width = 8,height = 5)
significant_snps <- res_single[res_single$p < 0.05, ] 

mr_scatter_plot(mr_results = mr(mydata,method_list = c("mr_egger_regression", "mr_ivw","mr_weighted_median")),mydata)

#########
mr_forest_plot(significant_snps)

#########
mr_funnel_plot(res_single)

#########
top_30 <- single[order(single$p), ][1:30, ]
mr_leaveoneout_plot(top_30)

dev.off()

###################
d <- subset(mydata, mr_keep)
index <- d$beta.exposure < 0
d$beta.exposure[index] <- d$beta.exposure[index] *-1
d$beta.outcome[index] <- d$beta.outcome[index] * -1

mrres <- subset(res)
mrres$a <- 0
d <- subset(d,d$id.exposure==mrres$id.exposure[1])
if ("MR Egger" %in% mrres$method) {
  temp <- mr_egger_regression(d$beta.exposure,
                              d$beta.outcome, d$se.exposure, d$se.outcome,
                              default_parameters())
  mrres$a[mrres$method == "MR Egger"] <- temp$b_i
}
if ("MR Egger (bootstrap)" %in% mrres$method) {
  temp <- mr_egger_regression_bootstrap(d$beta.exposure,
                                        d$beta.outcome, d$se.exposure, d$se.outcome,
                                        default_parameters())
  mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
}

ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure,
                                       y = beta.outcome)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome -
                                                                                                  se.outcome, ymax = beta.outcome + se.outcome),
                                                                                   colour = "grey", width = 0) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure -
                                                                                                                                                        se.exposure, xmax = beta.exposure + se.exposure),
                                                                                                                                         colour = "grey", height = 0) + ggplot2::geom_point() +
  ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a,
                                                  slope = b, colour = method), show.legend = TRUE) +
  ggplot2::scale_colour_manual(values = c("#a6cee3",
                                          "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
                                          "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                                          "#6a3d9a", "#ffff99", "#b15928")) + ggplot2::labs(colour = "MR Test",
                                                                                            x = paste("SNP effect on", d$exposure[1]), y = paste("SNP effect on",
                                                                                                                                                 d$outcome[1])) + ggplot2::theme(legend.position = "top",
                                                                                                                                                                                 legend.direction = "vertical") + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))


