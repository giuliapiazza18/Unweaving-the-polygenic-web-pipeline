# Script to calculate the polygenic score for ADHD (Demontis et al 2019) with LDPred2; based on https://privefl.github.io/bigsnpr/articles/LDpred2.html

library(data.table)
library(magrittr)
library(bigsnpr)
library(runonce)
library(dplyr)

# load GWAS data
sumstats <- bigreadr::fread2("ADHD_Demontis_2019.txt")
head(sumstats)
sumstats$logOR <- log(sumstats$OR)
sumstats$n_eff <- 4 / (1 / 20183 + 1 / 35191)
sumstats <- sumstats %>% select(CHR,SNP,BP,A1,A2,logOR,SE,P,n_eff)
names(sumstats) <- c("chr", "rsid", "pos", "a1", "a0", "beta", "beta_se", "p", "n_eff")


# filter out hapmap3 snps as suggested
info <- readRDS(runonce::download_file("https://ndownloader.figshare.com/files/25503788",fname = "map_hm3_ldpred2.rds"))
sumstats <- sumstats[sumstats$rsid%in% info$rsid,]


# load bed file, creaete .rds and attach to session
obj.bigSNP <- snp_attach("1000G_all_common_afterQC_children.rds")

# get useful aliases
G2   <- obj.bigSNP$genotypes
G <- snp_fastImputeSimple(G2)


CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
NCORES <- 10
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))

# match genotype to sumstats data (chr, pos, a1)
df_beta <- snp_match(sumstats, map, strand_flip = FALSE)
df_beta <- as_tibble(df_beta)

# compute LD reference from data with 3 cM window [genome-wide, on QCd data] ----
# - get genetic position in cM using function
POS2 <- snp_asGeneticPos(CHR, POS)

# - compute correlation between variants per chr
# - compute LD scores (not just corr matrix, need to compute scores)
# - convert matrix to be stored on disk to parallelize (more efficient; request like 60GB) as_SFBM

# open a temporary file
tmp <- tempfile(tmpdir = "tmp-data")

for (chr in 1:22) {
  # Extract SNPs that are included in the chromosome
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'G'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  # Calculate the LD
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                   infos.pos = POS2[ind.chr2], ncores = NCORES)
  
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}


# perform LD score regression ----
ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                               sample_size = n_eff, blocks = NULL))
h2_est <- ldsc[["h2"]]

# LDpred2-auto ----
# - h2 and p estimated from data (around 30 models)
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
                               allow_jump_sign = FALSE, shrink_corr = 0.95,
                               ncores = NCORES)


# # - keep prediction of chains that converged + average
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

pred_auto <- big_prodMat(G, beta_auto, ind.col = df_beta[["_NUM_ID_"]])

sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto <- rowMeans(beta_auto[, keep])
final_pred_auto <- rowMeans(pred_auto[, keep])

fam <- obj.bigSNP$fam$family.ID
fam_pred <- cbind(fam, final_pred_auto)


write.csv(fam_pred, file = "adhd_pgs.csv", row.names = F, quote = F)
