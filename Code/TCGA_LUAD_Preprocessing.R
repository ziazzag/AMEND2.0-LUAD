# TCGA-LUAD preprocessing
# data: http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUAD/20160128/

rm(list=ls())

path.data <- "Data/LUAD/"
path.results <- "Results/"
dir.create(path.results, showWarnings= FALSE)

library(data.table); library(limma); library(survival); library(missForest); library(dplyr); library(Matrix)

extract_string <- function(x, k, pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))

set.seed(351824)



#clinical
clin.dat <- data.table::fread(file= paste0(path.data, "All_CDEs.txt"), header= TRUE) |>
  as.data.frame()
clin.dat <- data.table::transpose(clin.dat)
colnames(clin.dat) <- as.character(clin.dat[1,])
clin.dat <- clin.dat[-1,]
clin.dat$patient_id <- toupper(clin.dat$patient_id)



#methylation
if(1){
  meth.dat <- data.table::fread(file= paste0(path.data, "LUAD.meth.by_mean.data.txt"), header= TRUE) |>
    as.data.frame()
  meth.dat <- na.omit(meth.dat)
  meth.dat <- meth.dat[-1,]
  rownames(meth.dat) <- meth.dat[,1]
  meth.dat <- meth.dat[,-1]

  tt.tmp <- substr(extract_string(colnames(meth.dat), "-", 4), 1, 2)
  meth.dat <- meth.dat[, tt.tmp %in% c("01", "11")]

  pt.ids <- extract_string(colnames(meth.dat), "-", 3)
  tt.tmp <- substr(extract_string(colnames(meth.dat), "-", 4), 1, 2)
  tissue.type <- ifelse(tt.tmp == "01", "tumor", ifelse(tt.tmp == "11", "normal", "other"))

  tmp <- matrix(0, nrow= nrow(meth.dat), ncol= ncol(meth.dat),
                dimnames= list(rownames(meth.dat), colnames(meth.dat)))
  for(i in seq_len(ncol(meth.dat))){
    tmp[,i] <- as.numeric(meth.dat[,i])
  }
  b.values <- tmp
  m.values <- log(b.values / (1 - b.values))

  multi.ids <- pt.ids[pt.ids %in% names(table(pt.ids))[table(pt.ids) > 1]]
  tt.tmp <- tissue.type[pt.ids %in% multi.ids]
  design <- model.matrix(~0 + tt.tmp)
  dat.tmp <- m.values[, pt.ids %in% multi.ids]
  corr.meth <- limma::duplicateCorrelation(dat.tmp, design, block=multi.ids)
  corr.meth <- corr.meth$consensus.correlation
  message(paste("Methylation duplicateCorrelation:", round(corr.meth, 7)))

  design <- model.matrix(~0 + tissue.type)
  fit <- lmFit(object= m.values, design= design, block= pt.ids, correlation= corr.meth)
  contrast_matrix <- makeContrasts("tissue.typetumor-tissue.typenormal", levels= colnames(design))
  fit <- contrasts.fit(fit, contrasts= contrast_matrix)
  fit_eb <- eBayes(fit, trend= TRUE)
  meth.res <- topTable(fit_eb, number= 100000, adjust.method= "BH")
}



#mrna
if(1){
  mrna.dat <- data.table::fread(file= paste0(path.data, "LUAD.uncv2.mRNAseq_RSEM_normalized_log2.txt"),
                                header= TRUE) |>
    as.data.frame()
  mrna.dat <- na.omit(mrna.dat)
  mrna.dat <- mrna.dat[!grepl("\\?", mrna.dat$gene),]

  symbol_entrez.map <- data.frame(
    Symbol= extract_string(mrna.dat$gene, "\\|", 1),
    Entrez= extract_string(mrna.dat$gene, "\\|", 2)
  )
  mrna.dat$gene <- extract_string(mrna.dat$gene, "\\|", 1)

  if(any(duplicated(mrna.dat$gene))){
    expr.cols <- setdiff(colnames(mrna.dat), "gene")
    row.means <- rowMeans(mrna.dat[, expr.cols], na.rm=TRUE)
    mrna.dat <- mrna.dat[order(mrna.dat$gene, -row.means), ]
    mrna.dat <- mrna.dat[!duplicated(mrna.dat$gene), ]
  }

  rownames(mrna.dat) <- mrna.dat$gene
  mrna.dat <- mrna.dat |> dplyr::select(-gene)
  mrna.dat <- mrna.dat[, extract_string(colnames(mrna.dat), "-", 4) %in% c("01", "11")]

  pt.ids <- extract_string(colnames(mrna.dat), "-", 3)
  tissue.type <- ifelse(extract_string(colnames(mrna.dat), "-", 4) == "01", "tumor",
                        ifelse(extract_string(colnames(mrna.dat), "-", 4) == "11", "normal", "other"))

  multi.ids <- pt.ids[pt.ids %in% names(table(pt.ids))[table(pt.ids) > 1]]
  tt.tmp <- tissue.type[pt.ids %in% multi.ids]
  design <- model.matrix(~0 + tt.tmp)
  dat.tmp <- mrna.dat[, pt.ids %in% multi.ids]
  corr <- limma::duplicateCorrelation(dat.tmp, design, block=multi.ids)
  corr$consensus.correlation

  design <- model.matrix(~0 + tissue.type)
  fit <- lmFit(object= mrna.dat, design= design, block= pt.ids, correlation= corr$consensus.correlation)
  contrast_matrix <- makeContrasts("tissue.typetumor-tissue.typenormal", levels= colnames(design))
  fit <- contrasts.fit(fit, contrasts= contrast_matrix)
  fit_eb <- eBayes(fit, trend= TRUE)
  mrna.res <- topTable(fit_eb, number= 100000, adjust.method= "BH")
}



#mirna
if(1){
  mirna.dat <- data.table::fread(file= paste0(path.data, "LUAD.miRseq_mature_RPM_log2.txt"),
                                 header= TRUE) |>
    as.data.frame()
  rownames(mirna.dat) <- mirna.dat[,1]
  mirna.dat <- mirna.dat[,-1]

  tt.tmp <- substr(extract_string(colnames(mirna.dat), "-", 4), 1, 2)
  mirna.dat <- mirna.dat[, tt.tmp %in% c("01", "11")]

  pt.ids <- extract_string(colnames(mirna.dat), "-", 3)
  tt.tmp <- substr(extract_string(colnames(mirna.dat), "-", 4), 1, 2)
  tissue.type <- ifelse(tt.tmp == "01", "tumor", ifelse(tt.tmp == "11", "normal", "other"))

  # remove mirnas with >25% NA
  na.counts <- apply(mirna.dat, 1, function(x) sum(is.na(x)))
  na.perc <- na.counts / ncol(mirna.dat)
  k <- 0.25
  mirna.dat <- mirna.dat[na.perc <= k,]

  # random forest imputation
  mirna.dat.t <- t(mirna.dat)
  mirna.dat.t.imp <- missForest(mirna.dat.t)
  mirna.dat.imp <- t(mirna.dat.t.imp$ximp)
  mirna.dat0 <- mirna.dat.imp

  multi.ids <- pt.ids[pt.ids %in% names(table(pt.ids))[table(pt.ids) > 1]]
  tt.tmp <- tissue.type[pt.ids %in% multi.ids]
  design <- model.matrix(~0 + tt.tmp)
  dat.tmp <- mirna.dat0[, pt.ids %in% multi.ids]
  corr <- limma::duplicateCorrelation(dat.tmp, design, block=multi.ids)
  corr$consensus.correlation

  design <- model.matrix(~0 + tissue.type)
  fit <- lmFit(object= mirna.dat0, design= design, block= pt.ids, correlation= corr$consensus.correlation)
  contrast_matrix <- makeContrasts("tissue.typetumor-tissue.typenormal", levels= colnames(design))
  fit <- contrasts.fit(fit, contrasts= contrast_matrix)
  fit_eb <- eBayes(fit, trend= TRUE)
  mirna.res <- topTable(fit_eb, number= 100000, adjust.method= "BH")
}



#save omics + DE
if(1){
  omics.dat <- list(mrna= mrna.dat, mirna= mirna.dat0, methyl= m.values)
  saveRDS(omics.dat, file= paste0(path.results, "LUAD_multi_omics.rds"))
  de.res <- list(mrna= mrna.res, mirna= mirna.res, methyl= meth.res)
  saveRDS(de.res, file= paste0(path.results, "LUAD_DE_results.rds"))
}




#survival analysis
# combined endpoint: days_to_death for deceased, last followup for censored
surv.dat <- clin.dat[, c("patient_id", "vital_status", "days_to_death", "days_to_last_followup")]
surv.dat$vital_status <- ifelse(surv.dat$vital_status == "alive", 0, 1)
surv.dat$days_to_death <- as.numeric(surv.dat$days_to_death)
surv.dat$days_to_last_followup <- as.numeric(surv.dat$days_to_last_followup)
surv.dat$days_to_event <- ifelse(surv.dat$vital_status == 1,
                                 surv.dat$days_to_death,
                                 surv.dat$days_to_last_followup)
surv.dat <- surv.dat[!is.na(surv.dat$days_to_event) & surv.dat$days_to_event > 0, ]

omics.dat <- lapply(omics.dat, function(x){
  tt.tmp <- substr(extract_string(colnames(x), "-", 4), 1, 2)
  x[, tt.tmp == "01"]
})

surv.dat.list <- vector("list", length(omics.dat))
names(surv.dat.list) <- names(omics.dat)
r.ids <- surv.dat$patient_id
for(i in seq_along(omics.dat)){
  message(i)
  omics.res <- matrix(0, nrow= length(r.ids), ncol= nrow(omics.dat[[i]]),
                      dimnames= list(r.ids, rownames(omics.dat[[i]])))
  pt.ids <- extract_string(colnames(omics.dat[[i]]), "-", 3)
  for(j in seq_along(r.ids)){
    if(r.ids[j] %in% pt.ids){
      omics.res[j,] <- omics.dat[[i]][, pt.ids == r.ids[j]]
    }else{
      omics.res[j,] <- rep(NA, ncol(omics.res))
    }
  }
  omics.res <- omics.res[match(r.ids, rownames(omics.res)),]
  surv.dat.list[[i]] <- cbind(surv.dat, omics.res)
  surv.dat.list[[i]] <- surv.dat.list[[i]][!is.na(omics.res[,1]),]
  omics.res <- as.matrix(surv.dat.list[[i]][, -c(1:ncol(surv.dat))])
  feat.sd <- apply(omics.res, 2, sd, na.rm=TRUE)
  feat.sd[feat.sd == 0 | is.na(feat.sd)] <- 1
  surv.dat.list[[i]][, -c(1:ncol(surv.dat))] <- omics.res %*% Matrix::Diagonal(x=1/feat.sd)
}

feature.names <- vector("list", length(surv.dat.list))
for(i in seq_along(feature.names)){
  res <- matrix("", nrow= ncol(surv.dat.list[[i]]) - ncol(surv.dat), ncol= 2,
                dimnames= list(1:(ncol(surv.dat.list[[i]]) - ncol(surv.dat)), c("Original", "Modified")))
  v <- colnames(surv.dat.list[[i]])[(ncol(surv.dat)+1):ncol(surv.dat.list[[i]])]
  res[,1] <- v
  res[,2] <- make.names(v, unique=TRUE)
  feature.names[[i]] <- res
  colnames(surv.dat.list[[i]])[(ncol(surv.dat)+1):ncol(surv.dat.list[[i]])] <- res[,2]
}




#cox models
cph.models <- vector("list", length(surv.dat.list))
names(cph.models) <- names(surv.dat.list)
for(i in seq_along(surv.dat.list)){
  cph.res <- vector("list", ncol(surv.dat.list[[i]]) - ncol(surv.dat))
  names(cph.res) <- colnames(surv.dat.list[[i]])[(ncol(surv.dat)+1):ncol(surv.dat.list[[i]])]
  message(length(cph.res))
  for(j in seq_along(cph.res)){
    if(j %% 500 == 0)
      message(paste0("i", i, ".j", j))
    v <- colnames(surv.dat.list[[i]])[j + ncol(surv.dat)]
    form <- as.formula(paste("Surv(days_to_event, vital_status)~", v))
    cph.res[[j]] <- tryCatch(
      coxph(form, data= surv.dat.list[[i]]),
      error= function(e) NULL
    )
  }
  cph.models[[i]] <- cph.res
}

for(i in seq_along(surv.dat.list)){
  colnames(surv.dat.list[[i]])[(ncol(surv.dat)+1):ncol(surv.dat.list[[i]])] <- feature.names[[i]][, "Original"]
  names(cph.models[[i]]) <- feature.names[[i]][, "Original"]
}

cph.models <- lapply(cph.models, function(x) Filter(Negate(is.null), x))

hr <- lapply(cph.models, function(x) lapply(x, function(y){
  tmp <- summary(y)$coefficients[, c(2, 5)]
  names(tmp) <- c("hr", "pval")
  tmp
}))

if(1){
  saveRDS(hr, file= paste0(path.results, "LUAD_hazard_ratios.rds"))
  saveRDS(surv.dat.list, file= paste0(path.results, "LUAD_surv_expr_data.rds"))
  saveRDS(de.res, file= paste0(path.results, "LUAD_DE_results.rds"))
  saveRDS(omics.dat, file= paste0(path.results, "LUAD_multi_omics.rds"))
}
