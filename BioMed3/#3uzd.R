#!/usr/bin/env Rscript

library(methods)

#Kaikuriuose palikau test if statementus kuriuos naudojau kai niekas nejo kad zinociau kurie failai man nesiskaito/kuriuose crashina

# Kad runninciau as laikiau dvorkin ir kita faila is skaidriu toje pacioje savo direktorijoje
rds_path <- 'dvorkin.rds'
obj <- readRDS(rds_path)
beta <- obj
annot <- attr(obj, '.annmatrix.rann')
samplekey <- attr(obj, '.annmatrix.cann')

# Amziu ir Cpg data veliasniems skaiciavimams
keep <- which(!is.na(samplekey$age))
beta_c <- beta[, keep, drop=FALSE]
samplekey_c <- as.data.frame(samplekey[keep, , drop=FALSE])
colnames(beta_c) <- rownames(samplekey_c)

sample_n <- 20000
if (!is.null(sample_n) && sample_n < nrow(beta_c)) {
  set.seed(12345)
  idx <- sort(sample(seq_len(nrow(beta_c)), sample_n))
  beta_sub <- beta_c[idx, , drop=FALSE]
} else beta_sub <- beta_c

# Limma testas
if (requireNamespace('limma', quietly=TRUE)) {
  library(limma)
  design <- model.matrix(~ samplekey_c$age, data=samplekey_c)
  fit <- lmFit(beta_sub, design); fit <- eBayes(fit)
  res <- topTable(fit, coef=2, number=Inf, sort.by='P')
  res$cpg <- rownames(res)
  res$P.Value <- res$P.Value
} else {
  res <- data.frame(cpg = rownames(beta_sub), P.Value = NA_real_, logFC = NA_real_, stringsAsFactors=FALSE)
  for (i in seq_len(nrow(beta_sub))) {
    y <- as.numeric(beta_sub[i, ])
    m <- lm(y ~ samplekey_c$age)
    res$logFC[i] <- coef(m)[2]
    res$P.Value[i] <- summary(m)$coefficients[2,4]
  }
}

# Cia sudaromi cpg resultatai kad veliau glaima skaiciavimams naudoti
if (exists('res') && is.data.frame(res)) {
  out_res <- res
  if (!('cpg' %in% colnames(out_res))) out_res$cpg <- rownames(out_res)
  if (!('logFC' %in% colnames(out_res))) out_res$logFC <- NA_real_
  if (!('P.Value' %in% colnames(out_res))) out_res$P.Value <- NA_real_
  write.csv(out_res[, c('cpg','logFC','P.Value')], 'amzius_lm.csv', row.names=FALSE)
  
  # 10 Cpgs svarbiausiu
  ord <- order(out_res$P.Value, na.last=TRUE)
  top_idx <- head(ord, 10)
  for (i in seq_along(top_idx)) {
    ri <- top_idx[i]
    if (is.na(ri) || ri < 1 || ri > nrow(out_res)) next
    cpgid <- out_res$cpg[ri]
    y <- NULL
    if (!is.null(beta_sub) && cpgid %in% rownames(beta_sub)) y <- as.numeric(beta_sub[cpgid, ])
    else if (!is.null(beta_c) && cpgid %in% rownames(beta_c)) y <- as.numeric(beta_c[cpgid, ])
    if (is.null(y)) next
    png(sprintf('top10_cpg_%02d.png', i), width=700, height=500)
    plot(samplekey_c$age, y, pch=16, xlab='Age', ylab='Methylation (beta)', main=cpgid)
    abline(stats::lm(y ~ samplekey_c$age), col='red')
    dev.off()
  }
  
  # P histograna tiesiog
  png('pvalue_histograma.png', width=900, height=600)
  hist(out_res$P.Value, breaks=seq(0,1,by=0.01), xlim=c(0,1), main='P-value histogram', xlab='P-value')
  dev.off()
  
  # Suannotuoju tuos CpGs kurie turi "significant value" (p <= 0.05)
  sig_thr <- 0.05
  sig_idx <- which(!is.na(out_res$P.Value) & out_res$P.Value <= sig_thr)
  sig_df <- out_res[sig_idx, , drop=FALSE]
  sig_annot <- data.frame(cpg = sig_df$cpg, logFC = sig_df$logFC, P.Value = sig_df$P.Value, stringsAsFactors=FALSE)
  if (!is.null(annot) && (is.data.frame(annot) || is.matrix(annot))) {
    ann_df <- as.data.frame(annot)
    if (is.null(rownames(ann_df))) rownames(ann_df) <- seq_len(nrow(ann_df))
    cn <- tolower(colnames(ann_df))
    chr_col <- if (length(grep('chr|chrom', cn))) colnames(ann_df)[grep('chr|chrom', cn)[1]] else NULL
    pos_col <- if (length(grep('pos|map|coordinate|mapinfo', cn))) colnames(ann_df)[grep('pos|map|coordinate|mapinfo', cn)[1]] else NULL
    gene_col <- if (length(grep('gene|symbol', cn))) colnames(ann_df)[grep('gene|symbol', cn)[1]] else NULL
    sig_annot$chr <- NA_character_
    sig_annot$pos <- NA_real_
    sig_annot$gene <- NA_character_
    rn <- rownames(ann_df)
    mm <- match(sig_annot$cpg, rn)
    if (!is.null(chr_col)) sig_annot$chr <- as.character(ann_df[[chr_col]])[mm]
    if (!is.null(pos_col)) sig_annot$pos <- as.numeric(as.character(ann_df[[pos_col]]))[mm]
    if (!is.null(gene_col)) sig_annot$gene <- as.character(ann_df[[gene_col]])[mm]
  }
  write.csv(sig_annot, 'sig_cpgs.csv', row.names=FALSE)
  
  up_count <- sum(sig_annot$logFC > 0, na.rm=TRUE)
  down_count <- sum(sig_annot$logFC < 0, na.rm=TRUE)
  message('Reiksmingi Cpgs (p<=', sig_thr, '): ', nrow(sig_annot), ' (up=', up_count, ', down=', down_count, ')')
  
  
  
  # Manhattano grafikiukas - rikiavau chromasomas pagal numerius
  if (exists('ann_df')) {
    man_df <- data.frame(cpg = out_res$cpg, P = out_res$P.Value, stringsAsFactors=FALSE)
    rn <- rownames(ann_df)
    mm2 <- match(man_df$cpg, rn)
    if (!is.null(chr_col)) man_df$chr <- as.character(ann_df[[chr_col]])[mm2]
    if (!is.null(pos_col)) man_df$pos <- as.numeric(as.character(ann_df[[pos_col]]))[mm2]
    man_df <- man_df[!is.na(man_df$chr) & !is.na(man_df$pos) & !is.na(man_df$P), , drop=FALSE]
    if (nrow(man_df) > 0) {
      man_df$chr2 <- as.character(man_df$chr)
      chr_max <- tapply(man_df$pos, man_df$chr2, max, na.rm=TRUE)
      numeric_names <- suppressWarnings(as.numeric(gsub('[^0-9]', '', names(chr_max))))
      if (all(is.na(numeric_names))) chr_order <- names(chr_max) else chr_order <- names(chr_max)[order(numeric_names, na.last=TRUE)]
      offsets <- setNames(c(0, cumsum(as.numeric(chr_max[chr_order]))[-length(chr_order)]), chr_order)
      man_df$cumpos <- man_df$pos + offsets[man_df$chr2]
      png('manhattan_all.png', width=1200, height=600)
      with(man_df, plot(cumpos, -log10(P), pch=20, cex=0.6, xaxt='n', xlab='Genomic position', ylab='-log10(P)', main='Manhattan plot'))
      ticks <- tapply(man_df$cumpos, man_df$chr2, mean)
      axis(1, at=ticks, labels=names(ticks), las=2, cex.axis=0.7)
      dev.off()
      
      # Tikais top chromasomos
      topchr <- names(sort(table(man_df$chr2[man_df$P < 0.05]), decreasing=TRUE))[1]
      if (!is.na(topchr) && nzchar(topchr)) {
        subdf <- man_df[man_df$chr2 == topchr, , drop=FALSE]
        if (nrow(subdf) > 0) {
          png('manhattan_topchr.png', width=1000, height=500)
          plot(subdf$pos, -log10(subdf$P), pch=20, cex=0.6, xlab=paste('Position, chr', topchr), ylab='-log10(P)', main=paste('Manhattan, chr', topchr))
          dev.off()
        }
      }
    }
  }
}

# Hovrath pragaro kodas 
coef_file_candidates <- c('13059_2013_3156_MOESM3_ESM.csv','horvath_coef.csv')
coef_file <- coef_file_candidates[file.exists(coef_file_candidates)][1]
if (is.na(coef_file) || is.null(coef_file)) {
  message('Hovrath failo nera :(')
} else {
  rl <- readLines(coef_file, n = 200)
  header_line <- NA_integer_
  for (i in seq_along(rl)) {
    line <- rl[i]
    if (!grepl(',', line)) next
    first_field <- tolower(gsub('^\\s*"|"\\s*$','', strsplit(line,',')[[1]][1], perl=TRUE))
    if (grepl('^cpg', first_field)) { header_line <- i; break }
  }
  if (!is.na(header_line)) {
    coef_df <- read.csv(coef_file, skip = header_line - 1, stringsAsFactors=FALSE, check.names=FALSE)
  } else {
    coef_df <- read.csv(coef_file, stringsAsFactors=FALSE, check.names=FALSE)
  }
  cn <- tolower(colnames(coef_df))
  idx_cpg <- grep('cpg|probeset|cg|cpgmarker', cn)
  idx_coef <- grep('coef|coefficient|coefficienttraining|training', cn)
  if (length(idx_cpg) == 0 || length(idx_coef) == 0) stop('Could not detect CpG/coefficient columns in Horvath file')
  cpg_col <- colnames(coef_df)[idx_cpg[1]]
  coef_col <- colnames(coef_df)[idx_coef[1]]
  message(sprintf('Horvath filas, kur CpG column "%s" ir Coef column "%s"', cpg_col, coef_col))
 
   # Numric koficientus sudarome
  raw_coefs <- coef_df[[coef_col]]
  coef_nums <- suppressWarnings(as.numeric(gsub(',', '.', as.character(raw_coefs))))
  coef_tab <- data.frame(CpG = as.character(coef_df[[cpg_col]]), Coef = coef_nums, stringsAsFactors=FALSE)
  na_coefs <- sum(is.na(coef_tab$Coef))
  if (na_coefs > 0) warning('NAs padarytas: ', na_coefs)
  intercept <- 0
  if (any(tolower(coef_tab$CpG) %in% c('(intercept)','intercept'))) {
    ii <- which(tolower(coef_tab$CpG) %in% c('(intercept)','intercept'))[1]
    if (!is.na(coef_tab$Coef[ii])) intercept <- coef_tab$Coef[ii]
    coef_tab <- coef_tab[!seq_len(nrow(coef_tab)) %in% ii, , drop=FALSE]
  }
  # Palikau tik tas dalis kurios turejo koficientus
  coef_tab <- coef_tab[!is.na(coef_tab$Coef) & coef_tab$CpG != '' & !is.na(coef_tab$CpG), , drop=FALSE]
  if (nrow(coef_tab) == 0) stop('Horvathe nera numeriniu cofficientu')
  common <- intersect(coef_tab$CpG, rownames(beta_c))
  message('Horvath: matched CpGs in data: ', length(common), ' / ', nrow(coef_tab))
  if (length(common) == 0) stop('Cpgs nesusiejami su duomenimis')
  if (length(common) < 10) warning('Tik keli Horvath cpgs susieti: ', length(common))
  cf <- coef_tab[match(common, coef_tab$CpG), , drop=FALSE]
  mat <- beta_c[cf$CpG, , drop=FALSE]
  if (any(is.na(mat))) {
    col_means <- colMeans(mat, na.rm=TRUE)
    for (k in seq_len(ncol(mat))) mat[is.na(mat[,k]), k] <- col_means[k]
  }
  # Vektoriu ir matricu uztikrinimas (kad poto nebreakingu del neteisingos rusies kodo - ilgai del to nes numeric nebuvo strigo)
  coef_vec <- as.numeric(cf$Coef)
  matm <- as.matrix(mat)
  pred_raw <- as.numeric(coef_vec %*% matm) + intercept
  n_finite_pred <- sum(is.finite(pred_raw))
  message('Horvatho RAW predicitons: ', n_finite_pred, ' / ', length(pred_raw))
  calib_idx <- which(is.finite(pred_raw) & !is.na(samplekey_c$age))
  message('Tinkami: ', length(calib_idx))
  if (length(calib_idx) < 5) {
    warning('Per mazai nuskaitoma tinkamu duomenu')
    horvath_age <- pred_raw
  } else {
    fit_cal <- lm(samplekey_c$age[calib_idx] ~ pred_raw[calib_idx])
    horvath_age <- as.numeric(predict(fit_cal, newdata = data.frame(pred_raw = pred_raw)))
    cc <- which(is.finite(horvath_age) & !is.na(samplekey_c$age))
    if (length(cc)>0) {
      mae <- mean(abs(horvath_age[cc] - samplekey_c$age[cc]), na.rm=TRUE)
      rmse <- sqrt(mean((horvath_age[cc] - samplekey_c$age[cc])^2, na.rm=TRUE))
      corr <- suppressWarnings(cor(horvath_age[cc], samplekey_c$age[cc], use='complete.obs'))
      message(sprintf('Horvath calib: n=%d, MAE=%.2f, RMSE=%.2f, Corr=%.3f', length(cc), mae, rmse, ifelse(is.na(corr), NA, corr)))
    }
  }
  # Outputai ir residualu grafikas
  write.csv(data.frame(sample=colnames(beta_c), horvath_age=horvath_age, chron_age=samplekey_c$age), 'horvath_comparison.csv', row.names=FALSE)
  if (exists('coef_tab') && is.data.frame(coef_tab)) write.csv(coef_tab, 'horvath_diag.csv', row.names=FALSE)
  if (exists('horvath_age')) {
    resid <- horvath_age - samplekey_c$age
    png('horvath_residuals.png', width=900, height=700)
    plot(samplekey_c$age, resid, pch=16, xlab='Chronological age', ylab='Horvath - Chron age', main='Horvath residuals vs age')
    abline(h=0, col='red', lty=2)
    dev.off()
    write.csv(data.frame(sample=colnames(beta_c), horvath_age=horvath_age, chron_age=samplekey_c$age, residual=as.numeric(resid)), 'horvath_with_residuals.csv', row.names=FALSE)
  }
  # Plotina grafika jeigu neloopina begalybe kartu - nebuvau tikras ar blogai data nuskaitinejo bet cia labai daznai pastrikdavo
  finite_idx <- which(is.finite(horvath_age) & is.finite(samplekey_c$age))
  if (length(finite_idx) > 1) {
    png('horvath_clock.png', width=900, height=700)
    plot(samplekey_c$age[finite_idx], horvath_age[finite_idx], pch=16, xlab='Chron age', ylab='Horvath age')
    abline(0,1,col='red',lty=2)
    dev.off()
  } else {
    message('neveikia horvath!!!.')
  }
}

cat('VEIKIA!!! RUNNINA KODAS YIPPIE!!!.\n')
# 