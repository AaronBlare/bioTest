rm(list=ls())

library(data.table)
library(sesame)
library(sesameData)
library(parallel)
library(ggplot2)

sesameDataCache()

sample_sheet_path <- "D:/Yandex.Disk/DNAm draft/Lesnoy_CVD/120/raw/samples_189.csv"
idat_dir <- dirname(sample_sheet_path)

ss <- fread(sample_sheet_path)
setDT(ss)

required_cols <- c(
  "Sample_ID", "sex", "age", "Sentrix_Position", "Sentrix_ID", "basename"
)

missing_cols <- setdiff(required_cols, names(ss))
if (length(missing_cols) > 0) {
  stop(
    "Sample sheet is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

ss[, Sample_ID := as.character(Sample_ID)]
ss[, sex := as.character(sex)]
ss[, age := as.numeric(age)]
ss[, basename := as.character(basename)]
ss[, Sentrix_ID := as.character(Sentrix_ID)]
ss[, Sentrix_Position := as.character(Sentrix_Position)]

if (anyDuplicated(ss$Sample_ID) > 0) {
  stop("Duplicated Sample_ID values found in sample sheet.")
}

if (anyDuplicated(ss$basename) > 0) {
  stop("Duplicated basename values found in sample sheet.")
}

idat_files <- list.files(
  path = idat_dir,
  pattern = "\\.idat(\\.gz)?$",
  full.names = TRUE
)

if (length(idat_files) == 0L) {
  stop("No IDAT files found in idat_dir: ", idat_dir)
}

idat_dt <- data.table(full_path = idat_files)
idat_dt[, fname := basename(full_path)]

idat_dt <- idat_dt[
  grepl("_(Grn|Red)\\.idat(\\.gz)?$", fname)
]

if (nrow(idat_dt) == 0L) {
  stop(
    "No *_Grn.idat(.gz) or *_Red.idat(.gz) files found in idat_dir: ",
    idat_dir
  )
}

## Shared file prefix used by readIDATpair/openSesame
idat_dt[, idat_prefix := sub("_(Grn|Red)\\.idat(\\.gz)?$", "", fname)]

## Check that each prefix has both colour channels
idat_channels <- idat_dt[, .(
  n_files = .N,
  has_grn = any(grepl("_Grn\\.idat(\\.gz)?$", fname)),
  has_red = any(grepl("_Red\\.idat(\\.gz)?$", fname))
), by = idat_prefix]

bad_prefixes <- idat_channels[!(has_grn & has_red)]
if (nrow(bad_prefixes) > 0L) {
  stop(
    "Some IDAT prefixes do not have both Grn and Red files. Examples: ",
    paste(head(bad_prefixes$idat_prefix, 10), collapse = ", ")
  )
}

idat_dt_unique <- unique(idat_dt[, .(idat_prefix)])

## NN sample sheet basename already equals the real IDAT prefix
ss_pref <- merge(
  ss,
  idat_dt_unique,
  by.x = "basename",
  by.y = "idat_prefix",
  all.x = TRUE,
  sort = FALSE
)

## Recreate explicit idat_prefix column for downstream chunks
ss_pref[, idat_prefix := basename]
ss_pref[, basename_matches_prefix := basename == idat_prefix]

message("Rows in sample sheet: ", nrow(ss))
message("Unique IDAT prefixes found: ", nrow(idat_dt_unique))
message("Samples with missing matched IDAT prefix: ", sum(is.na(ss_pref$idat_prefix)))
message("Duplicated Sample_ID after merge: ", sum(duplicated(ss_pref$Sample_ID)))
message(
  "basename == idat_prefix count: ",
  sum(ss_pref$basename_matches_prefix, na.rm = TRUE),
  " / ",
  nrow(ss_pref)
)

if (nrow(ss_pref) != nrow(ss)) {
  stop(
    "Merged sample sheet row count (", nrow(ss_pref),
    ") does not match original sample sheet row count (", nrow(ss), ")."
  )
}

if (sum(is.na(ss_pref$idat_prefix)) > 0) {
  stop(
    "Some samples could not be matched to IDAT prefixes via basename. ",
    "Inspect ss_pref[is.na(idat_prefix)]."
  )
}

if (sum(duplicated(ss_pref$Sample_ID)) > 0) {
  stop("Duplicated Sample_ID values detected after merging sample sheet with IDAT prefixes.")
}

setcolorder(ss_pref, c(
  "Sample_ID", "sex", "age",
  "Sentrix_ID", "Sentrix_Position",
  "basename", "idat_prefix", "basename_matches_prefix"
))

print(head(ss_pref))

rm(idat_files, idat_dt, idat_dt_unique, idat_channels, bad_prefixes)

mft <- sesameAnno_buildManifestGRanges(
        sesameAnno_download("EPICv2.hg38.manifest.tsv.gz"),
        columns = "nextBase"
    )
extR <- names(mft)[!is.na(mft$nextBase) & mft$nextBase == "R"]
extA <- names(mft)[!is.na(mft$nextBase) & mft$nextBase == "A"]

process_one_sample <- function(id_prefix, sample_id, idat_dir, verbose = TRUE) {
  start_time <- Sys.time()
  idat_path <- file.path(idat_dir, id_prefix)

  if (verbose) {
    message(sprintf("[%s] START sample=%s prefix=%s",
                    format(start_time, "%Y-%m-%d %H:%M:%S"), sample_id, id_prefix))
  }

  out <- tryCatch({

    if (verbose) {
      message(sprintf("[%s] %s : reading SigDF",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sample_id))
    }

    sdf_raw <- readIDATpair(idat_path)

    if (verbose) {
      message(sprintf("[%s] %s : computing sample-level QC stats",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sample_id))
    }

    qc_obj <- openSesame(
      idat_path,
      prep = "",
      func = sesameQC_calcStats
    )

    get_stat <- function(name) {
      out <- tryCatch(sesameQC_getStats(qc_obj, name), error = function(e) NA_real_)
      if (is.null(out) || length(out) == 0) return(NA_real_)
      out_num <- suppressWarnings(as.numeric(out))
      if (length(out_num) == 0 || is.na(out_num)) return(NA_real_)
      out_num
    }

    frac_dt <- get_stat("frac_dt")
    frac_dt_cg <- get_stat("frac_dt_cg")
    mean_intensity <- get_stat("mean_intensity")
    mean_intensity_total <- get_stat("mean_intensity_MU")
    RGratio <- get_stat("RGratio")
    RGdistort <- get_stat("RGdistort")

    bs_gct_score <- tryCatch({
      out <- as.numeric(bisConversionControl(sdf_raw, extR = extR, extA = extA))
      if (length(out) == 0 || is.na(out)) NA_real_ else out
    }, error = function(e) NA_real_)

    rank_fraction <- NA_real_
    rank_n_ref <- NA_real_
    qc_rank <- NULL

    if (verbose) {
      message(sprintf("[%s] %s : generating betas with openSesame(prep='CDB')",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sample_id))
    }

    betas_unmasked <- openSesame(idat_path, prep = "CDB")
    betas_masked <- betas_unmasked

    sex_call <- tryCatch({
      out <- inferSex(betas_unmasked)
      if (length(out) == 0 || is.null(out)) NA_character_ else as.character(out)
    }, error = function(e) NA_character_)

    end_time <- Sys.time()
    elapsed <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)

    if (verbose) {
      message(sprintf("[%s] DONE sample=%s elapsed=%ss frac_failed=%.4f mean_intensity=%.2f bs_gct=%.4f",
                      format(end_time, "%Y-%m-%d %H:%M:%S"), sample_id, elapsed,
                      if (is.na(frac_dt)) NA_real_ else 1 - frac_dt,
                      mean_intensity, bs_gct_score))
    }

    list(
      sample_id = as.character(sample_id),
      sdf_raw = sdf_raw,
      qc_obj = qc_obj,
      qc_rank = qc_rank,
      betas_unmasked = betas_unmasked,
      betas_masked = betas_masked,
      sex_inferred = sex_call,
      frac_dt = frac_dt,
      frac_dt_cg = frac_dt_cg,
      frac_failed = if (is.na(frac_dt)) NA_real_ else 1 - frac_dt,
      frac_failed_cg = if (is.na(frac_dt_cg)) NA_real_ else 1 - frac_dt_cg,
      mean_intensity = mean_intensity,
      mean_intensity_total = mean_intensity_total,
      RGratio = RGratio,
      RGdistort = RGdistort,
      bs_gct_score = bs_gct_score,
      rank_fraction = rank_fraction,
      rank_n_ref = rank_n_ref,
      status = "OK",
      error_message = NA_character_,
      start_time = start_time,
      end_time = end_time,
      elapsed_sec = elapsed
    )

  }, error = function(e) {
    end_time <- Sys.time()
    elapsed <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)

    message(sprintf("[%s] ERROR sample=%s elapsed=%ss message=%s",
                    format(end_time, "%Y-%m-%d %H:%M:%S"), sample_id, elapsed, conditionMessage(e)))

    list(
      sample_id = as.character(sample_id),
      sdf_raw = NULL,
      qc_obj = NULL,
      qc_rank = NULL,
      betas_unmasked = NULL,
      betas_masked = NULL,
      sex_inferred = NA_character_,
      frac_dt = NA_real_,
      frac_dt_cg = NA_real_,
      frac_failed = NA_real_,
      frac_failed_cg = NA_real_,
      mean_intensity = NA_real_,
      mean_intensity_total = NA_real_,
      RGratio = NA_real_,
      RGdistort = NA_real_,
      bs_gct_score = NA_real_,
      rank_fraction = NA_real_,
      rank_n_ref = NA_real_,
      status = "ERROR",
      error_message = as.character(conditionMessage(e)),
      start_time = start_time,
      end_time = end_time,
      elapsed_sec = elapsed
    )
  })

  out
}

library(parallel)
library(data.table)

n_cores <- 1

run_start <- Sys.time()
message(sprintf("[%s] START full processing: n_samples=%d, n_cores=%d",
                format(run_start, "%Y-%m-%d %H:%M:%S"), nrow(ss_pref), n_cores))

results_list <- mcmapply(
  FUN = process_one_sample,
  id_prefix = ss_pref$idat_prefix,
  sample_id = ss_pref$Sample_ID,
  MoreArgs = list(
    idat_dir = idat_dir,
    verbose = TRUE
  ),
  SIMPLIFY = FALSE,
  mc.cores = n_cores
)

names(results_list) <- as.character(ss_pref$Sample_ID)

run_end <- Sys.time()
message(sprintf("[%s] END full processing: elapsed=%.2f min",
                format(run_end, "%Y-%m-%d %H:%M:%S"),
                as.numeric(difftime(run_end, run_start, units = "mins"))))

message("Sanity checks on results_list:")
print(length(results_list))
print(table(vapply(results_list, is.null, logical(1)), useNA = "ifany"))
print(table(vapply(results_list, function(x) x$status, character(1)), useNA = "ifany"))
print(table(vapply(results_list, function(x) length(x$mean_intensity_total), integer(1)), useNA = "ifany"))

qc_dt <- rbindlist(lapply(results_list, function(x) {
  data.table(
    Sample_ID = as.character(x$sample_id),
    frac_dt = if (length(x$frac_dt) == 0) NA_real_ else as.numeric(x$frac_dt),
    frac_dt_cg = if (length(x$frac_dt_cg) == 0) NA_real_ else as.numeric(x$frac_dt_cg),
    frac_failed = if (length(x$frac_failed) == 0) NA_real_ else as.numeric(x$frac_failed),
    frac_failed_cg = if (length(x$frac_failed_cg) == 0) NA_real_ else as.numeric(x$frac_failed_cg),
    mean_intensity = if (length(x$mean_intensity) == 0) NA_real_ else as.numeric(x$mean_intensity),
    mean_intensity_total = if (length(x$mean_intensity_total) == 0) NA_real_ else as.numeric(x$mean_intensity_total),
    RGratio = if (length(x$RGratio) == 0) NA_real_ else as.numeric(x$RGratio),
    RGdistort = if (length(x$RGdistort) == 0) NA_real_ else as.numeric(x$RGdistort),
    bs_gct_score = if (length(x$bs_gct_score) == 0) NA_real_ else as.numeric(x$bs_gct_score),
    sex_inferred = if (length(x$sex_inferred) == 0) NA_character_ else as.character(x$sex_inferred),
    rank_fraction = if (length(x$rank_fraction) == 0) NA_real_ else as.numeric(x$rank_fraction),
    rank_n_ref = if (length(x$rank_n_ref) == 0) NA_real_ else as.numeric(x$rank_n_ref),
    status = if (length(x$status) == 0) NA_character_ else as.character(x$status),
    error_message = if (length(x$error_message) == 0) NA_character_ else as.character(x$error_message),
    start_time = if (length(x$start_time) == 0) NA_character_ else as.character(x$start_time),
    end_time = if (length(x$end_time) == 0) NA_character_ else as.character(x$end_time),
    elapsed_sec = if (length(x$elapsed_sec) == 0) NA_real_ else as.numeric(x$elapsed_sec)
  )
}), fill = TRUE)

ss_qc <- merge(ss_pref, qc_dt, by = "Sample_ID", all.x = TRUE)

message("QC table summaries:")
print(table(ss_qc$status, useNA = "ifany"))
print(summary(ss_qc$frac_failed))
print(summary(ss_qc$mean_intensity))
print(summary(ss_qc$bs_gct_score))

setDT(ss_qc)

extract_qc_stats <- function(x) {
  if (is.null(x) || is.null(x$qc_obj)) {
    return(data.table(
      Sample_ID = as.character(x$sample_id),
      mean_beta = NA_real_,
      median_beta = NA_real_,
      frac_unmeth = NA_real_,
      frac_meth = NA_real_,
      mean_intensity_MU = NA_real_,
      mean_ii = NA_real_,
      mean_oob_grn = NA_real_,
      mean_oob_red = NA_real_,
      num_na_cg = NA_real_,
      frac_na_cg = NA_real_,
      num_na_rs = NA_real_,
      frac_na_rs = NA_real_,
      InfI_switch_G2R = NA_real_,
      InfI_switch_R2G = NA_real_
    ))
  }

  s <- tryCatch(x$qc_obj@stat, error = function(e) NULL)

  pick_num <- function(nm) {
    if (is.null(s) || !(nm %in% names(s)) || is.null(s[[nm]]) || length(s[[nm]]) == 0) {
      return(NA_real_)
    }
    out <- suppressWarnings(as.numeric(s[[nm]])[1])
    if (length(out) == 0 || is.na(out)) return(NA_real_)
    out
  }

  data.table(
    Sample_ID = as.character(x$sample_id),
    mean_beta = pick_num("mean_beta"),
    median_beta = pick_num("median_beta"),
    frac_unmeth = pick_num("frac_unmeth"),
    frac_meth = pick_num("frac_meth"),
    mean_intensity_MU = pick_num("mean_intensity_MU"),
    mean_ii = pick_num("mean_ii"),
    mean_oob_grn = pick_num("mean_oob_grn"),
    mean_oob_red = pick_num("mean_oob_red"),
    num_na_cg = pick_num("num_na_cg"),
    frac_na_cg = pick_num("frac_na_cg"),
    num_na_rs = pick_num("num_na_rs"),
    frac_na_rs = pick_num("frac_na_rs"),
    InfI_switch_G2R = pick_num("InfI_switch_G2R"),
    InfI_switch_R2G = pick_num("InfI_switch_R2G")
  )
}

qc_extra <- rbindlist(lapply(results_list, extract_qc_stats), fill = TRUE)
ss_qc_full <- merge(ss_qc, qc_extra, by = "Sample_ID", all.x = TRUE)
setDT(ss_qc_full)

ss_qc_full[, sex_reported_std := toupper(substr(as.character(sex), 1, 1))]
ss_qc_full[, sex_inferred_std := toupper(substr(as.character(sex_inferred), 1, 1))]
ss_qc_full[, sex_mismatch := !is.na(sex_reported_std) &
                           !is.na(sex_inferred_std) &
                           sex_reported_std != sex_inferred_std]

frac_failed_cutoff_soft <- 0.05
frac_failed_cutoff_hard <- 0.15

mean_intensity_min <- quantile(ss_qc_full$mean_intensity, 0.05, na.rm = TRUE)
bs_gct_max <- quantile(ss_qc_full$bs_gct_score, 0.95, na.rm = TRUE)

message("QC thresholds:")
message("  frac_failed soft warning > ", frac_failed_cutoff_soft)
message("  frac_failed hard fail > ", frac_failed_cutoff_hard)
message("  mean_intensity 5th percentile cutoff = ", round(mean_intensity_min, 2))
message("  bs_gct_score 95th percentile cutoff = ", round(bs_gct_max, 4))

ss_qc_full[, fail_frac_soft := !is.na(frac_failed) & frac_failed > frac_failed_cutoff_soft]
ss_qc_full[, fail_frac_hard := !is.na(frac_failed) & frac_failed > frac_failed_cutoff_hard]
ss_qc_full[, warn_frac := fail_frac_soft & !fail_frac_hard]

ss_qc_full[, fail_int := !is.na(mean_intensity) & mean_intensity < mean_intensity_min]
ss_qc_full[, fail_bs := !is.na(bs_gct_score) & bs_gct_score > bs_gct_max]
ss_qc_full[, fail_status := is.na(status) | status != "OK"]
ss_qc_full[, fail_sex := !is.na(sex_mismatch) & sex_mismatch]

ss_qc_full[, qc_fail := fail_frac_hard | fail_int | fail_bs | fail_status | fail_sex]

setcolorder(ss_qc_full, c(
  "Sample_ID", "sex", "age",
  "Sentrix_ID", "Sentrix_Position", "basename", "idat_prefix", "basename_matches_prefix",
  "frac_dt", "frac_dt_cg", "frac_failed", "frac_failed_cg",
  "mean_intensity", "mean_intensity_total", "mean_intensity_MU", "mean_ii",
  "mean_oob_grn", "mean_oob_red", "RGratio", "RGdistort", "bs_gct_score",
  "mean_beta", "median_beta", "frac_unmeth", "frac_meth",
  "num_na_cg", "frac_na_cg", "num_na_rs", "frac_na_rs",
  "InfI_switch_G2R", "InfI_switch_R2G",
  "sex_inferred", "sex_reported_std", "sex_inferred_std", "sex_mismatch",
  "rank_fraction", "rank_n_ref",
  "fail_frac_soft", "fail_frac_hard", "warn_frac", "fail_int", "fail_bs", "fail_status", "fail_sex", "qc_fail",
  "status", "error_message", "start_time", "end_time", "elapsed_sec"
))

fwrite(ss_qc_full, "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_QC_table.csv")

message("QC fail table:")
print(table(ss_qc_full$qc_fail, useNA = "ifany"))

message("Status table:")
print(table(ss_qc_full$status, useNA = "ifany"))

message("Summary of frac_failed:")
print(summary(ss_qc_full$frac_failed))

message("Summary of mean_intensity:")
print(summary(ss_qc_full$mean_intensity))

message("Summary of bs_gct_score:")
print(summary(ss_qc_full$bs_gct_score))

qc_overall <- ss_qc_full[, .(
  n_total = .N,
  n_fail_frac_soft = sum(fail_frac_soft, na.rm = TRUE),
  n_fail_frac_hard = sum(fail_frac_hard, na.rm = TRUE),
  n_warn_frac = sum(warn_frac, na.rm = TRUE),
  n_fail_int = sum(fail_int, na.rm = TRUE),
  n_fail_bs = sum(fail_bs, na.rm = TRUE),
  n_fail_status = sum(fail_status, na.rm = TRUE),
  n_fail_sex = sum(fail_sex, na.rm = TRUE),
  n_fail_any = sum(qc_fail, na.rm = TRUE),
  n_pass_all = sum(qc_fail == FALSE, na.rm = TRUE)
)]

qc_overall[, `:=`(
  perc_fail_any = round(100 * n_fail_any / n_total, 1),
  perc_pass_all = round(100 * n_pass_all / n_total, 1),
  perc_fail_sex = round(100 * n_fail_sex / n_total, 1),
  perc_fail_hard = round(100 * n_fail_frac_hard / n_total, 1),
  perc_warn_frac = round(100 * n_warn_frac / n_total, 1),
  perc_fail_bs = round(100 * n_fail_bs / n_total, 1)
)]

fwrite(qc_overall, "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_QC_overall_summary.csv")

message("Overall QC summary:")
print(qc_overall)

library(data.table)
setDT(ss_qc_full)

# Define quantile cut‑offs for "high background"
oob_grn_p95 <- quantile(ss_qc_full$mean_oob_grn, 0.95, na.rm = TRUE)
oob_red_p95 <- quantile(ss_qc_full$mean_oob_red, 0.95, na.rm = TRUE)

message("OOB 95th percentile cutoffs: ",
        "mean_oob_grn = ", round(oob_grn_p95, 1),
        ", mean_oob_red = ", round(oob_red_p95, 1))

# Flags
ss_qc_full[, oob_grn_high := !is.na(mean_oob_grn) & mean_oob_grn > oob_grn_p95]
ss_qc_full[, oob_red_high := !is.na(mean_oob_red) & mean_oob_red > oob_red_p95]

# Combined background outlier flag
ss_qc_full[, oob_outlier := oob_grn_high | oob_red_high]

# Optional: a more conservative flag requiring both channels high
ss_qc_full[, oob_outlier_strict := oob_grn_high & oob_red_high]

message("Background outliers (any channel > 95th pct): ",
        sum(ss_qc_full$oob_outlier, na.rm = TRUE))
message("Background outliers (both channels > 95th pct): ",
        sum(ss_qc_full$oob_outlier_strict, na.rm = TRUE))

# If you want these in the saved QC CSV, re‑write it:
fwrite(ss_qc_full, "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_QC_table.csv")

library(ggplot2)
library(data.table)

keep_ids <- ss_qc_full[qc_fail == FALSE & status == "OK", Sample_ID]
keep_ids <- intersect(as.character(keep_ids), names(results_list))

message("Number of QC-passed samples available for density plots: ", length(keep_ids))

if (length(keep_ids) == 0) {
  stop("No QC-passed samples found for density plotting.")
}

set.seed(1)
n_cpg_sub <- 50000

beta_long_list <- lapply(keep_ids, function(sid) {
  x <- results_list[[sid]]
  if (is.null(x) || is.null(x$betas_masked)) return(NULL)

  betas <- x$betas_masked
  if (length(betas) == 0) return(NULL)

  cpg_ids <- names(betas)
  if (is.null(cpg_ids) || length(cpg_ids) == 0) return(NULL)

  good_idx <- !is.na(betas)
  if (!any(good_idx)) return(NULL)

  cpg_ids <- cpg_ids[good_idx]
  if (length(cpg_ids) == 0) return(NULL)

  if (length(cpg_ids) > n_cpg_sub) {
    cpg_ids <- sample(cpg_ids, n_cpg_sub)
  }

  data.table(
    Sample_ID = as.character(sid),
    Beta = as.numeric(betas[cpg_ids])
  )
})

beta_long <- rbindlist(beta_long_list, use.names = TRUE, fill = TRUE)

if (nrow(beta_long) == 0L) {
  stop("beta_long is empty: no usable masked betas for density plotting.")
}

beta_long <- merge(
  beta_long,
  ss_qc_full[, .(Sample_ID)],
  by = "Sample_ID",
  all.x = TRUE
)

beta_long <- beta_long[!is.na(Beta) & Beta >= 0 & Beta <= 1]

message("Rows in beta_long used for plotting: ", nrow(beta_long))
message("Unique samples represented in beta_long: ", uniqueN(beta_long$Sample_ID))

p_all <- ggplot(beta_long, aes(x = Beta, group = Sample_ID)) +
  geom_density(alpha = 0.2, linewidth = 0.4) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Per-sample Beta-value density (masked betas, QC-passed)",
    x = "Beta-value",
    y = "Density"
  ) +
  theme_bw(base_size = 12)

example_ids <- sample(unique(beta_long$Sample_ID), min(5, uniqueN(beta_long$Sample_ID)))
beta_example <- beta_long[Sample_ID %in% example_ids]

p_example <- ggplot(beta_example, aes(x = Beta, colour = Sample_ID)) +
  geom_density(linewidth = 0.6) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Beta-value densities for example samples (masked betas)",
    x = "Beta-value",
    y = "Density",
    colour = "Sample_ID"
  ) +
  theme_bw(base_size = 11)

print(p_all)
print(p_example)

library(ggplot2)
library(data.table)

plots_dir <- "C:/Git/bioTest/r_olga_script_test/plots"
if (!dir.exists(plots_dir)) dir.create(plots_dir)

p_frac_failed <- ggplot(ss_qc_full, aes(x = frac_failed)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "grey40") +
  geom_vline(xintercept = 0.05, colour = "orange", linetype = 2, linewidth = 0.7) +
  geom_vline(xintercept = 0.15, colour = "red", linetype = 2, linewidth = 0.7) +
  labs(
    title = "Fraction of failed CpGs per sample",
    x = "frac_failed",
    y = "Number of samples"
  ) +
  theme_bw(base_size = 11)

p_mean_intensity <- ggplot(ss_qc_full[!is.na(mean_intensity)], aes(x = mean_intensity)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "grey40") +
  geom_vline(
    xintercept = quantile(ss_qc_full$mean_intensity, 0.05, na.rm = TRUE),
    colour = "red", linetype = 2, linewidth = 0.7
  ) +
  labs(
    title = "Mean intensity per sample",
    x = "mean_intensity",
    y = "Number of samples"
  ) +
  theme_bw(base_size = 11)

p_bs <- ggplot(ss_qc_full[!is.na(bs_gct_score)], aes(x = bs_gct_score)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "grey40") +
  geom_vline(
    xintercept = quantile(ss_qc_full$bs_gct_score, 0.95, na.rm = TRUE),
    colour = "red", linetype = 2, linewidth = 0.7
  ) +
  labs(
    title = "Bisulfite conversion GCT score per sample",
    x = "bs_gct_score",
    y = "Number of samples"
  ) +
  theme_bw(base_size = 11)

p_rg_ratio <- ggplot(ss_qc_full[!is.na(RGratio)], aes(x = RGratio)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "grey40") +
  labs(
    title = "Red/Green ratio per sample",
    x = "RGratio",
    y = "Number of samples"
  ) +
  theme_bw(base_size = 11)

p_rg_distort <- ggplot(ss_qc_full[!is.na(RGdistort)], aes(x = RGdistort)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "grey40") +
  labs(
    title = "Red/Green distortion per sample",
    x = "RGdistort",
    y = "Number of samples"
  ) +
  theme_bw(base_size = 11)

p_mean_beta <- ggplot(ss_qc_full[!is.na(mean_beta)], aes(x = mean_beta)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "grey40") +
  labs(
    title = "Mean Beta per sample",
    x = "mean_beta",
    y = "Number of samples"
  ) +
  theme_bw(base_size = 11)

p_median_beta <- ggplot(ss_qc_full[!is.na(median_beta)], aes(x = median_beta)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "grey40") +
  labs(
    title = "Median Beta per sample",
    x = "median_beta",
    y = "Number of samples"
  ) +
  theme_bw(base_size = 11)

p_oob <- ggplot(
  ss_qc_full[!is.na(mean_oob_grn) & !is.na(mean_oob_red)],
  aes(x = mean_oob_grn, y = mean_oob_red, colour = qc_fail)
) +
  geom_point(alpha = 0.8) +
  labs(
    title = "Out-of-band background intensities",
    x = "mean_oob_grn",
    y = "mean_oob_red",
    colour = "qc_fail"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(plots_dir, "QC_frac_failed_hist.pdf"),
       plot = p_frac_failed, width = 6, height = 4)
ggsave(file.path(plots_dir, "QC_mean_intensity_hist.pdf"),
       plot = p_mean_intensity, width = 6, height = 4)
ggsave(file.path(plots_dir, "QC_bs_gct_hist.pdf"),
       plot = p_bs, width = 6, height = 4)
ggsave(file.path(plots_dir, "QC_RGratio_hist.pdf"),
       plot = p_rg_ratio, width = 6, height = 4)
ggsave(file.path(plots_dir, "QC_RGdistort_hist.pdf"),
       plot = p_rg_distort, width = 6, height = 4)
ggsave(file.path(plots_dir, "QC_mean_beta_hist.pdf"),
       plot = p_mean_beta, width = 6, height = 4)
ggsave(file.path(plots_dir, "QC_median_beta_hist.pdf"),
       plot = p_median_beta, width = 6, height = 4)
ggsave(file.path(plots_dir, "QC_oob_scatter.pdf"),
       plot = p_oob, width = 6, height = 5)

if (exists("p_all")) {
  ggsave(file.path(plots_dir, "DNAm_beta_density_all_samples.pdf"),
         plot = p_all, width = 7, height = 5)
  ggsave(file.path(plots_dir, "DNAm_beta_density_all_samples.png"),
         plot = p_all, width = 7, height = 5, dpi = 300)
}

if (exists("p_example")) {
  ggsave(file.path(plots_dir, "DNAm_beta_density_example_samples.pdf"),
         plot = p_example, width = 6, height = 4)
}

ggsave(file.path(plots_dir, "QC_oob_scatter.png"),
       plot = p_oob, width = 6, height = 5, dpi = 300)

message("All key QC plots saved in: ", plots_dir)

library(data.table)
library(ggplot2)

plots_dir <- "C:/Git/bioTest/r_olga_script_test/plots"
if (!dir.exists(plots_dir)) dir.create(plots_dir)

keep_ids_qc <- ss_qc_full[status == "OK", Sample_ID]
keep_ids_qc <- intersect(as.character(keep_ids_qc), names(results_list))

message("Number of status-OK samples entering chunk 9: ", length(keep_ids_qc))

if (length(keep_ids_qc) == 0) {
  stop("No status-OK samples available for SNP QC / SeSAMe diagnostic plots.")
}

## 1. Collect available masked betas and raw SigDFs
beta_list_ok <- lapply(keep_ids_qc, function(sid) {
  x <- results_list[[sid]]
  if (is.null(x) || is.null(x$betas_masked)) return(NULL)
  b <- x$betas_masked
  if (length(b) == 0 || is.null(names(b))) return(NULL)
  b
})
names(beta_list_ok) <- keep_ids_qc
beta_list_ok <- beta_list_ok[!vapply(beta_list_ok, is.null, logical(1))]

sdf_list_ok <- lapply(keep_ids_qc, function(sid) {
  x <- results_list[[sid]]
  if (is.null(x) || is.null(x$sdf_raw)) return(NULL)
  x$sdf_raw
})
names(sdf_list_ok) <- keep_ids_qc
sdf_list_ok <- sdf_list_ok[!vapply(sdf_list_ok, is.null, logical(1))]

qc_obj_list <- lapply(keep_ids_qc, function(sid) {
  x <- results_list[[sid]]
  if (is.null(x) || is.null(x$qc_obj)) return(NULL)
  x$qc_obj
})
names(qc_obj_list) <- keep_ids_qc
qc_obj_list <- qc_obj_list[!vapply(qc_obj_list, is.null, logical(1))]

message("Samples with non-null betas_masked: ", length(beta_list_ok))
message("Samples with non-null sdf_raw: ", length(sdf_list_ok))
message("Samples with non-null qc_obj: ", length(qc_obj_list))

## 2. Build SNP beta matrix from rs probes
if (length(beta_list_ok) < 2) {
  warning("Fewer than 2 samples with usable betas_masked; skipping SNP matrix.")
  snp_beta_mat <- NULL
  snp_qc_dt <- data.table()
} else {

  common_probes <- Reduce(intersect, lapply(beta_list_ok, names))
  rs_probes <- common_probes[grepl("^rs", common_probes, ignore.case = FALSE)]

  message("Common probes across samples: ", length(common_probes))
  message("Common rs probes across samples: ", length(rs_probes))

  if (length(rs_probes) < 10) {
    warning("Very few common rs probes detected; SNP matrix may not be informative.")
    snp_beta_mat <- NULL
    snp_qc_dt <- data.table()
  } else {

    snp_beta_mat <- sapply(names(beta_list_ok), function(sid) {
      as.numeric(beta_list_ok[[sid]][rs_probes])
    })

    rownames(snp_beta_mat) <- rs_probes
    colnames(snp_beta_mat) <- names(beta_list_ok)

    saveRDS(snp_beta_mat, file = "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_SNP_beta_matrix.rds")
    fwrite(
      data.table(Probe_ID = rownames(snp_beta_mat), snp_beta_mat),
      "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_SNP_beta_matrix.csv"
    )

    snp_qc_dt <- data.table(
      Sample_ID = colnames(snp_beta_mat),
      n_snp_probes = colSums(!is.na(snp_beta_mat))
    )

    fwrite(snp_qc_dt, "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_SNPcheck_table.csv")
  }
}

## 3. SeSAMe QC barplot across samples
if (length(qc_obj_list) > 0) {
  pdf(file.path(plots_dir, "QC_sesameQC_barplots.pdf"), width = 12, height = 8)
  tryCatch({
    sesame::sesameQC_plotBar(qc_obj_list)
  }, error = function(e) {
    message("sesameQC_plotBar failed: ", conditionMessage(e))
  })
  dev.off()
} else {
  message("Skipping sesameQC_plotBar: no usable qc_obj objects.")
}

## 4. Example single-sample diagnostic plots
first_sid <- keep_ids_qc[keep_ids_qc %in% names(sdf_list_ok)][1]

if (!is.na(first_sid) && length(first_sid) == 1 && !is.null(sdf_list_ok[[first_sid]])) {

  pdf(file.path(plots_dir, paste0("QC_redgrnQQ_", first_sid, ".pdf")), width = 6, height = 6)
  tryCatch({
    sesame::sesameQC_plotRedGrnQQ(sdf_list_ok[[first_sid]])
  }, error = function(e) {
    message("sesameQC_plotRedGrnQQ failed for ", first_sid, ": ", conditionMessage(e))
  })
  dev.off()

  pdf(file.path(plots_dir, paste0("QC_intensity_vs_beta_", first_sid, ".pdf")), width = 6, height = 6)
  tryCatch({
    sesame::sesameQC_plotIntensVsBetas(sdf_list_ok[[first_sid]])
  }, error = function(e) {
    message("sesameQC_plotIntensVsBetas failed for ", first_sid, ": ", conditionMessage(e))
  })
  dev.off()

} else {
  message("Skipping single-sample diagnostic plots: no usable sdf_raw object found.")
}

message("Chunk 9 completed: SNP matrix, SNP summary, and SeSAMe QC plots written.")

library(ggplot2)
library(data.table)

plots_dir <- "C:/Git/bioTest/r_olga_script_test/plots"
if (!dir.exists(plots_dir)) dir.create(plots_dir)

if (!exists("snp_beta_mat") || is.null(snp_beta_mat)) {
  stop("snp_beta_mat is not available; run chunk 9 before 9b.")
}

if (ncol(snp_beta_mat) < 2 || nrow(snp_beta_mat) < 10) {
  stop("snp_beta_mat has too few samples or SNP probes to build a meaningful heatmap.")
}

message("Chunk 9b: building SNP-based sample correlation heatmap.")

## 1. Sample–sample correlation matrix (using only non-missing pairs)
cor_mat <- suppressWarnings(
  stats::cor(snp_beta_mat, use = "pairwise.complete.obs", method = "pearson")
)

cor_mat[is.na(cor_mat)] <- 0  # ensure dist() works

## 2. Distance and hierarchical clustering
dist_mat <- as.dist(1 - cor_mat)
hc <- hclust(dist_mat, method = "average")

sample_order <- hc$labels[hc$order]
cor_mat_ord <- cor_mat[sample_order, sample_order, drop = FALSE]

## 3. Melt to long format and add cohort / region annotations
cor_dt <- as.data.table(cor_mat_ord, keep.rownames = "Sample1")
cor_dt <- melt(
  cor_dt,
  id.vars = "Sample1",
  variable.name = "Sample2",
  value.name = "Correlation"
)

annot <- ss_qc_full[, .(Sample_ID)]

cor_dt <- merge(
  cor_dt,
  annot,
  by.x = "Sample1",
  by.y = "Sample_ID",
  all.x = TRUE
)

cor_dt <- merge(
  cor_dt,
  annot,
  by.x = "Sample2",
  by.y = "Sample_ID",
  all.x = TRUE
)

## 4. Heatmap plot
p_heat <- ggplot(
  cor_dt,
  aes(
    x = factor(Sample1, levels = sample_order),
    y = factor(Sample2, levels = sample_order),
    fill = Correlation
  )
) +
  geom_raster() +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0.95,
    limits = c(0, 1),
    oob = scales::squish
  ) +
  labs(
    title = "Sample–sample correlation based on EPIC SNP probes",
    x = "Sample",
    y = "Sample",
    fill = "r (SNP betas)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )

## 5. Dendrogram plot
pdf(file.path(plots_dir, "SNP_sample_dendrogram.pdf"), width = 7, height = 5)
plot(
  as.dendrogram(hc),
  main = "Hierarchical clustering of samples based on EPIC SNP probes",
  ylab = "1 - Pearson correlation"
)
dev.off()

## 6. Save heatmap
ggsave(
  file.path(plots_dir, "SNP_sample_correlation_heatmap.pdf"),
  plot = p_heat,
  width = 7,
  height = 6
)
ggsave(
  file.path(plots_dir, "SNP_sample_correlation_heatmap.png"),
  plot = p_heat,
  width = 7,
  height = 6,
  dpi = 300
)

message("Chunk 9b completed: SNP correlation heatmap and dendrogram saved in 'plots/'.")

library(data.table)

if (!exists("snp_beta_mat") || is.null(snp_beta_mat)) {
  stop("snp_beta_mat is not available; run chunk 9 before 9c.")
}

if (ncol(snp_beta_mat) < 2 || nrow(snp_beta_mat) < 10) {
  stop("snp_beta_mat has too few samples or SNP probes to estimate relatedness.")
}

cat("\n=== SNP-based relatedness (pairwise correlations) ===\n")

# 1. Sample–sample Pearson correlation on SNP betas
cor_mat <- suppressWarnings(
  stats::cor(snp_beta_mat, use = "pairwise.complete.obs", method = "pearson")
)

# 2. Convert upper triangle to long format table
sample_ids <- colnames(snp_beta_mat)
pair_list <- which(upper.tri(cor_mat), arr.ind = TRUE)

related_dt <- data.table(
  Sample1     = sample_ids[pair_list[, 1]],
  Sample2     = sample_ids[pair_list[, 2]],
  SNP_corr    = cor_mat[upper.tri(cor_mat)]
)

# 3. Add cohort/region annotations if available
if (exists("ss_qc_full")) {
  annot <- ss_qc_full[, .(Sample_ID)]
  related_dt <- merge(
    related_dt, annot,
    by.x = "Sample1", by.y = "Sample_ID", all.x = TRUE
  )

  related_dt <- merge(
    related_dt, annot,
    by.x = "Sample2", by.y = "Sample_ID", all.x = TRUE
  )
}

# 4. Flag highly similar pairs (e.g. correlation >= 0.95)
high_thresh <- 0.95
related_dt[, high_related := SNP_corr >= high_thresh]

cat("Number of sample pairs with SNP_corr >=", high_thresh, ": ",
    sum(related_dt$high_related, na.rm = TRUE), "\n\n", sep = "")

if (sum(related_dt$high_related, na.rm = TRUE) > 0) {
  cat("Pairs with high SNP-based similarity (potential duplicates/relatives):\n")
  print(related_dt[high_related == TRUE][order(-SNP_corr)])
} else {
  cat("No sample pairs exceed the high-relatedness threshold.\n")
}

# 5. Save full relatedness table (optional but recommended)
fwrite(related_dt,
       file = "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_SNP_relatedness_table.csv")

cat("\nSNP-based relatedness table written to EPICv2_sesame_SNP_relatedness_table.csv\n")

library(data.table)

keep_ids_final <- ss_qc_full[qc_fail == FALSE & status == "OK", Sample_ID]
keep_ids_final <- intersect(as.character(keep_ids_final), names(results_list))

message("Number of QC-passed samples requested for final matrix: ", length(keep_ids_final))

if (length(keep_ids_final) == 0) {
  stop("No QC-passed samples available to build the final masked beta matrix.")
}

## 1. Collect betas_masked for each final sample, dropping any that are missing/empty
beta_list_final <- lapply(keep_ids_final, function(sid) {
  x <- results_list[[sid]]
  if (is.null(x) || is.null(x$betas_masked)) return(NULL)
  b <- x$betas_masked
  if (length(b) == 0 || is.null(names(b))) return(NULL)
  b
})
names(beta_list_final) <- keep_ids_final
beta_list_final <- beta_list_final[!vapply(beta_list_final, is.null, logical(1))]

message("Samples with non-null betas_masked in final set: ", length(beta_list_final))

if (length(beta_list_final) == 0) {
  stop("All QC-passed samples have NULL or empty betas_masked; cannot build final matrix.")
}

## 2. Intersect probes across the remaining samples
common_probes <- Reduce(intersect, lapply(beta_list_final, names))

message("Number of common probes across final samples: ", length(common_probes))

if (length(common_probes) == 0) {
  stop("No common probes across final samples; cannot build final matrix.")
}

## 3. Build the masked beta matrix
beta_mat_masked <- sapply(names(beta_list_final), function(sid) {
  as.numeric(beta_list_final[[sid]][common_probes])
})

rownames(beta_mat_masked) <- common_probes
colnames(beta_mat_masked) <- names(beta_list_final)

## 4. Save RDS and CSV
saveRDS(beta_mat_masked, file = "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_masked_beta_matrix.rds")

fwrite(
  data.table(Probe_ID = rownames(beta_mat_masked), beta_mat_masked),
  "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_masked_beta_matrix.csv"
)

## ------------ save all samples ------------
keep_ids_final_all <- ss_qc_full[, Sample_ID]
keep_ids_final_all <- intersect(as.character(keep_ids_final_all), names(results_list))

message("Number of all samples for final matrix: ", length(keep_ids_final_all))

## 1. Collect betas_masked for each final sample, dropping any that are missing/empty
beta_list_final_all <- lapply(keep_ids_final_all, function(sid) {
  x <- results_list[[sid]]
  if (is.null(x) || is.null(x$betas_masked)) return(NULL)
  b <- x$betas_masked
  if (length(b) == 0 || is.null(names(b))) return(NULL)
  b
})
names(beta_list_final_all) <- keep_ids_final_all
beta_list_final_all <- beta_list_final_all[!vapply(beta_list_final_all, is.null, logical(1))]

message("Samples with non-null betas_masked in final set: ", length(beta_list_final_all))

## 2. Intersect probes across the remaining samples
common_probes <- Reduce(intersect, lapply(beta_list_final_all, names))

message("Number of common probes across final samples: ", length(common_probes))

if (length(common_probes) == 0) {
  stop("No common probes across final samples; cannot build final matrix.")
}

## 3. Build the masked beta matrix
beta_mat_masked_all <- sapply(names(beta_list_final_all), function(sid) {
  as.numeric(beta_list_final_all[[sid]][common_probes])
})

rownames(beta_mat_masked_all) <- common_probes
colnames(beta_mat_masked_all) <- names(beta_list_final_all)

fwrite(
  data.table(Probe_ID = rownames(beta_mat_masked_all), beta_mat_masked_all),
  "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_masked_beta_matrix_all.csv"
)
## ------------------------

message(
  "Final masked beta matrix written: ",
  nrow(beta_mat_masked), " probes × ", ncol(beta_mat_masked), " samples.\n",
  "Files: EPICv2_sesame_masked_beta_matrix.rds and .csv"
)

library(data.table)

bad_probe_dir <- "C:/Git/bioTest/r_olga_script_test/probe_filters"
if (!dir.exists(bad_probe_dir)) dir.create(bad_probe_dir)

library(data.table)

mask_tsv_url <- "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.mask.tsv.gz"
local_tsv_gz <- "EPICv2.hg38.mask.tsv.gz"

if (!file.exists(local_tsv_gz)) {
  message("Скачиваем маску EPICv2...")
  tryCatch(
    download.file(mask_tsv_url, destfile = local_tsv_gz, mode = "wb", method = "auto", quiet = FALSE),
    error = function(e) {
      download.file(mask_tsv_url, destfile = local_tsv_gz, mode = "wb", method = "wininet", quiet = FALSE)
    }
  )
  if (!file.exists(local_tsv_gz) || file.size(local_tsv_gz) == 0) {
    stop("Не удалось скачать файл. Скачайте вручную: ", mask_tsv_url)
  }
}

mask_dt <- fread(local_tsv_gz, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

print(colnames(mask_dt))

if (is.logical(mask_dt$M_general)) {
  bad_ids <- mask_dt$Probe_ID[mask_dt$M_general == TRUE]
} else if (is.numeric(mask_dt$M_general)) {
  bad_ids <- mask_dt$Probe_ID[mask_dt$M_general == 1]
} else {
  bad_ids <- mask_dt$Probe_ID[as.logical(mask_dt$M_general)]
}

bad_probe_dt <- data.table(
  Probe_ID = bad_ids,
  bad_probe_reason = "masked_general"
)

bad_probe_ids_dt <- unique(bad_probe_dt[, .(Probe_ID)])

fwrite(bad_probe_dt, file.path(bad_probe_dir, "EPICv2_known_bad_probes_annotated.csv"))
fwrite(bad_probe_ids_dt, file.path(bad_probe_dir, "EPICv2_known_bad_probes_ids.csv"))

cat("\nKnown bad EPIC probe list created.\n")
cat("Annotated file: probe_filters/EPICv2_known_bad_probes_annotated.csv\n")
cat("Simple ID file: probe_filters/EPICv2_known_bad_probes_ids.csv\n")
cat("Number of unique known bad probes: ", nrow(bad_probe_ids_dt), "\n", sep = "")

library(data.table)

## Load final matrix if not already in memory
if (!exists("beta_mat_masked") || is.null(beta_mat_masked)) {
  if (!file.exists("C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_masked_beta_matrix.rds")) {
    stop("EPICv2_sesame_masked_beta_matrix.rds not found.")
  }
  beta_mat_masked <- readRDS("C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_masked_beta_matrix.rds")
}

if (is.null(rownames(beta_mat_masked))) {
  stop("beta_mat_masked must have probe IDs as rownames.")
}

## Load known bad probes
bad_probe_file <- file.path("C:/Git/bioTest/r_olga_script_test/probe_filters", "EPICv2_known_bad_probes_ids.csv")
if (!file.exists(bad_probe_file)) {
  stop("Known bad probe ID file not found: ", bad_probe_file)
}

bad_probe_ids <- fread(bad_probe_file)$Probe_ID
bad_probe_ids <- unique(as.character(bad_probe_ids))

## 1. Restrict to CpG probes only (drop rs and ch probes if present)
probe_ids <- rownames(beta_mat_masked)
is_cpg <- grepl("^cg", probe_ids, ignore.case = FALSE)

## 2. Missingness filter across samples
probe_missing_frac <- rowMeans(is.na(beta_mat_masked))
missing_thresh <- 0.05   # 5% missingness threshold

fail_missing <- probe_missing_frac > missing_thresh

## 3. Known bad probe filter
fail_known_bad <- probe_ids %in% bad_probe_ids

## 4. Final keep/drop flags
probe_qc_dt <- data.table(
  Probe_ID = probe_ids,
  is_cpg = is_cpg,
  missing_frac = probe_missing_frac,
  fail_missing = fail_missing,
  fail_known_bad = fail_known_bad
)

probe_qc_dt[, keep_probe := is_cpg & !fail_missing & !fail_known_bad]

## 5. Filter matrix
keep_probe_ids <- probe_qc_dt[keep_probe == TRUE, Probe_ID]
beta_mat_probe_filtered <- beta_mat_masked[keep_probe_ids, , drop = FALSE]

## 6. Save outputs
saveRDS(beta_mat_probe_filtered, "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_masked_beta_matrix_probe_filtered.rds")

fwrite(
  data.table(Probe_ID = rownames(beta_mat_probe_filtered), beta_mat_probe_filtered),
  "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_masked_beta_matrix_probe_filtered.csv"
)

fwrite(probe_qc_dt, "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_probe_QC_table.csv")

## 7. Report
cat("\n=== PROBE-LEVEL QC SUMMARY ===\n")
cat("Input probes:                         ", nrow(beta_mat_masked), "\n")
cat("CpG probes retained for consideration:", sum(probe_qc_dt$is_cpg), "\n")
cat("Removed for >5% missingness:          ", sum(probe_qc_dt$fail_missing), "\n")
cat("Removed as known bad probes:          ", sum(probe_qc_dt$fail_known_bad), "\n")
cat("Final probes retained:                ", nrow(beta_mat_probe_filtered), "\n")
cat("Files written:\n")
cat("  - EPICv2_sesame_masked_beta_matrix_probe_filtered.rds\n")
cat("  - EPICv2_sesame_masked_beta_matrix_probe_filtered.csv\n")
cat("  - EPICv2_sesame_probe_QC_table.csv\n")

library(data.table)

if (!exists("probe_qc_dt") || !exists("beta_mat_masked") || !exists("beta_mat_probe_filtered")) {
  cat("\n[probe_QC_overview] Required objects not found in memory.\n")
  cat("Run the probe-level QC chunk (10c) before this overview.\n\n")
} else {
  setDT(probe_qc_dt)

  n_input_all    <- nrow(beta_mat_masked)
  n_cpg          <- sum(probe_qc_dt$is_cpg)
  n_fail_miss    <- sum(probe_qc_dt$fail_missing)
  n_fail_bad     <- sum(probe_qc_dt$fail_known_bad)
  n_fail_both    <- sum(probe_qc_dt$fail_missing & probe_qc_dt$fail_known_bad)
  n_kept         <- nrow(beta_mat_probe_filtered)

  probe_QC_overview <- data.table(
    metric = c(
      "input_probes_all_types",
      "cpg_probes_considered",
      "fail_missing_only",
      "fail_known_bad_only",
      "fail_both_missing_and_known_bad",
      "final_probes_retained"
    ),
    count = c(
      n_input_all,
      n_cpg,
      n_fail_miss - n_fail_both,
      n_fail_bad - n_fail_both,
      n_fail_both,
      n_kept
    )
  )

  print(probe_QC_overview)

  fwrite(
    probe_QC_overview,
    file = "C:/Git/bioTest/r_olga_script_test/EPICv2_sesame_probe_QC_overview.csv"
  )

  cat("\n[probe_QC_overview] Summary written to EPICv2_sesame_probe_QC_overview.csv\n\n")
}

library(data.table)
setDT(ss_qc_full)

cat("========== DNAm QC REPORT ==========\n\n")
cat("Platform: Illumina MethylationEPIC v2.0; processing via SeSAMe;",
    "known-bad EPIC probes removed using Pidsley et al. supplementary filters.\n\n")

## 1. SAMPLE COUNTS
n_total <- nrow(ss_qc_full)
n_pass  <- sum(ss_qc_full$qc_fail == FALSE, na.rm = TRUE)
n_fail  <- sum(ss_qc_full$qc_fail == TRUE,  na.rm = TRUE)

cat("1) SAMPLE COUNTS\n")
cat("   Total samples in sample sheet:        ", n_total,
    "  (all samples with metadata and IDATs)\n", sep = "")
cat("   Samples passing main QC:              ", n_pass,
    "  (qc_fail == FALSE)\n", sep = "")
cat("   Samples failing main QC:              ", n_fail,
    "  (qc_fail == TRUE; excluded from main analyses)\n\n", sep = "")

## 2. REASONS FOR SAMPLE EXCLUSION
fail_dt <- ss_qc_full[qc_fail == TRUE]

cat("2) REASONS FOR SAMPLE EXCLUSION (among ", n_fail, " QC-failed samples)\n", sep = "")

if (nrow(fail_dt) == 0) {
  cat("   No samples flagged as QC failures.\n\n")
} else {
  n_fail_frac_hard <- sum(fail_dt$fail_frac_hard, na.rm = TRUE)
  n_fail_int       <- sum(fail_dt$fail_int,       na.rm = TRUE)
  n_fail_bs        <- sum(fail_dt$fail_bs,        na.rm = TRUE)
  n_fail_status    <- sum(fail_dt$fail_status,    na.rm = TRUE)
  n_fail_sex       <- sum(fail_dt$fail_sex,       na.rm = TRUE)

  cat("   Hard fraction-failed filter (fail_frac_hard): ", n_fail_frac_hard,
      "  [samples with frac_failed > hard cutoff]\n", sep = "")
  cat("   Low mean intensity (fail_int):               ", n_fail_int,
      "  [poor overall signal]\n", sep = "")
  cat("   High bisulfite GCT score (fail_bs):          ", n_fail_bs,
      "  [potential bisulfite conversion problems]\n", sep = "")
  cat("   Processing error (fail_status):              ", n_fail_status,
      "  [SeSAMe error or missing objects]\n", sep = "")
  cat("   Sex mismatch (fail_sex):                     ", n_fail_sex,
      "  [reported vs inferred sex disagrees]\n\n", sep = "")

  cat("   IDs of QC-failed samples (excluded):\n")
  cat("   ", paste(fail_dt$Sample_ID, collapse = ", "), "\n\n", sep = "")
}

## 3. SOFT WARNINGS AMONG QC-PASSED SAMPLES
pass_dt <- ss_qc_full[qc_fail == FALSE]

cat("3) SOFT WARNINGS AMONG QC-PASSED SAMPLES (", nrow(pass_dt),
    " samples kept but flagged for review)\n", sep = "")

if (nrow(pass_dt) == 0) {
  cat("   No QC-passed samples to inspect.\n\n")
} else {
  n_warn_frac <- sum(pass_dt$warn_frac, na.rm = TRUE)

  has_oob <- all(c("oob_outlier", "oob_outlier_strict") %in% names(ss_qc_full))
  if (!has_oob) {
    pass_dt[, `:=`(oob_outlier = FALSE, oob_outlier_strict = FALSE)]
  }

  n_oob_any    <- sum(pass_dt$oob_outlier,        na.rm = TRUE)
  n_oob_strict <- sum(pass_dt$oob_outlier_strict, na.rm = TRUE)

  cat("   warn_frac == TRUE:               ", n_warn_frac,
      "  [frac_failed above soft threshold but below hard]\n", sep = "")
  cat("   oob_outlier == TRUE:            ", n_oob_any,
      "  [OOB background above 95th percentile in at least one channel]\n", sep = "")
  cat("   oob_outlier_strict == TRUE:     ", n_oob_strict,
      "  [OOB background above 95th percentile in both channels]\n\n", sep = "")

  if (n_warn_frac > 0) {
    cat("   Samples with elevated frac_failed (warn_frac == TRUE, but still kept):\n")
    cat("   ", paste(pass_dt[warn_frac == TRUE, Sample_ID], collapse = ", "),
        "\n\n", sep = "")
  }

  if (n_oob_any > 0) {
    cat("   Samples with high OOB background (oob_outlier == TRUE; watch in downstream analyses):\n")
    cat("   ", paste(pass_dt[oob_outlier == TRUE, Sample_ID], collapse = ", "),
        "\n\n", sep = "")
  }
}

## 4. SNP PROBE COVERAGE
cat("4) SNP PROBE COVERAGE (59 EPIC SNP probes)\n")
if (exists("snp_qc_dt")) {
  setDT(snp_qc_dt)
  snp_min  <- min(snp_qc_dt$n_snp_probes, na.rm = TRUE)
  snp_med  <- median(snp_qc_dt$n_snp_probes, na.rm = TRUE)
  snp_mean <- mean(snp_qc_dt$n_snp_probes, na.rm = TRUE)
  snp_max  <- max(snp_qc_dt$n_snp_probes, na.rm = TRUE)

  cat("   n_snp_probes summary: min = ", snp_min,
      ", median = ", snp_med,
      ", mean = ", round(snp_mean, 2),
      ", max = ", snp_max,
      "  [ideal is 59 for all samples]\n", sep = "")

  low_snp_ids <- snp_qc_dt[n_snp_probes < 59, Sample_ID]
  if (length(low_snp_ids) == 0) {
    cat("   All samples have full (59-probe) SNP coverage.\n\n")
  } else {
    cat("   Samples with reduced SNP coverage (n_snp_probes < 59; identity QC less robust):\n")
    cat("   ", paste(low_snp_ids, collapse = ", "), "\n\n", sep = "")
  }
} else {
  cat("   SNP QC table not found in environment; skipping SNP coverage summary.\n\n")
}

## 5. PROBE-LEVEL QC AND MASKING SUMMARY
cat("5) PROBE-LEVEL QC AND MASKING SUMMARY\n")

if (!("num_na_cg" %in% names(ss_qc_full))) {
  cat("   Probe-level NA metrics (num_na_cg / frac_na_cg) not found; skipping.\n\n")
} else {
  na_cg_min  <- min(ss_qc_full$num_na_cg, na.rm = TRUE)
  na_cg_med  <- median(ss_qc_full$num_na_cg, na.rm = TRUE)
  na_cg_mean <- mean(ss_qc_full$num_na_cg, na.rm = TRUE)
  na_cg_max  <- max(ss_qc_full$num_na_cg, na.rm = TRUE)

  cat("   Number of CpGs set to NA per sample (num_na_cg):\n")
  cat("     min = ", na_cg_min,
      ", median = ", na_cg_med,
      ", mean = ", round(na_cg_mean, 1),
      ", max = ", na_cg_max,
      "  [higher values indicate more probe-level masking]\n", sep = "")

  na_cg_p95 <- quantile(ss_qc_full$frac_na_cg, 0.95, na.rm = TRUE)
  high_na_ids <- ss_qc_full[frac_na_cg > na_cg_p95 & !is.na(frac_na_cg), Sample_ID]

  cat("   95th percentile of frac_na_cg across samples: ",
      signif(na_cg_p95, 3),
      "  [fraction of CpGs masked/NA]\n", sep = "")

  if (length(high_na_ids) == 0) {
    cat("   No samples with exceptionally high masked-CpG fraction.\n\n")
  } else {
    cat("   Samples with high masked-CpG fraction (frac_na_cg > 95th percentile; inspect if kept):\n")
    cat("   ", paste(high_na_ids, collapse = ", "), "\n\n", sep = "")
  }

  if (all(c("InfI_switch_G2R", "InfI_switch_R2G") %in% names(ss_qc_full))) {
    infi_sw_mean <- mean(ss_qc_full$InfI_switch_G2R + ss_qc_full$InfI_switch_R2G, na.rm = TRUE)
    infi_sw_max  <- max(ss_qc_full$InfI_switch_G2R + ss_qc_full$InfI_switch_R2G, na.rm = TRUE)
    cat("   Type-I color-channel switch counts (InfI_switch_G2R + InfI_switch_R2G):\n")
    cat("     mean total switches per sample ~ ", round(infi_sw_mean, 1),
        ", max = ", infi_sw_max,
        "  [very high values could indicate technical artifacts]\n\n", sep = "")
  } else {
    cat("   Type-I switch metrics not found; skipping.\n\n")
  }
}

## 6. FINAL RECOMMENDED SAMPLE SET
keep_ids_final <- ss_qc_full[qc_fail == FALSE & status == "OK", Sample_ID]

cat("6) FINAL RECOMMENDED SAMPLE SET FOR DOWNSTREAM ANALYSES\n")
cat("   Definition: qc_fail == FALSE & status == 'OK'.\n")
cat("   Number of samples retained: ", length(keep_ids_final), "\n", sep = "")
cat("   Example retained IDs: ",
    paste(head(keep_ids_final, 10), collapse = ", "), "\n", sep = "")
cat("   (Optionally exclude 'watch list' samples with high OOB background\n",
    "    or high frac_na_cg in sensitivity analyses.)\n\n", sep = "")

## 7. SNP-BASED RELATEDNESS SUMMARY
cat("7) SNP-BASED RELATEDNESS (pairwise similarity from SNP betas)\n")
if (exists("related_dt")) {
  setDT(related_dt)
  high_thresh <- 0.95
  n_pairs_high <- sum(related_dt$SNP_corr >= high_thresh, na.rm = TRUE)

  cat("   Threshold for high relatedness: SNP_corr >=", high_thresh, "\n")
  cat("   Number of sample pairs above threshold: ", n_pairs_high, "\n", sep = "")

  if (n_pairs_high > 0) {
    top_pairs <- related_dt[SNP_corr >= high_thresh][order(-SNP_corr)]
    cat("   Highly related pairs (potential duplicates/relatives):\n")
    apply(head(top_pairs, 10), 1, function(row) {
      cat("     ", row[["Sample1"]], " – ", row[["Sample2"]],
          "  (SNP_corr = ", sprintf("%.3f", as.numeric(row[["SNP_corr"]])), ")\n", sep = "")
    })
    cat("\n")
  } else {
    cat("   No sample pairs exceed the high-relatedness threshold.\n\n")
  }
} else {
  cat("   Relatedness table (related_dt) not found; run chunk 9c to compute it.\n\n")
}

## 8) PROBE-LEVEL FILTERING SUMMARY (EPIC array)

cat("8) PROBE-LEVEL FILTERING SUMMARY (EPIC array)\n")

if (!exists("probe_qc_dt") || !exists("beta_mat_probe_filtered")) {
  cat("   Probe-level QC objects not found in memory.\n")
  cat("   (Run chunk 10c before generating this summary.)\n\n")
} else {
  n_input      <- nrow(beta_mat_masked)
  n_cpg        <- sum(probe_qc_dt$is_cpg)
  n_missing    <- sum(probe_qc_dt$fail_missing)
  n_known_bad  <- sum(probe_qc_dt$fail_known_bad)
  n_final      <- nrow(beta_mat_probe_filtered)

  cat("   Input probes (all types):                ", n_input,
      "  [rows in EPICv2_sesame_masked_beta_matrix]\n", sep = "")
  cat("   CpG probes considered:                   ", n_cpg,
      "  [probes with IDs starting 'cg']\n", sep = "")
  cat("   Removed for >5% missingness:             ", n_missing,
      "  [dataset-specific low-quality probes]\n", sep = "")
  cat("   Removed as known bad probes:             ", n_known_bad,
      "  [cross-reactive / SNP-affected / off-target]\n", sep = "")
  cat("   Final probes retained after filtering:   ", n_final,
      "  [used in EPICv2_sesame_masked_beta_matrix_probe_filtered]\n\n", sep = "")
}

cat("Full probe-level QC metrics are provided in: ",
    "EPICv2_sesame_probe_QC_table.csv\n\n")

cat("========== END OF DNAm QC REPORT ==========\n\n")

## 12. Final cleanup: drop auxiliary objects not needed downstream

# 1) Define which objects we want to KEEP in memory
keep_objects <- c(
  "ss_qc_full",                    # full QC table
  "beta_mat_masked",               # final masked beta matrix (if created in this session)
  "snp_beta_mat",                  # SNP beta matrix (optional, for interactive checks)
  "snp_qc_dt",                     # SNP coverage summary (optional)
  "related_dt"                     # SNP-relatedness table (optional)
)

# 2) Explicitly drop heavy components inside results_list
if (exists("results_list")) {
  for (sid in names(results_list)) {
    if (!is.null(results_list[[sid]]$sdf_raw)) {
      results_list[[sid]]$sdf_raw <- NULL
    }
    if (!is.null(results_list[[sid]]$qc_obj)) {
      results_list[[sid]]$qc_obj <- NULL
    }
    if (!is.null(results_list[[sid]]$betas_unmasked)) {
      results_list[[sid]]$betas_unmasked <- NULL
    }
    # betas_masked is already saved to disk; drop if you don't need it in memory:
    if (!is.null(results_list[[sid]]$betas_masked)) {
      results_list[[sid]]$betas_masked <- NULL
    }
  }
}

# 3) Remove obviously intermediate objects from the global environment
drop_candidates <- c(
  "ss", "ss_pref", "qc_dt", "qc_extra",
  "qc_overall", 
  "beta_long", "beta_long_list", "beta_list_final", "beta_list_ok",
  "sdf_list_ok", "qc_obj_list",
  "cor_mat", "cor_dt", "dist_mat", "hc", "sample_order",
  "p_all", "p_region", "p_example",
  "p_frac_failed", "p_mean_intensity", "p_bs",
  "p_rg_ratio", "p_rg_distort", "p_mean_beta", "p_median_beta", "p_oob",
  "oob_grn_p95", "oob_red_p95",
  "keep_ids", "keep_ids_qc", "keep_ids_final"
)

for (obj in drop_candidates) {
  if (exists(obj, envir = .GlobalEnv)) {
    rm(list = obj, envir = .GlobalEnv)
  }
}

# 4) Optionally drop results_list entirely if you don’t need it anymore
# (All key outputs are already saved to disk.)
if (exists("results_list")) {
  rm(results_list)
}

# 5) Final garbage collection
invisible(gc())

cat("\nFinal cleanup done. Kept objects in memory:\n")
cat("  ", paste(keep_objects[keep_objects %in% ls()], collapse = ", "), "\n\n", sep = "")
