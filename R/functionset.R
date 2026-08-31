
reverse <- function(x) {
  y <- switch(x,
              "G" = "C",
              "C" = "G",
              "T" = "A",
              "A" = "T",
              x)
  return(y)
}


matchAellel <- function(al_gwas, al_gtex) {
  p <- nrow(al_gwas)
  sign_factor <- rep(0, p)

  # matched major and minor allele
  case1 <- which((al_gwas[, 1] == al_gtex[, 1]) & (al_gwas[, 2] == al_gtex[, 2]))
  sign_factor[case1] <- 1

  # contrasted major and minor allele
  case2 <- which((al_gwas[, 1] == al_gtex[, 2]) & (al_gwas[, 2] == al_gtex[, 1]))
  sign_factor[case2] <- -1

  # the mis matched snps
  mis = (1:p)[-c(case1, case2)]
  for(i in mis)
  {
    ###############
    ## scenario 1 #
    ##       A1 A2#
    ## file1 T, G # | file1 A, G # | file1 A, G # | file1 A, G # | file1 A, C #
    ## file2 T, C # | file2 T, G # | file2 A, C # | file2 T, C # | file2 T, C #
    ###############
    if ((al_gwas[i, 1] == al_gtex[i, 1] | reverse(al_gwas[i, 1]) == al_gtex[i, 1]) &
        (al_gwas[i, 2] == al_gtex[i, 2] | reverse(al_gwas[i, 2]) == al_gtex[i, 2])) {
      sign_factor[i] <- 1
      ###############
      ## scenario 2 #
      ##       A1 A2#
      ## file1 T, G # | file1 A, G # | file1 A, G # | file1 A, G # | file1 A, C #
      ## file2 C, T # | file2 G, T # | file2 C, A # | file2 C, T # | file2 C, T #
      ###############
    } else if((al_gwas[i, 1] == al_gtex[i, 2] | reverse(al_gwas[i, 1]) == al_gtex[i, 2]) &
              (al_gwas[i, 2] == al_gtex[i, 1] | reverse(al_gwas[i, 2]) == al_gtex[i, 1])) {
      sign_factor[i] <- -1
    } else {
      print(cat(i,al_gwas[i, 1],al_gwas[i, 2],al_gtex[i, 1],al_gtex[i, 2]))
    }
  }
  return(sign_factor)
}


matchpanel <- function(ss, stringname3){

  message("Starting processing of file: ", basename(ss))
  start_time <- Sys.time()

  dir_path <- dirname(ss)
  file_name <- basename(ss)

  message("Step 1/6: Reading reference panel (bim file)...")
  bim <- data.table::fread(paste0(stringname3,".bim"), header = FALSE)
  bim <- data.frame(bim)

  message("Step 2/6: Reading GWAS summary statistics...")
  ss.gwas <- readr::read_delim(ss)
  ss.gwas <- data.frame(ss.gwas)
  cat("\n")

  message("Step 3/6: Preprocessing GWAS data...")
  ss.gwas <- na.omit(ss.gwas)
  ss.gwas <- ss.gwas[nchar(ss.gwas$A1)==1,]
  ss.gwas <- ss.gwas[nchar(ss.gwas$A2)==1,]
  ss.gwas$A1 <- toupper(ss.gwas$A1)
  ss.gwas$A2 <- toupper(ss.gwas$A2)
  ss.gwas <- ss.gwas[!duplicated(ss.gwas$SNP), ]
  ss.gwas <- ss.gwas[order(ss.gwas$SNP), ]
  ss.gwas <- ss.gwas[grep("^rs", ss.gwas$SNP), ]

  message("Step 4/6: Merging with reference panel...")
  ss.merge <- merge(ss.gwas, bim, by.x=c('SNP'), by.y=c('V2'))
  allele_trait <- cbind(ss.merge$A1, ss.merge$A2)
  allele_bim   <- cbind(ss.merge$V5, ss.merge$V6)

  message("Step 5/6: Matching allele directions...")
  sign_factor <- matchAellel(allele_trait, allele_bim)
  if (sum(sign_factor == 0) > 0) {
    message("Removing ", sum(sign_factor == 0), " SNPs with allele mismatch")
    ss.merge <- ss.merge[-which(sign_factor == 0),]
    sign_factor <- sign_factor[-which(sign_factor == 0)]
  }

  ss.match <- data.frame(SNP = ss.merge$SNP,
                         CHR = ss.merge$V1,
                         BP = ss.merge$V4,
                         A1 = ss.merge$V5,
                         A2 = ss.merge$V6,
                         BETA = ss.merge$BETA*sign_factor,
                         SE = ss.merge$SE,
                         P = ss.merge$P)

  new_file_name <- paste0("match.", file_name)
  new_file_path <- file.path(dir_path, new_file_name)

  message("Step 6/6: Saving matched data...")
  data.table::fwrite(ss.match, file = new_file_path, quote = FALSE, row.names = FALSE, sep = "\t")

  end_time <- Sys.time()
  duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 1)
  message("Successfully processed ", nrow(ss.match), " SNPs")
  message("Total processing time: ", duration, " seconds")
  message("Output saved to: ", new_file_path)

  return(list(data = ss.match, data_dir = new_file_path))
}


ivselect <- function(expgwas_dir, expgwis_dir = NULL, 
                     outgwas_dir, outgwis_dir = NULL,
                     stringname3, block_file, plink_dir = NULL,
                     pval_cutoff_gwas = 0.00000005, pval_cutoff_gwis = 0.00000005,
                     r2_cutoff = 0.01, kb_cutoff = 1024, maf_cutoff = 0.05,
                     lam = 0.1, coreNum = 1, intersect_mode = FALSE){
  
  start_time <- Sys.time()
  message("Starting MERLIN IV selection process...")
  
  # Step 1: Read Reference Panel
  message("Step 1/7: Reading reference panel (bim file)...")
  bim <- data.table::fread(paste0(stringname3,".bim"), header = FALSE)
  bim <- data.frame(bim)
  
  # Step 2: Read Summary Statistics Conditionally
  message("Step 2/7: Reading GWAS summary statistics...")
  
  expgwas <- data.frame(readr::read_delim(expgwas_dir, col_names=TRUE, show_col_types = FALSE))
  outgwas <- data.frame(readr::read_delim(outgwas_dir, col_names=TRUE, show_col_types = FALSE))
  
  expgwis <- NULL
  if (!is.null(expgwis_dir)) {
    message("  - Reading exposure GWIS...")
    expgwis <- data.frame(readr::read_delim(expgwis_dir, col_names=TRUE, show_col_types = FALSE))
  }
  
  outgwis <- NULL
  if (!is.null(outgwis_dir)) {
    message("  - Reading outcome GWIS...")
    outgwis <- data.frame(readr::read_delim(outgwis_dir, col_names=TRUE, show_col_types = FALSE))
  }
  
  # Step 3: PLINK Setup
  message("Step 3/7: Setting up PLINK...")
  if (is.null(plink_dir)) {
    plink_dir <- bigsnpr::download_plink()
  }
  message("  ✓ Using PLINK: ", plink_dir)
  
  # Step 4: LD Clumping for Exposure GWAS
  message("Step 4/7: Running LD clumping for exposure GWAS...")
  expgwas.LD.cmd <- paste(plink_dir, " --bfile ", stringname3,
                          " --clump-p1 ", pval_cutoff_gwas,
                          " --clump-r2 ", r2_cutoff,
                          " --clump-kb ", kb_cutoff,
                          " --maf ", maf_cutoff,
                          " --clump ", expgwas_dir,
                          " --clump-snp-field SNP --clump-field P --out ld_result1",
                          sep="")
  system(expgwas.LD.cmd, ignore.stdout = TRUE)
  
  if (!file.exists("ld_result1.clumped")) stop("PLINK clumping failed for exposure GWAS.")
  expgwas.LD <- data.table::fread("ld_result1.clumped")
  snp.expgwas.LD <- expgwas.LD$SNP
  file.remove("ld_result1.clumped", "ld_result1.log")
  
  # Step 5 & 6: LD Clumping for GWIS & Causal SNP Selection
  if (!is.null(expgwis_dir)) {
    message("Step 5/7: Running LD clumping for exposure GWIS...")
    expgwis.LD.cmd <- paste(plink_dir, " --bfile ", stringname3,
                            " --clump-p1 ", pval_cutoff_gwis,
                            " --clump-r2 ", r2_cutoff,
                            " --clump-kb ", kb_cutoff,
                            " --maf ", maf_cutoff,
                            " --clump ", expgwis_dir,
                            " --clump-snp-field SNP --clump-field P --out ld_result2",
                            sep="")
    system(expgwis.LD.cmd, ignore.stdout = TRUE)
    
    if (!file.exists("ld_result2.clumped")) stop("PLINK clumping failed for exposure GWIS.")
    expgwis.LD <- data.table::fread("ld_result2.clumped")
    snp.expgwis.LD <- expgwis.LD$SNP
    file.remove("ld_result2.clumped", "ld_result2.log")
    
    message("Step 6/7: Selecting causal SNPs...")
    message("  - Mode: ", ifelse(intersect_mode, "Intersection", "Union"))
    
    if (intersect_mode == TRUE) {
      snp.causal <- intersect(snp.expgwas.LD, snp.expgwis.LD)
    } else {
      union.snp <- union(snp.expgwas.LD, snp.expgwis.LD)
      union.snp.df <- data.frame(SNP = union.snp, P = 1)
      union.dir <- "union.snp.txt"
      write.table(union.snp.df, file = union.dir, quote = FALSE, row.names = FALSE)
      
      union.LD.cmd <- paste(plink_dir, " --bfile ", stringname3,
                            " --clump-p1 ", 1,
                            " --clump-r2 ", r2_cutoff,
                            " --clump-kb ", kb_cutoff,
                            " --maf ", maf_cutoff,
                            " --clump ", union.dir,
                            " --clump-snp-field SNP --clump-field P --out ld_result3",
                            sep="")
      system(union.LD.cmd, ignore.stdout = TRUE)
      union.LD <- data.table::fread("ld_result3.clumped")
      snp.causal <- union.LD$SNP
      file.remove(union.dir, "ld_result3.clumped", "ld_result3.log")
    }
  } else {
    message("Step 5 & 6/7: Exposure GWIS not provided. Using GWAS clumped SNPs as causal SNPs...")
    snp.causal <- snp.expgwas.LD
  }
  
  # Ensure causal SNPs exist in all provided datasets
  snp.causal <- intersect(snp.causal, expgwas$SNP)
  snp.causal <- intersect(snp.causal, outgwas$SNP)
  if (!is.null(expgwis)) snp.causal <- intersect(snp.causal, expgwis$SNP)
  if (!is.null(outgwis)) snp.causal <- intersect(snp.causal, outgwis$SNP)
  
  if (length(snp.causal) == 0) {
    stop("No SNPs remaining after intersecting across all provided datasets.")
  } else {
    message("  ✓ Final number of causal SNPs selected: ", length(snp.causal))
  }
  
  # Match with BIM file
  avbIndex <- match(snp.causal, bim$V2)
  avbIndex <- as.matrix(avbIndex[order(avbIndex)])
  snp.causal <- bim[avbIndex, ]$V2
  
  # Align all datasets
  expgwas.order <- expgwas[match(snp.causal, expgwas$SNP), ]
  outgwas.order <- outgwas[match(snp.causal, outgwas$SNP), ]
  
  bh11ld <- expgwas.order$BETA
  se11ld <- expgwas.order$SE
  bh21ld <- outgwas.order$BETA
  se21ld <- outgwas.order$SE
  
  # Extract conditionally
  bh12ld <- if (!is.null(expgwis)) expgwis[match(snp.causal, expgwis$SNP), ]$BETA else NULL
  se12ld <- if (!is.null(expgwis)) expgwis[match(snp.causal, expgwis$SNP), ]$SE else NULL
  
  bh22ld <- if (!is.null(outgwis)) outgwis[match(snp.causal, outgwis$SNP), ]$BETA else NULL
  se22ld <- if (!is.null(outgwis)) outgwis[match(snp.causal, outgwis$SNP), ]$SE else NULL
  
  bp <- expgwas.order$BP
  chr <- expgwas.order$CHR
  idx4panel <- matrix(numeric(0), nrow = 0, ncol = 1)
  
  # Step 7: Calculate R Matrix
  message("Step 7/7: Calculating correlation matrix...")
  Rblockres <- Cal_block_Rmatrix(bp, chr, avbIndex-1, idx4panel,
                                 block_file, stringname3, 1, coreNum, lam)
  R <- Rblockres$R; diag(R) <- 1
  
  end_time <- Sys.time()
  duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 1)
  message("\n", rep("-", 50))
  message("IV selection completed successfully!")
  message("Total SNPs selected: ", length(snp.causal))
  message("Total processing time: ", duration, " seconds")
  message(rep("-", 50))
  
  return(list(snp.causal = snp.causal, 
              gammah1 = bh11ld, gammah3 = bh12ld,
              Gammah1 = bh21ld, Gammah3 = bh22ld,
              se1 = se11ld, se2 = se12ld, se3 = se21ld, se4 = se22ld, 
              R = R))
}


MERLIN <- function(gammah1, gammah3 = NULL, Gammah1, Gammah3 = NULL,
                   se1 = NULL, se2 = NULL, se3 = NULL, se4 = NULL, 
                   R, rho_1 = NULL, rho_2 = NULL,
                   model = c("standard", "continuous_E", "discrete_E", "discrete_E_adj", "MO", "ME"),
                   p_common = NULL, p_exp = NULL, p_out = NULL,                   
                   maxIter = 12000, burnin = 5000, thin = 10, seed = 20262027) {
  
  # Match the model argument, defaulting to "standard"
  model <- match.arg(model)
  
  # Set random seed for reproducibility
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) stop("Error: 'seed' must be a single integer.")
    set.seed(as.integer(seed))
  }
  
  # Validate MCMC parameters
  if (!is.numeric(maxIter) || maxIter <= 0) stop("Error: 'maxIter' must be a positive integer.")
  if (!is.numeric(burnin) || burnin < 0) stop("Error: 'burnin' must be a non-negative integer.")
  if (!is.numeric(thin) || thin <= 0) stop("Error: 'thin' must be a positive integer.")
  if (maxIter %% thin != 0) stop("Error: 'maxIter' must be exactly divisible by 'thin'.")
  
  maxIter <- as.integer(maxIter)
  burnin <- as.integer(burnin)
  thin <- as.integer(thin)
  
  # Check core vectors required by all models
  core_vecs <- list(gammah1 = gammah1, Gammah1 = Gammah1)
  for (name in names(core_vecs)) {
    if (!is.numeric(core_vecs[[name]])) {
      stop(sprintf("Error: '%s' must be a numeric vector.", name))
    }
    if (anyNA(core_vecs[[name]]) || any(is.infinite(core_vecs[[name]]))) {
      stop(sprintf("Error: '%s' cannot contain NA, NaN, or Inf values.", name))
    }
  }
  
  # Check the LD matrix R
  if (!is.matrix(R) || !is.numeric(R)) stop("Error: 'R' must be a numeric matrix.")
  if (nrow(R) != ncol(R)) stop("Error: 'R' must be a square matrix.")
  if (anyNA(R) || any(is.infinite(R))) stop("Error: Matrix 'R' cannot contain NA, NaN, or Inf values.")
  
  p <- nrow(R)
  
  # Check dimension alignment between core vectors and matrix R
  for (name in names(core_vecs)) {
    if (length(core_vecs[[name]]) != p) {
      stop(sprintf("Error: The length of '%s' must match the dimension of matrix R.", name))
    }
  }
  
  start_time <- Sys.time()
  message(sprintf("Running MERLIN method (Model: %s)...", model))
  
  # Helper function to check model-specific required vectors
  check_required_vec <- function(vec, name) {
    if (is.null(vec)) stop(sprintf("Error: '%s' is required for the '%s' model.", name, model))
    if (!is.numeric(vec)) stop(sprintf("Error: '%s' must be a numeric vector.", name))
    if (length(vec) != p) stop(sprintf("Error: The length of '%s' must match the dimension of matrix R.", name))
    if (anyNA(vec) || any(is.infinite(vec))) stop(sprintf("Error: '%s' cannot contain NA, NaN, or Inf values.", name))
  }
  check_probability <- function(value, name) {
  if (is.null(value) || !is.numeric(value) || length(value) != 1 ||
      !is.finite(value) || value <= 0 || value >= 1) {
    stop(sprintf(
      "Error: '%s' must be a single numeric value strictly between 0 and 1, representing P(E = 1).",
      name
    ))
  }}
  
  # C++ Dispatch and Model-Specific Validation
  
  if (model %in% c("standard", "continuous_E", "discrete_E", "discrete_E_adj")) {
    
    check_required_vec(gammah3, "gammah3")
    check_required_vec(Gammah3, "Gammah3")
    check_required_vec(se1, "se1"); check_required_vec(se2, "se2")
    check_required_vec(se3, "se3"); check_required_vec(se4, "se4")
    
    if (any(se1 <= 0) || any(se2 <= 0) || any(se3 <= 0) || any(se4 <= 0)) {
      stop("Error: All standard errors (se1, se2, se3, se4) must be strictly greater than 0.")
    }
    if (is.null(rho_1) || is.null(rho_2) || !is.numeric(rho_1) || !is.numeric(rho_2)) {
      stop("Error: 'rho_1' and 'rho_2' must be single numeric values.")
    }
    if (abs(rho_1) >= 1 || abs(rho_2) >= 1) {
      stop("Error: 'rho_1' and 'rho_2' must be strictly between -1 and 1.")
    }
    
    if (model == "standard") {
      result <- MRGEI_Gam3seo_3to1(gammah1, gammah3, Gammah1, Gammah3, se1, se2, se3, se4, R, rho_1, rho_2, 0.5, maxIter, burnin, thin)
      res_final <- list(
        Beta1.hat  = result$Beta1.hat,
        Beta1.se   = result$Beta1.se,
        Beta1.pval = result$Beta1.pval,
        Beta4.hat  = result$Beta4.hat,
        Beta4.se   = result$Beta4.se,
        Beta4.pval = result$Beta4.pval)
    } else if (model == "continuous_E") {
      result <- MRGEI_Gam3seo_3to1(gammah1, gammah3, Gammah1, Gammah3, se1, se2, se3, se4, R, rho_1, rho_2, 0.5, maxIter, burnin, thin)
      res_final <- list(
        Beta1.hat  = result$Beta1.hat,
        Beta1.se   = result$Beta1.se,
        Beta1.pval = result$Beta1.pval,
        Beta4.hat  = result$Beta4.hat,
        Beta4.se   = result$Beta4.se,
        Beta4.pval = result$Beta4.pval)
    } else if (model == "discrete_E") {
      # discrete_E means that exposure and outcome have the same
      # P(E = 1). It does not necessarily mean p_common = 0.5.
      check_probability(p_common, "p_common")

      if (!is.null(p_exp) || !is.null(p_out)) {
        stop(
          "Error: For model = 'discrete_E', provide only 'p_common'. ",
          "Do not provide 'p_exp' or 'p_out'."
        )
      }

      p_target <- as.numeric(p_common)
      mu_target <- 2 * p_target - 1
      sd_target <- 2 * sqrt(p_target * (1 - p_target))

      message(sprintf(
        paste0(
          "\nBinary input scale check:\n",
          "  Model: discrete_E\n",
          "  Raw E coding: -1/+1\n",
          "  p_common = P(E=1) = %.4f\n",
          "  mean(E) = %.4f; SD(E) = %.4f\n",
          "  Exposure and outcome are assumed to have the same E proportion.\n",
          "  Input GWIS effects will be standardized automatically.\n",
          "  Returned beta estimates will be on the original E=-1/+1 scale."
        ),
        p_target, mu_target, sd_target
      ))

      # Marginal GWAS effects are already defined at the common mean of E.
      # Only GWIS BETA and SE need to be rescaled.
      gammah3_std <- sd_target * gammah3
      Gammah3_std <- sd_target * Gammah3
      se2_std <- sd_target * se2
      se4_std <- sd_target * se4

      result <- MRGEI_Gam3seo_3to1(
        gammah1, gammah3_std,
        Gammah1, Gammah3_std,
        se1, se2_std, se3, se4_std,
        R, rho_1, rho_2, p_target,
        maxIter, burnin, thin
      )

      # Back-transform paired MCMC draws to the original E=-1/+1 scale.
      beta1_std_draw <- as.numeric(result$Beta1res)
      beta4_std_draw <- as.numeric(result$Beta4res)

      beta4_raw_draw <- beta4_std_draw / sd_target
      beta1_raw_draw <- beta1_std_draw - (mu_target / sd_target) * beta4_std_draw

      beta1_hat <- mean(beta1_raw_draw)
      beta1_se <- sd(beta1_raw_draw)
      beta4_hat <- mean(beta4_raw_draw)
      beta4_se <- sd(beta4_raw_draw)

      res_final <- list(
        Beta1.hat  = beta1_hat,
        Beta1.se   = beta1_se,
        Beta1.pval = 2 * pnorm(abs(beta1_hat / beta1_se), lower.tail = FALSE),
        Beta4.hat  = beta4_hat,
        Beta4.se   = beta4_se,
        Beta4.pval = 2 * pnorm(abs(beta4_hat / beta4_se), lower.tail = FALSE)
      )
     } else if (model == "discrete_E_adj") {
       check_probability(p_exp, "p_exp")
       check_probability(p_out, "p_out")

       if (!is.null(p_common)) {
         stop(
           "Error: For model = 'discrete_E_adj', provide 'p_exp' and ",
           "'p_out'. Do not provide 'p_common'."
         )
       }

       p_exp <- as.numeric(p_exp)
       p_out <- as.numeric(p_out)

       if (isTRUE(all.equal(p_exp, p_out))) {
         warning(
           "'p_exp' and 'p_out' are equal. Consider using ",
           "model = 'discrete_E' with p_common = p_exp."
         )
       }
      mu_exp <- 2 * p_exp - 1
      mu_out <- 2 * p_out - 1
      sd_out <- 2 * sqrt(p_out * (1 - p_out))
      a_exp_to_out <- mu_out - mu_exp

      message(sprintf(
        paste0(
          "\nBinary input scale check:\n",
          "  Model: discrete_E_adj\n",
          "  Raw E coding: -1/+1\n",
          "  p_exp = P(E=1 in exposure) = %.4f\n",
          "  p_out = P(E=1 in outcome)  = %.4f\n",
          "  mean(E_exp) = %.4f; mean(E_out) = %.4f\n",
          "  Target analysis scale: outcome standardized E scale\n",
          "  Exposure GWAS will be projected to the outcome E proportion.\n",
          "  Exposure and outcome GWIS effects will be standardized using p_out.\n",
          "  Returned beta estimates will be on the original E=-1/+1 scale."
        ),
        p_exp, p_out, mu_exp, mu_out
      ))

      # Exposure marginal GWAS is moved from the exposure E mean
      # to the outcome E mean.
      gammah1_out <- gammah1 + a_exp_to_out * gammah3

      var_gammah1_out <- se1^2 + a_exp_to_out^2 * se2^2

      if (any(!is.finite(var_gammah1_out)) ||
          any(var_gammah1_out <= 0)) {
        stop(
          "Error: The transformed exposure GWAS variance is ",
          "non-positive. Check p_exp, p_out."
        )
      }

      se1_out <- sqrt(var_gammah1_out)

      # Convert both GWIS inputs to the outcome standardized E scale.
      gammah3_out <- sd_out * gammah3
      Gammah3_out <- sd_out * Gammah3
      se2_out <- sd_out * se2
      se4_out <- sd_out * se4

      # Outcome marginal GWAS is already at the outcome E mean.
      result <- MRGEI_Gam3seo_3to1(
        gammah1_out, gammah3_out,
        Gammah1, Gammah3_out,
        se1_out, se2_out, se3, se4_out,
        R, rho_1, rho_2, p_out,
        maxIter, burnin, thin
      )

      # Back-transform paired MCMC draws from the outcome standardized
      # E scale to the original E=-1/+1 scale.
      beta1_std_draw <- as.numeric(result$Beta1res)
      beta4_std_draw <- as.numeric(result$Beta4res)

      beta4_raw_draw <- beta4_std_draw / sd_out
      beta1_raw_draw <- beta1_std_draw - (mu_out / sd_out) * beta4_std_draw

      beta1_hat <- mean(beta1_raw_draw)
      beta1_se <- sd(beta1_raw_draw)
      beta4_hat <- mean(beta4_raw_draw)
      beta4_se <- sd(beta4_raw_draw)

      res_final <- list(
        Beta1.hat  = beta1_hat,
        Beta1.se   = beta1_se,
        Beta1.pval = 2 * pnorm(abs(beta1_hat / beta1_se), lower.tail = FALSE),
        Beta4.hat  = beta4_hat,
        Beta4.se   = beta4_se,
        Beta4.pval = 2 * pnorm(abs(beta4_hat / beta4_se), lower.tail = FALSE)
      )    
      }
    } else if (model == "MO") { 
    
    # The new MRGEI_Gamseo function acts as the MO model. 
    # It requires gammah3, se1, se2, se3, and rho_1 (passed as rho).
    check_required_vec(gammah3, "gammah3")
    check_required_vec(se1, "se1")
    check_required_vec(se2, "se2")
    check_required_vec(se3, "se3")
    
    if (any(se1 <= 0) || any(se2 <= 0) || any(se3 <= 0)) {
      stop("Error: Standard errors (se1, se2, se3) must be strictly greater than 0.")
    }
    if (is.null(rho_1) || !is.numeric(rho_1) || length(rho_1) != 1) {
      stop("Error: 'rho_1' must be a single numeric value.")
    }
    if (abs(rho_1) >= 1) {
      stop("Error: 'rho_1' must be strictly between -1 and 1.")
    }
    result <- MRGEI_Gamseo(gammah1, gammah3, Gammah1, se1, se2, se3, R, rho_1, maxIter, burnin, thin)
    res_final <- list(
      Beta1.hat  = result$Beta1.hat,
      Beta1.se   = result$Beta1.se,
      Beta1.pval = result$Beta1.pval,
      Beta4.hat  = result$Beta4.hat,
      Beta4.se   = result$Beta4.se,
      Beta4.pval = result$Beta4.pval)
    } else if (model == "ME") { 
    
    check_required_vec(Gammah3, "Gammah3")
    check_required_vec(se1, "se1")
    check_required_vec(se3, "se3")
    check_required_vec(se4, "se4")
    
    if (any(se1 <= 0) || any(se3 <= 0) || any(se4 <= 0)) {
      stop("Error: Standard errors (se1, se3, se4) must be strictly greater than 0.")
    }
    
    res.beta1 <- TwoSampleMR::mr_egger_regression(gammah1, Gammah1, se1, se3)
    res.beta4 <- TwoSampleMR::mr_egger_regression(gammah1, Gammah3, se1, se4)
    
    res_final <- list(
      Beta1.hat  = res.beta1$b,
      Beta1.se   = res.beta1$se,
      Beta1.pval = res.beta1$pval,
      
      Beta4.hat  = res.beta4$b,
      Beta4.se   = res.beta4$se,
      Beta4.pval = res.beta4$pval
    )
  }
  
  end_time <- Sys.time()
  duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 1)
  
  message("\n", rep("-", 50))
  message("MERLIN Analysis Completed Successfully!")
  message("Total processing time: ", duration, " seconds")
  message(rep("-", 50))
  
  return(res_final)
}



traceplot <- function(bhatpoint){
  y <- bhatpoint
  x <- 1:length(bhatpoint)
  da <- cbind(x, y)
  dat <- data.frame(da)

  p1 <- ggplot2::ggplot(data = dat, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = paste("Traceplot of", deparse(substitute(bhatpoint))),
      x = "GibbsSampleIndex",
      y = expression(hat(beta))) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12, face = "bold"))

  return(p1)
}


# copy from QingCheng0218/MR.CUE - GitHub
mhcstart = 28477797;
mhcend = 33448354;
summaryQC = function(mhcstart, mhcend, bh1, bh2, s12, s22, bp, chr,
                     rsname, avbIndex, idx4panel, xbound, ybound){
  # remove SNPs in MHC region:
  idxchr6 = which(chr==6);
  idxcut = idxchr6[which(bp[idxchr6]>=mhcstart & bp[idxchr6]<=mhcend)];
  pmhc = length(idxcut);

  if(pmhc!=0){
    bh1Rmhc = bh1[-idxcut];
    bh2Rmhc = bh2[-idxcut];
    s12Rmhc = s12[-idxcut];
    s22Rmhc = s22[-idxcut];
    bpRmhc = bp[-idxcut];
    chrRmhc = chr[-idxcut];
    rsnameRmhc = rsname[-idxcut];
    avbIndexRmhc = avbIndex[-idxcut];

    tmp0 = 1:length(bh1);
    tmp = tmp0[-idxcut];
    if(length(idx4panel)!=0){
      idx4panelRmhc = match(avbIndex[intersect((idx4panel + 1), tmp)], avbIndexRmhc) -1
    }else{
      idx4panelRmhc = idx4panel;
    }

  }else{
    bh1Rmhc = bh1;
    bh2Rmhc = bh2;
    s12Rmhc = s12;
    s22Rmhc = s22;
    bpRmhc = bp;
    chrRmhc = chr;
    rsnameRmhc = rsname;
    avbIndexRmhc = avbIndex;
    idx4panelRmhc = idx4panel;
  }


  # remove SNPs(exposure) with chi-square >80
  idx = which((bh1Rmhc/s12Rmhc)^2>xbound);
  px = length(idx);
  if(px!=0){
    bh1Rmhc_x = bh1Rmhc[-idx];
    bh2Rmhc_x = bh2Rmhc[-idx];
    s12Rmhc_x = s12Rmhc[-idx];
    s22Rmhc_x = s22Rmhc[-idx];
    bpRmhc_x = bpRmhc[-idx];
    chrRmhc_x = chrRmhc[-idx];
    rsnameRmhc_x = rsnameRmhc[-idx];
    avbIndexRmhc_x = avbIndexRmhc[-idx];
    # idx4panelRmhc_x = idx4panelRmhc[-idx];


    tmp0 = 1:length(bh1Rmhc);
    tmp = tmp0[-idx];
    if(length(idx4panel)!=0){
      idx4panelRmhc_x = match(avbIndexRmhc[intersect((idx4panelRmhc + 1), tmp)], avbIndexRmhc_x) -1;
    }else{
      idx4panelRmhc_x = idx4panel;
    }


  }else{
    bh1Rmhc_x = bh1Rmhc;
    bh2Rmhc_x = bh2Rmhc;
    s12Rmhc_x = s12Rmhc;
    s22Rmhc_x = s22Rmhc;
    bpRmhc_x = bpRmhc;
    chrRmhc_x = chrRmhc;
    rsnameRmhc_x = rsnameRmhc;
    avbIndexRmhc_x = avbIndexRmhc;
    idx4panelRmhc_x = idx4panelRmhc;
  }

  # remove SNPs(outcome) with chi-square >80
  idy = which((bh2Rmhc_x/s22Rmhc_x)^2>ybound);
  py = length(idy);
  if(py!=0){
    bh1Rmhc_xy = bh1Rmhc_x[-idy];
    bh2Rmhc_xy = bh2Rmhc_x[-idy];
    s12Rmhc_xy = s12Rmhc_x[-idy];
    s22Rmhc_xy = s22Rmhc_x[-idy];
    bpRmhc_xy = bpRmhc_x[-idy];
    chrRmhc_xy = chrRmhc_x[-idy];
    rsnameRmhc_xy = rsnameRmhc_x[-idy];
    avbIndexRmhc_xy = avbIndexRmhc_x[-idy];
    # idx4panelRmhc_xy = idx4panelRmhc_x[-idy];

    tmp0 = 1:length(bh1Rmhc_x);
    tmp = tmp0[-idx];

    if(length(idx4panel)!=0){
      idx4panelRmhc_xy = match(avbIndexRmhc_x[intersect((idx4panelRmhc_x + 1), tmp)], avbIndexRmhc_xy) -1;
    }else{
      idx4panelRmhc_xy = idx4panel;
    }

  }else{
    bh1Rmhc_xy = bh1Rmhc_x;
    bh2Rmhc_xy = bh2Rmhc_x;
    s12Rmhc_xy = s12Rmhc_x;
    s22Rmhc_xy = s22Rmhc_x;
    bpRmhc_xy = bpRmhc_x;
    chrRmhc_xy = chrRmhc_x;
    rsnameRmhc_xy = rsnameRmhc_x;
    avbIndexRmhc_xy = avbIndexRmhc_x;
    idx4panelRmhc_xy = idx4panelRmhc_x;
  }
  return(list(bh1new = bh1Rmhc_xy, bh2new = bh2Rmhc_xy, s12new = s12Rmhc_xy, s22new = s22Rmhc_xy,
              bpnew = bpRmhc_xy, chrnew = chrRmhc_xy, rsnamenew = rsnameRmhc_xy,
              avbIndexnew = avbIndexRmhc_xy, idx4panelnew = idx4panelRmhc_xy, pmhc = pmhc, px = px, py = py))
}


EstRhofun <- function(fileexposure, fileoutcome, stringname3,
                      ld_r2_thresh, lam, pth, coreNum = 1){

  # Estimate the rho
  res = matchsnp(fileexposure, fileoutcome, stringname3, FALSE);
  bh1 = as.numeric(res$bh1);
  bh2 = as.numeric(res$bh2);
  s12 = as.numeric(res$s12);
  s22 = as.numeric(res$s22);
  chr = as.numeric(res$chr);
  bp = res$bp;
  rsname = res$rsname
  avbIndex = res$idxin;
  idx4panel = res$idx4panel;
  QCresult = summaryQC(mhcstart, mhcend, bh1, bh2, s12, s22, bp,
                       chr, rsname, avbIndex, idx4panel, Inf, Inf);


  bh1new = QCresult$bh1new;
  bh2new = QCresult$bh2new;
  s12new = QCresult$s12new;
  s22new = QCresult$s22new;
  bpnew = QCresult$bpnew;
  chrnew = QCresult$chrnew;
  avbIndexnew = QCresult$avbIndexnew;
  rsnamenew = QCresult$rsnamenew;
  pmhc = QCresult$pmhc;
  px = QCresult$px;
  py = QCresult$py;

  max_cores <- parallel::detectCores()
  
  if (is.null(coreNum) || !is.numeric(coreNum) || coreNum < 1) {
    coreNum <- 1
    warning("Invalid 'coreNum' provided. Defaulting to 1 core.")
  } else if (coreNum > max_cores) {
    coreNum <- max_cores
    warning(sprintf("Requested 'coreNum' exceeds available cores. Adjusting to maximum available: %d cores.", max_cores))
  }
  
  coreNum <- as.integer(coreNum)

  IndSumRes = IndepSummary(bpnew, chrnew, avbIndexnew - 1, block_file, stringname3,
                           bh1new, bh2new, s12new, s22new, coreNum,
                           lam, ld_r2_thresh);
  bh1_ind = IndSumRes$bh1_ind;
  bh2_ind = IndSumRes$bh2_ind;
  se1_ind = IndSumRes$se1_ind;
  se2_ind = IndSumRes$se2_ind;

  z1_ind = bh1_ind / se1_ind;
  z2_ind = bh2_ind / se2_ind;

  # a = rep(-pth, 2);
  # b = rep(pth, 2);
  # z1_new = z1_ind [which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
  # z2_new = z2_ind[which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
  # rhores = truncEstfun(a, b, z1_new, z2_new, 4000, 1000, 10)
  # rhohat = mean(rhores);
  # p1 = length(z1_new);
  # pvalue = testR(rhohat, p1);
  #
  maxIter = 4000;
  thin = 10;
  burnin = 1000;

  nsave = maxIter / thin;

  if(length(pth)==1){
    Rhores = rep(NA, nsave);
    a = rep(-pth, 2);
    b = rep(pth, 2);
    z1_new = z1_ind [which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
    z2_new = z2_ind[which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
    rhores = truncEstfun(a, b, z1_new, z2_new, maxIter, burnin, thin)
    rhohat = mean(rhores);
    p1 = length(z1_new);
    pvalue = testR(rhohat, p1);
    Rhores = rhores;
    pres = p1;
  }else{
    rhohat = pvalue = rep(NA, length(pth));
    Rhores = matrix(NA, nrow = nsave, ncol = length(pth));
    pres = rep(NA, length(pth));
    for(i in 1:length(pth)){
      pth1 = pth[i];
      a = rep(-pth1, 2);
      b = rep(pth1, 2);
      z1_new = z1_ind [which(abs(z1_ind) < pth1&abs(z2_ind) < pth1)];
      z2_new = z2_ind[which(abs(z1_ind) < pth1&abs(z2_ind) < pth1)];
      rhores = truncEstfun(a, b, z1_new, z2_new, 4000, 1000, 10)
      rhohat[i] = mean(rhores);
      p1 = length(z1_new);
      pvalue[i] = testR(rhohat[i], p1);
      Rhores[, i] = rhores;
      pres[i] = p1;
    }
  }
  # ---------------------------------------------------------
  return(list(rhohat = rhohat, pvalue = pvalue, pres = pres, Rhores = Rhores))

}

