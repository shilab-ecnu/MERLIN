
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


ivselect <- function(expgwas_dir, expgwis_dir, outgwas_dir, outgwis_dir,
                     stringname3, block_file, plink_dir = NULL,
                     pval_cutoff_gwas = 0.00000005, pval_cutoff_gwis = 0.00000005,
                     r2_cutoff = 0.01, kb_cutoff = 1024, maf_cutoff = 0.05,
                     lam = 0.1, coreNum = 1, intersect_mode = FALSE){

  start_time <- Sys.time()
  message("Starting MERLIN IV selection process...")

  message("Step 1/7: Reading reference panel (bim file)...")
  bim <- data.table::fread(paste0(stringname3,".bim"), header = FALSE)
  bim <- data.frame(bim)

  message("Step 2/7: Reading GWAS summary statistics...")
  expgwas <- readr::read_delim(expgwas_dir, col_names=TRUE)
  expgwas <- data.frame(expgwas)
  expgwis <- readr::read_delim(expgwis_dir, col_names=TRUE)
  expgwis <- data.frame(expgwis)
  outgwas <- readr::read_delim(outgwas_dir, col_names=TRUE)
  outgwas <- data.frame(outgwas)
  outgwis <- readr::read_delim(outgwis_dir, col_names=TRUE)
  outgwis <- data.frame(outgwis)

  message("Step 3/7: Setting up PLINK...")
  if (is.null(plink_dir)) {
    plink_dir <- bigsnpr::download_plink()
  }
  message("  ✓ Using PLINK: ", plink_dir)

  message("Step 4/7: Running LD clumping for exposure GWAS...")
  expgwas.LD.cmd <- paste(plink_dir, " --bfile ", stringname3,
                          " --clump-p1 ", pval_cutoff_gwas,
                          " --clump-r2 ", r2_cutoff,
                          " --clump-kb ", kb_cutoff,
                          " --maf ", maf_cutoff,
                          " --clump ", expgwas_dir,
                          " --clump-snp-field SNP --clump-field P --out ld_result1",
                          sep="")
  system(expgwas.LD.cmd)
  expgwas.LD <- data.table::fread("ld_result1.clumped")
  snp.expgwas.LD <- expgwas.LD$SNP
  file.remove("ld_result1.clumped")
  file.remove("ld_result1.log")

  message("Step 5/7: Running LD clumping for exposure GWIS...")
  expgwis.LD.cmd <- paste(plink_dir, " --bfile ", stringname3,
                          " --clump-p1 ", pval_cutoff_gwis,
                          " --clump-r2 ", r2_cutoff,
                          " --clump-kb ", kb_cutoff,
                          " --maf ", maf_cutoff,
                          " --clump ", expgwis_dir,
                          " --clump-snp-field SNP --clump-field P --out ld_result2",
                          sep="")
  system(expgwis.LD.cmd)
  expgwis.LD <- data.table::fread("ld_result2.clumped")
  snp.expgwis.LD <- expgwis.LD$SNP
  file.remove("ld_result2.clumped")
  file.remove("ld_result2.log")

  message("Step 6/7: Selecting causal SNPs...")
  message("  - Mode: ", ifelse(intersect_mode, "Intersection", "Union"))
  if (intersect_mode == TRUE) {
    # intersect mode
    # appropriate when stringent quality criteria for instrumental variables are prioritized
    # potentially yielding fewer eligible variants.
    snp.causal <- intersect(snp.expgwas.LD, snp.expgwis.LD)

  } else {
    union.snp <- union(snp.expgwas.LD, snp.expgwis.LD)
    union.snp.df <- data.frame(SNP = union.snp, P = 1)
    union.dir <- "union.snp.txt"
    write.table(union.snp.df, file = union.dir,
                quote = FALSE, row.names = FALSE)

    union.LD.cmd <- paste(plink_dir, " --bfile ", stringname3,
                          " --clump-p1 ", 1,
                          " --clump-r2 ", r2_cutoff,
                          " --clump-kb ", kb_cutoff,
                          " --maf ", maf_cutoff,
                          " --clump ", union.dir,
                          " --clump-snp-field SNP --clump-field P --out ld_result3",
                          sep="")
    system(union.LD.cmd)
    union.LD <- data.table::fread("ld_result3.clumped")
    snp.causal <- union.LD$SNP
    file.remove(union.dir)
    file.remove("ld_result3.clumped")
    file.remove("ld_result3.log")
  }

  snp.causal <- intersect(snp.causal, expgwas$SNP)
  snp.causal <- intersect(snp.causal, expgwis$SNP)
  snp.causal <- intersect(snp.causal, outgwas$SNP)
  snp.causal <- intersect(snp.causal, outgwis$SNP)

  if (length(snp.causal) == 0) {
    stop("No SNPs remaining after filtering.")
  } else {
    message("Numbers of causal SNPs selected:",length(snp.causal))
  }

  avbIndex <- match(snp.causal, bim$V2)
  avbIndex <- as.matrix(avbIndex[order(avbIndex)])
  snp.causal <- bim[avbIndex, ]$V2

  expgwas.order <- expgwas[match(snp.causal, expgwas$SNP), ]
  expgwis.order <- expgwis[match(snp.causal, expgwis$SNP), ]
  outgwas.order <- outgwas[match(snp.causal, outgwas$SNP), ]
  outgwis.order <- outgwis[match(snp.causal, outgwis$SNP), ]

  bh11ld <- expgwas.order$BETA
  se11ld <- expgwas.order$SE
  bh12ld <- expgwis.order$BETA
  se12ld <- expgwis.order$SE
  bh21ld <- outgwas.order$BETA
  se21ld <- outgwas.order$SE
  bh22ld <- outgwis.order$BETA
  se22ld <- outgwis.order$SE

  bp <- expgwas.order$BP
  chr <- expgwas.order$CHR
  idx4panel <- matrix(numeric(0), nrow = 0, ncol = 1)

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

  return(list(snp.causal = snp.causal, gammah1 = bh11ld, gammah3 = bh12ld,
              Gammah1 = bh21ld, Gammah3 = bh22ld,
              se1 = se11ld, se2 = se12ld, se3 = se21ld, se4 = se22ld, R = R))

}


MERLIN <- function(gammah1, gammah3, Gammah1, Gammah3,
                  se1, se2, se3, se4, R, rho_1, rho_2) {

  if (!all.equal(
    length(gammah1), length(gammah3), length(Gammah1), length(Gammah3),
    length(se1), length(se2), length(se3), length(se4),
    nrow(R), ncol(R)
  )) {
    stop("All input vectors must have the same length as the LD matrix dimensions")
  }

  start_time <- Sys.time()
  message("Running MERLIN method...")

  result <- MRGEI_Gam3seo(
    gammah1 = gammah1,
    gammah3 = gammah3,
    Gammah1 = Gammah1,
    Gammah3 = Gammah3,
    se1 = se1,
    se2 = se2,
    se3 = se3,
    se4 = se4,
    R = R,
    rho_1 = rho_1,
    rho_2 = rho_2
  )

  end_time <- Sys.time()
  duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 1)

  message("\n", rep("-", 50))
  message("MERLIN Analysis Completed Successfully!")
  message("Total processing time: ", duration, " seconds")
  message(rep("-", 50))

  return(result)
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
                      ld_r2_thresh, lam, pth){

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

  coreNum = 24;

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
  return(list(rhohat = rhohat, pvalue = pvalue, pres = pres, Rhores = Rhores));

}



