#' Quickly run two-sample MR analyses of summary statistic data
#'
#' @param exposure A summary-statistic GWAS dataset, passed as a string, to be used as the exposure
#' @param outcome A summary-statistic GWAS dataset, passed as a string, to be used as the outcome
#' @param p p-value threshold for SNP retrieval from datasets
#'
#' @return A chart of MR results from various methods (an S4 object).
#' @export
#'
#' @examples
#' mr_analyze("prot-a-16", "finn-b-I9_CARDMYO", 1e-06)
#' Extracting data for 4 SNP(s) from 1 GWAS(s)
#' Harmonising  || id:prot-a-16 (prot-a-16) and Cardiomyopathy || id:finn-b-I9_CARDMYO 
#' (finn-b-I9_CARDMYO)
#'                     Method Estimate Std Error 95% CI         P-value
#'              Simple median    0.006     0.139  -0.268  0.279   0.968
#'            Weighted median   -0.013     0.139  -0.286  0.261   0.928
#'  Penalized weighted median   -0.013     0.139  -0.286  0.261   0.928
#'                                                                     
#'                        IVW   -0.096     0.113  -0.318  0.125   0.394
#'              Penalized IVW   -0.096     0.113  -0.318  0.125   0.394
#'                 Robust IVW   -0.043     0.338  -0.704  0.619   0.900
#'       Penalized robust IVW   -0.043     0.338  -0.704  0.619   0.900
#'                                                                     
#'                   MR-Egger    0.627     0.470  -0.294  1.547   0.182
#'                (intercept)   -0.229     0.145  -0.512  0.054   0.113
#'         Penalized MR-Egger    0.627     0.470  -0.294  1.547   0.182
#'                (intercept)   -0.229     0.145  -0.512  0.054   0.113
#'            Robust MR-Egger    0.628     0.158   0.319  0.937   0.000
#'                (intercept)   -0.230     0.045  -0.318 -0.141   0.000
#'  Penalized robust MR-Egger    0.628     0.158   0.319  0.937   0.000

mr_analyze <- function(exposure, outcome, p) {
    # extract exposures and outcomes from MRCIEU openGWAS
    exposures <- TwoSampleMR::extract_instruments(exposure, p1 = p) # nolint

    outcomes <- TwoSampleMR::extract_outcome_data(snps = exposures$SNP, outcomes = outcome) # nolint

    # harmonize phenotypes
    harmonized <- TwoSampleMR::harmonise_data(exposures, outcomes) # nolint

    # generate analyzable dataset
    rsid <- harmonized$SNP
    exposure_beta <- harmonized$beta.exposure
    exposure_se <- harmonized$se.exposure
    outcome_beta <- harmonized$beta.outcome
    outcome_se <- harmonized$se.outcome
    iv_snps <- cbind(rsid, exposure_beta, exposure_se, outcome_beta, outcome_se)
    iv_snps <- as.data.frame(iv_snps)

    colnames(iv_snps) <- c("rsid",
                        "exposure_beta",
                        "exposure_se",
                        "outcome_beta",
                        "outcome_se")

    bx <- as.numeric(exposure_beta)
    bxse <- as.numeric(exposure_se)
    by <- as.numeric(outcome_beta)
    byse <- as.numeric(outcome_se)
    snps <- rsid

    # Run MR
    input <- MendelianRandomization::mr_input(bx = bx, bxse = bxse, by = by, byse = byse, snps = snps) #nolint

    mr_all_result <- MendelianRandomization::mr_allmethods(input)
    print(mr_all_result)
}
