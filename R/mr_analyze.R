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
