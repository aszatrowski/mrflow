mr_analyze <- function(exposure, outcome, p, ivw_only) {
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
    input <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse, snps = snps) #nolint

    if (ivw_only == TRUE) {
        mr_ivw_result <- MendelianRandomization::mr_ivw(input)
        print(mr_ivw_result)
    } else {
        mr_all_result <- MendelianRandomization::mr_allmethods(input)
        print(mr_all_result)

    }
}