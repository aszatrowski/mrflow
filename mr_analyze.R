mr_analyze <- function(exposure, outcome, p, ivw_only) {
    exposures <- extract_instruments(exposure, p1 = p) # nolint

    outcomes <- extract_outcome_data(snps = exposures$SNP, outcomes = outcome) # nolint

    harmonized <- harmonise_data(exposures, outcomes) # nolint

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

    input <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse, snps = snps) #nolint
    mr_ivw_result <- MendelianRandomization::mr_ivw(input)
    mr_all_result <- MendelianRandomization::mr_allmethods(input)

    if (ivw_only == TRUE) {
        print(mr_ivw_result)
    }
    else {
        print(mr_all_result)
    }
}