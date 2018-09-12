#Which way to optimize based on flag


if (flagMidPercentile == TRUE) {
    corr.Cuts = seq(
    floor(sort(simMeasChange - 1)[[floor(((100 - midPercentile) / 200) * length(simMeasChange))]] * 100 / 5) * 5 / 100,
    ceil(sort(simMeasChange - 1)[[ceil(((midPercentile + ((100 - midPercentile) / 2)) / 100) * length(simMeasChange))]] * 100 / 5) * 5 / 100,
    by=0.05)
} else {
    if (flagCloseToClinical == TRUE) {
        if (flagCloseToClinicalNo0 == TRUE) {
            corr.Cuts = c(
            seq(-0.4, -0.2, by = 0.05),
            seq(0.1, 0.3, by = 0.05)
            )
        } else {
            corr.Cuts = seq(-0.5, 0.4, by = 0.05)
        }
    } else {
        corr.Cuts = seq(-1, 1, by = 0.05)
    }
}