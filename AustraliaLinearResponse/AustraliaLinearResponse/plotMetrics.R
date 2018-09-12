# plotting the max c-statistic and npv,ppv,spec,sens at the max c-statistic
tmpBaseDir = file.path('C:', 'Users', 'Justin', 'Documents', '_UC intern')
png(filename = file.path(tmpBaseDir,'allMetricsAtMaxC.png'))

plot(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2],xlim=c(50,500),
    ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
    main = "Values of Metrics for Different Numbers of Simulated Patients", pch = c(16))
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2])

points(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 3]))],
        na.omit(variedPatientNumbersData[, 3]), col = rgb(1, 0, 0), pch = c(16))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 3]))],
        na.omit(variedPatientNumbersData[, 3]), col = rgb(1, 0, 0))

points(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 4]))],
        na.omit(variedPatientNumbersData[, 4]), col = "blue", pch = c(16))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 4]))],
        na.omit(variedPatientNumbersData[, 4]), col = "blue")

points(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 5]))],
        na.omit(variedPatientNumbersData[, 5]), col = rgb(0, 0.7, 0), pch = c(16))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 5]))],
        na.omit(variedPatientNumbersData[, 5]), col = rgb(0, 0.7, 0))

points(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 6]))],
        na.omit(variedPatientNumbersData[, 6]), col = rgb(0.7, 0, 0.7), pch = c(16))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 6]))],
        na.omit(variedPatientNumbersData[, 6]), col = rgb(0.7, 0, 0.7))
legend("bottomleft", legend = c("Max c-statistic", "NPV at Max c", "PPV at Max c", "Specificity at Max c", "Sensitivity at Max c"), pch = c(16, 16,16, 16, 16),
       col = c("black", "red", "blue", rgb(0, 0.7, 0), rgb(0.7, 0, 0.7)), lty = c(1,1,1, 1, 1))
dev.off()

# plotting npv when max c and max npv
png(filename = file.path(tmpBaseDir, 'npv.png'))
plot(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 3]))],
        na.omit(variedPatientNumbersData[, 3]), col = "red", xlim = c(50, 500), ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
        main = "NPV at Maximum c-statistic and Maximum NPV", pch = 16,cex=1.2)
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 3]))],
        na.omit(variedPatientNumbersData[, 3]), col = "red")

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 7], col = "lightcoral", pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 7], col = "lightcoral", lty = 2)

legend("right", legend = c("NPV at Max c", "Max NPV"), pch = c(16, 17),
       col = c("red","lightcoral"), lty = c(1, 2), lwd = c(1, 1))
dev.off()

# plotting ppv when max c and max ppv
png(filename = file.path(tmpBaseDir, 'ppv.png'))
plot(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 4]))],
        na.omit(variedPatientNumbersData[, 4]), col = "blue", xlim = c(50, 500), ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
        main = "PPV at Maximum c-statistic and Maximum PPV", pch = 16, cex = 1.2)
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 4]))],
        na.omit(variedPatientNumbersData[, 4]), col = "blue", lwd = 1, lty = 1)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 8], col = rgb(0.1, 0.7, 1), pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 8], col = rgb(0.1, 0.7, 1), lty = 2)

legend("right", legend = c("PPV at Max c", "Max PPV"), pch = c(16, 17),
       col = c("blue", rgb(0.1, 0.7, 1)), lty = c(1, 2), lwd = c(1, 1))
dev.off()

# plotting spec when max c and max spec
png(filename = file.path(tmpBaseDir, 'spec.png'))
plot(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 5]))],
        na.omit(variedPatientNumbersData[, 5]), col = rgb(0, 0.7, 0), xlim = c(50, 500), ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
        main = "Specificity at Maximum c-statistic and Maximum Specificity", pch = 16, cex = 1.2)
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 5]))],
        na.omit(variedPatientNumbersData[, 5]), col = rgb(0, 0.7, 0), lwd = 1, lty = 1)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 9], col = rgb(0, 1, 0.5), pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 9], col = rgb(0, 1, 0.5), lty = 2)
legend("right", legend = c("Spec at Max c", "Max Spec"), pch = c(16, 17),
       col = c(col = rgb(0, 0.7, 0), col = rgb(0, 1, 0.5)), lty = c(1, 2), lwd = c(1, 1))
dev.off()

# plotting sens when max c and max sens
png(filename = file.path(tmpBaseDir, 'sens.png'))
plot(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 6]))],
        na.omit(variedPatientNumbersData[, 6]), col = rgb(0.7, 0, 0.7), xlim = c(50, 500), ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
        main = "Sensitivity at Maximum c-statistic and Maximum Sensitivity", pch = 16, cex = 1.2)
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 6]))],
        na.omit(variedPatientNumbersData[, 6]), col = rgb(0.7, 0, 0.7), lwd = 1, lty = 1)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 10], col =rgb(1,0.5,1),pch=17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 10], col = rgb(1, 0.5, 1), lty = 2)

legend("right", legend = c("Sens at Max c", "Max Sens"), pch = c(16, 17),
       col = c(col = rgb(0.7, 0, 0.7), col = rgb(1, 0.5, 1)), lty = c(1, 2), lwd = c(1, 1))
dev.off()

# plotting npvnonResp when max c and max npvnonResp
png(filename = file.path(tmpBaseDir, 'npvNonResp.png'))
plot(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 11]))],
        na.omit(variedPatientNumbersData[, 11]), col = rgb(1,0,0), xlim = c(50, 500), ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
        main = "NPV NonResp at Maximum c-statistic\nand Maximum NPV NonResp", pch = 16, cex = 1.2)
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 11]))],
        na.omit(variedPatientNumbersData[, 11]), col = rgb(1,0,0), lwd = 1, lty = 1)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 15], col = "lightcoral", pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 15], col = "lightcoral", lty = 2)

legend("right", legend = c("NPV NonResp at Max c", "Max NPV NonResp"), pch = c(16, 17),
       col = c(col = rgb(1, 0, 0), "lightcoral"), lty = c(1, 2), lwd = c(1, 1))
dev.off()

# plotting ppvProgFree when max c and max ppvProgFree
png(filename = file.path(tmpBaseDir, 'ppvProgFree.png'))
plot(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 12]))],
        na.omit(variedPatientNumbersData[, 12]), col = "blue", xlim = c(50, 500), ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
        main = "PPV ProgFree at Maximum c-statistic\nand Maximum PPV ProgFree", pch = 16, cex = 1.2)
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 12]))],
        na.omit(variedPatientNumbersData[, 12]), col = "blue", lwd = 1, lty = 1)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 16], col = rgb(0.1, 0.7, 1), pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 16], col = rgb(0.1, 0.7, 1), lty = 2)

legend("right", legend = c("PPV ProgFree at Max c", "Max PPV ProgFree"), pch = c(16, 17),
       col = c(col = "blue", col = rgb(0.1, 0.7, 1)), lty = c(1, 2), lwd = c(1, 1))
dev.off()

# plotting specNonResp when max c and max specNonResp
png(filename = file.path(tmpBaseDir, 'specNonResp.png'))
plot(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 13]))],
        na.omit(variedPatientNumbersData[, 13]), col = rgb(0, 0.7, 0), xlim = c(50, 500), ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
        main = "Spec NonResp at Maximum c-statistic\nand Maximum Spec NonResp", pch = 16, cex = 1.2)
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 12]))],
        na.omit(variedPatientNumbersData[, 13]), col = rgb(0, 0.7, 0), lwd = 1, lty = 1)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 17], col = rgb(0,1,0.5), pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 17], col = rgb(0, 1, 0.5), lty = 2)

legend("right", legend = c("Spec NonResp at Max c", "Max Spec NonResp"), pch = c(16, 17),
       col = c(col = rgb(0, 0.7, 0), col = rgb(0, 1, 0.5)), lty = c(1, 2), lwd = c(1, 1))
dev.off()

# plotting sensProgFree when max c and max sensProgFree
png(filename = file.path(tmpBaseDir, 'sensProgFree.png'))
plot(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 14]))],
        na.omit(variedPatientNumbersData[, 14]), col = rgb(0.7,0,0.7), xlim = c(50, 500), ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
        main = "Sens ProgFree at Maximum c-statistic\nand Maximum Sens ProgFree", pch = 16, cex = 1.2)
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 14]))],
        na.omit(variedPatientNumbersData[, 14]), col = rgb(0.7,0,0.7), lwd = 1, lty = 1)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 18], col = rgb(1, 0.5, 1), pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 18], col = rgb(1, 0.5, 1), lty = 2)

legend("right", legend = c("Sens ProgFree at Max c", "Max Sens ProgFree"), pch = c(16, 17),
       col = c(col = rgb(0.7, 0, 0.7), col = rgb(1, 0.5, 1)), lty = c(1, 2), lwd = c(1, 1))
dev.off()

# plotting Max NPV, c-stat, and dist between
png(filename = file.path(tmpBaseDir, 'Max NPV Nearest.png'))
plot(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2], xlim = c(50, 500),
    ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
    main = "Maximum c-statistic, Maximum NPV, \nand Distance Between", pch = c(16))
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2])

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 7], col = "lightcoral", pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 7], col = "lightcoral", lty = 2)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 19], pch = 16, col = rgb(0.5,0.5,0.5))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 19]))],
        na.omit(variedPatientNumbersData[, 19]), col = rgb(0.5,0.5,0.5), lwd = 1, lty =2)

legend("topright", legend = c("Max c-statistic", "Max NPV","Distance"), pch = c(16, 17),
       col = c("black", "lightcoral", rgb(0.5, 0.5, 0.5)), lty = c(1, 2, 2), lwd = c(1, 1, 1))
dev.off()

# plotting Max PPV, c-stat, and dist between
png(filename = file.path(tmpBaseDir, 'Max PPV Nearest.png'))
plot(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2], xlim = c(50, 500),
    ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
    main = "Maximum c-statistic, Maximum PPV, \nand Distance Between", pch = c(16))
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2])

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 8], col = rgb(0.1, 0.7, 1), pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 8], col = rgb(0.1, 0.7, 1), lty = 2)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 20], pch = 16, col = rgb(0.5, 0.5, 0.5))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 20]))],
        na.omit(variedPatientNumbersData[, 20]), col = rgb(0.5, 0.5, 0.5), lwd = 1, lty = 2)

legend("topright", legend = c("Max c-statistic", "Max PPV", "Distance"), pch = c(16, 17),
       col = c("black", rgb(0.1, 0.7, 1), rgb(0.5, 0.5, 0.5)), lty = c(1, 2, 2), lwd = c(1, 1, 1))
dev.off()

# plotting Max Spec, c-stat, and dist between
png(filename = file.path(tmpBaseDir, 'Max Spec Nearest.png'))
plot(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2], xlim = c(50, 500),
    ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
    main = "Maximum c-statistic, Maximum Spec, \nand Distance Between", pch = c(16))
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2])

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 9], col = rgb(0, 1, 0.5), pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 9], col = rgb(0, 1, 0.5), lty = 2)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 21], pch = 16, col = rgb(0.5, 0.5, 0.5))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 21]))],
        na.omit(variedPatientNumbersData[, 21]), col = rgb(0.5, 0.5, 0.5), lwd = 1, lty = 2)

legend("topright", legend = c("Max c-statistic", "Max Spec", "Distance"), pch = c(16, 17),
       col = c("black", rgb(0, 1, 0.5), rgb(0.5, 0.5, 0.5)), lty = c(1, 2, 2), lwd = c(1, 1, 1))
dev.off()

# plotting Max Sens, c-stat, and dist between
png(filename = file.path(tmpBaseDir, 'Max Sens Nearest.png'))
plot(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2], xlim = c(50, 500),
    ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
    main = "Maximum c-statistic, Maximum Sens, \nand Distance Between", pch = c(16))
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2])

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 10], col = rgb(1, 0.5, 1), pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 10], col = rgb(1, 0.5, 1), lty = 2)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 22], pch = 16, col = rgb(0.5, 0.5, 0.5))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 22]))],
        na.omit(variedPatientNumbersData[, 22]), col = rgb(0.5, 0.5, 0.5), lwd = 1, lty = 2)

legend("topright", legend = c("Max c-statistic", "Max Sens", "Distance"), pch = c(16, 17),
       col = c("black", rgb(1, 0.5, 1), rgb(0.5, 0.5, 0.5)), lty = c(1, 2, 2), lwd = c(1, 1, 1))
dev.off()

# plotting Max NPV NonResp, c-stat, and dist between
png(filename = file.path(tmpBaseDir, 'Max NPV NonResp Nearest.png'))
plot(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2], xlim = c(50, 500),
    ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
    main = "Maximum c-statistic, Maximum NPV NonResp, \nand Distance Between", pch = c(16))
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2])

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 15], col = "lightcoral", pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 15], col = "lightcoral", lty = 2)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 23], pch = 16, col = rgb(0.5, 0.5, 0.5))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 23]))],
        na.omit(variedPatientNumbersData[, 23]), col = rgb(0.5, 0.5, 0.5), lwd = 1, lty = 2)

legend("topright", legend = c("Max c-statistic", "Max NPV NonResp", "Distance"), pch = c(16, 17),
       col = c("black", col = "lightcoral", rgb(0.5, 0.5, 0.5)), lty = c(1, 2, 2), lwd = c(1, 1, 1))
dev.off()

# plotting Max PPV ProgFree, c-stat, and dist between
png(filename = file.path(tmpBaseDir, 'Max PPV ProgFree Nearest.png'))
plot(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2], xlim = c(50, 500),
    ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
    main = "Maximum c-statistic, Maximum PPV ProgFree, \nand Distance Between", pch = c(16))
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2])

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 16], col = rgb(0.1, 0.7, 1), pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 16], col = rgb(0.1, 0.7, 1), lty = 2)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 24], pch = 16, col = rgb(0.5, 0.5, 0.5))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 24]))],
        na.omit(variedPatientNumbersData[, 24]), col = rgb(0.5, 0.5, 0.5), lwd = 1, lty = 2)

legend("topright", legend = c("Max c-statistic", "Max PPV ProgFree", "Distance"), pch = c(16, 17),
       col = c("black", rgb(0.1, 0.7, 1), rgb(0.5, 0.5, 0.5)), lty = c(1, 2, 2), lwd = c(1, 1, 1))
dev.off()

# plotting Max Spec NonResp, c-stat, and dist between
png(filename = file.path(tmpBaseDir, 'Max Spec NonResp Nearest.png'))
plot(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2], xlim = c(50, 500),
    ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
    main = "Maximum c-statistic, Maximum Spec NonResp, \nand Distance Between", pch = c(16))
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2])

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 17], col = rgb(0, 1, 0.5), pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 17], col = rgb(0, 1, 0.5), lty = 2)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 25], pch = 16, col = rgb(0.5, 0.5, 0.5))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 25]))],
        na.omit(variedPatientNumbersData[, 25]), col = rgb(0.5, 0.5, 0.5), lwd = 1, lty = 2)

legend("topright", legend = c("Max c-statistic", "Max Spec NonResp", "Distance"), pch = c(16, 17),
       col = c("black", rgb(0, 1, 0.5), rgb(0.5, 0.5, 0.5)), lty = c(1, 2, 2), lwd = c(1, 1, 1))
dev.off()

# plotting Max Sens ProgFree, c-stat, and dist between
png(filename = file.path(tmpBaseDir, 'Max Sens ProgFree Nearest.png'))
plot(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2], xlim = c(50, 500),
    ylim = c(0, 1), xlab = "Number of Simulated Patients", ylab = "Value of Metric",
    main = "Maximum c-statistic, Maximum Sens ProgFree, \nand Distance Between", pch = c(16))
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 2])

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 18], col = rgb(1, 0.5, 1), pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 18], col = rgb(1, 0.5, 1), lty = 2)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 26], pch = 16, col = rgb(0.5, 0.5, 0.5))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 26]))],
        na.omit(variedPatientNumbersData[, 26]), col = rgb(0.5, 0.5, 0.5), lwd = 1, lty = 2)

legend("topright", legend = c("Max c-statistic", "Max Sens ProgFree", "Distance"), pch = c(16, 17),
       col = c("black", rgb(1, 0.5, 1), rgb(0.5, 0.5, 0.5)), lty = c(1, 2, 2), lwd = c(1, 1, 1))
dev.off()

# plotting Number of patients in each group PR,SD,PD
png(filename = file.path(tmpBaseDir, 'Group Dist.png'))
plot(variedPatientNumbersData[, 1], variedPatientNumbersData[, 27], xlim = c(50, 500),
    ylim = c(0, 500), xlab = "Number of Simulated Patients", ylab = "Number of Patients in Group",
    main = "Number of Patients in PR, SD, PD at the Maximum c-statistic", pch = c(16), col = rgb(0, 1, 0))
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 27], col = rgb(0, 1, 0))

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 28], pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 28], lty = 2)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 29], pch = 16, col = rgb(1, 0, 0))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 29]))],
        na.omit(variedPatientNumbersData[, 29]), col = rgb(1, 0, 0), lwd = 1, lty = 2)

legend("right", legend = c("PR", "SD", "PD"), pch = c(16, 17),
       col = c(rgb(0, 1, 0), "black", rgb(1, 0, 0)), lty = c(1, 2, 2), lwd = c(1, 1, 1))
dev.off()

# plotting Number of patients in each group PR,SD,PD scaled
png(filename = file.path(tmpBaseDir, 'Group Dist Scaled.png'))
plot(variedPatientNumbersData[, 1], variedPatientNumbersData[, 27] / variedPatientNumbersData[, 1], xlim = c(50, 500),
    ylim = c(0,1), xlab = "Number of Simulated Patients", ylab = "Percentage of Patients in Group",
    main = "Percentage of Patients in PR, SD, PD at the Maximum c-statistic", pch = c(16), col = rgb(0, 1, 0))
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 27]/variedPatientNumbersData[,1], col = rgb(0, 1, 0))

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 28] / variedPatientNumbersData[, 1], pch = 17)
lines(variedPatientNumbersData[, 1], variedPatientNumbersData[, 28] / variedPatientNumbersData[, 1], lty = 2)

points(variedPatientNumbersData[, 1], variedPatientNumbersData[, 29]/variedPatientNumbersData[,1], pch = 16, col = rgb(1, 0, 0))
lines(variedPatientNumbersData[, 1][which(complete.cases(variedPatientNumbersData[, 29] / variedPatientNumbersData[, 1]))],
        na.omit(variedPatientNumbersData[, 29] / variedPatientNumbersData[, 1]), col = rgb(1, 0, 0), lwd = 1, lty = 2)

legend("right", legend = c("PR", "SD", "PD"), pch = c(16, 17),
       col = c(rgb(0, 1, 0), "black", rgb(1, 0, 0)), lty = c(1, 2, 2), lwd = c(1, 1, 1))
dev.off()