# Nayak, S., Lee, D., Patel-Hett, S., Pittman, D., Martin, S., Heatherington, A., ¡¦ Hua, F. (2015). Using a Systems Pharmacology Model of the Blood Coagulation Network to Predict the Effects of Various Therapies on Biomarkers. CPT: Pharmacometrics & Systems Pharmacology, 4(7), 396???405. doi:10.1002/psp4.50

# Woo Seok Byun, Soo Yeun Lee, Hae Jin Yoon, Dong Jin Shin
# This project completed for 2019 Bioinformatics Class in Handong Global University.
# We modified vesion from MatLab to R to achieve more readability and writability.
# Unfortunately, compared to MatLab, R internals are single-threaded, making its performance much slower.
# We did not include the simulated annealing code block, as it was home-made and not optimal.

source("library.R")

options(scipen = 999)

#out <- data.frame(ode(y = state, times = times, func = Blood, parms = parameters))
#VIIa = c(4000, 2000, 1000, 500, 250, 125, 0.1)
#Xa = c(1, 0.5, 0.25, 0.13, 0.06, 0.03, 0)

draw("VIIa", c(4000, 2000, 1000, 500, 250, 125, 0.1))
draw("Xa", c(1, 0.5, 0.25, 0.13, 0.06, 0.03, 0))