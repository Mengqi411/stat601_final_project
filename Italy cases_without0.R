###################                           ######################
###################      GROUP 3 project      ######################
###################                           ######################
###################                           ######################

rm(list = ls())

library(stringr)
library(car)
library(tidyverse)
library(Hmisc)
library(gvlma)
library(xlsx)
# read data
case <- read.csv(file = "covid19_italy_region.csv", header = TRUE)
mobility <- read.csv(file = "2020_IT_Region_Mobility_Report.csv", header = TRUE)


# modify date col
case$Time <- gsub(x = case$Date, pattern = ".*T(.*)$", replacement = "\\1")
case$Date <- as.Date(gsub(x = case$Date, pattern = "(.*)T.*", replacement = "\\1"))
case <- relocate(.data = case, Time, .after = Date)
mobility$date <- as.Date(mobility$date)


# match the region names
case$RegionName[case$RegionName == "Valle d'Aosta"] <- "Aosta"
case$RegionName[case$RegionName == "Puglia"] <- "Apulia"
case$RegionName[case$RegionName == "Friuli Venezia Giulia"] <- "Friuli-Venezia Giulia"
case$RegionName[case$RegionName == "Lombardia"] <- "Lombardy"
case$RegionName[case$RegionName == "Toscana"] <- "Tuscany"
case$RegionName[case$RegionName == "Piemonte"] <- "Piedmont"
case$RegionName[case$RegionName == "Sardegna"] <- "Sardinia"
case$RegionName[case$RegionName == "Sicilia"] <- "Sicily"
case$RegionName[case$RegionName == "P.A. Trento"] <- "Trentino-South Tyrol"
case$NewPositiveCases[case$RegionName == "Trentino-South Tyrol"] <-
    case$NewPositiveCases[case$RegionName == "Trentino-South Tyrol"] +
    case$NewPositiveCases[case$RegionName == "P.A. Bolzano"]

# get region names
region_in_mobility <- unique(mobility$sub_region_1)[-1]
region_in_case <- unique(case$RegionName)
region_in_case <- region_in_case[order(region_in_case)]
region <- intersect(region_in_mobility, region_in_case)

result <- data.frame(region = region)

start <- as.Date("2020-03-20")
end <- as.Date("2020-11-20")

lambda <- c()
model <- list()
regression <- list()
vif <- data.frame(region = region)
corr <- list()
# layout(matrix(1:20, nrow = 5, ncol = 4, byrow = TRUE))
S.W_test <- c()
equal_variance <- c()
r_f <- data.frame(region = region, r.squared = NA, f.statistic.p.value = NA)
death_rate <- 0.034


# regression by regions
for (i in 1:length(region)) {
    mobility_by_region <- mobility[mobility$sub_region_1 == region[i] &
        mobility$sub_region_2 == "" &
        mobility$date >= start &
        mobility$date <= end, ]
    mobility_by_date <- mobility_by_region[order(mobility_by_region$date), ]
    case_by_region <- case[case$RegionName == region[i] & case$Date >= (start + 7) & case$Date <= (end + 14), ]
    case_by_date <- case_by_region[order(case_by_region$Date), ]
    for (j in 1:length(case_by_date$Date - 7)) {
        case_by_date$NewPositiveCases[j] <- mean(case_by_date$NewPositiveCases[j:(j + 6)])
    }
    case_by_date <- case_by_date[case_by_date$Date >= (start + 7) & 
        case_by_date$Date <= (end + 7), ]
    
    data_selected <- data.frame(
        NewPositiveCases = case_by_date$NewPositiveCases,
        retail_and_recreation = mobility_by_date$retail_and_recreation_percent_change_from_baseline,
        grocery_and_pharmacy = mobility_by_date$grocery_and_pharmacy_percent_change_from_baseline,
        parks = mobility_by_date$parks_percent_change_from_baseline,
        # transit_stations = mobility_by_date$transit_stations_percent_change_from_baseline,
        workplaces = mobility_by_date$workplaces_percent_change_from_baseline
        # residential = mobility_by_date$residential_percent_change_from_baseline
    )
    data_selected <- data_selected[data_selected$NewPositiveCases > 0, ] # remove non-positives
    bc <- boxCox(
        data = data_selected, lambda = seq(-5, 5, 0.01),
        data_selected$NewPositiveCases ~
            data_selected$retail_and_recreation +
            data_selected$grocery_and_pharmacy +
            data_selected$parks +
            data_selected$workplaces,
            plotit = FALSE,
            interp = TRUE
    )

    lambda[i] <- bc$x[which.max(bc$y)]
    names(lambda)[i] <- region[i]
    data_selected$NewPositiveCases_bc <- (data_selected$NewPositiveCases^lambda[i] - 1) / lambda[i]

    model[[region[i]]] <- lm(
        data = data_selected,
        NewPositiveCases_bc ~ retail_and_recreation + grocery_and_pharmacy + parks + workplaces
    )
    regression[[region[i]]] <- summary(model[[region[i]]])


    result[result$region == region[i], 2:length(model[[region[i]]]$coefficients)] <-
        model[[region[i]]]$coefficients[2:length(model[[region[i]]]$coefficients)]

    colnames(result)[2:length(model[[region[i]]]$coefficients)] <-
        names(model[[region[i]]]$coefficients)[2:length(model[[region[i]]]$coefficients)]

    vif[vif$region == region[i], 2:length(model[[region[i]]]$coefficients)] <-
        vif(model[[region[i]]])

    colnames(vif)[2:length(model[[region[i]]]$coefficients)] <-
        names(model[[region[i]]]$coefficients)[2:length(model[[region[i]]]$coefficients)]

    ############################ Diagnostics #######################################

    # Pearson's r
    corr[[region[i]]] <- rcorr(as.matrix(data_selected), type = "pearson")[[1]]

    # QQplot
    # plot(model[[region[i]]], which = 1)

    # Shapiro-Wilk normality test
    S.W_test[i] <- shapiro.test(model[[region[[i]]]]$residuals)$p.value
    names(S.W_test)[i] <- region[i]

    # equal variance
    if (gvlma(model[[region[i]]])[length(gvlma(model[[region[i]]]))]$GlobalTest$DirectionalStat4$Decision == 0) {
        equal_variance[i] <- FALSE
    } else {
        equal_variance[i] <- TRUE
    }
    names(equal_variance)[i] <- region[i]
    # r-squared and p-value
    r_f[r_f$region == region[[i]], 2] <- regression[[region[i]]]$r.squared
    f <- regression[[region[i]]]$fstatistic
    r_f[r_f$region == region[[i]], 3] <- pf(f[1], f[2], f[3], lower.tail = FALSE)
}
############################# Output ###########################################
# write.csv(x = result, file = "Italy_regression.csv")
# write.csv(x = vif, file = "Italy_vif.csv")

# cat(sep = "", "result:\n")
# print(result)
# cat("\n")

# cat(sep = "", "vif:\n")
# print(vif)
# cat("\n")

# # cat(sep = "", "corr:\n")
# # print(corr)
# # cat("\n")

# cat(sep = "", "S.W_test:\n")
# print(S.W_test)
# cat("\n")

# cat(sep = "", "Heteroscedasticity: TRUE indicates Heteroscedasticity exists\n")
# print(equal_variance)
# cat("\n")

# cat(sep = "", "r-squared and p-value:\n")
# print(r_f)
# cat("\n")

# write.xlsx(x = cbind(result,equal_variance), file = "03-20-11-20Italy_coefficient_without0.xlsx")

