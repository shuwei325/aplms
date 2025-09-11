#' @docType data
#' @title Global Annual Mean Surface Air Temperature Change
#' @usage data(temperature)
#' @description Land-ocean temperature index from 1880--2021 (with base period 1951-1980).
#' @format Time series data.
#' @source NASA/GISS/GISTEMP
#' @keywords datasets
#' @references
#' \href{https://data.giss.nasa.gov/gistemp/graphs/graph_data/Global_Mean_Estimates_based_on_Land_and_Ocean_Data/graph.txt}{Land-Ocean Temperature Index}.
#' @examples
#' data(temperature)
#'
"temperature"

#' @docType data
#' @title Respiratory diseases hospitalization Dataset
#' @description This dataset consists of respiratory diseases hospitalization in Sorocaba, São Paulo, Brazil.
#' The details of the statistical modeling using the APLMS-AR(p) approach can be found in Chou-Chen, et al. (2024). \doi{10.1007/s00362-024-01590-w}.
#' The hospitalization data of respiratory diseases in Sorocaba city, São Paulo, Brazil are obtained from the Hospital Information System of Brazil’s Unified National Health System (\href{https://datasus.saude.gov.br/}{SIH-SUS}), and the climatic and pollution data are provided by the \href{https://qualar.cetesb.sp.gov.br/qualar/}{QUALAR system}.
#' @format The "data" slot is a data frame with 932 weekly data on the following 29 variables.
#' \describe{
#' \item{date, year, epi.week, tdate}{Date, year, epidemiologic weeks, and time index.}
#' \item{y}{Respiratory hospitalization count.}
#' \item{MP10_max, MP10_min, MP10_avg}{Maximum, minimum and average of \eqn{MP10}.}
#' \item{NO_max,NO_min,NO_avg}{Maximum, minimum and average of \eqn{NO}.}
#' \item{NO2_max, NO2_min, NO2_avg}{Maximum, minimum and average of \eqn{NO_2}.}
#' \item{NOx_max, NOx_min, NOx_avg}{Maximum, minimum and average of \eqn{NO_x}.}
#' \item{O3_max, O3_min, O3_avg}{Maximum, minimum and average of \eqn{O_3}.}
#' \item{TEMP_max, TEMP_min, TEMP_avg}{Maximum, minimum and average of temperature.}
#' \item{RH_max, RH_min, RH_avg}{Maximum, minimum and average of relative humidity.}
#' \item{ampl_max, ampl_min, ampl_avg}{Maximum, minimum and average of daily temperature amplitude.}
#' }
#' @usage data(hospitalization)
#' @references Chou-Chen, S.W., Oliveira, R.A., Raicher, I. et al. (2024) Additive partial linear models with autoregressive symmetric errors and its application to the hospitalizations for respiratory diseases. Stat Papers 65, 5145–5166. \doi{10.1007/s00362-024-01590-w}
#' @keywords datasets
#' @examples
#' data(hospitalization)
#' head(hospitalization)
"hospitalization"


