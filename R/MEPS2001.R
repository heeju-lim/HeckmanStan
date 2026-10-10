#' @title MEPS 2001: Ambulatory Expenditures Data
#'
#' @description
#' This dataset is an extract from the 2001 Medical Expenditure Panel Survey (MEPS), providing information on ambulatory expenditures and various demographic and health-related variables. It has been used for illustrative examples by Cameron and Trivedi (2009, Chapter 16).
#'
#' @usage data(MEPS2001)
#' @format A data frame with 3,328 observations on the following 22 variables.
#' \describe{
#'   \item{educ}{Education status}
#'   \item{age}{Age}
#'   \item{income}{Income}
#'   \item{female}{Gender}
#'   \item{vgood}{Self-reported health status, very good}
#'   \item{good}{Self-reported health status, good}
#'   \item{hospexp}{Hospital expenditures}
#'   \item{totchr}{Total number of chronic diseases}
#'   \item{ffs}{Family support}
#'   \item{dhospexp}{Dummy variable for hospital expenditures}
#'   \item{age2}{Age squared}
#'   \item{agefem}{Interaction between age and gender}
#'   \item{fairpoor}{Self-reported health status, fair or poor}
#'   \item{year01}{Year of survey}
#'   \item{instype}{Type of insurance}
#'   \item{ambexp}{Ambulatory expenditures}
#'   \item{lambexp}{Log of ambulatory expenditures}
#'   \item{blhisp}{Ethnicity}
#'   \item{instype_s1}{Insurance type, version 1}
#'   \item{dambexp}{Dummy variable for ambulatory expenditures}
#'   \item{lnambx}{Log-transformed ambulatory expenditures}
#'   \item{ins}{Insurance status}
#' }
#'
#' @source
#' 2001 Medical Expenditure Panel Survey by the Agency for Healthcare Research and Quality.
#'
#' @references
#' Cameron, C.A. and Trivedi, P.K. (2009). *Microeconometrics Using Stata*. College Station, TX: Stata Press.
#' @docType data
#' @name MEPS2001
#' @examples
#'
#' data(MEPS2001)
#' head(MEPS2001)
#'
NULL
