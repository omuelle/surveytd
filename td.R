library(glmnet)
library(tempdisagg)
library(timeseriesdb)
library(data.table)
library(tstools)
library(purrr)
library(openxlsx)

source("utils.R")
source("td_func.R")

questions <- c(
  "q_ql_ass_order_blog",
  "q_ql_ass_stock_fin",
  "q_ql_ass_stock_intermed",
  "q_ql_chg_order_blog_pmppm",
  "q_ql_chg_order_in_pmppm",
  "q_ql_chg_order_in_pmpym",
  "q_ql_chg_prod_pmppm",
  "q_ql_chg_prod_pmpym",
  "q_ql_chg_stock_fin_pmppm",
  "q_ql_exp_chg_order_in_n3m",
  "q_ql_exp_chg_prod_n3m")

keys <- paste0("ch.kof.inu.ng08.fx.", questions, ".balance.d11")
con <- createConObj(dbname = "kofdb", dbhost = "archivedb.kof.ethz.ch", passwd = .rs.askForPassword("KOFDB password"))
tsl <- timeseriesdb::readTimeSeries(con, keys)
names(tsl) <- questions

# alternatively, load from csv
# tsl <- tstools::read_ts("tsl.csv", "csv")

tsl <- map(tsl, ~window(.x, start = 1990))
dates <- seq(ymDate(2000,1), ymDate(2019,4), by = "month")

t_start <- Sys.time()
# Real-time analysis.
res <- map(names(tsl), function(var) {
  res <- map(dates, function(date) {
    
    tslm <- map(tsl, ~window(.x, end=dateToTime(date)))
    tslq <- map(tslm, toQuarterlySeries)
    yq <- tslq[[var]]
    ym <- tslm[[var]]
    tslm <- tslm[setdiff(names(tsl), var)]
    tslq <- tslq[setdiff(names(tsl), var)]
    tslm_fullq <- map(tslm, ~window(.x, end=c(year(date), month(frequencyDate(date,4))+2), extend = T))
    tslm_fullq <- map(tslm_fullq, ~{ .x[is.na(.x)] <- 0; .x })
    
    # Disagg (comment in desired)
    
    # # Old chowlin (using tempdisagg package)
    # beta <- lassoCV(yq, tslq)
    # sel_vars <- names(tslq)[beta[2:length(beta)] != 0]
    # ytd <- tempDisagg(yq, tslm_fullq[sel_vars])
    # ytd <- window(ytd, end = dateToTime(date))
    # list(ytd=ytd, sel_vars=sel_vars)
    
    # # Standard chowlin
    # X <- do.call(cbind, tslq)
    # beta <- lassoBIC(yq, X)$beta
    # rho <- estimAR1(yq, X, beta)
    # ytd <- chowlin(yq, tslm_fullq, beta, rho)
    # ytd <- window(ytd, end = dateToTime(date))
    # list(ytd=ytd, beta=beta, rho=rho)
    
    # LassoAR1 chowlin
    res <- lassoAR1(yq, tslq)
    ytd <- chowlin(yq, tslm_fullq, res$beta, res$rho)
    ytd <- window(ytd, end = dateToTime(date))
    list(ytd=ytd, beta=res$beta, rho=res$rho, iter=res$iter)
    
    # # EM
    # ytd <- em(yq, tslm, nr_comp=1)
    # list(ytd=ytd)
    
    # Hot deck
    # ytd <- na.locf(window(toMonthlySeries(yq,T), end=dateToTime(date)))
    # list(ytd=ytd)
  })
  names(res) <- dates
  res
})
names(res) <- names(tsl)
Sys.time() - t_start

# Assign name and store
ver <- "ar1_chowlin"
assign(paste0("res_", ver), get("res"))
save(list = paste0("res_", ver), file = paste0("results/res_", ver, ".RData"))


# Aggregate real-time disaggregated series and compute stats.
computeRTStat <- function(res) {
  stat_rt <- map(names(res), function(var) {
    rt <- ts(map_dbl(res[[var]], ~tail(.x$ytd,1)), start = dateToTime(names(res[[var]])[1]), frequency = 12)
    list(rt = rt, rmse=rmse(rt, tsl[[var]]), cor=corTS(rt, tsl[[var]]))
  })
  names(stat_rt) <- names(res)
  stat_rt
}

# Create RRMSE table
dt <- data.table("Hot Deck"=map_dbl(computeRTStat(res_hd), ~.x$rmse),
                 "Standard Chowlin"=map_dbl(computeRTStat(res_std_chowlin), ~.x$rmse),
                 "LassoAR1 Chowlin"=map_dbl(computeRTStat(res_ar1_chowlin), ~.x$rmse),
                 "EM"=map_dbl(computeRTStat(res_em), ~.x$rmse))

dt <- do.call(data.table, c(list(Variable=names(res_hd)), map(dt, ~round(.x/dt[["Hot Deck"]],2))))
write.xlsx(dt, "rrmse.xlsx")



