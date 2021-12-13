library(glmnet)
library(tempdisagg)
library(timeseriesdb)
library(data.table)
library(tstools)
library(purrr)
library(openxlsx)
library(parallel)
library(nlme)

source("utils.R")
source("td_func.R")

monthly_questions <- c(
  "q_ql_ass_order_blog",
  "q_ql_ass_stock_fin",
  "q_ql_chg_order_blog_pmppm",
  "q_ql_chg_order_in_pmppm",
  "q_ql_chg_order_in_pmpym",
  "q_ql_chg_prod_pmppm",
  "q_ql_chg_prod_pmpym",
  "q_ql_chg_stock_fin_pmppm",
  "q_ql_exp_chg_order_in_n3m",
  "q_ql_exp_chg_prod_n3m",
  "q_ql_exp_chg_pur_intermed_n3m")
 
quarterly_questions <- c(
  "q_ql_chg_tech_cap_p3m",
  "q_ql_ass_tech_cap",
  "q_ql_chg_profit_p3m",
  "q_qn_capu_p3m",
  "ql_exp_chg_price_pur_n3m")

ql_items <- c("balance","share_neg","share_pos")
qn_items <- "value"

# reg <- paste0("^ch\\.kof\\.inu\\.ng08\\.fx\\.(sector_kof\\.[^.]*\\.)(", paste0(questions, collapse = "|"), ")\\.balance\\.d11$")
# keys <- dbGetQuery(con, sprintf("select ts_key from timeseries.timeseries_main where ts_key ~ '%s'", reg))$ts_key

questions <- monthly_questions
question_item <- c(kofbts::crossPaste(questions[grep("^q_ql", questions)], ql_items, sep = "."),
                   kofbts::crossPaste(questions[grep("^q_qn", questions)], qn_items, sep = "."))
keys <- paste0("ch.kof.inu.ng08.fx.", question_item, ".d11")

con <- createConObj(dbname = "kofdb", dbhost = "archivedb.kof.ethz.ch", passwd = .rs.askForPassword("KOFDB password"))
tsl <- timeseriesdb::readTimeSeries(con, keys)
names(tsl) <- question_item
tsl <- map(tsl, ~window(.x, end = c(2019,12)))

# alternatively, load from csv/RData
# tsl <- tstools::read_ts("tsl.csv", "csv")
# load("tsl.RData")

# # Add differences
# tsl_ind <- unlist(map(c(1,3,12), function(d) {
#   tsl <- map(tsl, ~diff(.x, d))
#   names(tsl) <- paste0(names(tsl), ".d", d, "m")
#   tsl
# }), recursive=F)
# tsl_ind <- c(tsl, tsl_ind)
# 
# # Add lags
# tsl_ind <- unlist(map(c(1,3,6,12), function(l) {
#   tsl <- map(tsl, ~lag(.x, -l))
#   names(tsl) <- paste0(names(tsl), ".l", l, "m")
#   tsl
# }), recursive=F)
# tsl_ind <- c(tsl, tsl_ind)
# 

# No transformations
tsl_ind <- tsl

# Impute quarterly questions
questions <- quarterly_questions
question_item <- c(kofbts::crossPaste(questions[grep("^q_ql", questions)], ql_items, sep = "."),
                   kofbts::crossPaste(questions[grep("^q_qn", questions)], qn_items, sep = "."))
keys <- paste0("ch.kof.inu.ng08.fx.", question_item, ".d11")
tsl <- timeseriesdb::readTimeSeries(con, keys)
names(tsl) <- question_item
tsl <- map(tsl, ~window(.x, end = c(2019,4)))

vars <- names(tsl)
end_date <- ymDate(2019,12)
mov_win <- 20*12
  
# Real-time analysis.
t_start <- Sys.time()
methods <- c("chowlin","lasso_chowlin","lassoar1_chowlin","em","locf","spline","ols")
# methods <- "em"
# methods <- c("locf", "em", "gls")
# methods <- "chowlin"
res <- map(methods, function(method) {
  res <- map(vars, function(var) {
    # var <- "q_ql_chg_stock_fin_pmppm.share_pos"
     
    cat(paste0(method, ": ", var), "\n")
    var_chunks <- strsplit(var, "\\.")[[1]]
    tslm <- tsl_ind[grepl(if(var_chunks[2]=="value") "balance" else var_chunks[2], names(tsl_ind)) & !grepl(var_chunks[1], names(tsl_ind))]
    tslq <- map(tslm, toQuarterlySeries)
    ym <- tsl[[var]]
    yq <- toQuarterlySeries(ym)
    min_date <- max(c(as.Date(map_chr(tslq, ~as.character(startDate(.x)))), startDate(yq)))
    dates <- seq(addMonths(min_date, mov_win-1), end_date, by = "month")
    mov_win <- mov_win
    
    f <- function(date) {
      # date <- as.Date("1999-01-01")
      end_time <- dateToTime(date)
      start_time <- dateToTime(frequencyDate(addMonths(date,-mov_win+1),4))
      tslm <- map(tslm, ~window(.x, start=start_time, end=end_time))
      tslq <- map(tslq, ~window(.x, start=start_time, end=end_time))
      ym <- window(ym, start=start_time, end=end_time)
      yq <- window(yq, start=start_time, end=end_time)
  
      # Make full quarter series    
      tslm_fullq <- map(tslm, ~window(.x, end=c(year(date), month(frequencyDate(date,4))+2), extend = T))
      tslm_fullq <- map(tslm_fullq, ~{ .x[is.na(.x)] <- 0; .x })
      
      if(method=="chowlin") {
        ytd <- tempDisagg(yq, tslm_fullq)
        ytd <- window(ytd, end = end_time)
        list(ytd=ytd)
      }
      else if(method=="lasso_chowlin") {
        beta <- lassoCV(yq, tslq)
        sel_vars <- names(tslq)[beta[2:length(beta)] != 0]
        ytd <- tempDisagg(yq, tslm_fullq[sel_vars])
        ytd <- window(ytd, end = end_time)
        list(ytd=ytd, sel_vars=sel_vars)
      }
      else if(method=="lassoar1_chowlin") {
        res <- lassoAR1(yq, tslq)
        ytd <- chowlin(yq, tslm_fullq, res$beta, res$rho)
        ytd <- window(ytd, end = end_time)
        list(ytd=ytd, beta=res$beta, rho=res$rho, iter=res$iter)
      }
      else if(method=="em") {
        ytd <- em(window(toMonthlySeries(yq,T), end=end_time), tslm, nr_comp=2)
        list(ytd=ytd)
      }
      else if(method=="lasso_em") {
        beta <- lassoBIC(yq, tslq)
        sel_vars <- names(tslq)[beta[2:length(beta)] != 0]
        ytd <- em(yq, tslm[sel_vars], nr_comp=2)
        list(ytd=ytd)
      }
      else if(method=="locf") {
        ytd <- na.locf(window(toMonthlySeries(yq,T), end=end_time))
        list(ytd=ytd)
      }
      else if(method=="spline") {
        ytd <- na.spline(window(toMonthlySeries(yq,T), end=end_time), method="natural")
        ytd <- ts(ytd, end=end_time, frequency = 12)
        list(ytd=ytd)
      }
      else if(method=="arma") {
        mod <- auto.arima(yq, 0, 0, ic = "bic")
        yq <- ts(c(yq, predict(mod, h=1)$pred), start=start(yq), frequency = 4)
        ytd <- window(na.approx(toMonthlySeries(yq,F)), end=end_time)
        list(ytd=ytd)
      }
      else if(method=="ols") {
        X <- do.call(cbind, tslq)
        mod <- gls(yq ~ X)
        X <- do.call(cbind, tslm)
        ym_estim <- ts(X %*% mod$coefficients[-1] + mod$coefficients[1], end=end_time, frequency = 12)
        ytd <- window(toMonthlySeries(yq,T), end=end_time)
        ytd[is.na(ytd)] <- ym_estim[is.na(ytd)]
        list(ytd=ytd)
      }
      else if(method=="step") {
        X <- do.call(cbind, tslq)
        mod <- step(lm(yq ~ X), trace = F, k=log(length(yq)))
        ytd <- ts(if(length(mod$coefficients) > 1) X %*% mod$coefficients[-1] + mod$coefficients[1] else mod$coefficients[1], end=end_time, frequency = 12)
        list(ytd=ytd)
      }
    }

    cluster <- makeCluster(4, type = "PSOCK")
    clusterEvalQ(cluster, source("source_all.R"))
    res <- clusterMap(cluster, f, dates, SIMPLIFY = F)
    stopCluster(cluster)
    
    # if(var=="q_ql_ass_order_blog.share_pos") browser()
    # dates <- tail(dates, 2)
    # res <- map(dates, f)
    
    names(res) <- dates
    res
  })
  names(res) <- vars
  res
})
names(res) <- methods
Sys.time() - t_start

res <- map2(res, names(res), function(ver, ver_name) {
  cat(ver_name)
  comp <- map(monthly_questions, function(question) {
    cat(question)
    bal <- ver[[grep(paste0(question, ".balance"),names(ver))]]
    neg <- ver[[grep(paste0(question, ".share_neg"),names(ver))]]
    pos <- ver[[grep(paste0(question, ".share_pos"),names(ver))]]
    if(any(map_lgl(pos, ~frequency(.x$ytd)==4))) browser()
    if(any(map_lgl(neg, ~frequency(.x$ytd)==4))) browser()
    diff <- map2(pos, neg, ~list(ytd = .x$ytd - .y$ytd))
    ll <- list(bal, diff, pos, neg)
    names(ll) <- paste0(question, c(".balance",".balance.composite",".share_pos",".share_neg"))
    ll
  })
  unlist(comp, recursive = F)
})

# Assign name and store
ver <- "may20q"
assign(paste0("res_", ver), get("res"))
save(list = paste0("res_", ver), file = paste0("results/res_", ver, ".RData"))

# Compute stats table for result
stats <- map(res, ~computeStat(.x, ymDate(1987,3)))
save(stats, file = "results/stats.RData")

library(openxlsx)
wb <- createWorkbook()
types <- c(rt_rmse_m23="RT RMSE M2+M3", rt_rmse_m2="RT RMSE M2", rt_rmse_m3="RT RMSE M3", rt_rrmse_m2_vs_m3="RT RRMSE M2 vs M3",
           post_rmse_m23="POST RMSE M2+M3", post_rmse_m2="POST RMSE M2", post_rmse_m3="POST RMSE M3", post_rrmse_m2_vs_m3="POST RRMSE M2 vs M3")
# export_vars <- vars
export_vars <- names(stats$locf)
export_methods <- c("locf","lassoar1_chowlin","lasso_chowlin","chowlin","em","ols","spline")
for(type in names(types)) {
  dt <- do.call(data.table, c(list(var = export_vars), map(stats[export_methods], ~map_dbl(.x[export_vars], ~.x[[type]] %>% round(2)))))
  sheet_name <- types[type]
  dtm <- addColumnAggregate(dt, function(x) round(mean(x),2), "Mean")
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = dtm)
  dt <- do.call(data.table, c(list(var = export_vars), map(stats[export_methods], ~map_dbl(.x[export_vars], ~.x[[type]]))))
  rel_to <- dt$locf
  for(col in setdiff(names(dt), "var")) dt[,(col) := (get(col) / rel_to) %>% round(2)]
  sheet_name <- sub("RMSE", "RRMSE", sheet_name)
  dtm <- addColumnAggregate(dt, function(x) round(mean(x),2), "Mean")
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = dtm)
}
saveWorkbook(wb, "exports/stats_1987-3.xlsx", overwrite = TRUE)

export_td(res, "may20q")
tstools::write_ts(tsl, fname = paste0("exports/may20q/actual"), format = "xlsx", wide = T)


# 
# dt_std <- dt
# dt_diffs <- dt
# dt_lags <- dt
# 
# names(dt_std) <- c("var","chowlin","lasso chowlin","lassoAR1 chowlin", "em", "locf")
# names(dt_diffs) <- c("var","lasso chowlin diffs","lassoAR1 chowlin diffs", "em diffs", "locf")
# names(dt_lags) <- c("var","lasso chowlin lags","lassoAR1 chowlin lags", "em lags", "locf")
# dt_std$locf <- NULL
# dt_diffs$locf <- NULL
# 
# xx <- dt_std[dt_diffs, on="var"]
# xx <- xx[dt_lags, on="var"]
# xx <- addColumnAggregate(xx, mean, "Mean")
# 
# apply(xx[,-"var"], 1, function(row) sort(row, index.return=T, decreasing = F)$ix) %>% t %>% data.table



