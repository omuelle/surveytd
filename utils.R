
addColumnAggregate <- function(table, func, row_name)
{
  total_row <- lapply(table, function(col) if(is.numeric(col)) func(col) else row_name)
  rbind(table, do.call("data.table", total_row))
}

roundTable <- function(table, digits=2)
{
  res <- do.call(data.table, lapply(table, function(col) if(is.numeric(col)) round(col, digits) else col))
  names(res) <- names(table)
  res
}

addMonths <- function(date, months) 
{
  m <- data.table::month(date)-1+months
  add_this <- m%%12+1
  
  d_str <- paste0(data.table::year(date)+m%/%12, "-",
                  ifelse(nchar(add_this) == 1, paste0("0", add_this), add_this), "-",
                  data.table::mday(date))
  
  tryCatch({as.Date(d_str)}, error = function(e) NA)
}

ymDate <- function(year, month)
{
  as.Date(paste0(year, "-", month, "-1"))
}

dateToTime <- function(date)
{
  data.table::year(date) + (data.table::month(date) - 1)/12
}

timeToDate <- function(time)
{
  time <- round(time, 3)
  ymDate(floor(time), round((time - floor(time)) / (1/12)) + 1)
}

startTime <- function(ts)
{
  time(ts)[1]
}

endTime <- function(ts)
{
  time(ts)[length(ts)]
}

startDate <- function(ts)
{
  timeToDate(startTime(ts))
}

endDate <- function(ts)
{
  timeToDate(endTime(ts))
}

frequencyDate <- function(date, frequency)
{
  ymDate(year(date), (ceiling(month(date) / (12/frequency)) - 1) * (12/frequency) + 1)
}

toQuarterlySeries <- function(series)
{
  dates <- zoo::as.Date(time(series))
  is_quarter <- month(dates) %in% c(1,4,7,10)
  ts(series[is_quarter], start=dateToTime(dates[is_quarter][1]), frequency=4)
}

toMonthlySeries <- function(x, fill_last_quarter = F)
{
  vals <- unlist(lapply(as.vector(x), function(val) c(val, NA, NA)))
  if(!fill_last_quarter) {
    vals <- vals[1:(length(vals)-2)]
  }
  ts(vals, start = startTime(x), frequency = 12)
}

stdize <- function(x)
{
  (x - mean(x, na.rm=T)) / sd(x, na.rm=T)
}

corTS <- function(a, b)
{
  from <- dateToTime(max(startDate(a), startDate(b)))
  to <- dateToTime(min(endDate(a), endDate(b)))
  cor(window(a, from, to), window(b, from, to))
}


extractMonths <- function(series, months)
{
  series[month(timeToDate(time(series))) %in% months]
}

mse <- function(a, b, months = NULL)
{
  diff <- a - b
  if(!is.null(months)) {
    diff <- extractMonths(diff, months)
  }
  mean(diff^2)
}

rmse <- function(a, b, months = NULL)
{
  sqrt(mse(a, b, months))
}

