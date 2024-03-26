#' @export
plot.backnow <- function(obj, plottype, ...){
  
  print("plotting backnow")
  
  stopifnot(plottype %in% c('est', 'rt'))
  
  if(plottype == 'est') {
    
    plot(x = obj$report_date, y = obj$report_cases,
         xlab = 'Date', ylab = 'N. Cases')
  
    newx <- obj$est_back_date
    lb <- obj$est_back[1, ]
    ub <- obj$est_back[3, ]
    
    polygon(c(rev(newx), newx), c(rev(ub), lb), 
            col = 'grey80', border = NA)
    
    lines(x = obj$est_back_date, y = obj$est_back[2, ], col = 'red', lt = '11')
    
  } else {
    
    plot(x = obj$est_rt_date, y = obj$est_rt[2,], col = 'white',
         xlab = 'Date', ylab = 'r(t)')
    
    newx <- obj$est_rt_date
    lb <- obj$est_rt[1, ]
    ub <- obj$est_rt[3, ]
    
    polygon(c(rev(newx), newx), c(rev(ub), lb), 
            col = 'grey80', border = NA)
    
    lines(x = obj$est_rt_date, y = obj$est_rt[2, ], col = 'red', lt = '11')
  }

  
}
