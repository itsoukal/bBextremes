#' Function to exctract block-maxima
#'
#' @param x a vector with the time series (not xts)
#' @param k length of block
#' @param n number of blocks
#' @param fun function to apply (default is max)
#'
#' @return
#' @export
#'
#' @examples
blockmaxxer=function(x, k, n, fun=max) {
  # x: a vector with the time series (not xts)
  # k: length of block
  # n: number of block

  if (length(x)!=(k*n))
  {
    print(paste0('The length of the time series is not equal with k*n'))
    break()
  } else {
    start=seq(from=1, to = (k*n), by = k)
    end=seq(from=k, to = (k*n), by = k)
    maxx=rep(NA, n)
    for (i in 1:n) {
      maxx[i]=fun(x[start[i]:end[i]])
    }
  }
  return(maxx)
}
