"RGBColVec" <-
function (nrgcols = 13) 
  {
    k <- trunc(nrgcols/2)
    if (2 * k == nrgcols) {
      r <- c(rev(seq(0, 1, length = k)), rep(0, k))
      g <- c(rep(0, k), seq(0, 1, length = k))
      colvec <- rgb(r, g, rep(0, 2 * k))
    }
    else {
      r <- c(rev(seq(1/(2 * k - 1), 1, length = k)), rep(0, 
                                         k + 1))
      g <- c(rep(0, k + 1), seq(1/(2 * k - 1), 1, length = k))
      colvec <- rgb(r, g, rep(0, 2 * k + 1))
    }
    colvec
  }

