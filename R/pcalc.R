#' Calculate PC reaction
#'
#' Returns a efficient protocol for pc reaction setup.
#'
#' @param v The final volume desired for each PC reaction. Numeric.
#' @param rxn The number of total reactions. Numeric.
#' @param initialConc Data Frame containing the initial concentrations of each component.
#' @param finalConc Data Frame containing the final concentrations of each component.
#' @param pError Pipetting error to consider for each reaction.
#' @importFrom dplyr bind_rows
#'
#' @export

pcalc <- function(v = 10, rxn = NULL, initialConc = NULL, finalConc = NULL, pError = 0){
  options(scipen=999)

if(length(initialConc) == 0){
  finalConc <- data.frame(
    buffer = 1,
    dNTP = 200,
    FPrimer = 10,
    RPrimer = 10,
    DNA = 10,
    Q5Poly = 0.02,
    optional = 0,

    row.names = "Q5"
  )

  initialConc <- data.frame(
    buffer = 5,
    dNTP = 10000,
    FPrimer = 100,
    RPrimer = 100,
    DNA = 100,
    Q5Poly = 2,
    optional = 0,

    row.names = "Q5"
  )
  rxn = 10
}

res <- list()
initialConc <- initialConc[,!initialConc == 0]
finalConc   <- finalConc  [,!finalConc   == 0]

for(i in 1:ncol(initialConc)){
  z <- initialConc[[i]]/finalConc[[i]]
  res[[i]] <- model.frame(y ~ I(x / z), data = data.frame(x = v, y = initialConc[[i]]))
}
res <- bind_rows(res)

res <- data.frame(initialConc = as.vector(t(initialConc)),
                  finalConc = as.vector(t(finalConc)),
                  rxnVols = res[["I(x/z)"]],
                  row.names = colnames(initialConc))

parc <- sum(res[["rxnVols"]])
water <- v - parc

if(parc > v){
  stop(
    paste(
      sprintf("Total volume cannot be bigger than %s", v),
      sprintf("Volume is at %s", parc),
      sep = "\n"))
}

pError <- 1 + pError

res["Water",] <- c(0,0, water)
res["Total",] <- c(0,0, water + parc)
res[, sprintf("rxnVols%s_pError%s", rxn, pError-1)] <- (res[,"rxnVols"] * rxn) * pError




#print(sprintf("The protocol for the reaction named: -- %s --", rownames(initialConc)))
return(res)
}

