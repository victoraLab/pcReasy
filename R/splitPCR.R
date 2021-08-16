#' Calculate PC reaction
#'
#' Returns a efficient protocol for pc reaction setup.
#'
#' @param nvars The number of variable. Numeric vector.
#' @param vars  Variable PCR components. Character vector.
#' @importFrom stringr str_extract
#'
#'
  splitPCR <- function(pcr = NULL, nVars=c(2,2,2), vars=c("MM", "DNA", "FPrimer", "RPrimer"), levels = c(0, 1, 2, 2)){
    options(scipen=999)

    #If empty, generate example from Q5 protocol
    if(length(initialConc) == 0){
      pcr <- pcalc(v = 25, pError = 0.1)
    }

    #Get the pipetting error rate used
    e <- as.numeric(str_extract(pattern = "[0-1].[0-9]$", string = colnames(pcr)[grep("_", colnames(pcr))]))
    e <- e + 1

    groups <- split(vars, levels)
    names(groups) <- sapply(names(groups), function(x){paste0(paste0(rep("s",as.numeric(x)), collapse = ""), "mm")})

    #Calculate the master mix (mm) volume
    total <- pcr["Total", grep("_", colnames(pcr))]
    varsTotal <- pcr[vars,grep("_", colnames(pcr))]

    mm <- total - sum(varsTotal, na.rm = T)
    varsTotal[1]  <- mm
    e.list <- list()
    #Calculate errors
    for( i in 1:length(groups)-1){
      e <- ((e-1)/2)+1
      e.list[i] <-  e
    }

    splitMM = (mm/e) * 1/nVars

    volSide <- (varsTotal/e) * 1/nVars



}
