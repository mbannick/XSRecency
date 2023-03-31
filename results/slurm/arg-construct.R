
get.epi.funcs <- function(itype, pt, bigT){
  if(itype == "constant"){
    inc.function <- XSRecency:::incidence.con
    infection.function <- infections.con
  } else if(itype == "linear"){
    inc.function <- XSRecency:::incidence.lin
    infection.function <- infections.lin
  } else if(itype == "exponential"){
    inc.function <- XSRecency:::incidence.exp
    infection.function <- infections.exp
  } else if(itype == "piecewise"){
    if(!pt){
      stop("Piecewise constant-linear incidence function only to be used
         with prior testing simulations.")
    }
    infection.function <- function(...) infections.lincon(bigT=bigT, ...)
  } else {
    stop("Unknown incidence function.")
  }
  return(list(
    infection.function=infection.function,
    inc.function=inc.function
  ))
}
