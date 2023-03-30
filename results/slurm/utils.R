replace.args <- function(main, defaults){
  for(arg in names(defaults)){
    if(!arg %in% names(main)){
      main[[arg]] <- defaults[arg]
    }
  }
}
