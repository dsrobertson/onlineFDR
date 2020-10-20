checkSTARdf <- function(d, version) {
  switch(as.character(version),
         #async
         "async" = {
           if(!("decision.times" %in% colnames(d))) {
             stop("Your dataset does not contain the column 'decision.times'.")
           }
           if(c("lags") %in% colnames(d)) {
             warning("Your dataset should not contain the columns 'lags'.")
           }
         },
         #dep
         "dep" = {
           if(!("lags" %in% colnames(d))) {
             stop("Your dataset does not contain the column 'lags'.")
           }
           if(c("decision.times") %in% colnames(d)) {
             warning("Your dataset should not contain the columns 'decision.times'.")
           }
         }
  )
}

