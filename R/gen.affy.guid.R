#' In the headers of affy mps/pgf/... files there is a line like:
#' #%guid=0000008635-1158704285-1969988809-0053097456-1768345064 This has 50
#' digits, in 5 blocks of 10; This function generates a random guid
#' 
#' @author Mark Cowley, 7/1/07
#' @export
gen.affy.guid <- function( ) {
    tmp <- ""
    for(i in 1:5) {
        if( i== 1 )
            tmp <- paste(round(runif(10,0,10)), collapse="")
        else
            tmp <- paste(c(tmp, paste(round(runif(10,0,10)), collapse="")), collapse="-")
    }

    tmp
}
#     tmp <- round(runif(50,0,10))
#     sprintf("%s%s%s%s%s%s%s%s%s%s-%s%s%s%s%s%s%s%s%s%s-%s%s%s%s%s%s%s%s%s%s-%s%s%s%s%s%s%s%s%s%s-%s%s%s%s%s%s%s%s%s%s",tmp)    
#     sprintf("%d%d%d%d%d%d%d%d%d%d-%d%d%d%d%d%d%d%d%d%d-%d%d%d%d%d%d%d%d%d%d-%d%d%d%d%d%d%d%d%d%d-%d%d%d%d%d%d%d%d%d%d",tmp)    
