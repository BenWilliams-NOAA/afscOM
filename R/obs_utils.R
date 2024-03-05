simulate_rpn <- function(q, naa, sel){
    total_naa <- apply(naa, 2, sum)
    return(q*sum(naa*sel))
}
