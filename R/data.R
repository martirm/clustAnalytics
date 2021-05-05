
#' Forex correlation network
#' 
#' Network built from correlations between time series of exchange rate returns.  
#' It was built from the 13 most traded currencies and with data of January 2009. 
#' It is a complete graph of 78 vertices (corresponding to pairs of currencies) 
#' and has edge weights bounded between 0 and 1.
#' 
#' @format An igraph object with 78 vertices and 3003 weighted edges
"g_forex"