# Tests of the rewireCpp function

# currently converting graphs to edgelists to check for equality, because it 
# seems some metadata of igraphs object prevents identical graphs created at 
# different times from being identical igraph objects. 
####### IMPORTANT: figure this out ############



test_that("trivial test. checks that the graph remains invariant when running 
      rewireCpp with 0 iterations",{
    data("karate", package = "igraphdata")
    expect_equal(igraph_to_edgelist(rewireCpp(karate, Q=0)), 
               igraph_to_edgelist(karate))
})


test_that("Check that total edge weight remains invariant with rewireCpp",{
    data("g_forex")
    total_weight <- function(g) {
        sum(as_adjacency_matrix(g, attr="weight"))
    }
    
    expect_equal(total_weight(karate), 
                 total_weight( rewireCpp(karate, Q=10, weight_sel="max_weight") ) )
    
    expect_equal(total_weight(g_forex),
                 total_weight( rewireCpp(g_forex, Q=10, weight_sel="const_var") ) )
        
}
)






