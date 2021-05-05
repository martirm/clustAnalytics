

library(igraphdata)
data("karate")
karate_gt_clustering <- c(1,1,1,1,1,1,1,1,1,2,1,1,1,1,2,2,1,1,
                          2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2)

test_that("test scoring functions", {
    expect_equal_to_reference(scoring_functions(karate, karate_gt_clustering, 
                                                no_clustering_coef = FALSE),
                              file="test_cache_scoring_functions_1")
    

})
