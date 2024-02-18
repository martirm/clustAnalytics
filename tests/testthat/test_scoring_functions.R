test_that("test scoring functions", {
    data(karate, package="igraphdata")
    expect_equal_to_reference(scoring_functions(karate, V(karate)$Faction, 
                                                no_clustering_coef = FALSE),
                              file="test_cache_scoring_functions_1")
    

})
