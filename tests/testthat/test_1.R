# This is an example test of the package

#library(rewire) #probably not needed (???)


# currently converting graphs to edgelists to check for equality, because it 
# seems some metadata of igraphs object prevents identical graphs created at 
# different times from being identical igraph objects. 
####### IMPORTANT: figure this out ############



test_that("trivial test. checks that the graph remains invariant when running 
          rewireCpp with 0 iterations",{
              expect_equal(igraph_to_edgelist(rewireCpp(karate, Q=0)), 
                           igraph_to_edgelist(karate))
              expect_equal_to_reference(igraph_to_edgelist(rewireCpp(karate, Q=0)), 
                                        file="test_cache1")
              #expect_equal_to_reference(seeded_test(), file="test_cache2")
          }
          )



test_that("runs the switching algorithm on several graphs using seeds, and
          checks that the results are consistent",{
              expect_equal_to_reference({
                  set.seed(0)
                  igraph_to_edgelist(rewireCpp(karate, Q=100))
              }, file="test_cache_karate_Q100")
              
              expect_equal_to_reference({
                  set.seed(0)
                  igraph_to_edgelist(rewireCpp(g_forex, Q=100, upper_bound=1))
              }, file="test_cache_forex_Q100")
              
          })




