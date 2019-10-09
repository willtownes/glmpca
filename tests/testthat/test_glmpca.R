context("Test glmpca and related functions.")
# this is where it would be helpful to have a very small example dataset so
# that we know the correct answer
# data("glmpcaMiniExample")
set.seed(1234)

test_that("glmpca works with intended input types", {
	Y <- matrix(rnbinom(50, 4, .8), ncol = 5)
	
	g1 <- glmpca(Y, L=2, fam = "poi")
	expect_is(g1, "list")
	expect_equal(dim(g1$factors), c(5,2))
	expect_equal(dim(g1$loadings), c(10,2))
	expect_is(g1$dev, "numeric")
	expect_is(g1$family, "glmpca_family")
	
	# try out different arguments, make sure it still runs
	g1 <- glmpca(Y, L=2, fam = "nb", verbose = TRUE, X = matrix(1:5, ncol=1),
				 Z = matrix(1:10, ncol=1))
	g1 <- glmpca(Y, L=2, fam = "mult", sz = 1:5)
	f0 <- matrix(rnorm(10)/10, ncol=2)
	l0 <- matrix(rnorm(20)/10, ncol=2)
	g1 <- glmpca(Y, L=2, init = list(factors = f0, loadings = l0))
	Y[Y > 0] <- 1
	g1 <- glmpca(Y, L=2, fam = "bern")
	
	# pathological cases
	Y[1,1] <- 2
	expect_error(glmpca(Y, L=2, fam = "bern"), "max(Y) <= 1 is not TRUE", fixed = TRUE)
	Y[2,2] <- -1
	expect_error(glmpca(Y, L=2, fam = "poi"), "min(Y) >= 0 is not TRUE", fixed = TRUE)
	Y[2,] <- 0
	expect_error(glmpca(Y, L=2, fam = "poi"), "Some rows were all zero,")
})


test_that("glmpca works with L=1", {
  counts <- matrix(rpois(n=100 * 5, lambda = 5), ncol = 5, nrow=100)
  res <- glmpca(counts, L=1)
  expect_equal(dim(res$factors), c(5,1))
  expect_equal(dim(res$loadings), c(100,1))
})
  