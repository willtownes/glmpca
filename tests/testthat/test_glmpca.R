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