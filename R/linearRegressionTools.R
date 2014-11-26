# linearRegressionTools.R


`lmSlopeDifference` <- function( lmAns1, lmAns2, coef1=2, coef2=2) {

	# extract the slope and standard error from each liner model
	lmSum1 <- summary( lmAns1)
	cm1 <- coef( lmSum1)
	lmSum2 <- summary( lmAns2)
	cm2 <- coef( lmSum2)

	# now extract the wanted coefficient from each
	slope1 <- cm1[ coef1, 1]
	se1 <- cm1[ coef1, 2]
	slope2 <- cm2[ coef2, 1]
	se2 <- cm2[ coef2, 2]

	# difference in slopes is 'model2 - model1'
	dslope <- slope2 - slope1
	t.value <- dslope / sqrt( se1*se1 + se2*se2)
	degFree <- lmSum1$df[2] + lmSum2$df[2]
	p.value <- pt( abs(t.value), df=degFree, lower.tail=FALSE) * 2

	return( list( "difference"=dslope, "t.value"=t.value, "p.value"=p.value))
}
