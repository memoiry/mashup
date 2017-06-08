#' compare betas and eigenvalues for networks
rm(list=ls())

require(R.matlab)
require(plotrix)
require(ggplot2)

# --------------------------------
# helper functions

# from https://susanejohnston.wordpress.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/ 
ggplotRegression <- function (fit) {

p <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 3),
                     "Intercept =",signif(fit$coef[[1]],2 ),"\n",
                     " Slope =",signif(fit$coef[[2]], 2),
                     " P =",signif(summary(fit)$coef[2,4], 2)))
p
}

# --------------------------------
# work begins
outDir <- "result/netdx/H_vs_beta"
if (!file.exists(outDir)) dir.create(outDir)

# load data
data <- readMat("result/netdx/eigenvalue_list.mat")
betas <- readMat("result/netdx/linear_regression_output.mat")
H <- t(data[[1]])
H <- cbind(H,betas[[1]])
colnames(H) <- c("age","grade","stage","beta")

# plot colored table
pdf(sprintf("%s/H_vs_beta_coloredMat.pdf",outDir),width=5,height=11)
plotrix::color2D.matplot(H,show.values=F,axes=F,extremes=c("yellow","red"),
border='white',las=1,ylab="K",xlab="N, beta")
axis(1,at=1:ncol(H)-0.5, labels=colnames(H),cex.axis=0.5)
dev.off()

# plot H vs beta
H <- as.data.frame(H)

pList <- list()
fit <- lm(beta ~ age, data = H); pList[[1]] <- ggplotRegression(fit)
fit <- lm(beta ~ grade, data = H); pList[[2]] <- ggplotRegression(fit)
fit <- lm(beta ~ stage, data = H); pList[[3]] <- ggplotRegression(fit)

# plot H vs beta - zoomed for smaller weights
require(ggplot2)
pList2 <- list()
pList2[[1]] <- ggplot(H,aes(age,beta))+geom_smooth(method='lm',formula=y~x)+geom_point() + xlim(0,0.5)+ggtitle("age")
pList2[[2]] <- ggplot(H,aes(stage,beta))+geom_smooth(method='lm',formula=y~x)+geom_point() + xlim(0,0.5)+ggtitle("stage")
pList2[[3]] <- ggplot(H,aes(grade,beta))+geom_smooth(method='lm',formula=y~x)+geom_point() + xlim(0,0.5)+ggtitle("grade")

source("multiplot.R")
pdf(sprintf("%s/H_vs_beta_corrDensity.pdf",outDir),width=8,height=8)
tryCatch({
	multiplot(plotlist=pList,layout=matrix(1:4,ncol=2))
	multiplot(plotlist=pList2,layout=matrix(1:4,ncol=2))

	# plot density distribution of values
	par(mfrow=c(2,2),las=1,bty='n')
	for (k in 1:4) {
		plot(density(H[,k]),main=colnames(H)[k])
	}
},error=function(ex) {
	print(ex)
},finally={
	dev.off()
})
