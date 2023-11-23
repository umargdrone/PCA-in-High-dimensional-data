#Computer Practical 3
#Principal component analysis (PCA)
################################################################################

#Applying PCA on trees, a low dimensional data set

# Load trees data
data(trees)

# Convert volume to factor
trees$vol_fac <- as.factor(ifelse(trees$Volume > 25, "large", "small"))

# Preview data
head(trees)

library(ggplot2)
trees_plot <- ggplot(trees) +
  geom_point(aes(Girth, Height, size=Volume, col=vol_fac),alpha=0.6) +
  labs(x = "Girth (inches)", y = "Height (inches)",
       color = "Volume class")
trees_plot +  ggtitle("Visualising the original trees data")

#PCA
# Create a matrix using Height and Girth variables
tree_mx <- cbind("height" = trees$Height, "girth" = trees$Girth)
tree_mx <- scale(tree_mx, center=TRUE, scale=TRUE)

tree_pca <- prcomp(tree_mx, center=FALSE, scale.=FALSE)#no need to standardise here as the variables are already scaled
#tree_pca <- prcomp(tree_mx)
summary(tree_pca)
tree_pca$rotation #loadings

pc_loadings <- t(tree_pca$rotation) * tree_pca$sdev

## Reuse plot from before and add PCs
library(tidyverse)
theme_set(theme_light())
trees_plot2 <- tree_mx %>%
  ## Convert to data.frame and re-add Volume columns for plotting
  data.frame(Volume = trees$Volume, vol_fac = trees$vol_fac) %>% 
  ggplot() +
  geom_point(aes(girth, height, size = Volume, col = vol_fac),
             alpha = 0.6) +
  labs(x = "Girth (inches), standardised", y = "Height (feet), standardised", 
       color = "")

trees_plot2+
  geom_segment(
    data = data.frame(pc_loadings),
    aes(x = 0, xend = girth, y = 0, yend = height),
    arrow = arrow(length=unit(0.2,"cm"))
  ) +
  geom_text(
    data = data.frame(pc_loadings),
    aes(x = girth, y = height, label = rownames(pc_loadings)),
    nudge_x = 0.1, nudge_y = 0.2
  ) +
  coord_equal(xlim = c(-2.3, 2.3), ylim = c(-2.1, 2.1)) +
  ggtitle("Trees data overlayed with PCs",
          subtitle = "Based on PCA on standardised data")


library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
ggbiplot(tree_pca, groups = trees$vol_fac, alpha = 0) +
  ## Add points layer to color and size points
  geom_point(aes(col = trees$vol_fac, size = trees$Volume), alpha = 0.6) +
  labs(size = "Volume", col = "") +
  theme(aspect.ratio = 0.6) +
  ggtitle("Biplot for the PCA on non-standardised data")


#Now to repeat the above calculations for PCA without scaling the variables 
#and by changing the unit of height from inches to millimeters
trees$vol_fac <- as.factor(ifelse(trees$Volume > 25, "large", "small"))
tree_mx <- cbind("height" = trees$Height, "girth" = trees$Girth)
trees$height <- trees$Height *25.4#change the unit from inches to mm
tree_mx <- cbind("height" = trees$height, "girth" = trees$Girth)

#PCA
tree_pca <- prcomp(tree_mx, center=FALSE, scale.=FALSE)
summary(tree_pca)
tree_pca$rotation
pc_loadings <- t(tree_pca$rotation) * tree_pca$sdev

ggbiplot(tree_pca, groups = trees$vol_fac, alpha = 0) +
  ## Add points layer to color and size points
  geom_point(aes(col = trees$vol_fac, size = trees$Volume), alpha = 0.6) +
  labs(size = "Volume", col = "") +
  theme(aspect.ratio = 0.6) +
  ggtitle("Biplot for the PCA on non-standardised data")


################################################################################

#Applying PCA on tumor_tissue, a high dimensional data set

load(file="tumor_tissue.rda")

str(tumor_tissue[, 1:10])

Y <- as.factor(tumor_tissue[, 1])
levels(Y) <- c("Normal", "Tissue")
X <- tumor_tissue[, -1]
X <- scale(X, center=TRUE, scale=TRUE)
length(Y)
dim(X)

#SVD 
svd_X <- svd(X)
V <- svd_X$v #eigenvectors of sample covariance X^T*X/n 
svd_X$d

Z <- svd_X$u %*% diag(svd_X$d) #calculate and store the scores

#PCA
pca_x <- prcomp(X, center=FALSE, scale.=FALSE)#no need to standardise here as the variables were already scaled
summary(pca_x)
pca_x$x #scores - same as Z in the above obtained directly using SVD
pca_x$rotation #loadings
View(pca_x$rotation-V) #to check the loadings obtained using prcomp and SVD are the same
princomp(X) #note that this R function is only applicable to low dimensional data when n>p

library(factoextra)

fviz_eig(pca_x, ncp=62, addlabels=TRUE)

fviz_pca_var(pca_x, col.var="contrib", gradient.cols=c("white", "blue", "red"))

ggbiplot(pca_x, groups = Y, alpha = 0) +
  ## Add points layer to color points
  geom_point(aes(col = Y), alpha = 0.6) +
  theme(aspect.ratio = 0.6) +
  ggtitle("Biplot for the first two PCs on tumor data")

#histograms
# First PC, first column of V
hist(V[, 1], breaks = 50, xlab = "PC 1 loadings", main = "")
# Add vertical line at 95% quantile
abline(v = quantile(V[, 1], 0.95), col = "red", lwd = 2)

# Second PC, second column of V
hist(V[, 2], breaks = 50, xlab = "PC 2 loadings", main = "")
abline(v = c(
  quantile(V[, 2], 0.05),
  quantile(V[, 2], 0.95)
), col = "red", lwd = 2)


#Sparse PCA
library(glmnet)

# For PC1
set.seed(1)
fit_loadings1 <- cv.glmnet(X, Z[, 1],
  alpha = 0.5, nfolds = 5
)
fit_loadings1

# For PC2
set.seed(1)
fit_loadings2 <- cv.glmnet(X, Z[, 2], alpha = 0.5, nfolds = 5)
fit_loadings2

par(mfrow = c(2, 1))
plot(fit_loadings1, main = "PC1")
plot(fit_loadings2, main = "PC2")

par(mfrow = c(2, 1))
plot(fit_loadings1$glmnet.fit, main = "PC1", xvar = "lambda")
abline(v = log(fit_loadings1$lambda.min), lty = 3)
abline(v = log(fit_loadings1$lambda.1se), lty = 3)
plot(fit_loadings2$glmnet.fit, main = "PC2", xvar = "lambda")
abline(v = log(fit_loadings2$lambda.min), lty = 3)
abline(v = log(fit_loadings2$lambda.1se), lty = 3)

################################################################################
