library(RStoolbox)
library(ggplot2)
library(raster)
data(rlogo)

## Original image a (+1 to prevent log(0))
img_a <-  rlogo + 1

## Degraded image b
img_b <- log(img_a)

## Cut-off half the image (just for better display)
img_b[, 1:50] <- NA

## Compare Images before histMatching
ggRGB(img_a,1,2,3) + 
  ggRGB(img_b, 1,2,3, ggLayer = TRUE, stretch = "lin", q = 0:1) + 
  geom_vline(aes(xintercept = 50)) + ggtitle("Img_a vs. Img_b")

## Do histogram matching
img_b_matched <- histMatch(img_b, img_a)

## Compare Images after histMatching
ggRGB(img_a, 1, 2, 3) + 
  ggRGB(img_b_matched, 1, 2, 3, ggLayer = TRUE, stretch = "lin", q = 0:1) + 
  geom_vline(aes(xintercept = 50)) + ggtitle("Img_a vs. Img_b_matched")

## Histogram comparison
opar <- par(mfrow = c(1, 3), no.readonly = TRUE)
img_a[, 1:50] <- NA
redLayers <- stack(img_a, img_b, img_b_matched)[[c(1,4,7)]]
names(redLayers) <- c("img_a", "img_b", "img_b_matched")

hist(redLayers)

## Reset par
par(opar)


#########################################################################################################
setwd("C:/Users/SRDhruba/Dropbox (Personal)/ResearchWork/Rtest/")

library(tiff)
library(RStoolbox)
library(raster)
library(ggplot2)
library(grid)
library(ggpubr)

# data(rlogo)

x1 <- readTIFF("tire.tif", native = FALSE)
x2 <- readTIFF("pout.tif", native = FALSE)

X1 <- raster(x1)
X2 <- raster(x2)

X1m <- histMatch(X1, ref = X2)
x1m <- as.matrix(X1m)


plot.new()
rasterImage(x1, 0, 0, 0.2, 0.8);    rasterImage(x2, 0.4, 0, 0.6, 0.8);    rasterImage(x1m, 0.8, 0, 1, 0.8)

op <- par(mfrow = c(1, 3))
hist(x1, main = "target");          hist(x2, main = "ref");               hist(x1m, main = "matched")

par(op)


####
library(RStoolbox)
library(raster)

hist.match <- function(src, ref, n.samp = 1e6){
  ## Estimate distributions...
  ss <- sample(src, size = n.samp, replace = TRUE)
  rr <- sample(ref, size = n.samp, replace = TRUE)
  ss.cdf <- ecdf(ss);     rr.cdf <- ecdf(rr)
  
  ## Calculate mapping...
  kn <- knots(rr.cdf);    kn.lims <- range(rr)
  fn <- rr.cdf(kn)
  rr.cdf.inv <- approxfun(x = fn, y = kn, yleft = kn.lims[1], yright = kn.lims[2], ties = "ordered")
  
  ## Return matched samples...
  src.matched <- rr.cdf.inv(ss.cdf(src))
  return(src.matched)
}

y1 <- runif(n = 100, min = -2, max = 1)
# y1 <- rnorm(n = 100, mean = 1.2, sd = 0.65)
y2 <- rnorm(n = 100, mean = 0, sd = 0.8)

op <- par(mfrow = c(1, 4))
hist(y1, freq = FALSE, main = "target")
hist(y2, freq = FALSE, main = "ref")

Y1 <- raster(as.matrix(y1))
Y2 <- raster(as.matrix(y2))

Y1m <- histMatch(Y1, ref = Y2)
y1m <- getValues(Y1m)

hist(y1m, freq = FALSE, main = "matched_1")

y1m2 <- hist.match(y1, ref = y2)

hist(y1m2, freq = FALSE, main = "matched_2")

par(op)


