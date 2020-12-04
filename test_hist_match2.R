# rm(list = ls())

## Set up system path... 
info <- Sys.getenv(c("USERNAME", "HOMEPATH"))
if (info["USERNAME"] == "SRDhruba"){
  info["DIRPATH"] <- sprintf("%s\\Dropbox%s\\ResearchWork\\Rtest\\", info["HOMEPATH"], " (Personal)")
} else {
  info["DIRPATH"] <- sprintf("%s\\Dropbox%s\\ResearchWork\\Rtest\\", info["HOMEPATH"], "")
}
setwd(info["DIRPATH"]);       cat("Current system path = ", getwd())


## Packages...
library(ggplot2)
library(ggpubr)
library(randomForest)
library(sROC)
library(locfit)
library(ks)


#### Functions...
printf <- function(..., end = "\n"){
  if ((nargs() > 1) & (grepl(list(... )[1], pattern = "%")))
    cat(sprintf(...), end)
  else
    cat(..., end)
}

dapply <- function(df, ...) as.data.frame(apply(df, ...))

norm01 <- function(x) (x - min(x)) / diff(range(x))

get.dist <- function(x, N = 1e6, dist.opt = "hist"){
  dist.opt <- tolower(dist.opt)
  
  xx <- sample(x, size = N, replace = TRUE)
  if (dist.opt == "hist") {
    dist.xx <- ecdf(xx)
  } else if (dist.opt == "kde") {
    dist.xx <- kCDF(xx, kernel = "normal", bw = bw.CDF.pi(xx, pilot = "nrd0"))
  }
}

calc.err <- function(y, y.pred, measure = "MSE"){
  measure <- toupper(measure)
  if (grepl(pattern = "MSE", measure)){
    err <- mean((y - y.pred)^2)
    if (measure == "RMSE"){
      err <- sqrt(err)
    } else if (measure == "NRMSE"){
      err <- sqrt(err / mean((y - mean(y))^2))
    }
  } else if (grepl(pattern = "MAE", measure)){
    err <- mean(abs(y - y.pred))
    if (measure == "NMAE"){
      err <- err / mean(abs(y - mean(y)))
    }
  }
  err
}

### Pick top genes...
get.top.genes <- function(ranks, m.top = 150, print.opt = FALSE) {
  
  ## Initialization...
  nI <- 0;        nGN <- 300
  gene.rank <- intersect(ranks[1:nGN, 1], ranks[1:nGN, 2])
  m  <- length(gene.rank);     m0 <- if (print.opt) m
  
  ## Run iterations...
  while(m < m.top) {
    nI <- nI + 1;       nGN <- nGN + 100
    gene.rank <- intersect(ranks[1:nGN, 1], ranks[1:nGN, 2])
    m  <- length(gene.rank)
  }
  gene.rank <- sort(gene.rank, decreasing = FALSE)    # Sort ranks
  
  ## Print results...
  if (print.opt)
    printf("#top genes chosen = %d (nGN = %d, nI = %d, m0 = %d)", m, nGN, nI, m0)
  
  gene.rank
}


#### Read tumor-cell line data...
Xdata1 <- read.table("Data/BRCA_gene_expression_METABRIC_26_Oct_2020.txt", sep = "\t", header = TRUE)
Xdata2 <- read.table("Data/BRCA_gene_expression_CCLE_26_Oct_2020.txt", sep = "\t", header = TRUE)
Xdata3 <- read.table("Data/BRCA_gene_expression_GDSC_26_Oct_2020.txt", sep = "\t", header = TRUE)

Ydata1 <- read.table("Data/BRCA_biomarker_expression_METABRIC_26_Oct_2020.txt", sep = "\t", header = TRUE)
Ydata2 <- read.table("Data/BRCA_biomarker_expression_CCLE_26_Oct_2020.txt", sep = "\t", header = TRUE)
Ydata3 <- read.table("Data/BRCA_biomarker_expression_GDSC_26_Oct_2020.txt", sep = "\t", header = TRUE)

rank1 <- read.table("Data/BRCA_biomarker_ranks_METABRIC_27_Oct_2020.txt", sep = "\t", header = TRUE)
rank2 <- read.table("Data/BRCA_biomarker_ranks_CCLE_27_Oct_2020.txt", sep = "\t", header = TRUE)
rank3 <- read.table("Data/BRCA_biomarker_ranks_GDSC_27_Oct_2020.txt", sep = "\t", header = TRUE)


#### Get data for a biomarker...
biomarkers <- colnames(Ydata1);       # printf("Biomarkers = ", biomarkers)
q <- length(biomarkers)

run <- function(k, s = 4321) {
# k <- 1;     s <- 4321 
bmChosen <- biomarkers[k];            printf("\nChosen biomarker = %s", bmChosen)
ranks <- cbind(rank1[, bmChosen], rank2[, bmChosen], rank3[, bmChosen])
gnRank <- get.top.genes(ranks[, 2:3], m.top = 150, print.opt = TRUE)
m <- length(gnRank)

X1 <- Xdata1[, gnRank];                   X2 <- rbind(Xdata2[, gnRank], Xdata3[, gnRank])
Y1 <- norm01(Ydata1[, bmChosen]);         Y2 <- norm01(c(Ydata2[, bmChosen], Ydata3[, bmChosen]))

X1n <- dapply(X1, MARGIN = 2, norm01);    X2n <- dapply(X2, MARGIN = 2, norm01)

grid.size <- 1e2
get.cdf <- function(z, ...) kcde(z, h = hscv(z, nstage = 2, binned = TRUE, bgridsize = grid.size*10), 
                                 binned = TRUE, bgridsize = grid.size, xmin = 0, xmax = 1,...)

X2n.map <- lapply(1:m, function(j) {
  x1j <- sample(X1n[, j], size = 1e4, replace = TRUE)
  x2j <- sample(X2n[, j], size = 1e4, replace = TRUE)
  
  # x2j.cdf <- kcde(x2j, binned = TRUE, bgridsize = grid.size, xmin = 0, xmax = 1)
  # map <- function(fn.xj) qkde(fhat = x2j.pdf, p = fn.xj)
  
  x2j.cdf <- get.cdf(x2j)
  map <- approxfun(x = x2j.cdf$estimate, y = x2j.cdf$eval.points, method = "linear", yleft = 0, yright = 1, ties = "ordered")
  
  # x1j.cdf <- kcde(x1j, binned = TRUE, bgridsize = grid.size, xmin = 0, xmax = 1, eval.points = X1n[, j])
  # fn.xj <- pkde(fhat = x1j.pdf, q = X1n[, j])
  
  x1j.cdf <- get.cdf(x1j, eval.points = X1n[, j])
  x2j.map <- map(x1j.cdf$estimate)
})
X2n.map <- as.data.frame(X2n.map, row.names = rownames(X1n), col.names = colnames(X1n))

## Build predictive model & predict response...
set.seed(s)
RF <- randomForest(x = X2n, y = Y2, ntree = 200, mtry = 5, replace = TRUE)
Y2.pred.map <- predict(RF, X2n.map)
Y2.pred.map[Y2.pred.map < 0] <- 0;    Y2.pred.map[Y2.pred.map > 1] <- 1


## Histogram matching for response...
y1 <- sample(Y1, size = 1e4, replace = TRUE)
y2 <- sample(Y2, size = 1e4, replace = TRUE)

# y1.cdf <- kcde(y1, binned = TRUE, bgridsize = grid.size, xmin = 0, xmax = 1)
# map.y <- function(fn.y) qkde(fhat = y1.pdf, p = fn.y)

y1.cdf <- get.cdf(y1)
map.y <- approxfun(x = y1.cdf$estimate, y = y1.cdf$eval.points, method = "linear", yleft = 0, yright = 1, ties = "ordered")

# y2.cdf <- kcde(y2, binned = TRUE, bgridsize = grid.size, xmin = 0, xmax = 1, eval.points = Y2.pred.map)
# fn.y <- pkde(fhat = y2.pdf, q = Y2.pred.map)

y2.cdf <- get.cdf(y2, eval.points = Y2.pred.map)
Y1.pred <- map.y(y2.cdf$estimate)
# Y1.pred <- map.y(fn.y)


#### Baseline model...
set.seed(s)
RF.base <- randomForest(x = X2n, y = Y2, ntree = 200, mtry = 5, replace = TRUE)
Y1.pred.base <- predict(RF.base, X1n)
Y1.pred.base[Y1.pred.base < 0] <- 0;    Y1.pred.base[Y1.pred.base > 1] <- 1


## Results table...
results <- data.frame("DMTL" = c(calc.err(Y1, Y1.pred, measure = "NRMSE"), calc.err(Y1, Y1.pred, measure = "NMAE"), 
                                 cor(Y1, Y1.pred, method = "spearman")), 
                      'BL' = c(calc.err(Y1, Y1.pred.base, measure = "NRMSE"), calc.err(Y1, Y1.pred.base, measure = "NMAE"), 
                               cor(Y1, Y1.pred.base, method = "spearman")), 
                      row.names = c("NRMSE", "NMAE", "SCC"))

printf("Results = \n");   print(results)

results
}

res <- sapply(1:q, run)

res2 <- list(DMTL = as.data.frame(t(sapply(1:q, function(i) res[1, ][[i]])), row.names = biomarkers), 
             BL = as.data.frame(t(sapply(1:q, function(i) res[2, ][[i]])), row.names = biomarkers))
res2 <- lapply(names(res2), function(k) { colnames(res2[[k]]) = c("NRMSE", "NMAE", "SCC");    res2[[k]] })
names(res2) <- c("DMTL", "BL")
res2$DMTL["Mean", ] = colMeans(res2$DMTL)
res2$BL["Mean", ] = colMeans(res2$BL)

rbind(DMTL = res2$DMTL["Mean", ], BL = res2$BL["Mean", ])



#############################################################################################
# h2j <- hscv(X2n[, j], nstage = 2, binned = TRUE, bgridsize = 1e4)
# x2j.cdf <- kcde(x2j, tail.flag = "lower.tail", h = h2j, xmin = 0, xmax = 1)
# 
# kn.j <- x2j.cdf$eval.points;     fn.j <- x2j.cdf$estimate
# x2j.cdf <- approxfun(x = kn.j, y = fn.j, method = "linear", ties = "ordered")
# 
# kn.j <- sort(sample(kn.j, size = 1e3, replace = TRUE));     fn.j <- x2j.cdf(kn.j)
# map <- approxfun(x = fn.j, y = kn.j, method = "linear", yleft = 0, yright = 1, ties = "ordered")
# 
# x2j.map <- map(x2j.cdf(X1n[, j]))

# h1 <- hscv(Y1, nstage = 2, binned = TRUE, bgridsize = 1e4)
# h2 <- hscv(Y2, nstage = 2, binned = TRUE, bgridsize = 1e4)
# Y1.cdf <- kcde(y1, tail.flag = "lower.tail", h = h1, xmin = 0, xmax = 1)
# Y2.cdf <- kcde(y2, tail.flag = "lower.tail", h = h2, xmin = 0, xmax = 1)
# 
# kn.y <- Y1.cdf$eval.points;     fn.y <- Y1.cdf$estimate
# Y1.cdf <- approxfun(x = kn.y, y = fn.y, method = "linear", ties = "ordered")
# 
# kn.y <- sort(sample(kn.y, size = 1e3, replace = TRUE));    fn.y <- Y1.cdf(kn.y)
# map.y <- approxfun(x = fn.y, y = kn.y, method = "linear", yleft = 0, yright = 1, ties = "ordered")
# 
# kn.y2 <- Y2.cdf$eval.points;     fn.y2 <- Y2.cdf$estimate
# Y2.cdf <- approxfun(x = kn.y2, y = fn.y2, method = "linear", ties = "ordered")
# 
# Y1.pred <- map.y(Y2.cdf(Y2.pred.map))


#############################################################################################
# X1.cdf <- lapply(1:m, function(j) {
#   # x1j <- sample(X1n[, j], size = 1e4, replace = TRUE)
#   # x2j <- sample(X2n[, j], size = 1e4, replace = TRUE)
# })

j <- 1
x1j <- sample(X1n[, j], size = 1e4, replace = TRUE)
x2j <- sample(X2n[, j], size = 1e4, replace = TRUE)

x2j.cdf.hist <- ecdf(x2j)
kn.2j <- knots(x2j.cdf.hist);   fn.2j <- x2j.cdf.hist(kn.2j)

# h.kde <- hlscv(x2j, binned = FALSE, deriv.order = 0)
# h.kde <- hns(x2j, deriv.order = 1)
h.kde <- hscv(x2j, nstage = 2, binned = TRUE, bgridsize = 1e4)
# h.kde <- hpi(x2j, binned = TRUE)
x2j.cdf.kde <- kcde(x2j, tail.flag = "lower.tail", h = h.kde, xmin = 0, xmax = 1)
# x2j.cdf.kde <- approxfun(x = x2j.cdf.kde$eval.points, y = x2j.cdf.kde$estimate, method = "linear", ties = "ordered")
# pt.2j <- sort(rnorm(n = 1e3, mean = mean(X2n[, j]), sd = sd(X2n[, j])))
# fv.2j <- x2j.cdf.kde(pt.2j)
pt.2j <- x2j.cdf.kde$eval.points;   fv.2j <- x2j.cdf.kde$estimate

plot.new()
op <- par(mfcol = c(1, 2))
barplot(kn.2j, fn.2j, width = 0.8)
barplot(pt.2j, fv.2j, width = 0.8)
par(op)

plot.new()
op <- par(mfcol = c(1, 2))
hist(kn.2j)
hist(pt.2j)
par(op)


# xx <- sample(X1n[, 1], size = 1e3, replace = TRUE)
# x.cdf <- ecdf(xx)
# plt.df <- data.frame(gene = seq(min(xx), max(xx), by = 1e-2))
# plt.df$dist <- x.cdf(plt.df$gene)
# 
# ggplot(data = plt.df, mapping = aes(x = gene, y = dist)) + 
#   geom_bar(stat = "identity", width = 0.2) + #, width = 0.8, fill = "cyan4", color = "gray2") + 
#   theme_light() + ggtitle(colnames(X1n[, 1]))



##################################################################################################

#### Get data for a biomarker...
biomarkers <- colnames(Ydata1);       # printf("Biomarkers = ", biomarkers)
q <- length(biomarkers)

bmChosen <- biomarkers[k];            printf("\nChosen biomarker = %s", bmChosen)
ranks <- cbind(rank1[, bmChosen], rank2[, bmChosen], rank3[, bmChosen])
  
## Top 'm' common genes...
m_opt <- 150;       nGN <- 300;               gnRank <- intersect(ranks[1:nGN, 2], ranks[1:nGN, 3])
nI <- 0;            m0 <- length(gnRank);     m <- m0
while(m < m_opt){
  nI <- nI + 1;     nGN <- nGN + 100
  gnRank <- intersect(ranks[1:nGN, 2], ranks[1:nGN, 3]);      m <- length(gnRank)
}
gnRank <- sort(gnRank, decreasing = FALSE)
printf("#top genes chosen = %d (nGN = %d, nI = %d, m0 = %d)", m, nGN, nI, m0)

gnRank <- get.top.genes(ranks[, 2:3], m.top = 150, print.opt = TRUE)

X1 <- Xdata1[, gnRank];                   X2 <- rbind(Xdata2[, gnRank], Xdata3[, gnRank])
Y1 <- norm01(Ydata1[, bmChosen]);         Y2 <- norm01(c(Ydata2[, bmChosen], Ydata3[, bmChosen]))

X1n <- dapply(X1, MARGIN = 2, norm01);    X2n <- dapply(X2, MARGIN = 2, norm01)
m <- if (ncol(X1n) == ncol(X2n)) ncol(X1n)


## Histogram matching for predictors...
X1n <- dapply(X1, MARGIN = 2, norm01);          X2n <- dapply(X2, MARGIN = 2, norm01)
X2n.mapped <- lapply(1:m, function(j){
  ## Calculate CDFs...
  X1j.cdf <- get.dist(X1n[, j], N = 1e3, dist.opt = "kde")
  X2j.cdf <- get.dist(X2n[, j], N = 1e3, dist.opt = "kde")
  
  ## Get mapping function...
  kn <- knots(X2j.cdf);     fn <- X2j.cdf(kn)
  map <- approxfun(x = fn, y = kn, method = "linear", yleft = 0, yright = 1, ties = "ordered")
  
  ## Mapped values...
  X2j.mapped <- map(X2j.cdf(X1n[, j]))
})
X2n.mapped <- as.data.frame(X2n.mapped, col.names = colnames(X1n), row.names = rownames(X1n))


## Build predictive model & predict response...
# set.seed(0)
RF <- randomForest(x = X2n, y = Y2, ntree = 200)
Y2.pred.mapped <- predict(RF, X2n.mapped)
Y2.pred.mapped[Y2.pred.mapped < 0] <- 0;    Y2.pred.mapped[Y2.pred.mapped > 1] <- 1


## Histogram matching for response...
Y1.cdf <- get.dist(Y1, N = 1e3, dist.opt = "kde")
Y2.cdf <- get.dist(Y2, N = 1e3, dist.opt = "kde")

kn.y <- knots(Y1.cdf);    fn.y <- Y1.cdf(kn.y)
map.y <- approxfun(x = fn.y, y = kn.y, method = "linear", yleft = 0, yright = 1, ties = "ordered")

Y1.pred <- map.y(Y2.cdf(Y2.pred.mapped))


#### Baseline model...
# set.seed(0)
RF.base <- randomForest(x = X2n, y = Y2, ntree = 200)
Y1.pred.base <- predict(RF.base, X1n)
Y1.pred.base[Y1.pred.base < 0] <- 0;    Y1.pred.base[Y1.pred.base > 1] <- 1


## Results table...
results <- data.frame("DMTL" = c(calc.err(Y1, Y1.pred, measure = "NRMSE"), calc.err(Y1, Y1.pred, measure = "NMAE"), 
                                 cor(Y1, Y1.pred, method = "spearman")), 
                      'BL' = c(calc.err(Y1, Y1.pred.base, measure = "NRMSE"), calc.err(Y1, Y1.pred.base, measure = "NMAE"), 
                               cor(Y1, Y1.pred.base, method = "spearman")), 
                      row.names = c("NRMSE", "NMAE", "SCC"))

printf("Results = \n");   print(results)


##
# gg.pp <- list()
# gg.pp[["Tumor"]] <- ggplot(MB.df, aes_string(x = bmChosen)) + geom_density(fill = "cyan4", alpha = 0.7) + theme_light() + 
#   xlab("") + ylab("") + ggtitle("Tumor") + theme(plot.title = element_text(hjust = 0.5))
# gg.pp[["Cell line"]] <- ggplot(CG.df, aes_string(x = bmChosen)) + geom_density(fill = "brown4", alpha = 0.7) + theme_light() + 
#   xlab("") + ylab("") + ggtitle("Cell line") + theme(plot.title = element_text(hjust = 0.5))
# gg.pp[["ncol"]] <- 2
# annotate_figure(do.call(ggarrange, gg.pp), top = text_grob(bmChosen, face = "bold"))
