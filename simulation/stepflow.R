library(gridExtra)
library(colorRamp2)
library(RColorBrewer)
library(grid)
library(ComplexHeatmap)
source("/Users/M216453/Documents/Mayo_Research/2022_08_04_pairwiseFDR/tdGBH-paper/simulation/func.R")
set.seed(123)
feature =20; outcome=20; signal.density= 0.1; affected.feature = 0.5; affected.outcome = 0.5
feature.no = feature * outcome; signals <- feature.no * signal.density
space <- affected.feature * affected.outcome * outcome * feature
affected.outcome.signal.density <- (feature.no* signal.density) / (affected.outcome * outcome) /feature #*affected.feature
outcome.signal.densities <- sample(c(rep(affected.outcome.signal.density, round(outcome * affected.outcome)), rep(0, round(outcome *(1-affected.outcome)))))
sim <- simulata.data(outcome =outcome, feature.no = feature.no,cor.struct ='None',signal.strength = 2,
                     feature.density=affected.feature, signal.densities = outcome.signal.densities)

p.mat <-  sim$p.mat
truth.mat <- sim$truth.mat
mean(truth.mat)
rownames(p.mat) <- paste0('f', 1:nrow(p.mat))
colnames(p.mat) <- paste0('o', 1:ncol(p.mat))
pdf("/Users/M216453/Documents/Mayo_Research/2022_08_04_pairwiseFDR/revision1/Figures/p.pdf", width = 4, height = 3.5)
Heatmap(p.mat, cluster_rows = FALSE, cluster_columns = FALSE,
        col = colorRamp2(c(0, median(p.mat, na.rm = TRUE), max(p.mat, na.rm = TRUE)), 
                         brewer.pal(9, 'Purples')[c(9,5,1)]),
        name = "pvalue",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(p.mat[i,j] <= 0.05) {
            grid.text('*', x = x, y = y, 
                      gp = gpar(col = "white"), # Set the color of the asterisk to white
                      just = "center") # Center the asterisk in the cell
          }
        })
dev.off()

pi0.method = 'lsl'; global.pi0.method = 'storey'; weight.method = 'new'; shrink = 0.1; renorm=F
pi0 <- estimate.pi0(as.vector(p.mat), method =  global.pi0.method)
pi0

apply(p.mat, 2, function(x) estimate.pi0(x, method = 'storey'))
apply(p.mat, 1, function(x) estimate.pi0(x, method = 'storey'))

pi0.o <- apply(p.mat, 2, function(x) estimate.pi0(x, method = pi0.method))
pi0.f <- apply(p.mat, 1, function(x) estimate.pi0(x, method = pi0.method))

pdf("/Users/M216453/Documents/Mayo_Research/2022_08_04_pairwiseFDR/revision1/Figures/pi0.o.pdf")
g <- tableGrob(cbind(pi0.o))
tg <- grid.draw(g)
dev.off()

pdf("/Users/M216453/Documents/Mayo_Research/2022_08_04_pairwiseFDR/revision1/Figures/pi0.g.pdf")
g <- tableGrob(cbind(pi0.f))
tg <- grid.draw(g)
dev.off()

pi0.o <- (1 - shrink) * pi0.o + shrink * pi0
pi0.f <- (1 - shrink)  * pi0.f + shrink * pi0

pdf("/Users/M216453/Documents/Mayo_Research/2022_08_04_pairwiseFDR/revision1/Figures/pi0.o.S.pdf")
g <- tableGrob(cbind(pi0.o))
tg <- grid.draw(g)
dev.off()

pdf("/Users/M216453/Documents/Mayo_Research/2022_08_04_pairwiseFDR/revision1/Figures/pi0.g.S.pdf")
g <- tableGrob(cbind(pi0.f))
tg <- grid.draw(g)
dev.off()


pi0.o.mat <- t(matrix(pi0.o, nrow = length(pi0.o), ncol = length(pi0.f)))
pi0.f.mat <- matrix(pi0.f, nrow = length(pi0.f), ncol = length(pi0.o))


if(weight.method == 'geo') pi0.og.mat <- sqrt(pi0.o.mat * pi0.g.mat)
if(weight.method == 'ari') pi0.og.mat <- (pi0.o.mat + pi0.g.mat) / 2
if(weight.method == 'new') {
  sd.r <- (sd(pi0.o) / sqrt(length(pi0.o)))  / (sd(pi0.o) / sqrt(length(pi0.o)) + sd(pi0.f) / sqrt(length(pi0.f)))
  pi0.og.mat <- sqrt((pi0.o.mat^(2 * sd.r)) * (pi0.f.mat^(2 * (1 - sd.r))))
}

pi0 <- mean(pi0.og.mat)
ws.og.mat <- (1 - pi0.og.mat) / pi0.og.mat
rownames(ws.og.mat) <- paste0('f', 1:nrow(ws.og.mat))
colnames(ws.og.mat) <- paste0('o', 1:ncol(ws.og.mat))
pdf("/Users/M216453/Documents/Mayo_Research/2022_08_04_pairwiseFDR/revision1/Figures/ws.og.mat.pdf", width = 4, height = 3.5)
Heatmap(ws.og.mat, cluster_rows = F, cluster_columns = F,
        col = colorRamp2(c(0, median(ws.og.mat), max(ws.og.mat)), brewer.pal(9, 'Blues')[c(1,5,9)]),
        rect_gp = gpar(col = "white", lwd = 1),name = "weights")
dev.off()


if (renorm == TRUE) {
  ws <- ws.og.mat / (1 - pi0)
  ws <- length(p.mat) / sum(ws) * ws
  p.ws.mat <- p.mat / ws
} else {
  p.ws.mat <- p.mat / ws.og.mat * (1 - pi0)
}

pdf("/Users/M216453/Documents/Mayo_Research/2022_08_04_pairwiseFDR/revision1/Figures/p.ws.mat.pdf", width = 4.5, height = 3.5)
Heatmap(p.ws.mat, cluster_rows = F, cluster_columns = F,
        col = colorRamp2(c(0, median(p.ws.mat), max(p.ws.mat)), brewer.pal(9, 'Purples')[c(1,5,9)]),
        rect_gp = gpar(col = "white", lwd = 1),name = "weighted pvalue")
dev.off()


p.adj <- matrix(p.adjust(as.vector(p.ws.mat), 'BH'), length(pi0.f), length(pi0.o), dimnames = list(rownames(p.mat),colnames(p.mat)))

pdf("/Users/M216453/Documents/Mayo_Research/2022_08_04_pairwiseFDR/revision1/Figures/padj.pdf", width = 4, height = 3.5)
Heatmap((p.adj), cluster_rows = F, cluster_columns = F,
        col = colorRamp2(c(0, median(p.adj), max(p.adj)), brewer.pal(9, 'Oranges')[c(9,4,1)]),
        name = "qvalue",cell_fun = function(j, i, x, y, width, height, fill) {
          if(p.adj[i,j] <= 0.05) {
            grid.text('*', x = x, y = y, r = unit(1/30,'cm'))
          }
        })
dev.off()

