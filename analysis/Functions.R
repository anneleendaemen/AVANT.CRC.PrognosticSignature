panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
	usr <- par("usr"); on.exit(par(usr))
	par(usr=c(0,1,0,1))
	r <- abs(cor(x,y, method="pearson"))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex.cor <- 1.1/strwidth(txt)
	text(0.5, 0.5, txt, cex=cex.cor*r)
}

panel.smooth<-function (x, y, col = "black", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "darkred", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}

cor.values <- function(tissue, method="pearson") {
	cor.genes <- vector("list", length(goi))
	names(cor.genes) <- goi
	for (g in goi) {
		cor.genes[[g]] <- cor(data.scaled[match(g, rownames(data.scaled)), tissues %in% tissue], avg.no[[g]][tissues %in% tissue], method=method)
	}
	cors <- unlist(cor.genes)
	return(cors)
}

cor.values.bysubtype <- function(tissue, subtype, method="pearson") {
	cor.genes <- vector("list", length(goi))
	names(cor.genes) <- goi
	for (g in goi) {
		cor.genes[[g]] <- cor(data.scaled[match(g, rownames(data.scaled)), tissues %in% tissue & CMS %in% subtype], avg.no[[g]][tissues %in% tissue & CMS %in% subtype], method=method)
	}
	cors <- unlist(cor.genes)
	return(cors)
}
