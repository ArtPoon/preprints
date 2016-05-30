require(ape, quietly=TRUE)
require(GapProcedure, quietly=TRUE)

# modify this as necessary
files <- Sys.glob('/Users/art/git/papers/clustrev/data/*.*_TRUE.fas')

res <- {}
for (file in files) {
	items <- strsplit(file, split='/')[[1]]
	filename <- items[length(items)]
	items <- strsplit(filename, split='\\.')[[1]]
	treatment <- items[1]
	rep <- as.integer(substr(items[2], 1, 1))
	print (treatment, rep)
	
	aln <- read.dna(file, format="fasta", comment.char=" ", as.matrix=FALSE)
	
	gp <- GapProcedure(as.matrix(aln), 'aK80')
	group.counts <- table(gp$class)
	labels <- gsub("[ ]+$", "", names(aln))
	
	# true condition, predicted condition
	tab <- table(grepl("_1_", names(aln)), group.counts[gp$class]>1)
	tn <- tab[1,1]
	fp <- tab[1,2]
	fn <- tab[2,1]
	tp <- tab[2,2]
	fpr <- fp / sum(tab[1,])
	tpr <- tp / sum(tab[2,])
	res <- rbind(res, c(treatment, rep, tn, fp, fn, tp, fpr, tpr))
}

df <- as.data.frame(res)
names(df) <- c('treatment', 'rep', 'true.negatives', 'false.positives', 'false.negatives', 'true.positives', 'FPR', 'TPR')
df$true.negatives <- as.integer(as.character(df$true.negatives))
df$false.positives <- as.integer(as.character(df$false.positives))
df$false.negatives <- as.integer(as.character(df$false.negatives))
df$true.positives <- as.integer(as.character(df$true.positives))

df$FPR <- as.double(as.character(df$FPR))
df$TPR <- as.double(as.character(df$TPR))

write.csv(df, file='~/git/papers/clustrev/results/gapprocedure.csv', quote=FALSE, row.names=FALSE)


#write.csv(res, file=output.csv, quote=FALSE, row.names=FALSE)
boxplot(split(df$TPR, df$treatment))
boxplot(split(df$FPR, df$treatment))

par(mar=c(5,5,2,2))
temp <- df[df$treatment=='Both',]
plot(temp$FPR, temp$TPR, xlim=c(0,1), ylim=c(0,1), col='red', pch='B', xlab='False positive rate', ylab='True positive rate', cex.lab=1.5, cex.axis=1.2)
idx <- chull(temp$FPR, temp$TPR)
polygon(temp$FPR[idx], temp$TPR[idx], border=NA, col=rgb(1,0,0,0.1))

abline(a=0, b=1, lty=2)

temp <- df[df$treatment=='fastTrans',]
points(temp$FPR, temp$TPR, col='blue', pch='T')
idx <- chull(temp$FPR, temp$TPR)
polygon(temp$FPR[idx], temp$TPR[idx], border=NA, col=rgb(0,0,1,0.1))

temp <- df[df$treatment=='fastSamp',]
points(temp$FPR, temp$TPR, col='forestgreen', pch='S')
idx <- chull(temp$FPR, temp$TPR)
polygon(temp$FPR[idx], temp$TPR[idx], border=NA, col=rgb(0,0.8,0,0.1))

legend(x=0.6, y=0.25, legend=c('Both', 'Transmission', 'Sampling'), pch=c('B', 'T', 'S'), col=c('red', 'blue', 'forestgreen'))
