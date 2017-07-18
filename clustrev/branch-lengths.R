nodes <- read.csv('~/git/papers/clustrev/scripts/extract-comments.csv', header=T)
nodes$is.internal <- (nodes$is.internal == 'True')

# normalize branch lengths by mean per tree
mean.bl <- sapply(split(nodes$branch.length, f=list(nodes$rep, nodes$treatment), drop=TRUE), mean)

n.bl <- sapply(split(nodes$branch.length, f=list(nodes$rep, nodes$treatment), drop=TRUE), length)

nodes$norm.bl <- nodes$branch.length# / rep(mean.bl, times=n.bl)



# evaluate differences between subpopulations by treatment and replicate
inner <- nodes[nodes$is.internal, ]
inner$treatment <- relevel(inner$treatment, 'Control')

fit <- glm(branch.length ~ location*treatment, data=inner, family='Gamma')


outer <- nodes[!nodes$is.internal, ]
outer$treatment <- relevel(outer$treatment, 'Control')
fit.2 <- glm(branch.length ~ location * treatment, data=outer, family='Gamma')

#require(lme4)
#fit <- lmer(branch.length ~ location * treatment + 1|rep, data=inner)
# fit.1 <- glmer(branch.length ~ location + (location|rep), data=inner, family='Gamma')

# calculate means
#temp <- nodes[nodes$treatment=='Control' & nodes$rep==0, ]
temp <- sapply(split(nodes$norm.bl[nodes$is.internal], f=list(nodes$location[nodes$is.internal], nodes$rep[nodes$is.internal], nodes$treatment[nodes$is.internal])), mean)
m <- matrix(unlist(strsplit(names(temp), split='\\.')), ncol=3, byrow=TRUE)
df <- data.frame(treatment=m[,3], rep=m[,2], location=m[,1], inner.bl=temp)
temp <- sapply(split(nodes$norm.bl[!nodes$is.internal], f=list(nodes$location[!nodes$is.internal], nodes$rep[!nodes$is.internal], nodes$treatment[!nodes$is.internal])), mean)
df$outer.bl <- temp


# generate plot
par(mar=c(5,5,2,2), cex=1)

temp <- df[df$treatment=='Both', ]
plot(temp$inner.bl, temp$outer.bl, pch=ifelse(temp$location == 0, 'B', 'b'), col='red', xlab='Internal branch lengths', ylab='Terminal branch lengths', cex.lab=1.5, cex.axis=1.2, xlim=range(df$inner.bl), ylim=range(df$outer.bl))
idx <- chull(temp$inner.bl, temp$outer.bl)
polygon(temp$inner.bl[idx], temp$outer.bl[idx], border=NA, col=rgb(1,0,0,0.1))

temp <- df[df$treatment=='Control', ]
points(temp$inner.bl, temp$outer.bl, pch=ifelse(temp$location==0, 'C', 'c'), col='grey60')
idx <- chull(temp$inner.bl, temp$outer.bl)
polygon(temp$inner.bl[idx], temp$outer.bl[idx], border=NA, col=rgb(0,0,0,0.05))

temp <- df[df$treatment=='fastTrans', ]
points(temp$inner.bl, temp$outer.bl, pch=ifelse(temp$location==0, 'T', 't'), col='blue')
idx <- chull(temp$inner.bl, temp$outer.bl)
polygon(temp$inner.bl[idx], temp$outer.bl[idx], border=NA, col=rgb(0,0,1,0.1))

temp <- df[df$treatment=='fastSamp', ]
points(temp$inner.bl, temp$outer.bl, pch=ifelse(temp$location==0, 'S', 's'), col='forestgreen')
idx <- chull(temp$inner.bl, temp$outer.bl)
polygon(temp$inner.bl[idx], temp$outer.bl[idx], border=NA, col=rgb(0,1,0,0.1))

legend(x=0.024, y=0.175, legend=c('Control', 'Transmission', 'Sampling', 'Both'), col=c('grey60', 'blue', 'forestgreen', 'red', 'black'), pch=c('C', 'T', 'S', 'B'), cex=0.9, bg='white')




