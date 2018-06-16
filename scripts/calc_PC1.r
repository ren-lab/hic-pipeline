# oe matrix
suppressMessages(require(data.table))
suppressMessages(require(GenomicRanges))
a=fread(commandArgs(trailing=T)[1],skip=1)
alist = as.list(a)

# correlation matrix
b=matrix(0,nrow=nrow(a),ncol=ncol(a))

rsums = rowSums(a,na.rm=T)
idx = which(rsums > 0)

for (i.idx in 1:length(idx)){
  for ( j.idx in i.idx:length(idx) ) {
  i = idx[i.idx]
  j = idx[j.idx]
#  print(c(i,j))
  b[i,j] = b[j,i] =  cor(alist[[i]],alist[[j]])

  }
  }

b2 = b[idx,idx]

pca2 = prcomp(b2)
pc1 = pca2$x[,1]
chr = commandArgs(trailing=T)[2]
start = (idx-1)*50000 
end = start + 50000
out = data.frame(chr,start,end,pc1)

b1=GRanges(out[which(out$pc1>0),1:3])
b2=GRanges(out[which(out$pc1<0),1:3])
tss = fread(commandArgs(trailing=T)[3])
colnames(tss) = c("chr","start","end","gene",".","strand")
tss2 = GRanges(tss[,1:3])

i1 = intersect(b1,tss2)
i2 = intersect(b2,tss2)
FE = sum(width(i1))/sum(width(b1)) / ( sum(width(i2))/sum(width(b2)) )

if (FE < 1) { out$pc1 = out$pc1 * -1 }

write.table(out,row.names=F,col.names=F,quote=F,sep='\t')

