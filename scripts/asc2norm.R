#Use Poisson regression remove systematic biases in Hi-C cis contact maps
#Ming Hu (minghu@fas.harvard.edu)
#Last update: 08.05.2012

args = commandArgs(trailingOnly=TRUE)

input=args[1]
feat=args[2]
chr=args[3]

output = sub("asc","norm.asc",input)

#read in input file
#u<-read.table("CAC_1.chr1.asc",skip=1) 
#print(c("input:", input));
#h <- readLines(input, n=1) 
u<-read.table(input,skip=1)
h=length(u)
if (h<2){
    h1=sprintf("%d %d 1 1",h,h)
    write.table(h1,output, row.names=F, col.names=F, sep=" ", quote=F)
    write.table(u, output, row.names=F, col.names=F, sep=" ", quote=F)
    quit()
}
#v<-read.table("mm10.AAGCTT.100K.chr1.txt",head=T) #user can change the name of this input file
#v<-read.table(feat,head=T) #user can change the name of this input file
#print(c("feature:", feat));
v<-read.table(feat,head=F) #user can change the name of this input file
colnames(v) <- c('chr', 'bin1', 'bin2', 'len', 'gcc', 'map')
v<-v[which(v$chr==chr),]
#head(v)
#remove bins with low effective length, low GC content and low mappability
#plot(density(v$len)) # len > 10000
#plot(density(v$gcc)) # gcc > 0.3
#plot(density(v$map)) # map >0.6

index<-I( v$len>10000 & v$gcc>0.3 & v$map>0.6 )
w<-which(v$len>10000 & v$gcc>0.3 & v$map>0.6)
inp_len=nrow(u)
#w
#inp_len

u<-u[ index==1, index==1 ]
v<-v[ index==1 ,]

#change matrix into vector
u<-as.matrix(u)
u_vec<-u[upper.tri(u,diag=F)]

#get cov matrix
len_m<-as.matrix(log(v[,4]%o%v[,4]))
gcc_m<-as.matrix(log(v[,5]%o%v[,5]))
map_m<-as.matrix(log(v[,6]%o%v[,6]))

#centralize cov matrix of enz, gcc
len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))

#change matrix into vector
len_vec<-len_m[upper.tri(len_m,diag=F)]
gcc_vec<-gcc_m[upper.tri(gcc_m,diag=F)]
map_vec<-map_m[upper.tri(map_m,diag=F)]

#fit Poisson regression: u~len+gcc+offset(map)
fit<-glm(u_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")

#user can use the following two lines to fit negative binomial regression: u~len+gcc+offset(map).
#library("MASS")
#fit<-glm.nb(u_vec~len_vec+gcc_vec+offset(map_vec))

#summary(fit)
coeff<-round(fit$coeff,4)
res<- round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)

#restore matrix
res2=mat.or.vec(inp_len, inp_len)
res2[w, w]=res

#output normalized cis contact map, user can change the name of this output file
#write.table(res2, file="output_normalized_hic_cis_contact_map.txt", row.names=F, col.names=F, sep="\t", quote=F)
h1=sprintf("%d %d 1 1",h,h)
write.table(h1, output, row.names=F, col.names=F, sep=" ", quote=F)
write.table(res2, output, row.names=F, col.names=F, sep=" ", quote=F)

