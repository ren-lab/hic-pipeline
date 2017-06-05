### R script to cpmpute insulatin scores from HiC matrix 
### Adapt ideas from https://github.com/blajoie/crane-nature-2015/blob/master/scripts/matrix2insulation.pl
### Written by Yunjiang Qiu
###########################
library("optparse")
options(scipen = 999)

option_list = list(
  make_option(c("-m", "--matrix"), type="character", default=NULL, 
              help="HiC matrix", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file name", metavar="character"),
  make_option(c("-w", "--win"), type="character", default=NULL, 
              help="window_size", metavar="character"),
  make_option(c("-c", "--chr"), type="character", default=NULL, 
              help="chr", metavar="character"),
  make_option(c("-b", "--bin"), type="character", default=NULL, 
              help="bin_size", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

bin <- gsub("Kb","000",opt$bin)
bin <- as.numeric(gsub("Mb","000000",bin))
win <- gsub("Kb","000",opt$win)
win <- as.numeric(gsub("Mb","000000",win))

mat <- opt$matrix
res <- win/bin

dat <- as.matrix(read.table(mat,head=F,skip=1))
## if length of dat is less than window size, then quit:
ins <- data.frame(chr = rep(opt$chr,nrow(dat)),
                  start = seq(0,bin*(nrow(dat)-1),bin),
                  end = seq(bin,bin*nrow(dat),bin),
                  ins = rep(NA,nrow(dat)))
if (nrow(dat)-res >= res+1) {
for (i in (res+1):(nrow(dat)-res))
{
  ins$ins[i] <- mean(dat[(i+1):(i+res),(i-res):(i-1)])
}
}
ins$ins <- log2(ins$ins/mean(na.omit(ins$ins)))

write.table(ins,opt$out,quote=F,row.name=F,col.name=F,sep = "\t")
