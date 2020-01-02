### asc2di.R
### Calculate Directionality Score. See Dixon et al. 
### Written by Yanxiao Zhang <zhangyx.shawn@gmail.com>
#############################
options(scipen = 999)

args=commandArgs(trailing=T)
input = args[1]
chr = args[2]
bin_size = as.integer(args[3])
bin_num = as.integer(args[4])
output = args[5]

## read the matrix
mat = read.table(input,skip=1)

if (nrow(mat)<bin_num) {
    system(paste0("touch ",output))
    quit()
    }

calc_DI = function(vec,bin_num,idx){
ll = ifelse(idx-bin_num>1,idx-bin_num,1)
lr = idx-1 
if (lr>=ll) {
A = sum(vec[ll:lr])
} else { A=0 }
rl = idx+1
rr = ifelse(idx+bin_num<=length(vec),idx+bin_num,length(vec))
if (rr>=rl){
B = sum(vec[rl:rr])
} else { B = 0 }

E = (A+B)/2
#DI score
DI = ifelse(E==0,0,sign(B-A)* ( (A-E)**2/E + (B-E)**2/E))
#print(c(A,B,E,DI))
DI
}

DIs = sapply(1:nrow(mat), function(x){calc_DI(mat[x,],bin_num,x)})

chrs = rep(chr,nrow(mat))
start = seq(from=0,by=bin_size,length.out=nrow(mat))
end = seq(from=bin_size, by=bin_size, length.out=nrow(mat))

table = data.frame(chrs,start,end, DIs)

write.table(table, output,col.names=F,row.names=F,quote=F,sep='\t')


