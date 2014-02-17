n=commandArgs(trailing=TRUE)

x=read.table(n[1],header=T)
z=which(x\$cprob >= 0.999)
y=x\$distance[z[1]]
write(y,n[2])
