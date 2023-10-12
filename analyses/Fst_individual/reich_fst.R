#############################################
fstHud <- function(g1,g2,ratio=F){
    n1 <- length(g1)*2
    n2 <- length(g2)*2
    p1 <- sum(g1)/n1
    p2 <- sum(g2)/n2
    n <- (p1-p2)^2-p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    d <- p1*(1-p1)+p2*(1-p2)
    if(ratio)
        return(n/d)
    else
        return(c(n,d))

  }
fstReich <- function(g1,g2,ratio=F){
    n1 <- length(g1)*2
    n2 <- length(g2)*2
    a1 <- sum(g1)
    a2 <- sum(g2)
    h1<- a1/n1*(n1-a1)/(n1-1)
    h2<- a2/n2*(n2-a2)/(n2-1)

    n <- (a1/n1-a2/n2)^2-h1/n1-h2/n2
    d <- n+h1+h2
    if(ratio)
        return(n/d)
    else
        return(c(n,d))

  }


num <- c()
denom <- c()
for(g1 in 0:2){
    for(g2 in 0:2){
        a=fstReich(g1, g2)
        num  <- c(num, a[1])
        denom <- c(denom, a[2])
    }
}

print(num)
print(denom)
l <- list()
n <-  as.character(read.table("2dsfs.names")[,1])
sfsS  <- read.table("2dsfs.data")

for(idx in 1:length(n)){
    vals <- sum(sfsS[idx, ] * num) / sum(sfsS[idx, ] * denom)
    l[idx] <- vals
}

res <- do.call(rbind, l)

master <- read.table("sampleinfo.txt", header=T, sep="\t")
#masterFile <- "sampleinfo.txt"
#master <- gsheet::gsheet2tbl(masterFile)
df <- as.data.frame(master)

inds <- df$Sample
locality <- df$Location
country <- df$Region

meta <- data.frame(inds, locality, country)
rownames(meta) <- inds
splitN <- do.call(rbind,lapply(n, function(x){unlist(strsplit(x, "_"))}))

res2 <- data.frame(splitN[,1], meta[splitN[,1], 2], splitN[,2], meta[splitN[,2], 2], res)
write.table(res2, file="reich.res", quote=F, col.names=F, row.names=F)
