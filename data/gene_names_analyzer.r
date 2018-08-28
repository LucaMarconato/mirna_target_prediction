a <- read.table("processed/sites_with_scored_interactions.tsv")
print(dim(a)[[1]])
b <- a[,c(1,2,3)]
c <- unique(b)
dim(c)
n <- dim(c)[[1]]
m1 <- dim(unique(c[1]))[[1]]
print(m1)
m2 <- dim(unique(c[2]))[[1]]
print(m2)
m3 <- dim(unique(c[3]))[[1]]
print(m3)
m <- min(m1, m2, m3)
if(m == n) {
    print("each of the first three columns can be taken as the key, at least considering this set of genes")
} else {
    print("find the rows occurring more than once")
}
