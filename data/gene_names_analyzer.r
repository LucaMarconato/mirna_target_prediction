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
    print("finding the rows where at least one of the elements in the first three columns are ambiguous")
}
print("it happens that only m2 != n, let us study m2")
d <- duplicated(c[2])
indexes <- which(d == TRUE)
elements <- c[[2]][indexes]
for(e in elements) {
    print(sprintf("rows with the value of the second column equal to %s", e))
    same_elements <- which(c[[2]] == e)    
    print(c[same_elements,])
}

## it appears that there are some errors in the data (see here: https://www.biostars.org/p/243381/
## pay attention with the genes printed by the above code
