# the following code is used to infer which is the correspondence between Seed.match and Site.Type

## # uncomment only the first time you run this code, since it is slow
## t0 <- read.table("raw/Conserved_Family_Info.txt", sep = "\t", header = T)
## t1 <- read.table("raw/Conserved_Site_Context_Scores.txt", sep = "\t", header = T)
## t2 <- read.table("raw/Nonconserved_Family_Info.txt", sep = "\t", header = T)
## t3 <- read.table("raw/Nonconserved_Site_Context_Scores.txt", sep = "\t", header = T)
## save(t0, t1, t2, t3, file = "site_type_analyzer.RData")
load("site_type_analyzer.RData")
print(dim(t0)[[1]])
print(dim(t1)[[1]])
print(dim(t2)[[1]])
print(dim(t3)[[1]])
t0 <- t0[t0$Species.ID == "9606",]
t1 <- t1[t1$Gene.Tax.ID == "9606",]
t2 <- t2[t2$Species.ID == "9606",]
t3 <- t3[t3$Gene.Tax.ID == "9606",]
print(dim(t0)[[1]])
print(dim(t1)[[1]])
print(dim(t2)[[1]])
print(dim(t3)[[1]])
seed_type0 <- t0$Seed.match
seed_type1 <- t1$Site.Type
seed_type2 <- t2$Seed.match
seed_type3 <- t3$Site.Type
table0 <- table(seed_type0)
table1 <- table(seed_type1)
table2 <- table(seed_type2)
table3 <- table(seed_type3)
seeds0 <- c("8mer", "7mer-m8", "7mer-a1")
sum0 <- table0[seeds0] + table2[seeds0]
seeds1 <- c("1","2","3")
sum1 <- table1[seeds1] + table3[seeds1]

par(mfrow = c(1,2))
barplot(sum0)
barplot(sum1)
par(mfrow = c(1,1))

# from the graph it seems that the correspondence is 
# 1 -> 7mer-a1
# 2 -> 7mer-m8
# 3 -> 8mer
# let us validate this conjecture by looking for corresponding rows in the data frames

t0 <- transform(t0, miR.Family = sprintf("hsa-%s", miR.Family)) 
t2 <- transform(t2, miR.Family = sprintf("hsa-%s", miR.Family))
df = merge(t2, t3, by.x=c("miR.Family", "Gene.ID", "Gene.Symbol", "Transcript.ID", "UTR.start", "UTR.end"), by.y=c("miRNA", "Gene.ID", "Gene.Symbol", "Transcript.ID", "UTR_start", "UTR.end"))
print(dim(df))
df[1:30,c("Seed.match","Site.Type")]
# this confirms the thesis

