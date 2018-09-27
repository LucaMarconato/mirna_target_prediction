source("my_device.r")

a <- read.table("processed/scored_interactions_processed.tsv", colClasses = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"), header = T)

new_maximized_device()
par(mfrow = c(2,3))

hist(2^(1 - a$context_score))

my_hist <- hist(a$context_score, plot = F)
my_hist$counts <- log10(my_hist$counts + 0.00000000001)
plot(my_hist, ylab = 'log10(counts)', ylim = c(0,7))

my_hist <- hist(2^(1 - a$context_score[2^(1 - a$context_score) > 0.9]), breaks = 100, plot = F)
my_hist$counts <- log10(my_hist$counts + 0.00000000001)
plot(my_hist, ylab = 'log10(counts)', ylim = c(0,7))

hist(2^(1 - a$weighted_context_score))

my_hist <- hist(a$weighted_context_score, plot = F)
my_hist$counts <- log10(my_hist$counts + 0.00000000001)
plot(my_hist, ylab = 'log10(counts)', ylim = c(0,7))

my_hist <- hist(2^(1 - a$weighted_context_score[2^(1 - a$weighted_context_score) > 0.9]), breaks = 100, plot = F)
my_hist$counts <- log10(my_hist$counts + 0.00000000001)
plot(my_hist, ylab = 'log10(counts)', ylim = c(0,7))
