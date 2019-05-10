# task: conduct an armitage test for every SNP (row) in the example dataset.
# get example data
armitage_example <- readRDS("R_dev/armitage_example.RDS")
head(genos)


# here's a function I wrote to do this for a single SNP, taking data in a 2x3 matrix:
# see Armitage (1955), Tests for Linear Trends in Proportions and Frequencies, Biometrics for equations

# table with cases on top, controls on bottom, then genotypes pp, pq, and qq in columns 1, 2, 3
m <- matrix(c(19, 29, 24, 497, 560, 269), 2, byrow = T) # example table, you probably won't need to remake this.
colnames(m) <- c("pp", "pq", "qq")
rownames(m) <- c("case", "control")
wv <- c(0,1,2) # these are weights. These should be 0 for all "pp", 1 for all "pq", and 2 for all "qq" genotypes


# the function
calc_single_amitage <- function(m, w){

  # define variables, names relate to Armitage 1955
  n <- m[1,] # cases
  N <- colSums(m) # Number of individuals with each genotype
  bT <- sum(m) # total number of individuals sequenced at this loci
  t <- sum(m[1,]) # total number of "case" individuals

  # get the sums that go into the equations, see the paper
  s1 <- sum(n*wv)
  s2 <- sum(N*wv)
  s3 <- sum(N*wv^2)

  # equation 5
  b <- (bT*s1 - t*s2)/(bT*s3 - (s2^2))

  # equation 6
  Vb <- (t*(bT - t))/(bT*(bT*s3 - s2^2))

  # equation 7
  chi <- (b^2)/Vb

  # use the pchisq function with 1 degree of freedom to get a p-value and return.
  return(pchisq(chi, 1, lower.tail = F))
}

# the output, returns a p-value
calc_single_amitage(m, wv)




# you'll need to expand this function to work on the whole cast_gs object!
# I haven't done this yet, but my tips are:
#    1) split the data into two sets, one with cases and one with controls
#    2) figure out weights for each column (pp = 0, pq = 1, qq = 2)
#    3) do the calculation. Remember that any column with a 0 in it won't effect the resulting p-values at all,
#       so there's no reason to remove them or worry about them for each row!


#----------------------------------------------------
#spliting case&control:
case <- armitage_example[,seq(from = 1, to = ncol(armitage_example), by = 2)]
control <- armitage_example[,2*c(1:10)]

# define variables
n <- case
N <- case + control
colnames(N) <- substr(colnames(N), start=1, stop=2)
bT <- rowSums(armitage_example)
t <- rowSums(n)

# select pp,pq and qq for case
homozygote <- case[, which(substr(colnames(case), start=1, stop=1) == substr(colnames(case), start=2, stop=2))]# adjust to pick out homozygotes by checking if the first allele matches the second allele
# case[,which(hom.status)]
pp.case <- matrixStats::rowMaxs(homozygote)
qq.case <- matrixStats::rowSums2(homozygote) - matrixStats::rowMaxs(homozygote)
pq.case <- t-pp.case-qq.case

n <- cbind(pp.case, pq.case, qq.case) # why use n again?


#select pp,pq and qq for all
homozygote <-N[, which(substr(colnames(case), start=1, stop=1) == substr(colnames(case), start=2, stop=2))]
pp <- matrixStats::rowMaxs(homozygote)
qq <- matrixStats::rowSums2(homozygote) - matrixStats::rowMaxs(homozygote)
pq <- rowSums(N)-pp-qq

m <- cbind(pp, pq, qq) #order of m is wrong


wv <- c(0,1,2)
sum1 <- t(n)*wv # transposed because matrices will add / multiply, whatever by column rather than by row
s1 <- colSums(sum1)
s2 <- colSums(t(m)*wv)


# get the other sums, then b, Vb, chi

# wrap in a function




# after this is done, we'll start work on the Excoffier method of estimating haplotype frequencies. It'll be tough!
# here a paper describing the method, I'll add an example file to work with soon.
# Maximum-Likelihood Estimation of Molecular Haplotype Frequencies in a Diploid Population
# Excoffier and Slatkin, 1995



test <- matrix(rnorm(1000000), ncol = 5)
test[sample(1000000, 500000, replace = F)] <- 0

# option 1
system.time(out <- qlcMatrix::rowMin(test, ignore.zero = TRUE))

# option 2
system.time(
  out <- matrixStats::rowSums2(test) - matrixStats::rowMaxs(test)
)

















# this just makes the imput data, don't worry about it. It's more for me for later.

# dat2 <- import.snpR.data(stickSNPs[,-c(1:3)],
#                          stickSNPs[,1:3],
#                          data.frame(pop = substr(colnames(stickSNPs)[-c(1:3)], 1, 3),
#                                     phenotype = rep(c("case", "control"), length = ncol(stickSNPs) - 3),
#                                     stringsAsFactors = F),
#                          "NN")
#
#
# dat2 <- add.facets.snpR.data(dat2, "phenotype")
#
# pheno.facet <- which(dat2@facet.meta == "phenotype")
#
# gs <- data.table::as.data.table(cbind(dat2@facet.meta[pheno.facet,], dat2@geno.tables$gs[pheno.facet,]))
# cast_gs <- data.table::dcast(data.table::setDT(gs), .snp.id ~ subfacet, value.var = colnames(dat2@geno.tables$gs))
# cast_gs <- cast_gs[,-1]
#
# saveRDS(as.matrix(cast_gs), "R_dev/armitage_example.RDS")
