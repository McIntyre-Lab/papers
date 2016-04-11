# Plot the quality score distribution of the variants in the filtered vcf file, 5.99mil variants


qual <- read.table('/home/fnew/dsim/vcf_filtered_qual_strings_sort.txt', header=TRUE)

quald <- density(qual$QUAL)


plot(qual$QUAL, type="b", xlim=c(0,6600))

hist(qual$QUAL, xlim=c(0,66000))

summary(qual$QUAL)
## Min.  1st Qu.   Median     Mean     3rd Qu.   Max.     ##
## 30.0  321.9     988.8      7195.0   7789.0    655500.0 ##

boxplot(qual$QUAL)


#het and qual data
hq <- read.table('/home/fnew/dsim/vcf/chr4_het_qual.txt', header=TRUE)
plot(y=hq$QUAL, x=hq$per_het, main="Percent Het vs Quality\nby position", xlab="Percent Heterozygosity", ylab="Variant Quality")


plot(y=hq$QUAL, x=hq$per_het, xlim=c(0,0.2), ylim=c(0,30000), main="Percent Het vs Quality\nby position\nZoomed in", xlab="Percent Heterozygosity", ylab="Variant Quality")

plot(y=hq$QUAL, x=hq$per_het, xlim=c(0,0.2), ylim=c(0,5000), main="Percent Het vs Quality\nby position\nZoomed in More", xlab="Percent Heterozygosity", ylab="Variant Quality")


