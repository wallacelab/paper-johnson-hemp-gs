library(ggplot2)
setwd('/home/jgwall/Projects/Hemp/CheckMattGenotyping/')

depths=read.delim('9_geno_depth.txt', check.names=F)  # VCFTools geno-depth report
depths = as.matrix(depths[, -c(1:2)])

total=rowSums(depths)
present=rowSums(depths>0)/ncol(depths)


qplot(x=log10(total), y=present, geom='bin2d')
ggsave("9_depth_plot.png")

# Restructure this to be able to play with it more
stats = data.frame(total_depth = total, 
                   mean_depth = total/ncol(depths),
                   sd_depth = apply(depths, FUN=sd, MARGIN=1),
                   fraction_present = present)
ggplot(stats, mapping=aes(x=log10(mean_depth), y=fraction_present)) +
    geom_bin2d() +
    scale_fill_gradient(trans='log')


substats = subset(stats, stats$mean_depth >= 1)
ggplot(substats, mapping=aes(x=log10(mean_depth), y=fraction_present)) +
    geom_bin2d()
ggsave("9_depth_plot.trimmed.png")
# The above indicates a concentration of things at about average depth 100 and 100% coverage.
# Are these paralogs or are these legitimate? I really need to know what my expected depth would be
# The only two concentrations I have are at ~1-2 reads per site, or 100 per individual. So either we
#    have way undersampling or the 100x is about correct.


# Looking at how hets are called

hetcalls = read.delim("tmp_hets.depth_only.txt", header=F)
names(hetcalls) = c("ref", "alt")
hetcalls$sum  =hetcalls$ref + hetcalls$alt
goodcalls = subset(hetcalls, hetcalls$sum > 50)
hist(goodcalls$alt / goodcalls$sum)  # Implies that at high depth, 50% is about average, but some spread
# Also, lowest have for any is 3 reads from either, which is about the filter I would put