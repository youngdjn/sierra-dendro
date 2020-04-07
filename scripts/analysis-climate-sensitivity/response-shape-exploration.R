#### Exploring the response shapes of tree growth to ppt variation ####



# It's kind of weird that a lot of studies use a linear growth~precip slope as an index of growth sensitivity. It's simple and easy to intepret, but it's probably quite wrong. It seems biologically likely that there's an underlying non-linear response function relating tree growth to moisture availability. This is presumably saturating. It may even appear unimodal, but where it does, this is probably because other factors associated with very high precipitation depress growth: increased snowpack, which shortens effective growing season, or decreased solar radiation during early parts of the growing season. 

# In the EDA in this script, I tried to figure out how strong this nonlinearity is, and whether the nonlinearity appears within individual trees or only across trees or plots. The short answer: individual trees' growth responses are generally nonlinear an saturating, and might be represented OK by a piecewise linear regression, but not a linear regression.

# The observed saturating response shape presumably varies among genotypes and even among phenotypes for the same genotype (or even individual over time). The hypotheses this study seeks to address can all be framed in terms of this response shape and its relationship to patterns of moisture availability across the region. 

# 1) As a species nears its range edge, we expect an increase in sensitivity. This prediction relies on the assumption that there is a relatively consistent response function for each tree species and that this response function steepens toward lower levels of moisture availability. In drier regions, the trees express greater sensitivity because they are exposed to lower rainfall. If these assumptions hold roughly true, we would expect species approaching their range limit to exhibit a more rapid increase in sensitivity, because their (species-level) response function drops more sharply as moisture declines. 

# Note on time scales. It's important to keep in mind that the response function we're looking at is relatively slow -- it's linked to an annual growth cycle, and because of exogenous and endogenous lags, its time scale is effectively longer than 1 year. So we are not talking about instantaneous rates of photosynthesis, or even daily averages, but rather the cumulative effects of weather and other conditions over time spans of more than a year. Why does this longer time-scale matter? Because it means that the "performance" metric reflects photosynthesis during early season when moisture is available, later season when it is low or even u avaibable, and during dormant winter months when metabolic costs and snow or frost damage imposes costs, and further the influence of all these factors in at least the previous year. 

# 2) as a species nears its range edge, we expect some degree of local adaptation (and/or acclimation). This should counteract the pattern predicted in (1). Effectively, this means that the response function will be shifted to the left, or lower moisture levels, for populations near the range limit. Therefore populations at the range edge should display, empirically, steeper responses, but when normalized to absolute moisture levels, they should exhibit shallower responses than central populations at equivalent moisture levels. 

# Idea: model individual tree responses as piecewise linear functions with a breakpoint. Then the location of the breakpoint as well as the steepness of the "low-moisture" slope may vary among populations. Alternative could be a quadratic -- the advantage of that parameterization is one fewer parameters; the disadvantage is that the interpretation of the 2 parameters is more intertwined. 

# To explore this idea, let's visualize the curves for populations and individual trees and see what the empirical responses look like, across the range of the explanatory variable. This seems like it requires working with the absolute ppt values for each time step and population. 

d <- read.csv("./working-data/tree-ring-data-long.csv")

# load and filter data
d <- dplyr::select(d, cluster.x, cluster.y, plot.id, tree.id, species, year, rwi, raw_width, ba, bai, voronoi.area, radius.external, ppt.norm, tmean.norm, rad.tot, ppt, tmean, starts_with("ppt.z"), starts_with("tmean.z"))
d <- filter(d, species=="PSME")

ggplot(d, aes(ppt, raw_width, color=cluster.y)) + geom_point() + theme_bw()
ggplot(d, aes(ppt, raw_width, color=ppt.norm)) + geom_point() + scale_color_gradient2(midpoint=mean(d$ppt.norm, na.rm=T), low="blue", mid="white", high="red", space ="Lab" ) # + theme_bw()

ggplot(d[d$cluster.y=="Yose",], aes(ppt, rwi, color=plot.id)) + geom_point()  + theme_bw()

ggplot(d[d$plot.id=="HARDIN",], aes(ppt, rwi, color=tree.id)) + geom_point()  + theme_bw()
ggplot(d[d$plot.id=="HARDIN B",], aes(ppt, rwi, color=tree.id)) + geom_point()  + theme_bw() + geom_smooth(method=lm, formula = y~poly(x, 2), se=FALSE, fullrange=FALSE)


ggplot(d[d$plot.id=="HARDIN B",], aes(ppt, rwi, color=tree.id)) + geom_point(cex=0.1)  + theme_bw() + geom_smooth(span=1, se=FALSE, fullrange=FALSE)

# smoothed response by plot within cluster
# rwi
ggplot(d, aes(ppt, rwi, color=plot.id))  + theme_bw() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE) + facet_wrap(~cluster.y) +  theme(legend.position="none")
# raw widths
ggplot(d, aes(ppt, raw_width, color=plot.id))  + theme_bw() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE) + facet_wrap(~cluster.y) +  theme(legend.position="none")

# smoothed response by cluster
#rwi
ggplot(d, aes(ppt, rwi, color=cluster.y))  + theme_bw() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE)
# raw widths
ggplot(d, aes(ppt, raw_width, color=cluster.y))  + theme_bw() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE) 
# what's with the weird "kink" they all have around 800-1200 mm? Must reflect a heterogeneity among plots? What? 


# Plot responses of trees within plots within clusters

# SIERRA
# quadratic
ggplot(d[d$cluster.y=="Sierra",], aes(ppt, rwi, group=tree.id)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  +  theme(legend.position="none")+ geom_smooth(method=lm, formula = y~poly(x, 2), se=FALSE, fullrange=FALSE) 
# loess
ggplot(d[d$cluster.y=="Sierra",], aes(ppt.z, rwi, group=tree.id, color=tree.id)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  +  theme(legend.position="none")+ geom_smooth(method=loess, se=FALSE, fullrange=FALSE) #+ geom_point(cex=0.5) + scale_color_gradient2(midpoint=mean(d$year, na.rm=T), low="blue", mid="white", high="red", space ="Lab" )
#linear
ggplot(d[d$cluster.y=="Sierra",], aes(ppt.z, rwi, group=tree.id, color=tree.id)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  +  theme(legend.position="none")+ geom_smooth(method=lm, se=FALSE, fullrange=FALSE)
# Interesting consistency here is apparent breakpoint for a lot of trees and plots around 1000 mm precip or a bit less, with some trees showing a second acceleration at high end -- this appears to be real, not driven by a single high-leverage year. Further, if we color the points by year, this reveals excellent interspersion of ppt and decade; this is a longer-term pattern not driven by a single extreme or recent change. Quite interesting - what accounts for the kink in some plots but not others? Are there differences in snowpack? 

# Does the pattern look qualitatively similar if we look at raw growth rates? 
ggplot(d[d$cluster.y=="Sierra",], aes(ppt, raw_width, group=tree.id, color=year)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  +  theme(legend.position="none")+ geom_smooth(method=loess, se=FALSE, fullrange=FALSE) 
# Yes pretty much. 

# Yose
# quadratic
ggplot(d[d$cluster.y=="Yose",], aes(ppt, rwi, group=tree.id)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  +  theme(legend.position="none")+ geom_smooth(method=lm, formula = y~poly(x, 2), se=FALSE, fullrange=FALSE)
# loess
ggplot(d[d$cluster.y=="Yose",], aes(ppt, rwi, group=tree.id, color=tree.id)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  +  theme(legend.position="none")+ geom_smooth(method=loess, se=FALSE, fullrange=FALSE) #+ geom_point(cex=0.5) + scale_color_gradient2(midpoint=mean(d$year, na.rm=T), low="blue", mid="white", high="red", space ="Lab" )
# Differences from Sierra: switchpoint seems slightly higher, like 1200 mm for most trees and plots. Consistent with local adaptation hypothesis, but could also reflect covariance of ppt and temperature, and/or snowpack? 

# Tahoe
# quadratic
ggplot(d[d$cluster.y=="Tahoe",], aes(ppt, rwi, group=tree.id)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  +  theme(legend.position="none")+ geom_smooth(method=lm, formula = y~poly(x, 2), se=FALSE, fullrange=FALSE)
# loess
ggplot(d[d$cluster.y=="Tahoe",], aes(ppt, rwi, group=tree.id, color=tree.id)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  +  theme(legend.position="none")+ geom_smooth(method=loess, se=FALSE, fullrange=FALSE) #+ geom_point(cex=0.5) + scale_color_gradient2(midpoint=mean(d$year, na.rm=T), low="blue", mid="white", high="red", space ="Lab" )
# Quite heterogeneous among plots, but fairly consistent among trees within plots. I think this supports that we can legitimately model a plot-level response shape. Hard to make any generalizations otherwise in this cluster, though linear and saturating responses seem to predominate, with downward slopes at highest precip for many plots. 

# Plumas
# quadratic
ggplot(d[d$cluster.y=="Plumas",], aes(ppt, rwi, group=tree.id)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  +  theme(legend.position="none")+ geom_smooth(method=lm, formula = y~poly(x, 2), se=FALSE, fullrange=FALSE)
# loess
ggplot(d[d$cluster.y=="Plumas",], aes(ppt, rwi, group=tree.id, color=tree.id)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  +  theme(legend.position="none")+ geom_smooth(method=loess, se=FALSE, fullrange=FALSE) #+ geom_point(cex=0.5)+ scale_color_gradient2(midpoint=mean(d$year, na.rm=T), low="blue", mid="white", high="red", space ="Lab" )
# Most linear and least sensitive region. But note, below 1000 mm they show sharp decline, like trees in other regions. This drought sensitivity signal is driven by 3 years early on in the record, and results in an inital regression breakpoint or saturation around 1200-1400mm (hard to tell by eye). There is not a consistent decline in RWI above 3000 mm, maybe because these wet plots are warmer than those of similar precip in the Tahoe area? Probably a snowpack effect, or growing degree days, or both ?


# Now let's look at Sierra cluster with the "big 3" species. 

d_all <- read.csv("./working-data/tree-ring-data-long.csv")

# load and filter data
d_all <- dplyr::select(d_all, cluster.x, cluster.y, plot.id, tree.id, species, year, rwi, raw_width, ba, bai, voronoi.area, radius.external, ppt.norm, tmean.norm, rad.tot, ppt, tmean, starts_with("ppt.z"), starts_with("tmean.z"))

ggplot(d_all[d_all$cluster.y=="Sierra",], aes(ppt, rwi, color=species)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  + geom_smooth(method=loess, se=T, fullrange=FALSE) # +  theme(legend.position="none")
ggplot(d_all[d_all$cluster.y=="Yose",], aes(ppt, rwi, color=species)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  + geom_smooth(method=loess, se=T, fullrange=FALSE)  +  theme(legend.position="none")
ggplot(d_all[d_all$cluster.y=="Tahoe",], aes(ppt, rwi, color=species)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  + geom_smooth(method=loess, se=T, fullrange=FALSE)  +  theme(legend.position="none")
ggplot(d_all[d_all$cluster.y=="Plumas",], aes(ppt, rwi, color=species)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  + geom_smooth(method=loess, se=T, fullrange=FALSE)  +  theme(legend.position="none")
ggplot(d_all[d_all$cluster.y=="Plumas",], aes(ppt, raw_width, color=species)) + facet_wrap(~plot.id, nrow=3)   + theme_bw()  + geom_smooth(method=loess, se=T, fullrange=FALSE)  +  theme(legend.position="none")
ggplot(d_all, aes(ppt, raw_width, color=species))   + theme_bw()  + geom_smooth(method=loess, se=T, fullrange=FALSE)  + facet_wrap(~cluster.y)

ggplot(d_all, aes(ppt, rwi, color=species))   + theme_bw()  + geom_smooth(method=loess, se=T, fullrange=FALSE)  + facet_wrap(~cluster.y)

ggplot(d, aes(ppt, raw_width, color=cluster.y))    + theme_bw()  + geom_smooth(method=loess, se=T, fullrange=FALSE)  
ggplot(d, aes(ppt, rwi, color=cluster.y))    + theme_bw()  + geom_smooth(method=loess, se=T, fullrange=FALSE)  
ggplot(d, aes(ppt, rwi, color=cluster.x))    + theme_bw()  + geom_smooth(method=loess, se=T, fullrange=FALSE)  

# compare to ppt z scores (local anomaly)
ggplot(d, aes(ppt.z, rwi, color=cluster.x))    + theme_bw()  + geom_smooth(method=loess, se=T, fullrange=FALSE) 
ggplot(d, aes(ppt.z, rwi, color=cluster.x))    + theme_bw()  + geom_smooth(method=lm, formula = y~poly(x, 1), se=T, fullrange=FALSE) 
# Steeper slope relative to local variation in lower-elevation plots. NL and SL are very similar, nearly identical. NH is lease sensitive, and even shows decline at high z-scores; SH is intermediate. Could this suggest that sensitivity increases to the same degree at the range limit in the south, as at the elevational local range limits? And further that adapation is not enough to buffer the response to drier extremes? 

# Conclusion on linearity for the purposes of model structure. We should not model the response shape as linear. It's not, even at the level of individual trees, ALTHOUGH sometimes it is -- apparently in low-elevation, wet sites where precipitation is less limiting and there is also little snowpack. But linearity at the tree level is the exception. To generalize: for most trees, there is a breakpoint below which sensitivity increases. Below this point, generally we don't have enough data to determine if this decline is accelerating, and mostly it seems like it could be reasonable to represent it as linear. Above this breakpoint is a shallower slope that's also possibly linear. In many trees, especially in higher-rainfall areas, but really in all the clusters, there's a second breakpoint beyond which the response of growth to ppt either re-accelerates, or alternatively becomes negative. The negative trend is most prevalent in Tahoe, though you also see it in Plumas. My guess is that this negative trend at high precipitation will be primarily in higher-elevation or colder sites, and is driven by snowpack. 

# Alternative modeling methods: a) piecewise linear with a breakpoint; b) quadradic; c) basis spline. I don't know how to do any of these in a hierarchical way but the piecewise linear and quadratic should be not too complicated. 
