ylistcpm <- ylist
ylistfpkm <- ylist

length(ylistcpm)
length(ylistfpkm)

summary(ylistcpm)
summary(ylistfpkm)

ylistcpm0 <- ylistcpm[ylistcpm != 0]
ylistfpkm0 <- ylistfpkm[ylistfpkm != 0]

length(ylistcpm0)
length(ylistfpkm0)

summary(ylistcpm0)
summary(ylistfpkm0)

head(y)
head(ynorm)
