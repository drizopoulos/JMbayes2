set1<-extractFrames_new(serBilir~ns(year,2)*drug+sex+(ns(year,2)|id),pbc2)
head(set1$X)
head(set1$Xhc)
head(set1$XhcC)
head(set1$XhcS)
head(set1$Z_)
head(set1$Xt)

set1<-extractFrames_new(serBilir~year+drug+sex+(year|id),pbc2)

set11<-extractFrames_new(serBilir~ns(year,2)+drug+sex+(ns(year,2)|id),pbc2)

head(set11$X)
head(set11$Xhc)
head(set11$XhcC)
head(set11$XhcS)
head(set11$Z_)
head(set11$means_X)
head(set11$means_Z)
head(set11$means_Xhc)
head(set11$SDs_X)
head(set11$SDs_Xhc)
head(set11$Xt)

set2<-extractFrames_new(serBilir~ns(year,2)*drug+sex+(1|id),pbc2)

head(set2$X)
head(set2$Xhc)
head(set2$XhcC)
head(set2$XhcS)
head(set2$Z_)
head(set2$means_X)
head(set2$means_Z)
head(set2$means_Xhc)
head(set2$SDs_X)
head(set2$SDs_Xhc)
head(set2$Xt)

set3<-extractFrames_new(serBilir~year+(1|id),pbc2)

set3a<-extractFrames_new(serBilir~sex+(1|id),pbc2)

set4<-extractFrames_new(serBilir~1+(1|id),pbc2)


set5<-extractFrames_new(serBilir~ns(year,2)*drug+sex+albumin+(1|id),pbc2)

set6<-extractFrames_new(serBilir~ns(year,2)*drug+sex+albumin+(albumin|id),pbc2)

set7<-extractFrames_new(serBilir~ns(year,2)*drug+sex+albumin+(ns(year,2)+albumin|id),pbc2)

set8<-extractFrames_new(serBilir~ns(year,2)*drug+sex+albumin+(ns(year,2)*drug|id),pbc2)

set9<-extractFrames_new(serBilir~year+drug+sex+albumin+(drug|id),pbc2)

