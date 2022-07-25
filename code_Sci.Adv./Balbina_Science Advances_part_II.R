# Palmeirim et al. 2022, Science Advances

# Summary of the script: 1 - Site-level modelling; 2 - Species-level modelling


# 1 - Site-level modelling

habitat<-read.table("habitat_variables.txt", header = T) 
area<-habitat[,1]
dist<-habitat[,2]
fire<-habitat[,7]
close.canopy<-habitat[,8]
lat<-habitat[,10]
long<-habitat[,11]

node_sites<-read.table("node_sites_metrics.txt", header = T) 

sites.sm<-node_sites[node_sites$Taxonomic_group == "sm",]
sites.liz<-node_sites[node_sites$Taxonomic_group == "liz",]
sites.mid.large<-node_sites[node_sites$Taxonomic_group == "mid.large",]
sites.birds<-node_sites[node_sites$Taxonomic_group == "birds",]
sites.frogs<-node_sites[node_sites$Taxonomic_group == "frogs",]
sites.beetles<-node_sites[node_sites$Taxonomic_group == "beetles",]
sites.bees<-node_sites[node_sites$Taxonomic_group == "bees",]
sites.trees<-node_sites[node_sites$Taxonomic_group == "trees",]
sites.all<-node_sites[node_sites$Taxonomic_group == "all_taxa",]

# 1.1  Normalised degree

norm.degree.sites.sm<-as.vector(sites.sm[,3])
norm.degree.sites.liz<-as.vector(sites.liz[,3])
norm.degree.sites.mid.large<-as.vector(sites.mid.large[,3])
norm.degree.sites.birds<-as.vector(sites.birds[,3])
norm.degree.sites.frogs<-as.vector(sites.frogs[,3])
norm.degree.sites.beetles<-as.vector(sites.beetles[,3])
norm.degree.sites.bees<-as.vector(sites.bees[,3])
norm.degree.sites.trees<-as.vector(sites.trees[,3])
norm.degree.sites.all<-as.vector(sites.all.taxa[,3])

# Small mammals
norm.degree.site.sm.data<-as.data.frame(cbind(norm.degree.sites.sm,area, dist, fire, close.canopy, lat, long))

m.degree.sm<-gls(norm.degree.sites.sm ~ log10(area) + dist + fire + close.canopy, data = norm.degree.site.sm.data, method = "REML") 
m.degree.sm.a<-update(m.degree.sm,correlation=corSpher(form=~lat + long))
#m.degree.sm.b<-update(m.degree.sm,correlation=corLin(form=~lat + long)) 
m.degree.sm.c<-update(m.degree.sm,correlation=corRatio(form=~lat + long))
m.degree.sm.d<-update(m.degree.sm,correlation=corGaus(form=~lat + long))
m.degree.sm.e<-update(m.degree.sm,correlation=corExp(form=~lat + long))
AIC(m.degree.sm,m.degree.sm.a,m.degree.sm.c,m.degree.sm.d,m.degree.sm.e) ### m.degree.sm best model

m.degree.sm2<-lm(norm.degree.sites.sm ~ scale(log10(area)) * scale(dist) + scale(fire)  + scale(close.canopy), data = norm.degree.site.sm.data)
summary(m.degree.sm2)
#plot(m.degree.sm2)
all.m.degree.sm2<-dredge(m.degree.sm2)
ave.all.m.degree.sm2<-model.avg(all.m.degree.sm2,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.degree.sm2) 
confint(ave.all.m.degree.sm2)
importance(ave.all.m.degree.sm2)

# Lizards
norm.degree.site.liz.data<-as.data.frame(cbind(norm.degree.sites.liz,area, dist, fire, close.canopy, lat, long))

m.degree.liz<-gls(norm.degree.sites.liz ~  log10(area) + dist +fire + close.canopy, data = norm.degree.site.liz.data, method = "REML")# a interacao faz com que o closed canopy seja importante!
m.degree.liz.a<-update(m.degree.liz,correlation=corSpher(form=~lat + long))
m.degree.liz.b<-update(m.degree.liz,correlation=corLin(form=~lat + long))
m.degree.liz.c<-update(m.degree.liz,correlation=corRatio(form=~lat + long))
m.degree.liz.d<-update(m.degree.liz,correlation=corGaus(form=~lat + long))
m.degree.liz.e<-update(m.degree.liz,correlation=corExp(form=~lat + long))
AIC(m.degree.liz,m.degree.liz.a,m.degree.liz.b,m.degree.liz.c,m.degree.liz.d,m.degree.liz.e) ### m.degree.liz best model

m.degree.liz2<-lm(norm.degree.sites.liz ~ scale(log10(area)) + scale(fire)  + scale(log10(area)) * scale(dist) +  scale(close.canopy), data = norm.degree.site.liz.data)# interaction is near significant
summary(m.degree.liz2)
#plot(m.degree.liz2)
all.m.degree.liz2<-dredge(m.degree.liz2)
ave.all.m.degree.liz2<-model.avg(all.m.degree.liz2,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.degree.liz2) 
confint(ave.all.m.degree.liz2)
importance(ave.all.m.degree.liz2)

# Mid-large mammals
norm.degree.site.mid.large.data <-as.data.frame(cbind(norm.degree.sites.mid.large,area, dist, fire, close.canopy, lat, long))
m.degree.mid.large<-gls(norm.degree.sites.mid.large ~ log10(area) + dist + fire + close.canopy, data = norm.degree.site.mid.large.data, method = "REML")# 
m.degree.mid.large.a<-update(m.degree.mid.large,correlation=corSpher(form=~lat + long))
#m.degree.mid.large.b<-update(m.degree.mid.large,correlation=corLin(form=~lat + long))
m.degree.mid.large.c<-update(m.degree.mid.large,correlation=corRatio(form=~lat + long))
m.degree.mid.large.d<-update(m.degree.mid.large,correlation=corGaus(form=~lat + long))
m.degree.mid.large.e<-update(m.degree.mid.large,correlation=corExp(form=~lat + long))
AIC(m.degree.mid.large,m.degree.mid.large.a,m.degree.mid.large.c,m.degree.mid.large.d,m.degree.mid.large.e) ### m.degree.mid.large.e best model -> if we do not account for the interaction between area and dist

m.degree.mid.large<-glm(norm.degree.sites.mid.large ~ scale(log10(area)) * scale(fire) + scale(dist) + scale(close.canopy), data = norm.degree.site.mid.large.data)# a interacao nao faz nada
summary(m.degree.mid.large)
#plot(m.degree.mid.large)
all.m.degree.mid.large<-dredge(m.degree.mid.large)
ave.all.m.degree.mid.large<-model.avg(all.m.degree.mid.large,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.degree.mid.large) 
confint(ave.all.m.degree.mid.large)
importance(ave.all.m.degree.mid.large)


# BIRDS
norm.degree.site.birds.data <-as.data.frame(cbind(norm.degree.sites.birds,area, dist, fire, close.canopy, lat, long))

m.degree.birds<-gls(norm.degree.sites.birds ~ log10(area) + dist + fire + close.canopy, data = norm.degree.site.birds.data, method = "REML")# 
m.degree.birds.a<-update(m.degree.birds,correlation=corSpher(form=~lat + long))
m.degree.birds.b<-update(m.degree.birds,correlation=corLin(form=~lat + long))
m.degree.birds.c<-update(m.degree.birds,correlation=corRatio(form=~lat + long))
m.degree.birds.d<-update(m.degree.birds,correlation=corGaus(form=~lat + long))
m.degree.birds.e<-update(m.degree.birds,correlation=corExp(form=~lat + long))
AIC(m.degree.birds,m.degree.birds.a,m.degree.birds.b, m.degree.birds.c,m.degree.birds.d,m.degree.birds.e) ### 

m.degree.birds2<-glm(norm.degree.sites.birds ~ scale(log10(area)) + scale(dist) + scale(fire) + scale(close.canopy),family = gaussian(link = "log"), data = norm.degree.site.birds.data) # a interacao nao faz nada
summary(m.degree.birds2)
#plot(m.degree.birds2)
all.m.degree.birds2<-dredge(m.degree.birds2)
ave.all.m.degree.birds2<-model.avg(all.m.degree.birds2,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.degree.birds2) 
confint(ave.all.m.degree.birds2)
importance(ave.all.m.degree.birds2)


# FROGS
norm.degree.site.frogs.data <-as.data.frame(cbind(norm.degree.sites.frogs,area, dist, fire, close.canopy, lat, long))

m.degree.frogs<-gls(norm.degree.sites.frogs ~ log10(area) + dist +  fire + close.canopy, data = norm.degree.site.frogs.data, method = "REML")# 
m.degree.frogs.a<-update(m.degree.frogs,correlation=corSpher(form=~lat + long))
#m.degree.frogs.b<-update(m.degree.frogs,correlation=corLin(form=~lat + long)) 
m.degree.frogs.c<-update(m.degree.frogs,correlation=corRatio(form=~lat + long))
m.degree.frogs.d<-update(m.degree.frogs,correlation=corGaus(form=~lat + long))
m.degree.frogs.e<-update(m.degree.frogs,correlation=corExp(form=~lat + long))
AIC(m.degree.frogs,m.degree.frogs.a,m.degree.frogs.c,m.degree.frogs.d,m.degree.frogs.e) ### 

m.degree.frogs<-lm(norm.degree.sites.frogs ~scale(log10(area)) + scale(dist) + scale(fire) + scale(close.canopy), data = norm.degree.site.frogs.data)# a interacao nao 'e significativa
summary(m.degree.frogs)
#plot(m.degree.frogs)
all.m.degree.frogs<-dredge(m.degree.frogs)
ave.all.m.degree.frogs<-model.avg(all.m.degree.frogs,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.degree.frogs) 
confint(ave.all.m.degree.frogs)
importance(ave.all.m.degree.frogs)

# Dung beetles
norm.degree.site.beetles.data <-as.data.frame(cbind(norm.degree.sites.beetles,area, dist, fire, close.canopy, lat, long))

m.degree.beetles<-gls(norm.degree.sites.beetles ~ log10(area) + dist + fire + close.canopy, data = norm.degree.site.beetles.data, method = "REML")# 
m.degree.beetles.a<-update(m.degree.beetles,correlation=corSpher(form=~lat + long))
m.degree.beetles.b<-update(m.degree.beetles,correlation=corLin(form=~lat + long)) 
m.degree.beetles.c<-update(m.degree.beetles,correlation=corRatio(form=~lat + long))
m.degree.beetles.d<-update(m.degree.beetles,correlation=corGaus(form=~lat + long))
m.degree.beetles.e<-update(m.degree.beetles,correlation=corExp(form=~lat + long))
AIC(m.degree.beetles,m.degree.beetles.a,m.degree.beetles.b,m.degree.beetles.c,m.degree.beetles.d,m.degree.beetles.e) ###

m.degree.beetles<-glm(log10(norm.degree.sites.beetles+1.1) ~ scale(log10(area)) * scale(log10(dist+1)) + scale(fire) + scale(close.canopy), family = gaussian(link = "log"), data = norm.degree.site.beetles.data) 
summary(m.degree.beetles)
#plot(m.degree.beetles2)
all.m.degree.beetles<-dredge(m.degree.beetles)
ave.all.m.degree.beetles<-model.avg(all.m.degree.beetles,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.degree.beetles) 
confint(ave.all.m.degree.beetles)
importance(ave.all.m.degree.beetles)

# Orchid bees
norm.degree.site.bees.data <-as.data.frame(cbind(norm.degree.sites.bees,area, dist, fire, close.canopy, lat, long))

m.degree.bees<-gls(norm.degree.sites.bees ~ log10(area) + dist + fire + close.canopy, data = norm.degree.site.bees.data, method = "REML")# 
m.degree.bees.a<-update(m.degree.bees,correlation=corSpher(form=~lat + long))
m.degree.bees.b<-update(m.degree.bees,correlation=corLin(form=~lat + long)) 
m.degree.bees.c<-update(m.degree.bees,correlation=corRatio(form=~lat + long))
m.degree.bees.d<-update(m.degree.bees,correlation=corGaus(form=~lat + long))
m.degree.bees.e<-update(m.degree.bees,correlation=corExp(form=~lat + long))
AIC(m.degree.bees,m.degree.bees.a,m.degree.bees.b,m.degree.bees.c,m.degree.bees.d,m.degree.bees.e) ### m.degree.beetles.b c and d are the best models
anova(m.degree.bees,m.degree.bees.c) # model c is better!

m1<-gls(norm.degree.sites.bees ~ scale(log10(area)) * scale(dist) + scale(fire)  + scale(close.canopy), correlation=corRatio(form=~lat + long), data= norm.degree.site.bees.data, method="ML")
summary(m1)
dd<-dredge(m1)
ave.all.m.degree.bees.c<-model.avg(dd,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.degree.bees.c, method = "REML") 
confint(ave.all.m.degree.bees.c)

# Trees
norm.degree.site.trees.data <-as.data.frame(cbind(norm.degree.sites.trees,area, dist, fire, close.canopy, lat, long))

m.degree.trees<-gls(norm.degree.sites.trees ~ log10(area) + dist + fire + close.canopy, data = norm.degree.site.trees.data, method = "REML")# 

m.degree.trees.a<-update(m.degree.trees,correlation=corSpher(form=~lat + long))
m.degree.trees.b<-update(m.degree.trees,correlation=corLin(form=~lat + long)) 
m.degree.trees.c<-update(m.degree.trees,correlation=corRatio(form=~lat + long))
m.degree.trees.d<-update(m.degree.trees,correlation=corGaus(form=~lat + long))
m.degree.trees.e<-update(m.degree.trees,correlation=corExp(form=~lat + long))
AIC(m.degree.trees,m.degree.trees.a,m.degree.trees.b,m.degree.trees.c,m.degree.trees.d,m.degree.trees.e) 
anova(m.degree.trees,m.degree.trees.a) # models are similar so I will keep the simplest (m.degree.trees)

m.degree.trees2<-lm(norm.degree.sites.trees ~  scale(log10(area)) *scale(fire) + scale(dist) + scale(close.canopy), data = norm.degree.site.trees.data)# a interacao nao muda nada
summary(m.degree.trees2)
#plot(m.degree.trees2)
all.m.degree.trees2<-dredge(m.degree.trees2)
ave.all.m.degree.trees2<-model.avg(all.m.degree.trees2,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.degree.trees2) 
confint(ave.all.m.degree.trees2)
importance(ave.all.m.degree.trees2)


# ALL TAXA
norm.degree.site.all.data <-as.data.frame(cbind(norm.degree.sites.all,area, dist, fire, close.canopy, lat, long))

m.degree.all<-gls(norm.degree.sites.all ~ log10(area) + dist + fire + close.canopy, data = norm.degree.site.all.data, method = "REML")# 
m.degree.all.a<-update(m.degree.all,correlation=corSpher(form=~lat + long))
m.degree.all.b<-update(m.degree.all,correlation=corLin(form=~lat + long)) 
m.degree.all.c<-update(m.degree.all,correlation=corRatio(form=~lat + long))
m.degree.all.d<-update(m.degree.all,correlation=corGaus(form=~lat + long))
m.degree.all.e<-update(m.degree.all,correlation=corExp(form=~lat + long))

AIC(m.degree.all,m.degree.all.a,m.degree.all.b,m.degree.all.c,m.degree.all.d,m.degree.all.e) 
anova(m.degree.all,m.degree.all.c)

m.degree.all2<-lm(norm.degree.sites.all ~  scale(log10(area)) * scale(fire) + scale(dist) + scale(close.canopy), data = norm.degree.site.all.data)# 
summary(m.degree.all2)
#plot(m.degree.all2)
all.m.degree.all2<-dredge(m.degree.all2)
ave.all.m.degree.all2<-model.avg(all.m.degree.all2,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.degree.all2) 
confint(ave.all.m.degree.all2)
importance(ave.all.m.degree.all2)



# 1.2  Nestedness contribution

nest.cont.sites.sm<-as.vector(sites.sm[,4])
nest.cont.sites.liz<-as.vector(sites.liz[,4])
nest.cont.sites.mid.large<-as.vector(sites.mid.large[,4])
nest.cont.sites.birds<-as.vector(sites.birds[,4])
nest.cont.sites.frogs<-as.vector(sites.frogs[,4])
nest.cont.sites.beetles<-as.vector(sites.beetles[,4])
nest.cont.sites.bees<-as.vector(sites.bees[,4])
nest.cont.sites.trees<-as.vector(sites.trees[,4])
nest.cont.sites.all<-as.vector(sites.all[,4])

# Small mammals
nest.sites.sm.data<-as.data.frame(cbind(nest.cont.sites.sm,area, dist, fire, close.canopy, lat, long))
m.sm<-gls(nest.cont.sites.sm ~ log10(area) + dist + fire + close.canopy, data = nest.sites.sm.data, method = "REML")# 
m.sm.a<-update(m.sm,correlation=corSpher(form=~lat + long))
#m.sm.b<-update(m.sm,correlation=corLin(form=~lat + long))
m.sm.c<-update(m.sm,correlation=corRatio(form=~lat + long))
m.sm.d<-update(m.sm,correlation=corGaus(form=~lat + long))
m.sm.e<-update(m.sm,correlation=corExp(form=~lat + long))
AIC(m.sm,m.sm.a,m.sm.c,m.sm.d,m.sm.e) ### m.sm best model

m.sm2<-lm(nest.cont.sites.sm ~ scale(log10(area)) +  scale(fire) + scale(dist) + scale(close.canopy), data = nest.sites.sm.data)  
summary(m.sm2)
#plot(m.sm)
all.m.sm2<-dredge(m.sm2)
ave.all.m.sm2<-model.avg(all.m.sm2,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.sm2) 
confint(ave.all.m.sm2)
importance(ave.all.m.sm2)

# Lizards
nest.sites.liz.data<-as.data.frame(cbind(nest.cont.sites.liz, area, dist, fire, close.canopy, lat, long))
m.liz<-gls(nest.cont.sites.liz~ log10(area) + dist + fire + close.canopy, data = nest.sites.liz.data, method = "REML")# 
m.liz.a<-update(m.liz,correlation=corSpher(form=~lat + long))
m.liz.b<-update(m.liz,correlation=corLin(form=~lat + long))
m.liz.c<-update(m.liz,correlation=corRatio(form=~lat + long))
m.liz.d<-update(m.liz,correlation=corGaus(form=~lat + long))
m.liz.e<-update(m.liz,correlation=corExp(form=~lat + long))
AIC(m.liz,m.liz.a,m.liz.c,m.liz.d,m.liz.e) ### m.liz best model

m.liz<-lm(nest.cont.sites.liz ~ scale(log10(area)) + scale(fire) + scale(dist) + scale(close.canopy), data = nest.sites.liz.data)# 
summary(m.liz)
#plot(m.liz)
all.m.liz<-dredge(m.liz)
ave.all.m.liz<-model.avg(all.m.liz,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.liz) 
confint(ave.all.m.liz)
importance(ave.all.m.liz)

# MID.LARGE MAMMALS
nest.sites.mid.large.data<-as.data.frame(cbind(nest.cont.sites.mid.large,area, dist, fire, close.canopy, lat, long))
m.mid.large<-gls(nest.cont.sites.mid.large ~ log10(area) + dist + fire + close.canopy, data = nest.sites.mid.large.data, method = "REML")# 
m.mid.large.a<-update(m.mid.large,correlation=corSpher(form=~lat + long))
m.mid.large.b<-update(m.mid.large,correlation=corLin(form=~lat + long)) # false convergence
m.mid.large.c<-update(m.mid.large,correlation=corRatio(form=~lat + long))
m.mid.large.d<-update(m.mid.large,correlation=corGaus(form=~lat + long))
m.mid.large.e<-update(m.mid.large,correlation=corExp(form=~lat + long))
AIC(m.mid.large,m.mid.large.a,m.mid.large.b,m.mid.large.c,m.mid.large.d,m.mid.large.e) ###

m.mid.large<-lm(nest.cont.sites.mid.large ~ scale(log10(area)) * scale(fire) + scale(dist) + scale(close.canopy), data = nest.sites.mid.large.data)# 
summary(m.mid.large)
#plot(m.mid.large)
all.m.mid.large<-dredge(m.mid.large)
ave.all.m.mid.large<-model.avg(all.m.mid.large,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.mid.large) 
confint(ave.all.m.mid.large)
importance(ave.all.m.mid.large)

# BIRDS 
nest.sites.birds.data<-as.data.frame(cbind(nest.cont.sites.birds,area, dist, fire, close.canopy, lat, long))
m.birds<-gls(log10(nest.cont.sites.birds) ~ log10(area) + dist + fire + close.canopy, data = nest.sites.birds.data, method = "REML")# 
m.birds.a<-update(m.birds,correlation=corSpher(form=~lat + long))
m.birds.b<-update(m.birds,correlation=corLin(form=~lat + long))
m.birds.c<-update(m.birds,correlation=corRatio(form=~lat + long))
m.birds.d<-update(m.birds,correlation=corGaus(form=~lat + long))
m.birds.e<-update(m.birds,correlation=corExp(form=~lat + long))
AIC(m.birds,m.birds.a,m.birds.b,m.birds.c,m.birds.d,m.birds.e) ### m.birds best model
anova(m.birds,m.birds.b) # MODEL B is better than without geographic coordinates!

m1<-gls(log10(nest.cont.sites.birds) ~ scale(log10(area)) + scale(fire) + scale(dist) + scale(close.canopy), correlation=corLin(form=~lat + long), data= nest.sites.birds.data, method="ML")
summary(m1)
dd<-dredge(m1)
ave.all.m.degree.bees.c<-model.avg(dd,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.degree.bees.c, method = "REML") 
confint(ave.all.m.degree.bees.c)

# FROGS
nest.sites.frogs.data<-as.data.frame(cbind(nest.cont.sites.frogs,area, dist, fire, close.canopy, lat, long))
m.frogs<-gls(nest.cont.sites.frogs ~ log10(area) + dist + fire + close.canopy, data = nest.sites.frogs.data, method = "REML")# 
m.frogs.a<-update(m.frogs,correlation=corSpher(form=~lat + long))
#m.frogs.b<-update(m.frogs,correlation=corLin(form=~lat + long))
m.frogs.c<-update(m.frogs,correlation=corRatio(form=~lat + long))
m.frogs.d<-update(m.frogs,correlation=corGaus(form=~lat + long))
m.frogs.e<-update(m.frogs,correlation=corExp(form=~lat + long))
AIC(m.frogs,m.frogs.a,m.frogs.c,m.frogs.d,m.frogs.e) 

m.frogs<-lm(nest.cont.sites.frogs ~ scale(log10(area)) + scale(fire) + scale(dist) + scale(close.canopy), data = nest.sites.frogs.data)# 
summary(m.frogs)
#plot(m.frogs)
all.m.frogs<-dredge(m.frogs)
ave.all.m.frogs<-model.avg(all.m.frogs,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.frogs) 
confint(ave.all.m.frogs)
#importance(ave.all.m.frogs)

# Dung-beetles
nest.sites.beetles.data<-as.data.frame(cbind(nest.cont.sites.beetles,area, dist, fire, close.canopy, lat, long))
m.beetles<-gls(nest.cont.sites.beetles ~ log10(area) + dist + fire + close.canopy, data = nest.sites.beetles.data, method = "REML")# 
m.beetles.a<-update(m.beetles,correlation=corSpher(form=~lat + long))
m.beetles.b<-update(m.beetles,correlation=corLin(form=~lat + long)) 
m.beetles.c<-update(m.beetles,correlation=corRatio(form=~lat + long))
m.beetles.d<-update(m.beetles,correlation=corGaus(form=~lat + long))
m.beetles.e<-update(m.beetles,correlation=corExp(form=~lat + long))
AIC(m.beetles,m.beetles.a,m.beetles.b,m.beetles.c,m.beetles.d,m.beetles.e) ### m.beetles best model

m.beetles<-lm(nest.cont.sites.beetles ~ scale(log10(area)) * scale(fire) + scale(dist) + scale(close.canopy), data = nest.sites.beetles.data)# 
summary(m.beetles)
#plot(m.beetles)
all.m.beetles<-dredge(m.beetles)
ave.all.m.beetles<-model.avg(all.m.beetles,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.beetles) 
confint(ave.all.m.beetles)
importance(ave.all.m.beetles)

# BEES
nest.sites.bees.data<-as.data.frame(cbind(nest.cont.sites.bees,area, dist, fire, close.canopy, lat, long))
m.bees<-gls(nest.cont.sites.bees ~ log10(area) + dist + fire + close.canopy, data = nest.sites.bees.data, method = "REML")# 
m.bees.a<-update(m.bees,correlation=corSpher(form=~lat + long))
#m.bees.b<-update(m.bees,correlation=corLin(form=~lat + long)) # false convergence
m.bees.c<-update(m.bees,correlation=corRatio(form=~lat + long))
m.bees.d<-update(m.bees,correlation=corGaus(form=~lat + long))
m.bees.e<-update(m.bees,correlation=corExp(form=~lat + long))
AIC(m.bees,m.bees.a,m.bees.c,m.bees.d,m.bees.e) ### m.bees best model

m.bees<-glm(nest.cont.sites.bees ~ scale(log10(area)) + scale(fire) + scale(dist) + scale(close.canopy), data = nest.sites.bees.data)# 
summary(m.bees)
#plot(m.bees)
all.m.bees<-dredge(m.bees)
ave.all.m.bees<-model.avg(all.m.bees,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.bees) 
confint(ave.all.m.bees)
importance(ave.all.m.bees)

# Trees
nest.sites.trees.data<-as.data.frame(cbind(nest.cont.sites.trees,area, dist, fire, close.canopy, lat, long))
m.trees<-gls(nest.cont.sites.trees ~ log10(area) + dist + fire + close.canopy, data = nest.sites.trees.data, method = "REML")# 
m.trees.a<-update(m.trees,correlation=corSpher(form=~lat + long))
#m.trees.b<-update(m.trees,correlation=corLin(form=~lat + long)) 
m.trees.c<-update(m.trees,correlation=corRatio(form=~lat + long))
m.trees.d<-update(m.trees,correlation=corGaus(form=~lat + long))
m.trees.e<-update(m.trees,correlation=corExp(form=~lat + long))
AIC(m.trees,m.trees.a,m.trees.c,m.trees.d,m.trees.e) ### m.trees best model

m.trees<-lm(nest.cont.sites.trees ~ scale(log10(area)) * scale(fire) + scale(dist) + scale(close.canopy), data = nest.sites.trees.data)# 
summary(m.trees)
#plot(m.trees)
all.m.trees<-dredge(m.trees)
ave.all.m.trees<-model.avg(all.m.trees,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.trees) 
confint(ave.all.m.trees)
importance(ave.all.m.trees)

# All taxa
nest.sites.all.data<-as.data.frame(cbind(nest.cont.sites.all, area, dist, fire, close.canopy, lat, long))
m.all<-gls(log10(nest.cont.sites.all) ~log10(area) + dist + fire + close.canopy, data = nest.sites.all.data, method = "REML")# 
m.all.a<-update(m.all,correlation=corSpher(form=~lat + long))
#m.all.b<-update(m.all,correlation=corLin(form=~lat + long))
m.all.c<-update(m.all,correlation=corRatio(form=~lat + long))
m.all.d<-update(m.all,correlation=corGaus(form=~lat + long))
m.all.e<-update(m.all,correlation=corExp(form=~lat + long))
AIC(m.all,m.all.a,m.all.c,m.all.d,m.all.e) ### 

m.all<-glm(log10(nest.cont.sites.all) ~ scale(log10(area)) + scale(fire) + scale(dist) + scale(close.canopy), family = gaussian(link = "log"), data = nest.sites.all.data) # 
summary(m.all)
#plot(m.all)
all.m.all<-dredge(m.all)
ave.all.m.all<-model.avg(all.m.all,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.all) 
confint(ave.all.m.all)
importance(ave.all.m.all)

# 1.3  Figure 5

# 1 - Species-level modelling
node_species<-read.table("node_species_metrics.txt", header = T) 

# 1.1  Normalised degree

# Small mammals
sm.traits<-read.table("sm.traits.txt", header = T) 
body.size<-sm.traits[,1]
diet<-sm.traits[,2]
locomotion<-sm.traits[,3]
matrix.tol<-sm.traits[,4]

species.sm<-node_species[node_species$taxa == "sm",]
degree.species.sm<-as.vector(species.sm[,5])
degree.species.sm.data<-as.data.frame(cbind(degree.species.sm, body.size, diet, locomotion, matrix.tol))

m.traits.sm<-lm(degree.species.sm ~ scale(log10(body.size)) + scale(diet) + scale(locomotion) + scale(log10(matrix.tol+1)), data = degree.species.sm.data)# 
summary(m.traits.sm)
#plot(m.traits.sm)
all.m.traits.sm<-dredge(m.traits.sm)
ave.all.m.traits.sm<-model.avg(all.m.traits.sm,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.sm) 
confint(ave.all.m.traits.sm)
importance(ave.all.m.traits.sm)


# Lizards
liz.traits<-read.table("liz.traits.txt", header = T) 
cor(liz.traits) # habitat and thermo are highly correlated -> I did not keep habitat
thermo<-liz.traits[,1]
habitat<-liz.traits[,2]
svl<-liz.traits[,3]
prey<-liz.traits[,4]
species.liz<-node_species[node_species$taxa == "liz",]
degree.species.liz<-as.vector(species.liz[,5])
degree.species.liz.data<-as.data.frame(cbind(degree.species.liz, thermo, svl, prey))

m.traits.liz<-lm(degree.species.liz ~ scale(thermo) + scale(log10(svl)) + scale(prey), data = degree.species.liz.data)# 
summary(m.traits.liz)
#plot(m.traits.liz)
all.m.traits.liz<-dredge(m.traits.liz)
ave.all.m.traits.liz<-model.avg(all.m.traits.liz,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.liz) 
confint(ave.all.m.traits.liz)
importance(ave.all.m.traits.liz)


# Mid-large mammals
mid.large.traits<-read.table("mid.large.traits.txt", header = T) 
body_mass<-mid.large.traits[,1]
group_size<-mid.large.traits[,2]
home_range<-mid.large.traits[,3]
diet<-mid.large.traits[,4]
species.mid.large<-node_species[node_species$taxa == "mid.large",]
degree.species.mid.large<-as.vector(species.mid.large[,5])
degree.species.mid.large.data<-as.data.frame(cbind(degree.species.mid.large,body_mass, group_size, home_range, diet))

m.traits.degree.mid.large<-lm(degree.species.mid.large ~ scale(log10(body_mass)) + scale(log10(group_size)) + scale(log10(home_range)) + scale(diet), data = degree.species.mid.large.data)# 
summary(m.traits.degree.mid.large)
#plot(m.traits.degree.mid.large)
all.m.traits.degree.mid.large<-dredge(m.traits.degree.mid.large)
ave.all.m.traits.degree.mid.large<-model.avg(all.m.traits.degree.mid.large,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.degree.mid.large) 
confint(ave.all.m.traits.degree.mid.large)
importance(ave.all.m.traits.degree.mid.large)

# Understorey birds
setwd("~/R/redes/data")
bird.traits<-read.table("birds.traits.txt", header = T) 
forest.dep<-bird.traits[,1]; hist(forest.dep)
trophic.level<-bird.traits[,2]; hist(trophic.level)
body.mass<-bird.traits[,3]; hist(log10(body.mass)) # usar log
species.birds<-node_species[node_species$taxa == "birds",]
degree.species.birds<-as.vector(species.birds[,5])
degree.species.birds.data<-as.data.frame(cbind(degree.species.birds,forest.dep, trophic.level,body.mass))

m.traits.degree.birds<-lm(log10(degree.species.birds) ~ scale(forest.dep) + scale(trophic.level) + scale(log10(body.mass)), data = degree.species.birds.data)# 
summary(m.traits.degree.birds)
#plot(m.traits.degree.birds)
all.m.traits.degree.birds<-dredge(m.traits.degree.birds)
ave.all.m.traits.degree.birds<-model.avg(all.m.traits.degree.birds,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.degree.birds) 
confint(ave.all.m.traits.degree.birds)
importance(ave.all.m.traits.degree.birds)

# Frogs
frog.traits<-read.table("node_species_metrics_frogs.txt", header = T)
rep.mode<-frog.traits[,8]
habitat.div<-frog.traits[,9]
body.size<-as.numeric(frog.traits[,10])
degree.species.frogs<-as.vector(frog.traits[,3])
degree.species.frogs.data<-as.data.frame(cbind(degree.species.frogs,forest.dep, trophic.level,body.mass))

m.traits.degree.frogs<-glm(degree.species.frogs ~ scale(rep.mode) + scale(habitat.div) + scale(log10(body.size)))# 
summary(m.traits.degree.frogs)
#plot(m.traits.degree.frogs)
all.m.traits.degree.frogs<-dredge(m.traits.degree.frogs)
ave.all.m.traits.degree.frogs<-model.avg(all.m.traits.degree.frogs,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.degree.frogs) 
confint(ave.all.m.traits.degree.frogs)
importance(ave.all.m.traits.degree.frogs)

# BEETLES
beetles.traits<-read.table("dung-beetles.traits.txt", header = T) 
body.size<-beetles.traits[,1]
relocation<-beetles.traits[,2]
diet<-beetles.traits[,3]
species.beetles<-node_species[node_species$taxa == "beetles",]
degree.species.beetles<-as.vector(species.beetles[,5])
degree.species.beetles.data0<-as.data.frame(cbind(degree.species.beetles,body.size,relocation,diet))
degree.species.beetles.data<-degree.species.beetles.data0[-c(32),]
beetles<-as.numeric(log10(degree.species.beetles[-c(32)]))
body.size2<-log10(body.size[-c(32)])
relocation2<-relocation[-c(32)]
diet2<-diet[-c(32)]

m.traits.degree.beetles<-lm(beetles ~ scale(body.size2) + relocation2 + diet2, data = degree.species.beetles.data) 
summary(m.traits.degree.beetles)
#plot(m.traits.degree.beetles)
all.m.traits.degree.beetles<-dredge(m.traits.degree.beetles)
ave.all.m.traits.degree.beetles<-model.avg(all.m.traits.degree.beetles,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.degree.beetles) 
confint(ave.all.m.traits.degree.beetles)
importance(ave.all.m.traits.degree.beetles)

# BEES
bees.traits<-read.table("bees.traits.txt", header = T) 
cor(bees.traits) # all the measures are highly correlated
body.length<-bees.traits[,1]
species.bees<-node_species[node_species$taxa == "bees",]
degree.species.bees<-as.vector(species.bees[,5]) 
degree.species.bees.data<-as.data.frame(cbind(degree.species.bees,body.length))

m.traits.degree.bees<-lm(log10(degree.species.bees) ~ body.length, data = degree.species.bees.data)# 
summary(m.traits.degree.bees)
plot(m.traits.degree.bees)
confint(m.traits.degree.bees)


# TREES
trees.traits<-read.table("trees.traits.txt", header = T) 
wood_density<-trees.traits[,1]
seed_mass<-trees.traits[,2]
vertical_strat<-trees.traits[,3]
dispersal_mode<-trees.traits[,4]
regeneration<-trees.traits[,5]

trees.traits2<-read.table("trees.traits_numeric.txt", header = T) 
cor(trees.traits2[,1:5])

species.trees<-node_species[node_species$taxa == "trees",]
#degree.species.trees<-as.vector(species.trees[,4])
degree.species.trees<-as.numeric(trees.traits[,7])
degree.species.trees.data<-as.data.frame(cbind(degree.species.trees,wood_density,seed_mass,vertical_strat,dispersal_mode,regeneration))

wood_density2<-scale(wood_density)
seed_mass2<-scale(seed_mass)
degree.species.trees2<-log10(degree.species.trees)

m.traits.degree.trees<-lm(degree.species.trees2 ~ wood_density2 + seed_mass2 + vertical_strat + dispersal_mode + regeneration, data = degree.species.trees.data)# 
summary(m.traits.degree.trees)
#plot(m.traits.degree.trees)
all.m.traits.degree.trees<-dredge(m.traits.degree.trees)
ave.all.m.traits.degree.trees<-model.avg(all.m.traits.degree.trees,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.degree.trees) 
confint(ave.all.m.traits.degree.trees)
importance(ave.all.m.traits.degree.trees)

# 1.2  Nestedness contribution

nest.cont.species.sm<-as.vector(species.sm[,8])
nest.cont.species.liz<-as.vector(species.liz[,8])
nest.cont.species.mid.large<-as.vector(species.mid.large[,8])
nest.cont.species.birds<-as.vector(species.birds[,8])
nest.cont.species.frogs<-as.vector(species.frogs[,8])
nest.cont.species.beetles<-as.vector(species.beetles[,8])
nest.cont.species.bees<-as.vector(species.bees[,8])
nest.cont.species.trees<-as.vector(species.trees[,8])

# Small mammals
sm.traits<-read.table("sm.traits.txt", header = T) 
body.size<-sm.traits[,1]
diet<-sm.traits[,2]
locomotion<-sm.traits[,3]
matrix.tol<-sm.traits[,4]
nest.species.sm.data<-as.data.frame(cbind(nest.cont.species.sm, log10(body.size), diet, locomotion, log10(matrix.tol+1)))

m.traits.sm<-lm(nest.cont.species.sm ~ scale(log10(body.size)) + scale(diet) + scale(locomotion) + scale(log10(matrix.tol+1)), data = nest.species.sm.data)# 
summary(m.traits.sm)
#plot(m.traits.sm)
all.m.traits.sm<-dredge(m.traits.sm)
ave.all.m.traits.sm<-model.avg(all.m.traits.sm,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.sm) 
confint(ave.all.m.traits.sm)
importance(ave.all.m.traits.sm)

# Lizards
liz.traits<-read.table("liz.traits.txt", header = T) 
cor(liz.traits) # habitat and thermo are highly correlated -> I did not keep habitat
thermo<-liz.traits[,1]
habitat<-liz.traits[,2]
svl<-liz.traits[,3]
prey<-liz.traits[,4]
nest.species.liz.data<-as.data.frame(cbind(nest.cont.species.liz,thermo, log10(svl), prey))

m.traits.liz<-lm(log10(nest.cont.species.liz) ~ scale(thermo) + scale(log10(svl)) + scale(prey), data = nest.species.sm.data)# 
summary(m.traits.liz)
#plot(m.traits.liz)
all.m.traits.liz<-dredge(m.traits.liz)
ave.all.m.traits.liz<-model.avg(all.m.traits.liz,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.liz) 
confint(ave.all.m.traits.liz)
importance(ave.all.m.traits.liz)

# Mid-large mammals
mid.large.traits<-read.table("mid.large.traits.txt", header = T) 
body_mass<-mid.large.traits[,1]
group_size<-mid.large.traits[,2]
home_range<-mid.large.traits[,3]
diet<-mid.large.traits[,4]
nest.species.mid.large.data<-as.data.frame(cbind(nest.cont.species.mid.large,body_mass, group_size, home_range, diet))

m.traits.mid.large<-glm(nest.cont.species.mid.large+1 ~ scale(log10(body_mass)) + scale(log10(group_size)) + scale(log10(home_range)) + scale(diet), family=gaussian(link = "log"), data = nest.species.mid.large.data)# 
summary(m.traits.mid.large)
#plot(m.traits.mid.large) # 
all.m.traits.mid.large<-dredge(m.traits.mid.large)
ave.all.m.traits.mid.large<-model.avg(all.m.traits.mid.large,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.mid.large) 
confint(ave.all.m.traits.mid.large)
importance(ave.all.m.traits.mid.large)

# Birds
bird.traits<-read.table("birds.traits.txt", header = T) 
forest.dep<-bird.traits[,1]
trophic.level<-bird.traits[,2]
body.mass<-bird.traits[,3]
nest.cont.species.birds.data<-as.data.frame(cbind(nest.cont.species.birds,forest.dep, trophic.level,body.mass))

m.traits.nest.cont.birds<-lm(nest.cont.species.birds ~ scale(forest.dep) + scale(trophic.level) + scale(log10(body.mass)), data = nest.cont.species.birds.data)# 
summary(m.traits.nest.cont.birds)
#plot(m.traits.nest.cont.birds)
all.m.traits.nest.cont.birds<-dredge(m.traits.nest.cont.birds)
ave.all.m.traits.nest.cont.birds<-model.avg(all.m.traits.nest.cont.birds,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.nest.cont.birds) 
confint(ave.all.m.traits.nest.cont.birds)
importance(ave.all.m.traits.nest.cont.birds)

# Frogs
frog.traits<-read.table("node_species_metrics_frogs.txt", header = T)
rep.mode<-frog.traits[,8]
habitat.div<-frog.traits[,9]
body.size<-as.numeric(frog.traits[,10])
nest.cont.species.frogs<-as.vector(frog.traits$nested_contribution)
nest.cont.species.frogs.data<-as.data.frame(cbind(nest.cont.species.frogs,rep.mode, habitat.div, body.size))

m.traits.nest.cont.frogs<-glm(nest.cont.species.frogs+3 ~ scale(rep.mode) + scale(habitat.div) + scale(log10(body.size)), family = gaussian(link = "log"), data = nest.cont.species.frogs.data)# 
summary(m.traits.nest.cont.frogs)
#plot(m.traits.nest.cont.frogs)
all.m.traits.nest.cont.frogs<-dredge(m.traits.nest.cont.frogs)
ave.all.m.traits.nest.cont.frogs<-model.avg(all.m.traits.nest.cont.frogs,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.nest.cont.frogs) 
confint(ave.all.m.traits.nest.cont.frogs)
importance(ave.all.m.traits.nest.cont.frogs)


# Dung beetles
beetles.traits<-read.table("dung-beetles.traits.txt", header = T) 
body.size<-beetles.traits[,1]
relocation<-beetles.traits[,2]
diet<-beetles.traits[,3]
nest.cont.species.beetles.data<-as.data.frame(cbind(nest.cont.species.beetles,body.size,relocation,diet))

body.size2<-log10(body.size[-c(32)] ) # excluding one species for which we do not have trait data
nest.cont.species.beetles2<-nest.cont.species.beetles[-c(32)] 
relocation2<-relocation[-c(32)] 
diet2<-diet[-c(32)]
nest.cont.species.beetles.data2<-nest.cont.species.beetles.data[-c(32),]

m.traits.nest.cont.beetles<-glm(nest.cont.species.beetles2 ~ scale(body.size2) + relocation2 + diet2, data = nest.cont.species.beetles.data2)
summary(m.traits.nest.cont.beetles)
#plot(m.traits.nest.cont.beetles)
all.m.traits.nest.cont.beetles<-dredge(m.traits.nest.cont.beetles)
ave.all.m.traits.nest.cont.beetles<-model.avg(all.m.traits.nest.cont.beetles,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.nest.cont.beetles) 
confint(ave.all.m.traits.nest.cont.beetles)
importance(ave.all.m.traits.nest.cont.beetles)

# BEES
setwd("~/R/redes/data")
bees.traits<-read.table("bees.traits.txt", header = T) 

cor(bees.traits) # all the measures are highly correlated
body.length<-bees.traits[,1]; hist(body.length)
hist(log10(body.length))

nest.cont.species.bees.data<-as.data.frame(cbind(nest.cont.species.bees,body.length))
nest.cont.species.bees<-bees$nested_contribution

hist(nest.cont.species.bees); descdist(nest.cont.species.bees, discrete = FALSE) 

m.traits.nest.cont.bees<-lm(nest.cont.species.bees~ body.length)# 
summary(m.traits.nest.cont.bees)
#plot(m.traits.nest.cont.bees)
confint(m.traits.nest.cont.bees)

m.traits.nest.cont.bees<-lm(nest.cont.species.bees ~ log10(body.length))# 
summary(m.traits.nest.cont.bees)

plot(body.length,nest.cont.species.bees)
abline(glm(nest.cont.species.bees[-c(22)]~body.length[-c(22)]))

# TREES
trees.traits<-read.table("trees.traits.txt", header = T) 
wood_density<-trees.traits[,1]
seed_mass<-trees.traits[,2]
vertical_strat<-trees.traits[,3]
dispersal_mode<-trees.traits[,4]
regeneration<-trees.traits[,5]
nest.cont.species.trees.data<-as.data.frame(cbind(nest.cont.species.trees,wood_density,seed_mass,vertical_strat,dispersal_mode,regeneration))

wood_density2<-scale(wood_density)
seed_mass2<-scale(seed_mass)

m.traits.nest.cont.trees<-lm(nest.cont.species.trees ~ wood_density2 + seed_mass2 + vertical_strat + dispersal_mode + regeneration, data = nest.cont.species.trees.data)# 
summary(m.traits.nest.cont.trees)
#plot(m.traits.nest.cont.trees)
all.m.traits.nest.cont.trees<-dredge(m.traits.nest.cont.trees)
ave.all.m.traits.nest.cont.trees<-model.avg(all.m.traits.nest.cont.trees,rank = "AICc", revised.var = TRUE)
summary(ave.all.m.traits.nest.cont.trees) 
confint(ave.all.m.traits.nest.cont.trees)
importance(ave.all.m.traits.nest.cont.trees)


