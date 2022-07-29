############################################
## Purpose: trying to align transmission (incidence) between leapfrog and Spectrum
## Requires pulling dev_add_hivp_entrant from leapfrog
############################################


library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-no-hiv-deaths_spectrum-v6.13_2022-02-12.pjnz"
#pjnz <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art_spectrum-v6.13_2022-02-12.pjnz"


## Check that prevalence, deaths and incidence  matches between
## the two models
pjnz1 <- test_path(pjnz)

demp <- prepare_leapfrog_demp(pjnz1)
hivp <- prepare_leapfrog_projp(pjnz1)
hivp$incrr_sex = round(hivp$incrr_sex, 2)
hivp$incrr_age = round(hivp$incrr_age, 4)

if(pjnz1 == "tests/testthat/../testdata/spectrum/v6.13/bwa_aim-adult-no-art-no-hiv-deaths_spectrum-v6.13_2022-02-12.pjnz"){
  hivp$cd4_mort_full[hivp$cd4_mort_full > 0] <- 0
}

# ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
# ## in EPP-ASM preparation
# demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
# demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)

lmod <- leapfrogR(demp, hivp)

specres <- eppasm::read_hivproj_output(pjnz1)


## Prevalence

prev.l = lmod$hivpop1[16:81,,] 
prev.s = specres$hivpop[16:81,,]

prev.l = as.data.frame.table(prev.l, responseName = "lfrog")
prev.l = left_join(prev.l, data.frame(Var1 = unique(prev.l$Var1), age = 15:80))
prev.l = left_join(prev.l, data.frame(Var2 = unique(prev.l$Var2), sex = c("male", "female")))
prev.l = left_join(prev.l, data.frame(Var3 = unique(prev.l$Var3), year = 1970:2030))
prev.l = select(prev.l, c("age", "sex", "year", "lfrog"))

prev.s = as.data.frame.table(prev.s, responseName = "spec")
colnames(prev.s)[1:3] = c("age", "sex", "year")
prev.s$age = as.integer(as.character(prev.s$age))
prev.s$sex = (as.character(prev.s$sex))
prev.s$year = as.integer(as.character(prev.s$year))

prev <- full_join(prev.l, prev.s)

prev <- prev %>% summarise(age, sex, year, lfrog, spec, diff = lfrog - spec)


##shows where differences between lfrog and spectrum are in prevalence
ggplot(prev, aes(year, age, fill = diff)) + geom_tile() + ggtitle("Prevalence") + scale_fill_gradientn(colours = c("red", "white", "green"),
                                                                                                                    values = rescale(c(-225,0,160)),
                                                                                                                    guide = "colorbar", limits=c(-225,160)) +
  labs(caption = "Red means Spectrum larger than lfrog and green is v.v.") + facet_wrap(~sex,)



prev <- prev %>% select(-diff) %>% gather(variable, value, -age, -sex, -year)
prev <- prev %>% summarise(age = floor(age / 5) * 5, sex, year, variable, value)
prev <- prev %>% filter(!is.nan(value))
prev <- prev %>% group_by(age, sex, year, variable) %>% summarise(value = sum(value))
prev <- unique(prev)
ggplot(prev, aes(year, value, col = as.factor(variable))) + geom_line(aes(lty= as.factor(sex)), lwd = 1.1) + facet_wrap(~age, scales = "free") + theme_bw() +
  ggtitle("Prevalence")




## Incidence, only going to compare the single year age groups bc thats what doesn't align
inc.l <- lmod$infections[16:81,,-1]
inc.l = as.data.frame.table(inc.l, responseName = "lfrog")
inc.l = left_join(inc.l, data.frame(Var1 = unique(inc.l$Var1), age = 15:80))
inc.l = left_join(inc.l, data.frame(Var2 = unique(inc.l$Var2), sex = c("male", "female")))
inc.l = left_join(inc.l, data.frame(Var3 = unique(inc.l$Var3), year = 1971:2030))
inc.l = select(inc.l, c("age", "sex", "year", "lfrog"))
inc.l$age = as.factor(inc.l$age)
inc.l$sex = as.factor(inc.l$sex)
inc.l$year = as.factor(inc.l$year)


inc.s <- specres$infections[16:81,,-1]
inc.s = as.data.frame.table(inc.s, responseName = "spec")

inc <- full_join(inc.l, inc.s)
inc <- inc %>% summarise(age = as.integer(as.character(age)), sex, year = as.integer(as.character(year)), lfrog, spec, diff = lfrog - spec)

##shows where differences between lfrog and spectrum are in prevalence
ggplot(inc, aes(year, age, fill = diff)) + geom_tile() + ggtitle('Incidence') + scale_fill_gradientn(colours = c("red", "white", "green"),
                                                                                                                   values = rescale(c(-25,0,25)),
                                                                                                                   guide = "colorbar", limits=c(-25,25)) +
  labs(caption = "Red means Spectrum larger than lfrog and green is v.v.") + facet_wrap(~sex)

##HIV DEATHS: Note this is set to zero so not super informative here 
hivdeaths.l <- lmod$hivdeaths[16:81,,-1]
hivdeaths.l = as.data.frame.table(hivdeaths.l, responseName = "lfrog")
hivdeaths.l = left_join(hivdeaths.l, data.frame(Var1 = unique(hivdeaths.l$Var1), age = 15:80))
hivdeaths.l = left_join(hivdeaths.l, data.frame(Var2 = unique(hivdeaths.l$Var2), sex = c("male", "female")))
hivdeaths.l = left_join(hivdeaths.l, data.frame(Var3 = unique(hivdeaths.l$Var3), year = 1971:2030))
hivdeaths.l = select(hivdeaths.l, c("age", "sex", "year", "lfrog"))
hivdeaths.l$age = as.factor(hivdeaths.l$age)
hivdeaths.l$sex = as.factor(hivdeaths.l$sex)
hivdeaths.l$year = as.factor(hivdeaths.l$year)


hivdeaths.s <- specres$hivdeaths[16:81,,-1]
hivdeaths.s = as.data.frame.table(hivdeaths.s, responseName = "spec")

hivdeaths <- full_join(hivdeaths.l, hivdeaths.s)
hivdeaths <- hivdeaths %>% summarise(age = as.integer(as.character(age)), sex, year = as.integer(as.character(year)), lfrog, spec, diff = lfrog - spec)

##shows where differences between lfrog and spectrum are in prevalence
ggplot(hivdeaths, aes(year, age, fill = diff)) + geom_tile() + ggtitle("Single year age groups") + scale_fill_gradientn(colours = c("red", "white", "green"),
                                                                                                                  values = rescale(c(-12,0,10)),
                                                                                                                  guide = "colorbar", limits=c(-12,10)) +
  labs(caption = "Red means Spectrum larger than lfrog and green is v.v.", subtitle = 'Deaths')
hivdeaths <- hivdeaths %>% select(-diff) %>% gather(variable, value, -age, -sex, -year)

# ggplot(hivdeaths, aes(year, value, col = as.factor(variable))) + geom_line(aes(lty= as.factor(sex))) + facet_wrap(~age, scales = "free") + theme_bw() + 
#   ggtitle("Single year age groups")
#


##Manually calculating infections
t = 7

hivn_ag <- lmod$totpop1[,,t-1] - lmod$hivpop1[,,t-1]
hivn_ag = hivn_ag[16:81,]

Xhivn <- apply(hivn_ag, 2, sum)
Xhivn_incagerr <- hivp$incrr_age[,,t] * hivn_ag

incrate_g <- hivp$incidinput[t] * hivp$incrr_sex[t] * (Xhivn[1] + Xhivn[2]) / (Xhivn[1] + hivp$incrr_sex[t] * Xhivn[2])

infections <- c()
for(age in 1:66){
  x = hivn_ag[age, 2] * incrate_g * hivp$incrr_age[age, 2, t] * Xhivn[2] / sum(Xhivn_incagerr[,2])
  infections <- c(infections, x)
}

lmod$infections[16,2,t]
infections[1]
infections-lmod$hivdeaths[16:81,2,t]-lmod$natdeaths_hivpop[16:81,2,t]

specres$infections[16,2,t]




lmod$hivpop[16,,25] * (1 - demp$Sx[17,,26])
lmod$natdeaths_hivpop[17,,26]

lmod$infections[16,,26] - lmod$hivdeaths[16,,26]
lmod$hivpop[16,,26] 

specres$natdeaths_hivpop <- specres$hivpop[16,,25]  * (1 - demp$Sx[17,,26])

specres$infections[16,,26] - lmod$hivdeaths[16,,26] - specres$natdeaths_hivpop 
specres$hivpop[16,,26]








