

setwd('~/Dropbox/biophys_tools')
library(NicheMapR)

rm(list=ls())

load('./micro_ast.RData')

# extract full sun conditions
metout <- as.data.frame(micro$metout)
soil <- as.data.frame(micro$soil)

# get required inputs
TAs <- metout$TALOC
TGRDs <- soil$D0cm
TSKYs <- metout$TSKYC
VELs <- metout$VLOC
RHs <- metout$RHLOC
QSOLRs <- metout$SOLR
Zs <- metout$ZEN
K_subs <- micro$tcond[, 3]
# micro$tcond[, 3] <- 4; K_subs <- micro$tcond[, 3]

# use ectoR_devel to compute body temperature in open without respiratory heat loss
# or metabolic heat gain
TC <- unlist(lapply(1:length(TAs), function(x){ectoR_devel(
  Ww_g = 30, # wet weight, g
  shape = 3, # using frog geometry
  pct_wet = 80,
  alpha = 0.85, # solar absorptivity
  M_1 = 0, # turn of metabolic heat
  postur = 0, # average posture, half way between normal and parallel to sun
  pantmax = 0, # turn off respiratory heat exchange
  pct_cond = 30, 
  K_sub = K_subs[x],
  #K_sub = 0.1,
  alpha_sub = (1 - micro$REF), # substrate absorptivity used in microclimate model
  elev = micro$elev, # elevation from microclimate model
  TA = TAs[x], # air temperature at frog height from microclimate model, deg C
  TGRD = TGRDs[x], # ground temperature from microclimate model, deg C
  TSUBST = TGRDs[x], # ground temperature from microclimate model, deg C
  TSKY = TSKYs[x], # sky temperature from microclimate model, deg C
  VEL = VELs[x], # wind speed from microclimate model, m/s
  RH = RHs[x], # relative humidity from microclimate model, %
  QSOLR = QSOLRs[x], # total horizontal plane solar radiation from microclimate model, W/m2
  Z = Zs[x] # solar zenith angle from microclimate model, degrees
)$TC})) # run ectoR_devel across environments

# run ectotherm model for a non-behaving animal without respiratory heat loss
# or metabolic heat gain
# micro$tcond[, 3] <- 0.1 # set K_sub (equally, you could make the input for ectoR_devel come from micro$tcond)
ecto <- ectotherm(
  Ww_g = 30, # wet weight, g
  shape = 3, # using frog geometry
  pct_wet = 0.1,
  alpha_min = 0.85, # minimum solar absorptivity
  alpha_max = 0.85, # maximum solar absorptivity
  M_1 = 0, # turn of metabolic heat
  postur = 0, # average posture, half way between normal and parallel to sun
  pantmax = 0, # turn off respiratory heat exchange
  pct_cond = 40, 
  live = 0
)

# extract results
environ <- as.data.frame(ecto$environ)
enbal <- as.data.frame(ecto$enbal)
masbal <- as.data.frame(ecto$masbal)
TC_ectotherm <- environ$TC

# compare
time <- micro$dates
plot(time, TC, type = 'l', ylim = c(0, 40), ylab = 'body temperature, deg C', xlab = 'month of year')
points(time, TC_ectotherm, type = 'l', col = 2)

# write.csv(TC, file="TCnmr.csv",row.names=F)
