

setwd('~/Dropbox/biophys_tools')

library(NicheMapR)

# micro <- micro_global(loc = cbind(-5.41, 43.4), runmoist = T, snowmodel = T, writecsv = T)
# save(micro, file='./micro_ast.RData')

load('micro_ast.RData')
metout <- data.frame(micro$metout)
head(metout)

plot(metout$TALOC, type='l')
