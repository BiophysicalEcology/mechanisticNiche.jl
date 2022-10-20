

#### SKIN EVAPORATION

function skinevap(;
  # Tc = 25,
  Tskin=25.1,
  Gevap=1.177235e-09,
  ψbody=-7.07 * 100,
  pct_wet=0.1, # fraction of the skin that is wet
  Aeff=1.192505e-05, # effective area acting as a free-water exchanger
  Atot=0.01325006,
  hd=0.02522706,
  pct_eyes=0.03 / 100,
  Ta=20,
  rh=50,
  vel=0.1,
  bp=101325
)

  MW = 0.018 #! molar mass of water, kg/mol
  RGC = 8.314 #! gas constant, J/mol/K
  rh_skin = exp(ψbody / (RGC / MW * (Tskin + 273.15))) * 100 #

  wb = 0
  dp = 999

  wairSurf = wetair(db=Tskin, wb=wb, rh=rh_skin, dp=dp, bp=bp)
  ρsurf = wairSurf.vd


  wair = wetair(db=Ta, wb=wb, rh=rh, dp=dp, bp=bp)
  ρair = wair.vd

  # peyes = pct_eyes / 100
  peyes = pct_eyes
  Weyes = hd * peyes * Atot * (ρsurf - ρair)

  Wresp = Gevap / 1000

  if (Weyes > 0)
    Wcut = (Aeff - peyes * Atot * pct_wet) * hd * (ρsurf - ρair)
  else
    Wcut = Aeff * hd * (ρsurf - ρair)
  end

  water = Weyes + Wresp + Wcut

  HTOVPR = 2.5012E+06 - 2.3787E+03 * Ta  # FROM DRYAIR: LATENT HEAT OF VAPORIZATION
  Qsevap = (Weyes + Wcut) * HTOVPR

  # KG/S TO G/S
  Weyes = Weyes * 1000
  Wresp = Wresp * 1000
  Wcut = Wcut * 1000
  Wevap = water * 1000

  return (QSEVAP=Qsevap, WEVAP=Wevap, WRESP=Wresp, WCUT=Wcut, WEYES=Weyes)
end



function skinevap(o::fixedParams,
  e::fixedEnvParams,
  v::envVars,
  s::stateVariables,
  Gevap=1.177235e-9,
  hd=0.02522706
)

  @unpack ψbody, pct_wet, pct_eyes = o
  @unpack bp = e
  @unpack Ta, rh, vel = v
  @unpack Tskin = s

  Ta = Ta[1]
  rh = rh[1]
  vel = vel[1]

  # or maybe in a later version include a geom object as input (only call geom ~once per time-step)
  geom = geometry(o)
  Atot = geom.AREA
  Aeff = geom.AEFF

  MW = 0.018 #! molar mass of water, kg/mol
  RGC = 8.314 #! gas constant, J/mol/K
  rh_skin = exp(ψbody / (RGC / MW * (Tskin + 273.15))) * 100 #

  wb = 0
  dp = 999

  wairSurf = wetair(db=Tskin, wb=wb, rh=rh_skin, dp=dp, bp=bp)
  ρsurf = wairSurf.vd


  wair = wetair(db=Ta, wb=wb, rh=rh, dp=dp, bp=bp)
  ρair = wair.vd

  peyes = pct_eyes / 100
  Weyes = hd * peyes * Atot * (ρsurf - ρair)

  Wresp = Gevap / 1000

  if (Weyes > 0)
    Wcut = (Aeff - peyes * Atot * pct_wet) * hd * (ρsurf - ρair)
  else
    Wcut = Aeff * hd * (ρsurf - ρair)
  end

  water = Weyes + Wresp + Wcut

  HTOVPR = 2.5012E+06 - 2.3787E+03 * Ta  # FROM DRYAIR: LATENT HEAT OF VAPORIZATION
  Qsevap = (Weyes + Wcut) * HTOVPR

  # KG/S TO G/S
  Weyes = Weyes * 1000
  Wresp = Wresp * 1000
  Wcut = Wcut * 1000
  Wevap = water * 1000

  return (QSEVAP=Qsevap, WEVAP=Wevap, WRESP=Wresp, WCUT=Wcut, WEYES=Weyes)
end


# for zero solver
function skinevap(Tx,
  o::fixedParams,
  e::fixedEnvParams,
  v::envVars,
  Gevap=1.177235e-9,
  hd=0.02522706
)

  Tskin = Tx
  
  @unpack ψbody, pct_wet, pct_eyes = o
  @unpack bp = e
  @unpack Ta, rh, vel = v

  Ta = Ta[1]
  rh = rh[1]
  vel = vel[1]

  # or maybe in a later version include a geom object as input (only call geom ~once per time-step)
  geom = geometry(o)
  Atot = geom.AREA
  Aeff = geom.AEFF

  MW = 0.018 #! molar mass of water, kg/mol
  RGC = 8.314 #! gas constant, J/mol/K
  rh_skin = exp(ψbody / (RGC / MW * (Tskin + 273.15))) * 100 #

  wb = 0
  dp = 999

  wairSurf = wetair(db=Tskin, wb=wb, rh=rh_skin, dp=dp, bp=bp)
  ρsurf = wairSurf.vd


  wair = wetair(db=Ta, wb=wb, rh=rh, dp=dp, bp=bp)
  ρair = wair.vd

  # peyes = pct_eyes / 100
  peyes = pct_eyes
  Weyes = hd * peyes * Atot * (ρsurf - ρair)

  Wresp = Gevap / 1000

  if (Weyes > 0)
    Wcut = (Aeff - peyes * Atot * pct_wet) * hd * (ρsurf - ρair)
  else
    Wcut = Aeff * hd * (ρsurf - ρair)
  end

  water = Weyes + Wresp + Wcut

  HTOVPR = 2.5012E+06 - 2.3787E+03 * Ta  # FROM DRYAIR: LATENT HEAT OF VAPORIZATION
  Qsevap = (Weyes + Wcut) * HTOVPR

  # KG/S TO G/S
  Weyes = Weyes * 1000
  Wresp = Wresp * 1000
  Wcut = Wcut * 1000
  Wevap = water * 1000

  return (QSEVAP=Qsevap, WEVAP=Wevap, WRESP=Wresp, WCUT=Wcut, WEYES=Weyes)
end


# @show skinevap()