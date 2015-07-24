# load the necessary functions
## -- by Martin Lercher; June 2014 --
rm(list=ls());
## -- 
library(sybil)
library(glpkAPI)

######################################################################################################
#
# Questions:
#
#   Ratio of nucleotide investment / amino acid investment (in biomass)
#   % of energy used for nucleotide production relativ to total cellular investment
#   % of energy used for amino acid production relativ to total cellular investment
#   ATP required for each nucleotide
#   variance of energy investment into nucleotides
#   variance of energy investment into amino acids
#
######################################################################################################

######################################################################################################
# The modified E. coli model
######################################################################################################

# load E. coli model
m <- readTSVmod(prefix="iJO1366", quoteChar = "\"")
# remove the requirement to produce maintenance ATP:
m <- changeBounds(m, "ATPM", lb=0)

# What's the environment?
ex <- findExchReact(m)
# set glucose uptake limit to 1.0:
m <- changeBounds(m, ex["EX_glc(e)"], lb = -1)

upt <- uptReact(ex)
ex[upt]
react_name(m)[76] # check if cobalamin can be used as energy source; if yes switch it off

## !!! remove oxygen from the environment? !!! 
## -- uncomment the following line if you prefer anaerobic enviroment --
## m <- changeBounds(m, ex["EX_o2(e)"], lb = 0)  # 

# biomass reaction:
shrinkMatrix(m, j="Ec_biomass_iJO1366_core_53p95M")

######################################################################################################
# Energy production 
######################################################################################################

# how much ADP -> ATP can E. coli produce from the glucose?
m2 <- addReact(m, "burnATP", c("atp[c]", "h2o[c]", "adp[c]", "pi[c]", "h[c]"), c(-1, -1, 1, 1, 1))   # same as ATPM
m2 <- changeObjFunc(m2, react = "burnATP")  
optimizeProb(m2)  # 23.500
# alternatively:
# m3 <- changeObjFunc(m, react = "ATPM")  
# optimizeProb(m3)

# how much GDP -> GTP can E. coli produce from the glucose? The same as ATP, of course!
m3 <- addReact(m, "burnGTP", c("gtp[c]", "h2o[c]", "gdp[c]", "pi[c]", "h[c]"), c(-1, -1, 1, 1, 1))   # similar to ATPM
m3 <- changeObjFunc(m3, react = "burnGTP")  
optimizeProb(m3) # 23.500

# how much of ATP de novo synthesis is possible?
m4 <- addReact(m, "atp_tr", c("atp[c]", "atp[e]"), c(-1, 1))   # transport reaction to extracellular
m4 <- addReact(m4, "EX_atp(e)", "atp[e]", -1)   # exchange reaction extracellular -> NIL
m4 <- changeObjFunc(m4, react = "EX_atp(e)")
optimizeProb(m4)  # 0.530

######################################################################################################
# Nucleotide production
######################################################################################################
# Exchange reactions for nucleotides:
# "EX_thym(e)", A: "EX_ade(e)", C: "EX_csn(e)", U: "EX_ura(e)", G: "EX_gua(e)"

# minimize the amount of glucose taken up:
m <- changeBounds(m, ex["EX_glc(e)"], lb = -1000)    # allow unlimited uptake of glucose
mg <- changeObjFunc(m, react = "EX_glc(e)", obj_coef = 1)   # but minimize this uptake (given what needs to be produced)
glc.per.mol <- c()  # use this to store the results

# how much of glucose is needed per Guanine (gua) de novo synthesis?
# Guanine exchange: EX_gua(e)
m.nt <- addReact(mg, "gua_tr", c("gua[c]", "gua[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- changeBounds(m.nt, "EX_gua(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
o5 <- optimizeProb(m.nt)
glc.per.mol$G <- -lp_obj(o5)

# how much of glucose is needed per Adenine (ade) de novo synthesis?
# Adenine exchange: EX_ade(e)
m.nt <- addReact(mg, "ala_tr", c("ade[c]", "ade[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- changeBounds(m.nt, "EX_ade(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
o5 <- optimizeProb(m.nt)
glc.per.mol$A <- -lp_obj(o5)

# how much of glucose is needed per Thymine (thym) de novo synthesis?
# Adenine exchange: EX_thym(e)
m.nt <- addReact(mg, "thym_tr", c("thym[c]", "thym[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- changeBounds(m.nt, "EX_thym(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
o5 <- optimizeProb(m.nt)
glc.per.mol$T <- -lp_obj(o5)

# how much of glucose is needed per Uracil (ura) de novo synthesis?
# Adenine exchange: EX_ura(e)
m.nt <- addReact(mg, "ura_tr", c("ura[c]", "ura[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- changeBounds(m.nt, "EX_ura(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
o5 <- optimizeProb(m.nt)
glc.per.mol$U <- -lp_obj(o5)

# how much of glucose is needed per Cytosine (csn) de novo synthesis?
# Adenine exchange: EX_ura(e)
m.nt <- addReact(mg, "csn_tr", c("csn[c]", "csn[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- changeBounds(m.nt, "EX_csn(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
o5 <- optimizeProb(m.nt)
glc.per.mol$C <- -lp_obj(o5)

glc.per.mol
# The amount of glucose needed to make one molecule of each nucleotide (this includes both carbon and energy!):
# G 0.8665298; A 0.899384; T 0.8971963; U 0.5700935; C 0.6714579
# => U < C < G < T ~ A; A+U < G+C

