## -- for testing purpose only --

library(sybil)
library(sybilSBML);
library(glpkAPI)

## load an xml model --
m <- readSBMLmod(filename="~/05_skews_bacteria/7_math_models/Bacillus_subtilis.xml");
# remove the requirement to produce maintenance ATP:
m <- changeBounds(m, "ATPM", lb=0)

# set glucose uptake limit to 1.0:
ex <- findExchReact(m)
m <- changeBounds(m, ex["EX_glc(e)"], lb = -1);

m2 <- addReact(m, "burnATP", c("atp[c]", "h2o[c]", "adp[c]", "pi[c]", "h[c]"), c(-1, -1, 1, 1, 1))   # same as ATPM
m2 <- changeObjFunc(m2, react = "burnATP")  
atp.per.glc <- lp_obj( optimizeProb(m2) )  # 23.500
atp.per.glc

upt <- uptReact(ex)
ex[upt]

mg <- changeBounds(m, ex["EX_glc(e)"], lb = -1000)    # allow unlimited uptake of glucose
mg <- changeObjFunc(mg, react = "EX_glc(e)", obj_coef = 1)   # but minimize this uptake (given what needs to be produced)
glc.per.mol.RNA <- c()  # use this to store the results

# how much of glucose is needed per Guanine (gua) de novo synthesis etc.?
# add missing (hypothetical) transport reactions:
m.nt <- addReact(mg, "ade_tr", c("ade[c]", "ade[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(mg, "gua_tr", c("gua[c]", "gua[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(mg, "ura_tr", c("ura[c]", "ura[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(mg, "csn_tr", c("csn[c]", "csn[e]"), c(-1, 1))   # transport reaction to extracellular

m.nt.o <- changeBounds(m.nt, "EX_ade(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
glc.per.mol.RNA["A"] <- -lp_obj(optimizeProb(m.nt.o))
m.nt.o <- changeBounds(m.nt, "EX_csn(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
glc.per.mol.RNA["C"] <- -lp_obj(optimizeProb(m.nt.o))
m.nt.o <- changeBounds(m.nt, "EX_gua(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
glc.per.mol.RNA["G"] <- -lp_obj(optimizeProb(m.nt.o))
m.nt.o <- changeBounds(m.nt, "EX_ura(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
glc.per.mol.RNA["U"] <- -lp_obj(optimizeProb(m.nt.o))

atp.per.RNA <- glc.per.mol.RNA * atp.per.glc
atp.per.RNA



m2 <- readTSVmod(prefix="iJO1366", quoteChar = "\"")
m2 <- changeBounds(m2, "ATPM", lb=0);
