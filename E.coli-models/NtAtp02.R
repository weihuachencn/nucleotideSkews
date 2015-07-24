# load the necessary functions
## -- by Martin Lercher: June 2014 --
rm(list=ls());

library(sybil)
library(glpkAPI)

######################################################################################################
#
# Questions:
#
# + Ratio of nucleotide investment / amino acid investment (in biomass)
# +  % of energy used for nucleotide production relativ to total cellular investment
# +  % of energy used for amino acid production relativ to total cellular investment
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

ex <- findExchReact(m)

# set glucose uptake limit to 1.0:
m <- changeBounds(m, ex["EX_glc(e)"], lb = -1)

# What's the environment?
upt <- uptReact(ex)
ex[upt]
react_name(m)[76] # check if cobalamin can be used as energy source; if yes switch it off

# remove oxygen from the environment?
# m <- changeBounds(m, ex["EX_o2(e)"], lb = 0)  # 

# biomass reaction:
shrinkMatrix(m, j="Ec_biomass_iJO1366_core_53p95M")

######################################################################################################
# Energy production 
######################################################################################################

# how much ADP -> ATP can E. coli produce from the glucose?
m <- changeBounds(m, ex["EX_glc(e)"], lb = -1)  # set glucose uptake limit to 1.0
m2 <- addReact(m, "burnATP", c("atp[c]", "h2o[c]", "adp[c]", "pi[c]", "h[c]"), c(-1, -1, 1, 1, 1))   # same as ATPM
m2 <- changeObjFunc(m2, react = "burnATP")  
atp.per.glc <- lp_obj( optimizeProb(m2) )  # 23.500
atp.per.glc
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
# RNA nucleotide production (see Suppl. Table 6 of Orth et al.)
######################################################################################################
# Exchange reactions for RNA nucleotides:
# A: "EX_ade(e)", C: "EX_csn(e)", U: "EX_ura(e)", G: "EX_gua(e)"

# minimize the amount of glucose taken up:
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
# with oxygen: U<C<G<A, U+A<G+C
# A        C        G        U 
# 21.13552 15.77926 20.36345 13.39720 
# without oxygen: U<C<A<G, U+A<G+C 
#  A     C     G     U 
# 14.00  8.50 15.25  6.00 

# The amount of glucose needed to make one molecule of each RNA nucleotide (this includes both carbon and energy!):
# G 0.8665298; A 0.899384; U 0.5700935; C 0.6714579
# => U < C < G < A; A+U < G+C


# average difference between production costs:
avg.diff <- 0
for (i in 1:3) {
  for (j in (i+1):4) {
   avg.diff <- avg.diff + abs( atp.per.RNA[i] - atp.per.RNA[j] )
  }
}
avg.diff <- avg.diff / 6  
avg.diff  # 4.633195 = average difference between RNA molecules measured in ATP


######################################################################################################
# DNA nucleotide production (see Suppl. Table 6 of Orth et al.)
######################################################################################################
# Exchange reactions for DNA nucleotides do not yet exist
glc.per.mol.DNA <- c()  # use this to store the results

# Add exchange reactions for the four DNA nucleotides:
m.nt <- mg
m.nt <- addReact(m.nt, "datp_tr", c("datp[c]", "datp[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(m.nt, "dctp_tr", c("dctp[c]", "dctp[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(m.nt, "dgtp_tr", c("dgtp[c]", "dgtp[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(m.nt, "dttp_tr", c("dttp[c]", "dttp[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(m.nt, "EX_datp(e)", c("datp[e]"), c(-1))   # transport reaction to NIL
m.nt <- addReact(m.nt, "EX_dgtp(e)", c("dgtp[e]"), c(-1))   # transport reaction to NIL
m.nt <- addReact(m.nt, "EX_dctp(e)", c("dctp[e]"), c(-1))   # transport reaction to NIL
m.nt <- addReact(m.nt, "EX_dttp(e)", c("dttp[e]"), c(-1))   # transport reaction to NIL

# how much of glucose is needed per Guanine (dgtp) etc. de novo synthesis?
m.nt.2 <- changeBounds(m.nt, "EX_datp(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
glc.per.mol.DNA["A"] <- -lp_obj(optimizeProb(m.nt.2))
m.nt.2 <- changeBounds(m.nt, "EX_dctp(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
glc.per.mol.DNA["C"] <- -lp_obj(optimizeProb(m.nt.2))
m.nt.2 <- changeBounds(m.nt, "EX_dgtp(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
glc.per.mol.DNA["G"] <- -lp_obj(optimizeProb(m.nt.2))
m.nt.2 <- changeBounds(m.nt, "EX_dttp(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
glc.per.mol.DNA["T"] <- -lp_obj(optimizeProb(m.nt.2))
glc.per.mol.DNA

atp.per.DNA <- glc.per.mol.DNA * atp.per.glc
atp.per.DNA
# A        C        G        T 
# 46.12748 40.40561 45.08045 45.72897 

# The number of glucose molecules needed to make one molecule of each DNA nucleotide (this includes both carbon and energy!):
# A 1.962871; C 1.719388; G: 1.918317; T 1.945913
# => C < G < T < A; G+C < A+T
# without oxygen: T<C<A<G, T+A<G+C 
# A        C        G        T 
# 7.733333 5.433333 8.133333 5.233333 

# average difference between production costs:
avg.diff <- 0
for (i in 1:3) {
  for (j in (i+1):4) {
    avg.diff <- avg.diff + abs( atp.per.DNA[i] - atp.per.DNA[j] )
  }
}
avg.diff <- avg.diff / 6  
avg.diff  # 2.969018 = average difference between DNA molecules measured in ATP

######################################################################################################
# How much glucose is needed to make different parts of the biomass? (see Suppl. Table 6 of Orth et al.)
######################################################################################################

# Glucose needed for 1 unit total biomass:
m.nt <- changeBounds(mg, "Ec_biomass_iJO1366_core_53p95M", lb=1) # require at least 1 biomass unit to be produced
m.nt <- changeBounds(m.nt, ex["EX_glc(e)"], lb = -1000)
m.nt <- changeObjFunc(m.nt, react = "EX_glc(e)", obj_coef = 1)  
g.tot <- lp_obj(optimizeProb(m.nt)) # EX_glc(e) = -10.045840

# Biomass reaction restricted to RNA (assuming all de-novo synthesized ATP is for RNA):
m.nt <- addReact(mg, "biomass.RNA", c("atp[c]", "ctp[c]", "gtp[c]", "utp[c]"), c(-0.169975, -0.129799, -0.209121, -0.140101))
# Produce RNA for one unit of biomass and estimate glucose required:
m.nt <- changeBounds(m.nt, "biomass.RNA", lb=1) # require at least 1 biomass unit to be produced
# minimize glucose production for this:
m.nt <- changeBounds(m.nt, ex["EX_glc(e)"], lb = -1000)
m.nt <- changeObjFunc(m.nt, react = "EX_glc(e)", obj_coef = 1)  
g.RNA <- lp_obj( optimizeProb(m.nt) ) # EX_glc(e) = -1.137113
g.RNA / g.tot   # 0.1131924

# Biomass reaction restricted to DNA:
m.nt <- addReact(mg, "biomass.DNA", c("datp[c]", "dctp[c]", "dgtp[c]", "dttp[c]"), c(-0.024805, -0.025612, -0.025612, -0.024805))
# Produce DNA for one unit of biomass and estimate glucose required:
m.nt <- changeBounds(m.nt, "biomass.DNA", lb=1) # require at least 1 biomass unit to be produced
# minimize glucose production for this:
m.nt <- changeBounds(m.nt, ex["EX_glc(e)"], lb = -1000)
m.nt <- changeObjFunc(m.nt, react = "EX_glc(e)", obj_coef = 1)  
g.DNA <- lp_obj( optimizeProb(m.nt) ) # EX_glc(e) = -0.189771 
g.DNA / g.tot  #  0.01889048

# Biomass reaction restricted to Amino Acids:
m.aa <- addReact(mg, "biomass.aa", c("ala_DASH_L[c]", "arg_DASH_L[c]", "asn_DASH_L[c]", "asp_DASH_L[c]", "cys_DASH_L[c]", "gln_DASH_L[c]", "glu_DASH_L[c]", "gly[c]", "his_DASH_L[c]", "ile_DASH_L[c]", "leu_DASH_L[c]", "lys_DASH_L[c]", "met_DASH_L[c]", "phe_DASH_L[c]", "pro_DASH_L[c]", "ser_DASH_L[c]", "thr_DASH_L[c]", "trp_DASH_L[c]", "tyr_DASH_L[c]", "val_DASH_L[c]"), 
                 c(-0.499149, -0.28742, -0.234232, -0.234232, -0.088988, -0.255712, -0.255712, -0.595297, -0.092056,  -0.282306, -0.437778, -0.333448, -0.149336, -0.180021, -0.214798, -0.209684, -0.246506, -0.055234, -0.133993, -0.411184) )
# Produce DNA for one unit of biomass and estimate glucose required:
m.aa <- changeBounds(m.aa, "biomass.aa", lb=1) # require at least 1 biomass unit to be produced
# minimize glucose production for this:
m.aa <- changeBounds(m.aa, ex["EX_glc(e)"], lb = -1000)
m.aa <- changeObjFunc(m.aa, react = "EX_glc(e)", obj_coef = 1)  
g.AA <- lp_obj( optimizeProb(m.aa) ) # EX_glc(e) = -4.753151
g.AA / g.tot   # 0.4731462

# relative investments:
g.info <- g.DNA + g.RNA + g.AA  # all information molecules (DNA, RNA, AA) together
g.info / g.tot   # 0.6052291 of total glucose investment
(g.DNA + g.RNA) / g.tot  # 0.1320829 of total glucose investment
g.RNA / g.DNA  # 5.992033 glucose investment into RNA vs. DNA


######################################################################################################
# Average difference in AA production costs
######################################################################################################

# add (hypothetical) transport reactions:
m.aa <- mg
m.aa <- addReact(m.aa, "ala_DASH_L_tr", c("ala_DASH_L[c]", "ala_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_ala_DASH_L(e)", c("ala_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "arg_DASH_L_tr", c("arg_DASH_L[c]", "arg_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_arg_DASH_L(e)", c("arg_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "asn_DASH_L_tr", c("asn_DASH_L[c]", "asn_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_asn_DASH_L(e)", c("asn_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "asp_DASH_L_tr", c("asp_DASH_L[c]", "asp_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_asp_DASH_L(e)", c("asp_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "cys_DASH_L_tr", c("cys_DASH_L[c]", "cys_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_cys_DASH_L(e)", c("cys_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "gln_DASH_L_tr", c("gln_DASH_L[c]", "gln_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_gln_DASH_L(e)", c("gln_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "glu_DASH_L_tr", c("glu_DASH_L[c]", "glu_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_glu_DASH_L(e)", c("glu_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "gly_tr", c("gly[c]", "gly[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_gly(e)", c("gly[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "his_DASH_L_tr", c("his_DASH_L[c]", "his_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_his_DASH_L(e)", c("his_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "his_DASH_L_tr", c("ile_DASH_L[c]", "ile_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_ile_DASH_L(e)", c("ile_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "leu_DASH_L_tr", c("leu_DASH_L[c]", "leu_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_leu_DASH_L(e)", c("leu_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "lys_DASH_L_tr", c("lys_DASH_L[c]", "lys_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_lys_DASH_L(e)", c("lys_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "met_DASH_L_tr", c("met_DASH_L[c]", "met_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_met_DASH_L(e)", c("met_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "phe_DASH_L_tr", c("phe_DASH_L[c]", "phe_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_phe_DASH_L(e)", c("phe_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "pro_DASH_L_tr", c("pro_DASH_L[c]", "pro_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_pro_DASH_L(e)", c("pro_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "ser_DASH_L_tr", c("ser_DASH_L[c]", "ser_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_ser_DASH_L(e)", c("ser_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "thr_DASH_L_tr", c("thr_DASH_L[c]", "thr_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_thr_DASH_L(e)", c("thr_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "trp_DASH_L_tr", c("trp_DASH_L[c]", "trp_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_trp_DASH_L(e)", c("trp_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "tyr_DASH_L_tr", c("tyr_DASH_L[c]", "tyr_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_tyr_DASH_L(e)", c("tyr_DASH_L[e]"), c(-1))   # transport reaction to NIL
m.aa <- addReact(m.aa, "val_DASH_L_tr", c("val_DASH_L[c]", "val_DASH_L[e]"), c(-1, 1))   # transport reaction to extracellular
m.aa <- addReact(m.aa, "EX_val_DASH_L(e)", c("val_DASH_L[e]"), c(-1))   # transport reaction to NIL

aa.transports <- c("EX_ala_DASH_L(e)", "EX_arg_DASH_L(e)", "EX_asn_DASH_L(e)", "EX_asp_DASH_L(e)", "EX_cys_DASH_L(e)", "EX_gln_DASH_L(e)", "EX_glu_DASH_L(e)", 
                    "EX_gly(e)", "EX_his_DASH_L(e)", "EX_ile_DASH_L(e)", "EX_leu_DASH_L(e)", "EX_lys_DASH_L(e)", "EX_met_DASH_L(e)", "EX_phe_DASH_L(e)", 
                    "EX_pro_DASH_L(e)", "EX_ser_DASH_L(e)", "EX_thr_DASH_L(e)", "EX_trp_DASH_L(e)", "EX_tyr_DASH_L(e)", "EX_val_DASH_L(e)")

# find minimal glucose for each: 
glc.per.aa <- c()
for (aa in aa.transports) {
  m.aa.o <- changeBounds(m.aa, aa, lb=1) # require at least 1 mol of this AA to be produced
  glc.per.aa[aa] <- -lp_obj(optimizeProb(m.aa.o))
}
summary( glc.per.aa)
atp.per.aa <- glc.per.aa * atp.per.glc
summary( atp.per.aa )

save.image(file = "working_NtAtp02.RData");

# average difference between production costs:
avg.diff <- 0
for (i in 1:19) {
  for (j in (i+1):20) {
    avg.diff <- avg.diff + abs( atp.per.aa[i] - atp.per.aa[j] )
  }
}
avg.diff <- avg.diff / 190 
avg.diff  # 13.18006 = average difference between AA molecules measured in ATP

