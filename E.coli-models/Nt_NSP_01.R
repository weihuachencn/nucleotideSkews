# load the necessary functions

library(sybil)
library(glpkAPI)

######################################################################################################
#
# Questions:
#
# + Ratio of nucleotide investment / amino acid investment (in biomass)
# +  % of NSP used for nucleotide production relativ to total cellular investment
# +  % of NSP used for amino acid production relativ to total cellular investment
#   NSP required for each nucleotide
#   variance of NSP investment into nucleotides
#   variance of NSP investment into amino acids
#
######################################################################################################

######################################################################################################
# The modified E. coli model
######################################################################################################

# load E. coli model
m <- readTSVmod(prefix="iJO1366", quoteChar = "\"")
# remove the requirement to produce maintenance ATP:
m <- changeBounds(m, "ATPM", lb=0)
# set glucose uptake limit to 1.0:
ex <- findExchReact(m)
m <- changeBounds(m, ex["EX_glc(e)"], lb = -1000)

# What's the environment?
upt <- uptReact(ex)
ex[upt]
react_name(m)[76] # check if cobalamin can be used as energy source; if yes switch it off

# remove oxygen from the environment?
# m <- changeBounds(m, ex["EX_o2(e)"], lb = 0)  # 

# biomass reaction:
# shrinkMatrix(m, j="Ec_biomass_iJO1366_core_53p95M")
shrinkMatrix(m, j="Ec_biomass_iJO1366_WT_53p95M")

######################################################################################################
# SET THE OBJECTIVE FUNCTION TO A COMBINATION OF NITROGEN, SULPHUR, AND PHOSPHOR
######################################################################################################

# minimize the sum of N + S + P taken up:
m2 <- changeObjFunc(m, react = c("EX_nh4(e)", "EX_so4(e)", "EX_pi(e)"), obj_coef = c(1, 1, 1))   # minimize this uptake (given what needs to be produced)
nsp.per.mol.RNA <- c()  # use this to store the results

######################################################################################################
# RNA nucleotide production (see Suppl. Table 6 of Orth et al.)
######################################################################################################
# Exchange reactions for RNA nucleotides:
# A: "EX_ade(e)", C: "EX_csn(e)", U: "EX_ura(e)", G: "EX_gua(e)"

# how much of NSP is needed per Guanine (gua) de novo synthesis etc.?
# add missing (hypothetical) transport reactions:
m.nt <- m2
m.nt <- addReact(m.nt, "ade_tr", c("ade[c]", "ade[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(m.nt, "gua_tr", c("gua[c]", "gua[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(m.nt, "ura_tr", c("ura[c]", "ura[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(m.nt, "csn_tr", c("csn[c]", "csn[e]"), c(-1, 1))   # transport reaction to extracellular

m.nt.o <- changeBounds(m.nt, "EX_ade(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
nsp.per.mol.RNA["A"] <- -lp_obj(optimizeProb(m.nt.o))
m.nt.o <- changeBounds(m.nt, "EX_csn(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
nsp.per.mol.RNA["C"] <- -lp_obj(optimizeProb(m.nt.o))
m.nt.o <- changeBounds(m.nt, "EX_gua(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
nsp.per.mol.RNA["G"] <- -lp_obj(optimizeProb(m.nt.o))
m.nt.o <- changeBounds(m.nt, "EX_ura(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
nsp.per.mol.RNA["U"] <- -lp_obj(optimizeProb(m.nt.o))

# result:
nsp.per.mol.RNA
# with oxygen: 
# A C G U 
# 5 3 5 2 
# ie., U < C < A = G
# i.e., U < A, C < G, U+A < G+C (just as in the energy case)
# without oxygen: 
# same!

# average difference between production costs:
avg.diff <- 0
for (i in 1:3) {
  for (j in (i+1):4) {
   avg.diff <- avg.diff + abs( nsp.per.mol.RNA[i] - nsp.per.mol.RNA[j] )
  }
}
avg.diff <- avg.diff / 6  
avg.diff  #  1.833333 = average difference between RNA molecules measured in NSP




##############################################################
### how many Ns required for the synthesis of each RNA --
# minimize the sum of N :
m2 <- changeObjFunc(m, react = c("EX_nh4(e)"), obj_coef = c(1))   # minimize this uptake (given what needs to be produced)
n.per.mol.RNA <- c()  # use this to store the results


######################################################################################################
# RNA nucleotide production (see Suppl. Table 6 of Orth et al.)
######################################################################################################
# Exchange reactions for RNA nucleotides:
# A: "EX_ade(e)", C: "EX_csn(e)", U: "EX_ura(e)", G: "EX_gua(e)"

# how much of n is needed per Guanine (gua) de novo synthesis etc.?
# add missing (hypothetical) tranort reactions:
m.nt <- m2
m.nt <- addReact(m.nt, "ade_tr", c("ade[c]", "ade[e]"), c(-1, 1))   # tranort reaction to extracellular
m.nt <- addReact(m.nt, "gua_tr", c("gua[c]", "gua[e]"), c(-1, 1))   # tranort reaction to extracellular
m.nt <- addReact(m.nt, "ura_tr", c("ura[c]", "ura[e]"), c(-1, 1))   # tranort reaction to extracellular
m.nt <- addReact(m.nt, "csn_tr", c("csn[c]", "csn[e]"), c(-1, 1))   # tranort reaction to extracellular

m.nt.o <- changeBounds(m.nt, "EX_ade(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
n.per.mol.RNA["A"] <- -lp_obj(optimizeProb(m.nt.o))
m.nt.o <- changeBounds(m.nt, "EX_csn(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
n.per.mol.RNA["C"] <- -lp_obj(optimizeProb(m.nt.o))
m.nt.o <- changeBounds(m.nt, "EX_gua(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
n.per.mol.RNA["G"] <- -lp_obj(optimizeProb(m.nt.o))
m.nt.o <- changeBounds(m.nt, "EX_ura(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
n.per.mol.RNA["U"] <- -lp_obj(optimizeProb(m.nt.o))

n.per.mol;


######################################################################################################
# DNA nucleotide production (see Suppl. Table 6 of Orth et al.)
######################################################################################################
# Exchange reactions for DNA nucleotides do not yet exist
nsp.per.mol.DNA <- c()  # use this to store the results

# Add exchange reactions for the four DNA nucleotides:
m.nt <- m2
m.nt <- addReact(m.nt, "datp_tr", c("datp[c]", "datp[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(m.nt, "dctp_tr", c("dctp[c]", "dctp[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(m.nt, "dgtp_tr", c("dgtp[c]", "dgtp[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(m.nt, "dttp_tr", c("dttp[c]", "dttp[e]"), c(-1, 1))   # transport reaction to extracellular
m.nt <- addReact(m.nt, "EX_datp(e)", c("datp[e]"), c(-1))   # transport reaction to NIL
m.nt <- addReact(m.nt, "EX_dgtp(e)", c("dgtp[e]"), c(-1))   # transport reaction to NIL
m.nt <- addReact(m.nt, "EX_dctp(e)", c("dctp[e]"), c(-1))   # transport reaction to NIL
m.nt <- addReact(m.nt, "EX_dttp(e)", c("dttp[e]"), c(-1))   # transport reaction to NIL

# how much of NSP is needed per Guanine (dgtp) etc. de novo synthesis?
m.nt.2 <- changeBounds(m.nt, "EX_datp(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
nsp.per.mol.DNA["A"] <- -lp_obj(optimizeProb(m.nt.2))
m.nt.2 <- changeBounds(m.nt, "EX_dctp(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
nsp.per.mol.DNA["C"] <- -lp_obj(optimizeProb(m.nt.2))
m.nt.2 <- changeBounds(m.nt, "EX_dgtp(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
nsp.per.mol.DNA["G"] <- -lp_obj(optimizeProb(m.nt.2))
m.nt.2 <- changeBounds(m.nt, "EX_dttp(e)", lb=1) # require at least 1 mol of this nucleotide to be produced
nsp.per.mol.DNA["T"] <- -lp_obj(optimizeProb(m.nt.2))

# results:
nsp.per.mol.DNA
# with oxygen:
# A C G T 
# 8 6 8 5   = RNA + 3 (i.e., order as in RNA case)
# ie., T < C < A = G
# i.e., T < A, C < G, T+A < G+C (just as in the energy case)
# without oxygen: 
# same!

# average difference between production costs:
avg.diff <- 0
for (i in 1:3) {
  for (j in (i+1):4) {
    avg.diff <- avg.diff + abs( nsp.per.mol.DNA[i] - nsp.per.mol.DNA[j] )
  }
}
avg.diff <- avg.diff / 6  
avg.diff  # 1.833333 = average difference between DNA molecules measured in NSP


######################################################################################################
# How much NSP is needed to make different parts of the biomass? (see Suppl. Table 6 of Orth et al.)
######################################################################################################

# NSP needed for 1 unit total biomass:
m.nt <- changeBounds(m2, "Ec_biomass_iJO1366_WT_53p95M", lb=1) # require at least 1 biomass unit to be produced
nsp.tot <- lp_obj(optimizeProb(m.nt)) # N+S+P = 11.68151
nsp.tot

# Biomass reaction restricted to RNA:
m.nt <- addReact(m2, "biomass.RNA", c("atp[c]", "ctp[c]", "gtp[c]", "utp[c]"), c(-0.169975, -0.129799, -0.209121, -0.140101))
# Produce RNA for one unit of biomass and estimate NSP required:
m.nt <- changeBounds(m.nt, "biomass.RNA", lb=1) # require at least 1 biomass unit to be produced
# minimize NSP uptake for this:
nsp.RNA <- lp_obj( optimizeProb(m.nt) ) # N+S+P=4.512067
nsp.RNA / nsp.tot   # 0.3862572 !!!!

# Biomass reaction restricted to DNA:
m.nt <- addReact(m2, "biomass.DNA", c("datp[c]", "dctp[c]", "dgtp[c]", "dttp[c]"), c(-0.024805, -0.025612, -0.025612, -0.024805))
# Produce DNA for one unit of biomass and estimate NSP required:
m.nt <- changeBounds(m.nt, "biomass.DNA", lb=1) # require at least 1 biomass unit to be produced
nsp.DNA <- lp_obj( optimizeProb(m.nt) ) # 0.681033 
nsp.DNA / nsp.tot  #  0.05830009

# Biomass reaction restricted to Amino Acids:
m.aa <- addReact(m2, "biomass.aa", c("ala_DASH_L[c]", "arg_DASH_L[c]", "asn_DASH_L[c]", "asp_DASH_L[c]", "cys_DASH_L[c]", 
                                     "gln_DASH_L[c]", "glu_DASH_L[c]", "gly[c]", "his_DASH_L[c]", "ile_DASH_L[c]", 
                                     "leu_DASH_L[c]", "lys_DASH_L[c]", "met_DASH_L[c]", "phe_DASH_L[c]", "pro_DASH_L[c]",
                                     "ser_DASH_L[c]", "thr_DASH_L[c]", "trp_DASH_L[c]", "tyr_DASH_L[c]", "val_DASH_L[c]"), 
                 c(-0.499149, -0.28742, -0.234232, -0.234232, -0.088988, 
                   -0.255712, -0.255712, -0.595297, -0.092056, -0.282306, 
                   -0.437778, -0.333448, -0.149336, -0.180021, -0.214798, 
                   -0.209684, -0.246506, -0.055234, -0.133993, -0.411184) )
# Produce AAs for one unit of biomass and estimate NSP required:
m.aa <- changeBounds(m.aa, "biomass.aa", lb=1) # require at least 1 biomass unit to be produced
# minimize NSP uptake for this:
nsp.AA <- lp_obj( optimizeProb(m.aa) ) # -7.360408
nsp.AA / nsp.tot   # 0.6130002

# relative investments:
nsp.info <- nsp.DNA + nsp.RNA + nsp.AA  # all information molecules (DNA, RNA, AA) together
nsp.info / nsp.tot   # 1.074648 of total NSP investment - more than 100% because of pi left over from ATP usage???
(nsp.DNA + nsp.RNA) / nsp.tot  # 0.4445573 of total NSP investment
nsp.RNA / nsp.DNA  # 6.625328 NSP investment into RNA vs. DNA
nsp.RNA / nsp.AA

######################################################################################################
# Average difference in AA production costs
######################################################################################################

# add (hypothetical) transport reactions:
m.aa <- m2
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

# find minimal NSP for each: 
nsp.per.aa <- c()
for (aa in aa.transports) {
  m.aa.o <- changeBounds(m.aa, aa, lb=1) # require at least 1 mol of this AA to be produced
  nsp.per.aa[aa] <- -lp_obj(optimizeProb(m.aa.o))
}
summary( nsp.per.aa) # median = 1, max=4

# average difference between production costs:
avg.diff <- 0
for (i in 1:19) {
  for (j in (i+1):20) {
    avg.diff <- avg.diff + abs( nsp.per.aa[i] - nsp.per.aa[j] )
  }
}
avg.diff <- avg.diff / 190 
avg.diff  # 0.7947368 = average difference between AA molecules measured in NSP

