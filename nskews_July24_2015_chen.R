## -- Authors: Wei-Hua Chen (chenwh550@gmail.com), Martin Lercher (...  ...);
## -- please send us email if you have any questions ...
## -- last modified May 18, 2014 --
## -- last modified May 23, 2014 --
##      added figure 4 and sup figure 2 --
## -- last modified June 30, 2014 --
##      added sup figure 3 (observed aa cost) --
## -- last modified : Sep 4, 2014 --
##    --> (alternatively) limit our analyses to 364 nonredundant bacterial species (data from Wu et al, 2012; PMID: 22230424) --
## -- last modified : Nov 3, 2014 --
##    --> added analyses for the impacts of life-styles on skews --
## -- last modified: Mar 11, 2015 -- MJL
##    --> changed formulae to reflect small N limit (following Bulmer)
## -- Last modified : March 20, 2015 -- MJL; current version of the model is 8.2 and the funciton script is 'MJL06_funcs.r'
## -- last modified : july 24, 2015;
## -- last modified : Sep 3, 2015;
##   added fig.s10, fig.s11; renamed s.11 to s13; removed original s.9

## --
# save.image(file = "working.RData");

## -- load data, libraries and functions --
rm(list=ls()); ## -- clean up --

## -- load necessary libraries; will install them if haven't been done before ... --
if (!require("ggplot2")) {
    install.packages("ggplot2", dependencies = TRUE);
    library(ggplot2);
}
## -- load the latest version of related R functions --
#setwd("~/1.Projects/Chen/ATskew/Nature2revision/data/nucleotideSkews2014Rproject") ;
source("data_and_other_R_scripts/MJL06_funcs.R");
source("data_and_other_R_scripts/multiplot.R"); ## multiplot for ggplot2 --


##############################################################################################################
### NOTE ###
### dMu.io.M6 : contains nucleotide skews attributed to replication (mutation)
### z.4s.M6 : contains nucleotide skews attributed to transcription (selection)
### zAA.4s.M6 : contains nucleotide skews attributed to translation at nonsynonymous sites (selection)
##############################################################################################################


## -- global parameters for plots  --
{
    color.transcription <- "#22AD5C";
    color.rep <- "#2E6ED5";
    color.ns <- "#F93E4F";
    color.2s <- "#BE58D3";
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"));

    ## -- others --
    boxplot.width <- 0.6;

    ## -- figure suffix --
    fig.suffix <- ""; ## by default

}

## ========================= load all data ========================
## >>>>>> Sep 11, 2014; also added a subset of 366 / 344 manually curated non-redundant bacterial genomes
## >>>>>> Sep 19, 2014; added ecosystem infor obtained from NCBI bioproject (ftp://ftp.ncbi.nlm.nih.gov/bioproject/);
##        NOTE the ecosystem / enviroment data could be imcomplete
{
    ## -- table S2 : all data --
    ## -- values in the following columns are unique: accession, NCBI_taxonomy_id, NCBI_bioProjectAccession, and NCBI_bioProjectID
    dat <- read.delim(file = "data_and_other_R_scripts/alldata_skews_by_genomic_regions_or_codon_sites_sep2014.txt", header = T, as.is = T, row.names = 1, sep = "\t");

    ## -- table S5 : skews by tai groups --
    tai5 <- read.delim(file = "data_and_other_R_scripts/skews_by_tai.txt", header = T, as.is = T, sep  = "\t");

    ## -- table S7 : the average gc contents and cost per coded amino acids by tAI groups --
    costs.quantiles <- read.delim( file = "data_and_other_R_scripts/gc_and_aacost_by_tai.txt", header = T, as.is = T, sep = "\t" );

    ## -- the random data for Figure 3 --
    sim2 <- read.delim(file = "data_and_other_R_scripts/random_data.txt", header = T, as.is = T, sep = "\t");

    ## -- taboe S6 : the tradeoff --
    tradeoff <- read.table(file = "data_and_other_R_scripts/the_tradeoff.txt", header = T, as.is = T, sep = "\t");

    ## -- load extra data : observed AA costs for selected bacterial genomes --
    aaCostsObserved <- read.delim(file = "data_and_other_R_scripts/aaCosts_observed_june2014.txt", header = T, sep = "\t", as.is = T, row.names = 1);

    ## -- data.frame contains 364 nonredundant bacterial species --
    nr364 <- read.delim(file = "data_and_other_R_scripts/nonredundant_364_eubacteria_Wu_et_al_2012_PMID22230424.txt", sep = "\t", header = T, as.is = T );

    ## -- intracellular organisms; Table S13 ...
    intra <- read.delim(file = "data_and_other_R_scripts/intracellular_nov3_2014_TableS13.txt", header = T, as.is = T, sep = "\t");
    dat[, "intra"] <- "N";

    dat[ row.names(dat) %in% intra$Genbank.Accession , "intra"] <- "Y";
    table(dat$intra);

}

## -- !!!! important !!!!
## -- July 24, 2015 ...
## -- section: how to generate extended data figures 3~5
## --          1. unmask the codes in the brackets --
## --          2. run the codes to plot figures 1,2,3
## --          3. now: figure1 == figure s3, figure2 == figure s4, figure 3 == figure s5  !!!!
## -- sep 4, 2014--
## -- run the codes of this section if you want to limit the analyses to the 364 nonredundant bacterial species --
## -- by default, the following section is masked and will not run --
{
    ## --
    #dat <- dat[ rownames(dat) %in% nr364$Acc, ];
    #tai5 <- tai5[ tai5$acc %in% nr364$Acc,  ];
    #fig.suffix <- "_nr364_S";

}

## -- calculate dMus --
## ================= estimate replication from interopnic regions, 4s, 2s and ns sites ==================
## -- structure of resulting lists --
## -- list[["AT"]] = dMu for AT skews, [["GC"]] = dMu for GC skews,
## --    the following are only availabe if type != "i"
## --     [["AT.fac"]], [["GC.fac"]] = factor for AT/GC skews in comparison with interoperonic regions
## --     [["raw.GC.le"]], [["raw.GC.la"]]  = raw gc skews at leading and lagging strands respectively --
## --     [["raw.AT.le"]], [["raw.AT.la"]]  = raw at skews at leading and lagging strands respectively --
{

    ## -- calculate (Model08) --
    dMu.io.M6 <- dMuCalc("i", model="Model08");
    dMu.4s.M6 <- dMuCalc("4", model="Model08");
    dMu.2s.M6 <- dMuCalc("2", model="Model08");
    dMu.1s.M6 <- dMuCalc("1", model="Model08");

}

## ================= calculate relative strength of selection ... ==================
## -- structure of resulting lists --
## -- list [["GC.le"]], [["GC.la"]] -- z score of GC skews at leading and lagging strands respectively --
## --      [["AT.le"]], [["AT.la"]] -- z score of AT skews at leading and lagging strands respectively --
{
    ## -- calculate z --
    z.4s.M6 <- zCalc( dMu.4s.M6, model="Model08", dMuIO= dMu.io.M6 );
    z.2s.M6 <- zCalc( dMu.2s.M6, model="Model08", dMuIO= dMu.io.M6 );
    z.1s.M6 <- zCalc( dMu.1s.M6, model="Model08", dMuIO= dMu.io.M6 );

    ### correlations between leading and lagging strand:
    #cor.test(z.4s[["AT.le"]], z.4s[["AT.la"]])
    cor.test(z.4s.M6[["AT.le"]], z.4s.M6[["AT.la"]]);   # 0.771 $$$
    # Figure 1C: $$$
    #plot(z.4s.M6[["AT.le"]], z.4s.M6[["AT.la"]], cex=0.3)
    #abline(0, 1, col="grey")

    #cor.test(z.4s[["GC.le"]], z.4s[["GC.la"]]);   # 0.844
    cor.test(z.4s.M6[["GC.le"]], z.4s.M6[["GC.la"]]);    # 0.842 $$$
    # Figure 1D: $$$
    #plot(z.4s.M6[["GC.le"]], z.4s.M6[["GC.la"]], cex=0.3)
    #abline(0, 1, col="grey");
}

## =============  calculate factors for transcription and zAA ==============
## -- structure of resulting lists --
## -- [["GC.le"]], [["GC.la"]] -- zAA for GC at leading and lagging strand respectively --
##    [["AT.le"]], [["AT.la"]] -- zAA for AT at leading and lagging strand respectively --
##    [["fac.AT"]], [["fac.GC"]] -- estimated factor x of transcription (z4s) for AT and GC respectively --
##    [["fac.raw"]] -- an dataframe shows how fac.AT and fac.GC are determined --
{

    zAA.1s.M6 <- zAAcalcM08( dMu.xs= dMu.1s.M6, dMu.4s= dMu.4s.M6, dMu.io= dMu.io.M6, z.4s = z.4s.M6 );
    zAA.2s.M6 <- zAAcalcM08( dMu.xs= dMu.2s.M6, dMu.4s= dMu.4s.M6, dMu.io= dMu.io.M6, z.4s = z.4s.M6 );

}

## ==== process skews_by_tai file (table S2) ====
{
    ## -- do stats --
    ## --
    tai5.stats <- taiStats( data = tai5, dat= dat, dMuIO=dMu.io.M6, z.4s.M6,  model = "Model08" ); ## --

    ## -- get data --
    tmp5.combined <- getDataFromTaiStats1( tai.stats= tai5.stats, strand="combined" );

    df.all <- cbind(
        dat[, c("genomicGC", "genomicLen")],
        data.frame( dMuioAT = dMu.io.M6[["AT"]], dMuioGC = dMu.io.M6[["GC"]],
                    dMu4sAT = dMu.4s.M6[["AT"]], dMu4sGC = dMu.4s.M6[["GC"]],
                    zRNA.AT.le = z.4s.M6[["AT.le"]], zRNA.AT.la = z.4s.M6[["AT.la"]],
                    zRNA.GC.le = z.4s.M6[["GC.le"]], zRNA.GC.la = z.4s.M6[["GC.la"]],
                    zAA.AT.le = zAA.1s.M6[["AT.le"]], zAA.AT.la = zAA.1s.M6[["AT.la"]],
                    zAA.GC.le = zAA.1s.M6[["GC.le"]], zAA.GC.la = zAA.1s.M6[["GC.la"]] )
        );

    #  write.table(df.all, file = "skews_nov5_2014.txt", row.names = T, col.names = T, sep = "\t", quote = F);

}


## ===================================
## === use ggplot2 to replot everything ===
## -- main figures; figures 1,2,3,4, --
## -- and sup figure 2 --
{
    ## -- fig 1, density plot --
    {
        df.tais <- rbind(
            ## -- transcription / dMu.io --
            data.frame( cat = "rep.at", group = "AT", dat = dMu.io.M6[["AT"]] ),
            data.frame( cat = "rep.gc", group = "GC", dat = dMu.io.M6[["GC"]] ),

            ## -- z.4s --
            data.frame( cat = "trans.at", group = "AT", dat = (z.4s.M6[["AT.le"]] + z.4s.M6[["AT.la"]])/2),
            data.frame( cat = "trans.gc", group = "GC", dat = (z.4s.M6[["GC.le"]] + z.4s.M6[["GC.la"]])/2 ),

            ## -- zAA.1s --
            data.frame( cat = "AA1s.at", group = "AT", dat = ( zAA.1s.M6[["AT.le"]] + zAA.1s.M6[["AT.la"]] )/2 ),
            data.frame( cat = "AA1s.gc", group = "GC", dat = ( zAA.1s.M6[["GC.le"]] + zAA.1s.M6[["GC.la"]] )/2 )
        );

        ## -- overall characteristics: leading; some photoshopping is needed to make the plot publishable ... 汗!!  --
        pdf(file = paste("fig.1", fig.suffix, ".pdf", sep  = ""), width = 5, height = 5);
        ggplot(df.tais, aes( x=dat, group = cat )) + geom_density( aes(colour = factor(cat), linetype = group )  ) +
            scale_colour_manual( values = c( color.rep, color.rep, color.transcription, color.transcription, color.ns, color.ns ),
                                 name = "Skews attributed to\nvarious factors",
                                 labels = c( expression(paste("mutational bias ", italic(mu^{io}), " for AT attributed to replication")),
                                             expression(paste("mutational bias ", italic(mu^{io}), " for GC attributed to replication")),
                                             ## -- transcription --
                                             expression(paste("selection coefficient ", italic(S[RNA]), " for AT attributed to transcription")),
                                             expression(paste("selection coefficient ", italic(S[RNA]), " for GC attributed to transcription")),
                                             ## -- translation --
                                             expression(paste("selection coefficient ", italic(S[AA]), " for AT attributed to translation")),
                                             expression(paste("selection coefficient ", italic(S[AA]), " for GC attributed to translation")) )
            ) +
            theme(legend.justification=c(1,1), legend.position=c(1,1), legend.background = element_rect(colour = NA, fill = NA), legend.text = element_text( size = 5 )  ) +
            xlim( c( -1, 2.5 ) ) +
            #theme( legend.position = "none" ) +
            xlab( "" ) + ylab( "density" );
        dev.off();

    }
    ## -- fig 2. skews attributed to three fundamental celular processes --
    ## -- in total six panels --
    {
        ## -- six panels --
        ## panel A, rep AT --
        type <- "AT"; dMu <- dMu.4s.M6;
        fac <- dMu[[  paste(type, ".fac", sep = "")  ]];
        tmp <- data.frame( x = dMu.io.M6[[ type ]], y = dMu[[ type ]], gc = dat$genomicGC );
        tmp <- tmp[ complete.cases(tmp), ];
        t <- bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( tmp$x, tmp$y )$estimate , digits=3) , sep = "")) ) );
        t2 <- parse(text=paste( "slope== ", round( fac , digits=3) , sep = ""));
        p.rep.at <-
            ggplot(tmp, aes( x, y ) ) +
            geom_point( aes(colour = gc ), size = 0.8 ) +
            scale_colour_gradientn( colours = jet.colors(7), name = "genomic GC%" ) +
            #" for AT attributed to replication \nderived from interoperonic regions"
            xlab( expression(paste(italic(mu^{io})," for AT derived from interoperonic regions")) ) +
            ylab( expression(paste(italic(mu^{io})," for AT derived from 4s sites"))  ) +
            ggtitle( t ) +
            geom_abline( slope = fac, linetype = 2, colour = "red", size = 0.5 ) +
            annotate("text", x = max( tmp[,"x"] ), y = max(tmp[,"y"]), label = paste("slope=", round( dMu.4s.M6[["AT.fac"]], 3), sep = "" ), hjust = 1.4, yjust = 2);

        ## panel B, rep GC --
        type <- "GC"; dMu <- dMu.4s.M6;
        fac <- dMu[[  paste(type, ".fac", sep = "")  ]];
        tmp <- data.frame( x = dMu.io.M6[[ type ]], y = dMu[[ type ]], gc = dat$genomicGC );
        tmp <- tmp[ complete.cases(tmp), ];
        t <- bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( tmp$x, tmp$y )$estimate , digits=3) , sep = "")) ) );
        t2 <- parse(text=paste( "slope== ", round( fac , digits=3) , sep = ""));
        p.rep.gc <-
            ggplot(tmp, aes( x, y ) ) +
            geom_point( aes(colour = gc ), size = 0.8 ) +
            scale_colour_gradientn( colours = jet.colors(7), name = "genomic GC%" ) +
            xlab( expression(paste(italic(mu^{io})," for GC derived from interoperonic regions")) )  +
            ylab( expression(paste(italic(mu^{io})," for GC derived from 4s sites"))  ) +
            ggtitle( t ) +
            geom_abline( slope = fac, linetype = 2, colour = "red", size = 0.5 ) +
            annotate("text", x = max( tmp[,"x"] ), y = max(tmp[,"y"]), label = paste("slope=", round( dMu.4s.M6[["GC.fac"]], 3), sep = "" ), hjust = 1.4, yjust = 2);

        ## panel C&D, trans AT --
        t <- bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( z.4s.M6[["AT.le"]], z.4s.M6[["AT.la"]] )$estimate , digits=3), sep = "" ))) );
        pz.at <-
            ggplot( data.frame( x = z.4s.M6[["AT.le"]], y = z.4s.M6[["AT.la"]], gc = dat$genomicGC ), aes( x,y ) ) +
            geom_point( aes(colour = gc ), size = 0.8 ) +
            scale_colour_gradientn( colours = jet.colors(7), name = "genomic GC%" ) +
            xlab( expression(paste(italic(S[RNA])," for AT, leading strand")) )  +
            ylab( expression(paste(italic(S[RNA])," for AT, lagging strand"))  ) +
            ggtitle( t ) +
            geom_abline( slope = 1, linetype = 2, colour = "darkgray", size = 0.5 );

        t <- bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( z.4s.M6[["GC.le"]], z.4s.M6[["GC.la"]] )$estimate , digits=3), sep = "" ))) );
        pz.gc <-
            ggplot( data.frame( x = z.4s.M6[["GC.le"]], y = z.4s.M6[["GC.la"]], gc = dat$genomicGC ), aes( x, y ) ) +
            geom_point( aes(colour = gc ), size = 0.8 ) +
            scale_colour_gradientn( colours = jet.colors(7), name = "genomic GC%" ) +
            xlab( expression(paste(italic(S[RNA])," for GC, leading strand")) )  +
            ylab( expression(paste(italic(S[RNA])," for GC, lagging strand"))  ) +
            ggtitle( t ) +
            geom_abline( slope = 1, linetype = 2, colour = "darkgray", size = 0.5 );

        ## -- panel E&F, AA at translation, AT and GC --
        t <- bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( zAA.1s.M6[["AT.le"]], zAA.1s.M6[["AT.la"]] )$estimate , digits=3), sep = "" ))) );
        pAAns.at <-
            ggplot( data.frame( x = zAA.1s.M6[["AT.le"]], y = zAA.1s.M6[["AT.la"]], gc = dat$genomicGC ), aes( x,y ) ) +
            geom_point( aes(colour = gc ), size = 0.8 ) +
            scale_colour_gradientn( colours = jet.colors(7), name = "genomic GC%" ) +
            xlab( expression(paste(italic(S[AA])," for AT, leading strand")) )  +
            ylab( expression(paste(italic(S[AA])," for AT, lagging strand"))  ) +
            ggtitle( t ) +
            geom_abline( slope = 1, linetype = 2, colour = "darkgray", size = 0.5 );

        t <- bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( zAA.1s.M6[["GC.le"]], zAA.1s.M6[["GC.la"]] )$estimate , digits=3), sep = "" ))) );
        pAAns.gc <-
            ggplot( data.frame( x = zAA.1s.M6[["GC.le"]], y = zAA.1s.M6[["GC.la"]], gc = dat$genomicGC ), aes( x, y ) ) +
            geom_point( aes(colour = gc ), size = 0.8 ) +
            scale_colour_gradientn( colours = jet.colors(7), name = "genomic GC%" ) +
            xlab( expression(paste(italic(S[AA])," for GC, leading strand")) )  +
            ylab( expression(paste(italic(S[AA])," for GC, lagging strand"))  ) +
            ggtitle( t ) +
            geom_abline( slope = 1, linetype = 2, colour = "darkgray", size = 0.5 );

        ## -- plot --
        pdf(file = paste("fig.2", fig.suffix, ".pdf", sep = ""), width = 13, height = 15);
        multiplot(p.rep.at, pz.at, pAAns.at, p.rep.gc, pz.gc,  pAAns.gc, cols= 2);
        dev.off();
    }

    ## -- figure 3 / fig 3. skews by tai groups --
    ## -- in total four panels --
    {
        ## -- panel A --
        ## -- assemble data; get data ready --
        percentages <- paste( round( c(
            sum(tmp5.combined$trans.at[ tmp5.combined$tai == 5 ] > tmp5.combined$trans.at[ tmp5.combined$tai == 1 ] ),
            sum(tmp5.combined$trans.gc[ tmp5.combined$tai == 5 ] > tmp5.combined$trans.gc[ tmp5.combined$tai == 1 ] ),

            sum(tmp5.combined$zAA2.at[ tmp5.combined$tai == 5 ] > tmp5.combined$zAA2.at[ tmp5.combined$tai == 1 ] ),
            sum(tmp5.combined$zAA2.gc[ tmp5.combined$tai == 5 ] > tmp5.combined$zAA2.gc[ tmp5.combined$tai == 1 ] ),

            sum(tmp5.combined$zAA1.at[ tmp5.combined$tai == 5 ] > tmp5.combined$zAA1.at[ tmp5.combined$tai == 1 ] ),
            sum(tmp5.combined$zAA1.gc[ tmp5.combined$tai == 5 ] > tmp5.combined$zAA1.gc[ tmp5.combined$tai == 1 ] )
        ) / dim(dat)[1] * 100, digits = 1), "%", sep = "");

        quantiles = 5;
        labels <- paste( "Q", 1:quantiles, sep = "" );
        labels.with.Random <- c("R*", paste( "Q", 1:quantiles, sep = "" ) );


        lty.gc <- 2;
        a <-
            ggplot( tmp5.combined, aes( factor(tai), trans.at ) ) + geom_boxplot( fill = color.transcription, linetype = 1 ,outlier.size = 1, width = boxplot.width) +
            xlab( "tAI group" ) + ylab( expression( paste( italic(S[RNA]) ) ) ) +
            theme( axis.text = element_text(size = 10) , plot.title = element_text(size = rel(0.8), colour = "darkblue"))+
            ggtitle( paste("AT(", percentages[1], ")", sep = "") ) + ylim( range( c( tmp5.combined$trans.at, tmp5.combined$trans.gc ) ) ) +
            scale_x_discrete(breaks=c(1:quantiles), labels=labels);

        b <- ggplot( tmp5.combined, aes( factor(tai), trans.gc ) ) + geom_boxplot( fill = color.transcription, linetype = lty.gc ,outlier.size = 1, width = boxplot.width) +
            xlab( "tAI group" ) + ylab( expression( paste( italic(S[RNA]) ) ) ) +
            theme( axis.text = element_text(size = 10) , plot.title = element_text(size = rel(0.8), colour = "darkblue"))+
            ggtitle( paste("GC(", percentages[2], ")", sep = "") ) + ylim( range( c( tmp5.combined$trans.at, tmp5.combined$trans.gc ) ) ) +
            scale_x_discrete(breaks=c(1:quantiles), labels=labels);

        ## -- panel B --
        data.aa1 <- rbind( tmp5.combined[ , c("tai", "acc", "zAA1.at", "zAA1.gc") ], data.frame( tai = 0, acc = "acc", zAA1.at = sim2$skewAT1, zAA1.gc = sim2$skewGC1 ) );

        e <- ggplot( data.aa1, aes( factor(tai), zAA1.at ) ) + geom_boxplot( fill = c( "white", rep(color.ns, 5)), linetype = 1 ,outlier.size = 1, width = boxplot.width) +
            xlab( "tAI group" ) + ylab( expression( paste( italic(S[AA]) ) ) ) +
            theme( axis.text = element_text(size = 10), plot.title = element_text(size = rel(0.8), colour = "darkblue") )+
            ggtitle( paste("AT(", percentages[5], ")", sep = "") )+ ylim( range( c( data.aa1$zAA1.at, data.aa1$zAA1.gc ) ) ) +
            scale_x_discrete(breaks=c(0:quantiles), labels=labels.with.Random);

        f <- ggplot( data.aa1, aes( factor(tai), zAA1.gc ) ) + geom_boxplot( fill = c( "white", rep(color.ns, 5)), linetype = lty.gc ,outlier.size = 1, width = boxplot.width) +
            xlab( "tAI group" ) + ylab( expression( paste( italic(S[AA]) ) ) ) +
            theme( axis.text = element_text(size = 10), plot.title = element_text(size = rel(0.8), colour = "darkblue") )+
            ggtitle( paste("GC(", percentages[6], ")", sep = "") )+ ylim( range( c( data.aa1$zAA1.at, data.aa1$zAA1.gc ) ) ) +
            scale_x_discrete(breaks=c(0:quantiles), labels=labels.with.Random);

        ## -- panel c : GC content --
        gc <- ggplot( costs.quantiles, aes( factor(quantile), (GCLe + GCLa) / 2 * 100) ) +
            geom_boxplot(outlier.size = 1, width = boxplot.width, fill = "#6baed6", colour = "#3182bd" ) +
            xlab( "tAI group" ) + theme( axis.text = element_text(size = 10), plot.title = element_text(size = rel(0.8), colour = "darkblue") )+
            ylab( "GC%" ) + ggtitle( "GC%" );

        aa <- ggplot( costs.quantiles, aes( factor(quantile), (costLe + costLa )/2 ) ) +
            geom_boxplot(outlier.size = 1, width = boxplot.width, fill = "#74c476", colour = "#31a354" ) +
            xlab( "tAI group" ) + theme( axis.text = element_text(size = 10), plot.title = element_text(size = rel(0.8), colour = "darkblue") )+
            ylab( "average costs (~P) in amino acid syntheses" ) + ggtitle("AA costs");


        ## -- PS is needed to make the resulting figure publishable...
        pdf(file = paste("fig.3", fig.suffix, ".pdf", sep = ""), height = 5, width = 13);
        multiplot( a, b, e, f, gc, aa, cols= 6 );
        dev.off();

        ## ===========================================================
        ## -- the following is for testing porpuses only --
        data.aa2 <- rbind( tmp5.combined[ , c("tai", "acc", "zAA2.at", "zAA2.gc") ], data.frame( tai = 0, acc = "acc", zAA2.at = sim2$skewAT2, zAA2.gc = sim2$skewGC2 ) );
        c <- ggplot( data.aa2, aes( factor(tai), zAA2.at ) ) + geom_boxplot( fill = c( "white", rep(color.2s, 5)), linetype = 1 ,outlier.size = 1, width = boxplot.width) +
            xlab( "tAI group" ) + ylab( "AT skews attributed to AA selection at ns sites" ) +
            theme( axis.text = element_text(size = 10), plot.title = element_text(size = rel(0.8), colour = "darkblue") )+
            ggtitle( paste("AT(", percentages[3], ")", sep = "") )+ ylim( range( c( data.aa2$zAA2.at, data.aa2$zAA2.gc ) ) ) +
            scale_x_discrete(breaks=c(0:quantiles), labels=labels.with.Random);

        d <- ggplot( data.aa2, aes( factor(tai), zAA2.gc ) ) + geom_boxplot( fill = c( "white", rep(color.2s, 5)), linetype = lty.gc ,outlier.size = 1, width = boxplot.width) +
            xlab( "tAI group" ) + ylab( "GC skews attributed to AA selection at ns sites" ) +
            theme( axis.text = element_text(size = 10), plot.title = element_text(size = rel(0.8), colour = "darkblue") )+
            ggtitle( paste("GC(", percentages[4], ")", sep = "") )+ ylim( range( c( data.aa2$zAA2.at, data.aa2$zAA2.gc ) ) ) +
            scale_x_discrete(breaks=c(0:quantiles), labels=labels.with.Random);

        pdf(file= paste("fig.zaa2", fig.suffix, ".pdf", sep = ""), height = 5, width = 5);
        multiplot(c,d, cols=2);
        dev.off();
    }

    ## -- figure 4, the tradeoff --
    {
        ## -- ## --
        pdf(file = "fig.4.pdf", height = 5, width = 6);
        ggplot( data = tradeoff, aes( x = overallSkews, y = aaAvgCost, group = GC, colour = GC ) ) +
            geom_line(  ) +
            scale_colour_gradientn( colours = jet.colors(7), name = "coding GC%" ) +
            xlab("overall AT / GC skews") + ylab("mean cost per amino acid");
        dev.off();
    }

    ## -- sup figures --
    ## -- last modified : july 24, 2015 ...
    ## -- last modified : sep 3, 2015 ...
    ## -- fig.s9. two panles; both are similar to figure 4,
    ##             but in one GC skews increase but AT skews remain unchanged,
    ##             in the other AT skews increase but GC skews remain unchanged
    {
        ## -- load data --
        tradeoffs1 <- read.delim(file = "data_and_other_R_scripts/sim5_lite_ATonly_noGCskews_march14_2014.txt", header = T, as.is = T);
        tradeoffs2 <- read.delim(file = "data_and_other_R_scripts/sim5_lite_GConly_noATskews_march14_2014.txt", header = T, as.is = T);

        ## -- add a new column: skew --
        x.lvl <- seq( 10, 90, 5 );
        y <- seq( -0.8, 0.8, 0.1 );
        tradeoffs1[, "skew"] <- y[ match( tradeoffs1$overallATGCskew, x.lvl ) ];
        tradeoffs2[, "skew"] <- y[ match( tradeoffs2$overallATGCskew, x.lvl ) ];

        ## -- plot --
        a <- ggplot( data = tradeoffs1, aes( x = skew, y = aaAvgCost, group = GC, colour = GC ) ) +
            geom_line(  ) + ggtitle("AT skews only") +
            scale_colour_gradientn( colours = jet.colors(7), name = "coding GC%" ) +
            xlab("overall AT skews") + ylab("mean cost per amino acid");

        b <- ggplot( data = tradeoffs2, aes( x = skew, y = aaAvgCost, group = GC, colour = GC ) ) +
            geom_line(  ) + ggtitle("GC skews only") +
            scale_colour_gradientn( colours = jet.colors(7), name = "coding GC%" ) +
            xlab("overall GC skews") + ylab("mean cost per amino acid");

        ## -- save to output file --
        pdf(file = "fig.s9_tradeoffs_AT_or_GC_skews_only.pdf", height = 5, width = 12);
        multiplot(a,b, cols=2);
        dev.off();
    }

    ## --------------------------------------
    ## --- fig S10 / fig.s10 ...
    ## -- last modified : Sep 3, 2015 ; new Figure S10: boxplot of raw skews as a function of tais --
    {
        library(reshape2);

        pRaw.ns.at <-
            ggplot( melt( as.data.frame(data.matrix(tai5[, c("quantile", "ATle1", "ATla1")])), id.vars = "quantile" ),
                    aes( x = factor(quantile), y = value, fill = variable ) ) +
            geom_boxplot( ) + xlab( "tAI group" ) + ylab( expression( paste( "raw skew (", italic(gamma), ")" ) ) ) +
            scale_fill_discrete(name="Strand",
                                breaks=c("ATle1", "ATla1"),
                                labels=c("leading", "lagging")) +
            ggtitle( "AT at ns site" )

        pRaw.ns.gc <-
            ggplot( melt( as.data.frame(data.matrix(tai5[, c("quantile", "GCle1", "GCla1")])), id.vars = "quantile" ),
                    aes( x = factor(quantile), y = value, fill = variable ) ) +
            geom_boxplot( ) + xlab( "tAI group" ) + ylab( expression( paste( "raw skew (", italic(gamma), ")" ) ) ) +
            scale_fill_discrete(name="Strand",
                                breaks=c("GCle1", "GCla1"),
                                labels=c("leading", "lagging")) +
            ggtitle( "GC at ns site" )

        pRaw.4s.at <-
            ggplot( melt( as.data.frame(data.matrix(tai5[, c("quantile", "ATle4", "ATla4")])), id.vars = "quantile" ),
                    aes( x = factor(quantile), y = value, fill = variable ) ) +
            geom_boxplot( ) + xlab( "tAI group" ) + ylab( expression( paste( "raw skew (", italic(gamma), ")" ) ) ) +
            scale_fill_discrete(name="Strand",
                                breaks=c("ATle4", "ATla4"),
                                labels=c("leading", "lagging")) +
            ggtitle( "AT at 4s site" )

        pRaw.4s.gc <-
            ggplot( melt( as.data.frame(data.matrix(tai5[, c("quantile", "GCle4", "GCla4")])), id.vars = "quantile" ),
                    aes( x = factor(quantile), y = value, fill = variable ) ) +
            geom_boxplot( ) + xlab( "tAI group" ) + ylab( expression( paste( "raw skew (", italic(gamma), ")" ) ) ) +
            scale_fill_discrete(name="Strand",
                                breaks=c("GCle4", "GCla4"),
                                labels=c("leading", "lagging")) +
            ggtitle( "GC at 4s site" )


        pdf(file = "fig.s10_rawSkews_ns_4s_stratified_by_tai.pdf", height = 10, width = 8);
        multiplot(pRaw.ns.at, pRaw.4s.at, pRaw.ns.gc, pRaw.4s.gc,  cols= 2);
        dev.off();

    }

    ## --------------------------------------
    ## -- last modified : Sep 3, 2015;
    ## -- fig.s11 (new since sep 3, 2015) : S_rc_and_nrc_flanking_nucleotides.ver.sep2015 --
    {
        ## - load data --
        datRCNRC <- read.delim(file = "data_and_other_R_scripts/skews_4s_io_for_rc_and_nonrc_flanking_nucs_july_27_2015.txt", header = T, as.is = T, sep = "\t")
        ## -- NOTE: here RC == triplets with reverse-complementary flanking nucleotides --
        ##  NRC == triplets with non-reverse-complemnetary flanking nucleotides--

        atY <- datRCNRC[ datRCNRC$cat == "AT" & datRCNRC$rcfl == "Y", ];
        sRC.at <- log( ( (1 + atY$le4s) / (1 - atY$le4s) ) * ( (1 + atY$la4s) / (1-atY$la4s) ) )/2;

        atN <- datRCNRC[ datRCNRC$cat == "AT" & datRCNRC$rcfl == "N", ];
        sNRC.at <- log( ( (1 + atN$le4s) / (1 - atN$le4s) ) * ( (1 + atN$la4s) / (1-atN$la4s) ) )/2;

        atA <- datRCNRC[ datRCNRC$cat == "AT" & datRCNRC$rcfl == "A", ];
        sA.at <- log( ( (1 + atA$le4s) / (1 - atA$le4s) ) * ( (1 + atA$la4s) / (1-atA$la4s) ) )/2;

        ## -- cor --
        cor.test(sA.at, sRC.at);
        cor.test(sA.at, sNRC.at);
        cor.test(sRC.at, sNRC.at);


        ## =========================================
        ## -- GC --
        gcY <- datRCNRC[ datRCNRC$cat == "GC" & datRCNRC$rcfl == "Y", ];
        sRC.gc <- log( ( (1 + gcY$le4s) / (1 - gcY$le4s) ) * ( (1 + gcY$la4s) / (1-gcY$la4s) ) )/2;

        gcN <- datRCNRC[ datRCNRC$cat == "GC" & datRCNRC$rcfl == "N", ];
        sNRC.gc <- log( ( (1 + gcN$le4s) / (1 - gcN$le4s) ) * ( (1 + gcN$la4s) / (1-gcN$la4s) ) )/2;

        gcA <- datRCNRC[ datRCNRC$cat == "GC" & datRCNRC$rcfl == "A", ];
        sA.gc <- log( ( (1 + gcA$le4s) / (1 - gcA$le4s) ) * ( (1 + gcA$la4s) / (1-gcA$la4s) ) )/2;

        ## -- cor --
        cor.test(sA.gc, sRC.gc);
        cor.test(sA.gc, sNRC.gc);
        cor.test(sRC.gc, sNRC.gc);


        ## -- plot a figure as supplementary figure ... --
        df.gc <- rbind(
            data.frame( cat = "All", datRCNRC = sA.gc ),
            data.frame( cat = "rc", datRCNRC = sRC.gc),
            data.frame( cat = "nrc", datRCNRC = sNRC.gc)
        );

        df.at <- rbind(
            data.frame( cat = "All", datRCNRC = sA.at ),
            data.frame( cat = "rc", datRCNRC = sRC.at),
            data.frame( cat = "nrc", datRCNRC = sNRC.at)
        );

        p.dens.gc <-
            ggplot(df.gc, aes( x=datRCNRC, group = cat )) + geom_density( aes(colour = factor(cat), linetype = cat )  ) +
            scale_colour_manual( values = c( "black", "green", "red" ) ) +
            geom_vline( xintercept = 0, linetype = 1, size = 0.5 ) +
            xlab( expression(paste(italic(S[RNA])," for GC")) )  +
            theme(legend.position="none");

        p.dens.at <-
            ggplot(df.at, aes( x=datRCNRC, group = cat )) + geom_density( aes(colour = factor(cat), linetype = cat )  ) +
            scale_colour_manual( values = c( "black", "green", "red" ) ) +
            geom_vline( xintercept = 0, linetype = 1,  size = 0.5 )+
            xlab( expression(paste(italic(S[RNA])," for AT")) )  +
            theme(legend.position="none");

        p.points.gc <-
            ggplot( data.frame( x = sA.gc, y = sRC.gc ), aes(x,y) ) + geom_point() +
            geom_abline(slope = 1, linetype = 2, colour = "darkgray", size = 0.5) +
            ggtitle( bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( sA.gc, sRC.gc )$estimate , digits=3), sep = "" ))) ) ) +
            xlab( expression(paste(italic(S[RNA])," for GC, all")) )  +
            ylab( expression(paste(italic(S[RNA])," for GC, sites with reverse-complementary flanking nucleotides"))  )


        p.points.at <-
            ggplot( data.frame( x = sA.at, y = sRC.at ), aes(x,y) ) + geom_point() +
            geom_abline(slope = 1, linetype = 2, colour = "darkgray", size = 0.5) +
            ggtitle( bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( sA.at, sRC.at )$estimate , digits=3), sep = "" ))) ) ) +
            xlab( expression(paste(italic(S[RNA])," for AT, all")) )  +
            ylab( expression(paste(italic(S[RNA])," for AT, sites with reverse-complementary flanking nucleotides"))  )

#         p.points2.gc <-
#             ggplot( data.frame( x = sA.gc, y = sNRC.gc ), aes(x,y) ) + geom_point() +
#             geom_abline(slope = 1, linetype = 2, colour = "darkgray", size = 0.5) +
#             ggtitle( bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( sA.gc, sNRC.gc )$estimate , digits=3), sep = "" ))) ) ) +
#             xlab( expression(paste(italic(S[RNA])," for GC, all")) )  +
#             ylab( expression(paste(italic(S[RNA])," for GC, sites with other flanking nucleotides"))  )
#
#
#         p.points2.at <-
#             ggplot( data.frame( x = sA.at, y = sNRC.at ), aes(x,y) ) + geom_point() +
#             geom_abline(slope = 1, linetype = 2, colour = "darkgray", size = 0.5) +
#             ggtitle( bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( sA.at, sNRC.at )$estimate , digits=3), sep = "" ))) ) ) +
#             xlab( expression(paste(italic(S[RNA])," for AT, all")) )  +
#             ylab( expression(paste(italic(S[RNA])," for AT, sites with other flanking nucleotides"))  )

        pdf(file = "fig.s11_rc_and_nrc_flanking_nucleotides_sep2015.pdf", height = 10, width = 10.2);
        multiplot(p.dens.at, p.points.at, p.dens.gc, p.points.gc, cols= 2);
        dev.off();
    }

    ## -- last modified : july 24, 2015 ...
    ## -- last modified : sep 3, 2015 ...
    ## -- fig.s13; similar to fig.4, but use alternative source of aa cost i.e. aerobic heterotrophe  --
    {
        ## -- ## --
        pdf(file = "fig.s13_tradefoff_alternative_source_of_aaCosts.pdf", height = 5, width = 6);
        ggplot( data = tradeoff, aes( x = overallSkews, y = aaAvgCost_Aer_Het, group = GC, colour = GC ) ) +
            geom_line(  ) +
            scale_colour_gradientn( colours = jet.colors(7), name = "coding GC%" ) +
            xlab("overall AT / GC skews") + ylab("mean cost per amino acid");
        dev.off();
    }


    ## -- figs. s3~5
    ## -- see section: '## -- section: how to generate extended data figures 3~5' for the instructions for generating figures 3~5
}

## ===================   Sup figure 2, simulation ; july 24, 2015 ====================================
{
    ## -- Figure S2, this is to illustrate the expected dependence of γ on Ne, s, and μ.
    ## -- Y-axis shows the skews (γ values), and X-axis shows the S values (here S = 2Nes).
    mu <- seq(-1, 1, 0.01);
    mu <- mu[ 2:(length(mu)-1) ]; ## it cannot be - 1
    #s <- 10 ^ seq(-2, 1, 0.001 );
    s <- seq(0.01, 10, 0.1);
    simu <- data.frame( mu = rep(mu, length(s)), S = rep(s, each = length( mu )) );

    ## -- simulation: mu ranges from -1 to 1, S ranges from 10-3 to 10^3 --
    ## -- here lamda = (1+u) / (1 - u) --
    lamda <- ( 1 + simu[, "mu"] ) / ( 1 - simu[, "mu"] );
    simu[, "skew" ] <- ( lamda - exp( - simu[, "S"]) ) / ( lamda + exp( - simu[, "S"] ) );

    ## -- plot --
    pdf(file = "fig.s2_simulation.pdf", height = 4, width = 5.5);
    library("ggplot2");
    ggplot( data = simu, aes( x = S, y = skew, group = mu, colour = mu ) ) +
        geom_line(size=0.3) +
        scale_x_log10( breaks = c(0.01, 0.1,  1,  10) ) +
        scale_colour_gradientn( colours = jet.colors(7), name = expression(mu) ) +
        xlab( expression( italic(S) ) ) + ylab( expression(paste( "Skew (",  gamma,")")) );
    dev.off();
} ## -- end --


#### ============= statistical tests ===============
{
    zAA.at <- ( zAA.1s.M6[["AT.le"]] + zAA.1s.M6[["AT.la"]])/2;
    r <- sum(zAA.at > 0, na.rm=TRUE) # 1542 -> 99.5%
    binom.test(r, dim(dat)[1], p = 0.5, alternative = "two.sided")

    zAA.gc <- ( zAA.1s.M6[["GC.le"]] + zAA.1s.M6[["GC.la"]])/2;
    r <- sum(zAA.gc > 0, na.rm=TRUE) # 1542 -> 99.5%
    binom.test(r, dim(dat)[1], p = 0.5, alternative = "two.sided")


    zRNA.at <- (z.4s.M6[["AT.le"]] + z.4s.M6[["AT.la"]])/2;
    r <- sum(zRNA.at < 0, na.rm=TRUE) # 559 -> 95.3%
    binom.test(r, dim(dat)[1], p = 0.5, alternative = "two.sided")

    zRNA.gc <- (z.4s.M6[["GC.le"]] + z.4s.M6[["GC.la"]])/2;
    r <- sum(zRNA.gc < 0, na.rm=TRUE) # 559 -> 95.3%
    binom.test(r, dim(dat)[1], p = 0.5, alternative = "two.sided")


    dMu.at <- dMu.io.M6[["AT"]];
    r <- sum(dMu.at < 0, na.rm=TRUE) # 1423 -> 66.8
    binom.test(r, dim(dat)[1], p = 0.5, alternative = "two.sided")
}

###########################################################################
##  how could codon usage bias at 4s sites be driven by nucleotide skews?
##  i.e. among preferred codons of all bacterial,
##    are there more T-ending than A-ending ones, and more C-ending than G-ending ones?
###########################################################################
{
    ## -- load codons and corresponding AAs --
    codons <- read.delim(file = "data_and_other_R_scripts/codons.txt", header = T, as.is = T, sep = "\t");
    ## -- load preferred codons in tested bacteria --
    preferredCodons <- read.delim(file = "data_and_other_R_scripts/bactarial_preferred_codons_april_2014.txt", header = T, as.is = T, na.strings = c("none") );
    preferredCodons <- preferredCodons[ rownames(preferredCodons) %in% rownames(dat),  ];
    preferredCodons[, "gc"] <- dat[ rownames(preferredCodons), "genomicGC"];

    ## -- analyze --
    four <- foursAnalysis( df = preferredCodons, codons= codons);
}


#### ============= correlations ===============
{
    ## -- correlation of replication skews calculated from io and 4s sites --
    cor.test( dMu.io.M6[["AT"]], dMu.4s.M6[["AT"]] ); # R=0.780
    cor.test( dMu.io.M6[["GC"]], dMu.4s.M6[["GC"]] ); # R=0.899

    ## -- transcription skews from 4s sites leading and lagging strands --
    cor.test( z.4s.M6[["AT.le"]], z.4s.M6[["AT.la"]] ); # R=0.771
    cor.test( z.4s.M6[["GC.le"]], z.4s.M6[["GC.la"]] ); # R=0.842

    ## -- translation skews from 2s sites leading and lagging strands --
    cor.test( zAA.2s.M6[["AT.le"]],  zAA.2s.M6[["AT.la"]] );
    cor.test( zAA.2s.M6[["GC.le"]], zAA.2s.M6[["GC.la"]] );
    # plot(zAA.2s.M6[["AT.le"]],  zAA.2s.M6[["AT.la"]], type = "p");

    ## -- translation skews from ns sites leading and lagging strands --
    cor.test( zAA.1s.M6[["AT.le"]],  zAA.1s.M6[["AT.la"]] );
    cor.test( zAA.1s.M6[["GC.le"]], zAA.1s.M6[["GC.la"]] );
}


## == June 8, 2014 ==
{
    ## -- correlation of skews attributed to replication, transcription and translation with genome GC content --
    cor.test(dMu.io.M6[["GC"]], dat$genomicGC); ## == -0.54
    cor.test(dMu.io.M6[["AT"]], dat$genomicGC); ## == -0.31

    cor.test(z.4s.M6[["GC.le"]], dat$genomicGC); ## == 0.053
    cor.test(z.4s.M6[["GC.la"]], dat$genomicGC); ## == 0.074
    cor.test(z.4s.M6[["AT.le"]], dat$genomicGC); ## == -0.31
    cor.test(z.4s.M6[["AT.la"]], dat$genomicGC); ## == -0.16

    cor.test(zAA.1s.M6[["GC.le"]], dat$genomicGC); ## == -0.28
    cor.test(zAA.1s.M6[["GC.la"]], dat$genomicGC); ## == -0.40
    cor.test(zAA.1s.M6[["AT.le"]], dat$genomicGC); ## == -0.60
    cor.test(zAA.1s.M6[["AT.la"]], dat$genomicGC); ## == -0.77
}

#################################################################
### Sep 19 ~ nov 4, 2014
###    does life-style contribute to the biases?? --
#################################################################
## =============================================================================================
## ---- we don't estimate zRNA for leading and lagging separately --
## ---- & we don't calculate tau ----
##   u4s.leading = zRNA + tau * io
##   u4s.lagging = zRNA - tau * io
##      then zRNA = ( u4s.leading + urs.lagging ) / 2
##   --> similarly --
##   u.ns.leading = zAA + zRNA + tau * io
##   u.ns.lagging = zAA + zRNA - tau * io
##      then zAA = ( u.ns.leading + u.ns.lagging ) / 2 - zRNA --
##
##   -- let's do it !! --
{
    ## -- March 22, 2015 --
    ## -------------- last modified nov 3, 2014  ---------------
    ## -- last modified : March 22, 2015 by WHC --
    df.lifestyle <- cbind(
        dat[ ,  c( "genomicGC", "genomicLen", "habitatNCBI", "intra", "habitatImg", "habitatGold" )],
        data.frame(
            zRNA.at = ( df.all[["zRNA.AT.le"]] + df.all[["zRNA.AT.la"]] ) / 2,
            zRNA.gc = ( df.all[["zRNA.GC.le"]] + df.all[["zRNA.GC.la"]] ) / 2,
            zAA.at = ( df.all[["zAA.AT.le"]] + df.all[["zAA.AT.la"]] ) / 2,
            zAA.gc = ( df.all[["zAA.GC.le"]] + df.all[["zAA.GC.la"]] ) / 2
        ) );

        table(df.lifestyle$intra);

    ## -- combine all data together --
    df.lifestyle[, "habitat"] <-  rowSums( data.frame( df.lifestyle$habitatNCBI != "HostAssociated",
                                                       df.lifestyle$habitatImg == "Environmental",
                                                       df.lifestyle$habitatGold == "Environmental"), na.rm = T );

    # ----
    df.size3 <- df.lifestyle;
    df.size3 <- df.size3[ df.size3$intra == "Y" | df.size3$habitat >= 3, ];
    table(df.size3$intra);
}


##########################################################################################
# use independent contrasts: caper package, brunch function
# see http://cran.r-project.org/web/packages/caper/vignettes/caper.pdf
# each row of data at the tips is used only once. This follows Burt (1989):
#  contrasts whose paths do not meet or cross at any point will be phylogenetically independent.
##########################################################################################
{
    ## -- load caper, install it if not yet ...
    if (!require("caper")){
        install.packages("caper", dependencies = TRUE);
        library(caper);
    }

    # reduce data matrix to what we'll really need:
    x <- df.size3
    df3 <- data.frame(acc=row.names(x), L=x$genomicLen/1e6, GC=x$genomicGC, intra=x$intra, zRNA.gc=x$zRNA.gc, zRNA.at=x$zRNA.at, zAA.gc=x$zAA.gc, zAA.at=x$zAA.at, row.names=row.names(x), stringsAsFactors=TRUE)

    # read tree (Newick format, rooted bewteen bacteria and archaea):
    tree.raw <- read.tree("data_and_other_R_scripts/RAxML_bestTree.multioriginal_rerooted.tre")
    is.rooted(tree.raw)
    # set branches of length 0 to 1e-6:
    tree.e6 <- tree.raw
    tree.e6$edge.length[tree.e6$edge.length<1e-6] <- 1e-6
    summary(tree.e6$edge.length)

    # make the tree ultrametric (unless we did this before):
    ultrametric.tree.exists <- TRUE
    if (!ultrametric.tree.exists) {
        tree.e6.chronos <- chronos(tree.e6)
        is.ultrametric(tree.e6.chronos)
        is.binary.tree(tree.e6.chronos)
        min(tree.e6.chronos$edge.length)
        write.tree(tree.e6.chronos, file="RAxML_bestTree.multioriginal_rerooted_chronos.tre")
    }
    tree <- read.tree("data_and_other_R_scripts/RAxML_bestTree.multioriginal_rerooted_chronos.tre") # re-reading this gives a proper phylo object!
}

#####################################################################################################################
#
# work with pgls function of caper, with an ULTRAMETRIC tree:
#
#####################################################################################################################
{
    # restrict tree to those included in the df3 dataset:
    unused.taxa <- setdiff(tree$tip.label, df3$acc)
    tree.df3 <- drop.tip(tree, unused.taxa);

    # make comparative data object (pgls doesn't like polytomies or branch lenght=0):
    comp <- comparative.data(phy=tree.df3, data=df3, names.col="acc", vcv=TRUE, vcv.dim=3, warn.dropped=TRUE, na.omit=TRUE)



    ################################################################
    # zRNA.gc
    ################################################################

    #crunch.zRNA.gc.2 <- crunch(zRNA.gc ~ intra * L * GC, comp, ref.var=intra, factor.action="allow")
    #summary(crunch.zRNA.gc.2)

    crunch.zRNA.gc <- crunch(zRNA.gc ~ intra * L * GC - GC , comp, ref.var=intra, factor.action="allow")
    summary(crunch.zRNA.gc)
    #par(mfrow=c(3,2));  caic.diagnostics(crunch.zRNA.gc)

    ################################################################
    # zRNA.at
    ################################################################

    #crunch.zRNA.at <- crunch(zRNA.at ~ intra * L * GC, comp, ref.var=intra, factor.action="allow")
    #summary(crunch.zRNA.at)
    #par(mfrow=c(3,2));  caic.diagnostics(crunch.zRNA.at)

    crunch.zRNA.at.2 <- crunch(zRNA.at ~ intra * L * GC - L - intra:GC:L - intra:GC, comp, ref.var=intra, factor.action="allow")
    summary(crunch.zRNA.at.2)
    par(mfrow=c(3,2));  caic.diagnostics(crunch.zRNA.at)

    AIC(crunch.zRNA.at, crunch.zRNA.at.2)

    ################################################################
    # zAA.gc
    ################################################################

    #crunch.zAA.gc <- crunch(zAA.gc ~ intra * L * GC, comp, ref.var=intra, factor.action="allow")
    #summary(crunch.zAA.gc)
    #par(mfrow=c(3,2));  caic.diagnostics(crunch.zRNA.gc)

    crunch.zAA.gc.2 <- crunch(zAA.gc ~ intra * L * GC  - L - GC:L, comp, ref.var=intra, factor.action="allow")
    summary(crunch.zAA.gc.2)
    #par(mfrow=c(3,2));  caic.diagnostics(crunch.zRNA.gc)

    AIC(crunch.zAA.gc, crunch.zAA.gc.2)

    ################################################################
    # zAA.at
    ################################################################

    #crunch.zAA.at <- crunch(zAA.at ~ intra * L * GC, comp, ref.var=intra, factor.action="allow")
    #summary(crunch.zAA.at)
    #par(mfrow=c(3,2));  caic.diagnostics(crunch.zAA.at)

    crunch.zAA.at.3 <- crunch(zAA.at ~ intra * L * GC - intra:GC , comp, ref.var=intra, factor.action="allow")
    summary(crunch.zAA.at.3)

    AIC(crunch.zAA.at, crunch.zAA.at.2);

}

## -- recalculate the correlations between leading and lagging strands with phylogenetic relateness is considered for ...  --
{
    # restrict tree to those included in the df3 dataset:
    df.all$acc = rownames(df.all);
    unused.taxa <- setdiff(tree$tip.label, df.all$acc)
    tree.dfall <- drop.tip(tree, unused.taxa)

    # make comparative data object (pgls doesn't like polytomies or branch lenght=0):
    comp.dfall <- comparative.data(phy=tree.dfall, data=df.all, names.col="acc", vcv=TRUE, vcv.dim=3, warn.dropped=TRUE, na.omit=TRUE)

    mod.mu.at <- pgls( dMuioAT ~ dMu4sAT , comp.dfall ); #
    mod.mu.gc <- pgls( dMuioGC ~ dMu4sGC , comp.dfall ); #

    mod.zRNA.at <- pgls( zRNA.AT.le ~ zRNA.AT.la , comp.dfall ); #
    mod.zRNA.gc <- pgls( zRNA.GC.le ~ zRNA.GC.la , comp.dfall ); #

    mod.zAA.at <- pgls( zAA.AT.le ~ zAA.AT.la , comp.dfall ); #
    mod.zAA.gc <- pgls( zAA.GC.le ~ zAA.GC.la , comp.dfall ); #


    ## -- correlation --
    ## == Mu: fitted ==
    cor.test( mod.mu.at[["fitted"]][, 1], df.all[ rownames( mod.mu.at[["fitted"]]), "dMuioAT" ] );
    cor.test( mod.mu.gc[["fitted"]][, 1], df.all[ rownames(mod.mu.gc[["fitted"]]), "dMuioGC" ] );

    ## -- original --
    cor.test( df.all$dMuioAT, df.all$dMu4sAT );
    cor.test( df.all$dMuioGC, df.all$dMu4sGC );

    ## == zRNA : fitted ==
    cor.test( mod.zRNA.at[["fitted"]][, 1], df.all[ rownames( mod.zRNA.at[["fitted"]]), "zRNA.AT.le" ] );
    cor.test( mod.zRNA.gc[["fitted"]][, 1], df.all[ rownames( mod.zRNA.gc[["fitted"]]), "zRNA.GC.le" ] );

    ## -- original --
    cor.test( df.all$zRNA.AT.le, df.all$zRNA.AT.la );
    cor.test( df.all$zRNA.GC.le, df.all$zRNA.GC.la );

    ## == zAA : fitted ==
    cor.test( mod.zAA.at[["fitted"]][, 1], df.all[ rownames( mod.zAA.at[["fitted"]]), "zAA.AT.le" ] );
    cor.test( mod.zAA.gc[["fitted"]][, 1], df.all[ rownames( mod.zAA.gc[["fitted"]]), "zAA.GC.le" ] );

    ## -- original --
    cor.test( df.all$zAA.AT.le, df.all$zAA.AT.la );
    cor.test( df.all$zAA.GC.le, df.all$zAA.GC.la );
}


#############################################
### march 30, 2015 ; some additional plots ...
### july 24, 2015; supplementary figures 6~8 ....
#############################################
{

    ############################
    ## -- for GC skews ...
    raw.gc <- data.frame( io = dMu.io.M6[["GC"]],
                          gc.4s.le = dMu.4s.M6[["raw.GC.le"]], gc.4s.la = dMu.4s.M6[["raw.GC.la"]],
                          gc.1s.le = dMu.1s.M6[["raw.GC.le"]], gc.1s.la = dMu.1s.M6[["raw.GC.la"]]);

    raw.gc.melted <- melt(raw.gc);
    pdf(file = "figN1.1.pdf", width = 5, height = 5);
    ## -- 1.	densities of all 5 together (similar to Figure 1, but only for GC)\
    p1.gc <- ggplot(raw.gc.melted, aes( x=value, group = variable )) + geom_density( aes(colour = factor(variable), linetype = factor(variable) )  ) +
        scale_colour_manual( values = c( color.rep, color.transcription, color.transcription, color.ns, color.ns ),
                             name = "raw skews at various sites",
                             labels = c( "interoperonic region", "4s sites on leading strand", "4s sites on lagging strand", "ns sits on leading strand", "ns sites on lagging strand" )

        ) +
        theme(legend.justification=c(1,1), legend.position=c(1,1), legend.background = element_rect(colour = NA, fill = NA) ) +
        #theme( legend.position = "none" ) +
        xlim( c(-0.65, 1.2) ) +
        xlab( "skews" ) + ylab( "density" );

    ## -- 2. scatter plot for 4s sites versus io: GC at leading and lagging respectively ...
    p2.1.gc <- assemData4FigN1("GC", "le", dMu.io.M6, dMu.4s.M6, dat$genomicGC, "raw GC skew of interoperonic region", "raw GC skew of 4s sites \nof the leading strand" );
    p2.2.gc <- assemData4FigN1("GC", "la", dMu.io.M6, dMu.4s.M6, dat$genomicGC, "raw GC skew of interoperonic region", "raw GC skew of 4s sites \nof the lagging strand" );

    ## -- 3. for 1s sites; GC at leading and lagging respecitvely ...
    p3.1.gc <- assemData4FigN1("GC", "le", dMu.io.M6, dMu.1s.M6, dat$genomicGC, "raw GC skew of interoperonic region", "raw GC skew of 1s sites \nof the leading strand" );
    p3.2.gc <- assemData4FigN1("GC", "la", dMu.io.M6, dMu.1s.M6, dat$genomicGC, "raw GC skew of interoperonic region", "raw GC skew of 1s sites \nof the lagging strand" );

    ## -- save to pdf ...
    pdf(file="fig.s6.panelA.pdf", height = 5, width = 5);
    p1.gc;
    dev.off();

    pdf(file = "Fig.s7.pdf", width = 13, height = 12);
    multiplot(p2.1.gc, p2.2.gc, p3.1.gc, p3.2.gc, cols= 2);
    dev.off();

    #############################
    ## for AT skews --
    raw.at <- data.frame( io = dMu.io.M6[["AT"]],
                          at.4s.le = dMu.4s.M6[["raw.AT.le"]], at.4s.la = dMu.4s.M6[["raw.AT.la"]],
                          at.1s.le = dMu.1s.M6[["raw.AT.le"]], at.1s.la = dMu.1s.M6[["raw.AT.la"]]);

    raw.at.melted <- melt(raw.at);
    ## -- 1.	densities of all 5 together (similar to Figure 1, but only for GC)\
    p1.at <- ggplot(raw.at.melted, aes( x=value, group = variable )) + geom_density( aes(colour = factor(variable), linetype = factor(variable) )  ) +
        scale_colour_manual( values = c( color.rep, color.transcription, color.transcription, color.ns, color.ns ),
                             name = "raw skews at various sites",
                             labels = c( "interoperonic", "4s sites leading strand", "4s sites lagging strand", "ns sits leading strand", "ns sites lagging strand" )

        ) +
        theme(legend.justification=c(1,1), legend.position=c(1,1), legend.background = element_rect(colour = NA, fill = NA) ) +
        #theme( legend.position = "none" ) +
        xlim( c(-0.65, 1.2) ) +
        xlab( "skews" ) + ylab( "density" );

    ## -- 2. scatter plot for 4s sites versus io: at at leading and lagging respectively ...
    p2.1.at <- assemData4FigN1("AT", "le", dMu.io.M6, dMu.4s.M6, dat$genomicGC, "raw AT skew of interoperonic region", "raw AT skew of 4s sites \nof the leading strand" );
    p2.2.at <- assemData4FigN1("AT", "la", dMu.io.M6, dMu.4s.M6, dat$genomicGC, "raw AT skew of interoperonic region", "raw AT skew of 4s sites \nof the lagging strand" );

    ## -- 3. for 1s sites; at at leading and lagging respecitvely ...
    p3.1.at <- assemData4FigN1("AT", "le", dMu.io.M6, dMu.1s.M6, dat$genomicGC, "raw AT skew of interoperonic region", "raw AT skew of 1s sites \nof the leading strand" );
    p3.2.at <- assemData4FigN1("AT", "la", dMu.io.M6, dMu.1s.M6, dat$genomicGC, "raw AT skew of interoperonic region", "raw AT skew of 1s sites \nof the lagging strand" );

    ## -- plot
    pdf(file="fig.s6.panelB.pdf", height = 5, width = 5);
    p1.at;
    dev.off();

    pdf(file = "fig.s8.pdf", width = 13, height = 12);
    multiplot(p2.1.at, p2.2.at, p3.1.at, p3.2.at, cols= 2);
    dev.off();


    #### figN2; not used in the final version; noted July 24, 2015 ... ####
    tmp <- data.frame( x = dMu.1s.M6[["AT"]], y = dMu.1s.M6[["AT.fac"]] * dMu.io[["AT"]], gc = dat$genomicGC );
    tmp <- tmp[ complete.cases(tmp), ];
    t <- bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( tmp$x, tmp$y )$estimate , digits=3) , sep = "")) ) );

    p.at <- ggplot(tmp, aes( x, y ) ) +
        geom_point( aes(colour = gc ), size = 0.8 ) +
        scale_colour_gradientn( colours = jet.colors(7), name = "genomic GC%" ) +
        xlab( expression(paste(italic(mu)," for AT derived from ns sites")) ) +
        ylab( expression(paste(italic(mu)," for AT derived from interoperonic sites") %*% tau)   ) +
        ggtitle( t ) +
        geom_abline( slope = 1, linetype = 2, colour = "darkgray", size = 0.5 );

    tmp <- data.frame( x = dMu.1s.M6[["GC"]], y = dMu.1s.M6[["GC.fac"]] * dMu.io[["GC"]], gc = dat$genomicGC );
    tmp <- tmp[ complete.cases(tmp), ];
    t <- bquote(.(parse(text=paste( "italic(R)== ", round( cor.test( tmp$x, tmp$y )$estimate , digits=3) , sep = "")) ) );

    p.gc <- ggplot(tmp, aes( x, y ) ) +
        geom_point( aes(colour = gc ), size = 0.8 ) +
        scale_colour_gradientn( colours = jet.colors(7), name = "genomic GC%" ) +
        xlab( expression(paste(italic(mu)," for GC derived from ns sites")) ) +
        ylab( expression(paste(italic(mu)," for GC derived from interoperonic sites") %*% tau)   ) +
        ggtitle( t ) +
        geom_abline( slope = 1, linetype = 2, colour = "darkgray", size = 0.5 );

    pdf(file = "FigN2.pdf", width = 13, height = 6);
    multiplot(p.at, p.gc, cols= 2);
    dev.off();
}

###############################################
### how long are the interopoeronic regions??
{
    ## mean = 8.8% and median = 7.9% ...
    summary( rowSums( dat[, c("ioAle", "ioTle", "ioGle", "ioCle")] ) / dat$genomicLen * 100 );
}


###############################################
### == correlations of triplet distributions from leading and lagging starands
### == 4s sites plus flanking nucleotides ==
{
    ## -- load data --
    codons <- read.table(file="~/Dropbox/perl_scripts/codon_analysis/valid_codons_and_revcoms_april2015.txt", header = T, as.is = T);
    files <- read.table(file = "~/databases/bacteria/analysis_bacterial_2012_2013/genomeAcc2Info/acc2flanking_nucleotides_4s_april_2015.txt", header = F, as.is = T);

    results <- c();
    for( i in 1:dim(files)[1]){
        acc <- files[i, "V1"];
        tmp <- read.table( file = files[i, "V2"], row.names = 1, header = T, as.is = T );

        tmp <- tmp[codons$codon, ];

        ## -- cor 1 --
        r1 <- cor.test( tmp$leading, tmp$lagging )$estimate;

        ## -- cor 2 --
        df <- data.frame( leading = tmp$leading, lagging = tmp[ as.character(codons$revcomp), "lagging" ] );
        df <- df[ complete.cases(df), ];
        r2 <- cor.test( df$leading, df$lagging )$estimate;

        ## -- leading versus leading rev --
        df <- data.frame( a = tmp$leading, b = tmp[ as.character(codons$revcomp), "leading" ] );
        df <- df[ complete.cases(df), ];
        r3 <- cor.test( df$a, df$b )$estimate;

        ## -- lagging versus lagging rev --
        df <- data.frame( a = tmp$lagging, b = tmp[ as.character(codons$revcomp), "lagging" ] );
        df <- df[ complete.cases(df), ];
        r4 <- cor.test( df$a, df$b )$estimate;


        results <- rbind(results, c( r1, r2, r3, r4));

        if( i %% 100 == 0 ){
            print(i);
        }
    }
    rownames(results) <- files[, "V1"];
    colnames(results) <- c("ce1", "ce2", "ce3", "ce4");

    results <- as.data.frame(results);

    r3 <- results;
    r3[, "gc"] <- dat[ rownames(r3), "genomicGC"]

    ## -- correlations between the correlation coefficient values and gc
    cor.test(r3$gc, r3$ce1);
    cor.test(r3$gc, r3$ce2);
    cor.test(r3$gc, r3$ce3);
    cor.test(r3$gc, r3$ce4);

    #r3 <- r3[ complete.cases(r3),  ];
    #heatmap(cor(r3));

    r2 <- results;
    r2[, "acc"] <- files[, "V1"];

    ## -- plot --
    library("ggplot2");
    library("reshape2");
    rmelt <- melt(r2, id.vars = "acc");
    rmelt[, "gc"] <- dat[ as.character(rmelt$acc), "genomicGC"];
    rmelt <- rmelt[ complete.cases(rmelt),  ];

    p1 <- ggplot( rmelt, aes(x = value, group = variable, colour = factor(variable)) ) + geom_density() +
        scale_colour_discrete( "Correlation coefficient",
                            breaks=c("ce1", "ce2", "ce3", "ce4"),
                            labels=c("leading versus lagging", "leading versus revcomp lagging",
                                     "leading versus revcomp leading", "lagging versus revcomp lagging")  );

    p2 <- ggplot(rmelt, aes(x = gc, y = value, group = variable, colour = factor(variable))) + geom_point() +
        scale_colour_discrete( "Correlation coefficient",
                               breaks=c("ce1", "ce2", "ce3", "ce4"),
                               labels=c("leading versus lagging", "leading versus revcomp lagging",
                                        "leading versus revcomp leading", "lagging versus revcomp lagging")  );

    pdf(file = "FigNx_ver2.pdf", width = 15, height = 6);
    multiplot(p1, p2, cols= 2);
    dev.off();


}

###########################################
## -- NOTE: not used; Sep 3, 2015 --
## ============== Sup figure 9, comparing Srna obtained from Eq3&4 and Eq4 of the supplemntary ... =======
## -- July 24, 2015 ...
{
    ## -- calculate raw skews --
    skew.le.4s.at <- (dat$A4le - dat$T4le) / (dat$A4le + dat$T4le);
    skew.la.4s.at <- (dat$A4la - dat$T4la) / (dat$A4la + dat$T4la);

    skew.le.4s.gc <- (dat$G4le - dat$C4le) / (dat$G4le + dat$C4le);
    skew.la.4s.gc <- (dat$G4la - dat$C4la) / (dat$G4la + dat$C4la);

    ## -- sRNAs using eq4 of the supplementary ... --
    sRNA.at <- log( ( (1 + skew.le.4s.at) / (1 - skew.le.4s.at) ) * ( (1 + skew.la.4s.at) / (1-skew.la.4s.at) ) )/2;
    sRNA.gc <- log( ( (1 + skew.le.4s.gc) / (1 - skew.le.4s.gc) ) * ( (1 + skew.la.4s.gc) / (1-skew.la.4s.gc) ) )/2;

    ## -- combine everything together ...
    df.all.2 <- cbind(
        dat[, c("genomicGC", "genomicLen")],
        data.frame( dMuioAT = dMu.io.M6[["AT"]], dMuioGC = dMu.io.M6[["GC"]],
                    dMu4sAT = dMu.4s.M6[["AT"]], dMu4sGC = dMu.4s.M6[["GC"]],
                    zRNA.AT.le = z.4s.M6[["AT.le"]], zRNA.AT.la = z.4s.M6[["AT.la"]],
                    zRNA.GC.le = z.4s.M6[["GC.le"]], zRNA.GC.la = z.4s.M6[["GC.la"]],
                    zAA.AT.le = zAA.1s.M6[["AT.le"]], zAA.AT.la = zAA.1s.M6[["AT.la"]],
                    zAA.GC.le = zAA.1s.M6[["GC.le"]], zAA.GC.la = zAA.1s.M6[["GC.la"]],
                    zRNA.AT.eq4 = sRNA.at, zRNA.GC.eq4 = sRNA.gc)
    );

    ## -- plot --
    pa <- plotPanelEq4Eq3( data.frame( x = df.all.2$zRNA.AT.eq4, y = df.all.2$zRNA.AT.le, gc = df.all.2$genomicGC ),
                           xlab = expression(paste(italic(S[RNA])," for AT, derived using EQ4")),
                           ylab =  expression(paste(italic(S[RNA])," for AT, derived using EQ3a, leading strand")) );


    pb <- plotPanelEq4Eq3( data.frame( x = df.all.2$zRNA.AT.eq4, y =df.all.2$zRNA.AT.la , gc = df.all.2$genomicGC ),
                           xlab = expression(paste(italic(S[RNA])," for AT, derived using EQ4")),
                           ylab =  expression(paste(italic(S[RNA])," for AT, derived using EQ3a, lagging strand")) );


    pc <- plotPanelEq4Eq3( data.frame( x = df.all.2$zRNA.GC.eq4, y = df.all.2$zRNA.GC.le, gc = df.all.2$genomicGC ),
                           xlab = expression(paste(italic(S[RNA])," for GC, derived using EQ4")),
                           ylab =  expression(paste(italic(S[RNA])," for GC, derived using EQ3a, leading strand")) );


    pd <- plotPanelEq4Eq3( data.frame( x = df.all.2$zRNA.GC.eq4, y =df.all.2$zRNA.GC.la , gc = df.all.2$genomicGC ),
                           xlab = expression(paste(italic(S[RNA])," for GC, derived using EQ4")),
                           ylab =  expression(paste(italic(S[RNA])," for GC, derived using EQ3a, lagging strand")) );


    pdf(file = "fig.s9_sRNA_eq34_and_eq4.pdf", height = 10, width = 12);
    multiplot(pa, pc, pb, pd,  cols= 2);
    dev.off();
}

