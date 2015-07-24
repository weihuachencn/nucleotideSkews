## -- created on April 29, 2015;
## -- for flanking nucleotide analysis

## -- input : type = i or 1,2,4
##            data = list of two data.frames
##                  each df contains five columns, acc, leading, io, cat, lagging --
## -- output :
dMuCalcFlanking <- function( type, data = dat, flankingNuc = "TxG"  ){
    results <- list();

    at <- as.data.frame( data[["at"]] ); colnames( at )<- c("acc", "le", "io", "cat", "la")  ;
    gc <- as.data.frame( data[["gc"]] ); colnames( gc )<- c("acc", "le", "io", "cat", "la")  ;

    ## -- limit to a certain type of flanking nucleotides --
    at <- at[ at[, "cat"] == flankingNuc,  ];
    gc <- gc[ gc[, "cat"] == flankingNuc,  ];

    ## -- interoperonic -- UNCHANGED --
    dMu.GC <- as.numeric(as.vector( gc[, "io" ] ));
    dMu.AT <- as.numeric(as.vector( at[, "io" ] ));

    if( type != "i" ){
        ## -- gc --
        gc.le <- as.numeric(as.vector( gc[, "le" ] ));
        gc.la <- as.numeric(as.vector( gc[, "la" ] ));

        # CALCULATE MUTATIONAL CONTRIBUTION:
        # dMu.gc.coding = \Delta \mu' = \tau \Delta \mu
        paste("Using Model08\n") ;
        temp <- sqrt( (1 + gc.le) / (1 - gc.le) * (1 - gc.la) / (1 + gc.la) )
        dMu.gc.coding <- (temp - 1) / (temp + 1) ;

        ## -- at --
        at.le <- as.numeric(as.vector( at[, "le" ] ));
        at.la <- as.numeric(as.vector( at[, "la" ] ));

        # dMu.at.coding = \Delta \mu' = \tau \Delta \mu
        # -- Model08 --
        paste("Using Model08\n") ;
        temp <- sqrt( (1 + at.le) / (1 - at.le) * (1 - at.la) / (1 + at.la) )
        dMu.at.coding <- (temp - 1) / (temp + 1) ;

        ## -- factors --
        factor.gc <- lm(  dMu.gc.coding ~ dMu.GC + 0 )$coefficients;  # = \tau
        factor.at <- lm(  dMu.at.coding ~ dMu.AT + 0 )$coefficients;  # = \tau

        results[["GC.fac.ci"]] <- confint( lm(  dMu.gc.coding ~ dMu.GC + 0 ) );
        results[["AT.fac.ci"]] <- confint( lm(  dMu.at.coding ~ dMu.AT + 0 ) );

        ## -- save --
        results[["AT"]] <- dMu.at.coding; # = \Delta \mu' = \tau \Delta \mu
        results[["GC"]] <- dMu.gc.coding; # = \Delta \mu' = \tau \Delta \mu
        results[["AT.fac"]] <- factor.at; # = \tau
        results[["GC.fac"]] <- factor.gc; # = \tau

        ## -- also save raw data --
        results[["raw.GC.le"]] <- gc.le;  # G - C
        results[["raw.GC.la"]] <- gc.la;  # G^- - C^-
        results[["raw.AT.le"]] <- at.le;  # A - T
        results[["raw.AT.la"]] <- at.la;  # A^- - T^-

    } else {
        results[["AT"]] <- dMu.AT;  # = \Delta \mu
        results[["GC"]] <- dMu.GC;  # = \Delta \mu
    }

    return( results );
}

## -- last modified: Dec 25, 2013 --
## -- Jan 20, 2014: simplest model --
## -- March 10, 2015: small-N limit (Bulmer) in Model08
##
## -- input :
##            dMu  = delta Mu at 1,2 or 4 sites ...
##            dMu.io = delta Mu derived from interoperonic regions, default = dMu.io --
## -- output : list of :
##             [["GC.le"]], [["GC.la"]] = calculated from dMu and raw skews ...
##             [["AT.le"]], [["AT.la"]] = ...
zRNACalcdMuCalcFlanking <- function( dMu, dMuIO = dMu.io ){
    results <- list();

    dMu.GC <- dMuIO[["GC"]] * dMu[["GC.fac"]]
    dMu.AT <- dMuIO[["AT"]] * dMu[["AT.fac"]]

    # -- Model08 --
    gc.le <- -log( (1 - dMu[[ "raw.GC.le" ]]) / (1 + dMu[[ "raw.GC.le" ]]) * (1 + dMu.GC) / (1 - dMu.GC) ) ;
    gc.la <- -log( (1 - dMu[[ "raw.GC.la" ]]) / (1 + dMu[[ "raw.GC.la" ]]) * (1 - dMu.GC) / (1 + dMu.GC) ) ;
    at.le <- -log( (1 - dMu[[ "raw.AT.le" ]]) / (1 + dMu[[ "raw.AT.le" ]]) * (1 + dMu.AT) / (1 - dMu.AT) ) ;
    at.la <- -log( (1 - dMu[[ "raw.AT.la" ]]) / (1 + dMu[[ "raw.AT.la" ]]) * (1 - dMu.AT) / (1 + dMu.AT) ) ;

    ## -- remove NaNs --
    gc.le[ is.nan(gc.le) ] <- 0;
    gc.la[ is.nan(gc.la) ] <- 0;
    at.le[ is.nan(at.le) ] <- 0;
    at.la[ is.nan(at.la) ] <- 0;

    ## --
    results[[ "GC.le" ]] <- gc.le;
    results[[ "GC.la" ]] <- gc.la;

    results[[ "AT.le" ]] <- at.le;
    results[[ "AT.la" ]] <- at.la;

    ## -- return -
    return( results );
}
