## -- created in Dec 2013 (estimated);  
## -- last modified : May 18, 2014 .. 
## -- last modified : June 9, 2014 ..

## -- fourfold synonymous codon analyses --


foursAnalysis <- function( df, codons = codons ){
    gc <- mean( df$gc );
    fours <- table( unlist(df[, unique( codons$AA[ codons$Three == 4 ] ) ]) );
    fours <- fours[ names( fours ) %in% unique( codons$Codon[ codons$Three == 4 ] ) ];
    
    a.4s <- sum( fours[ names( fours ) %in% unique( codons$Codon[ codons$Three == 4 & codons$L3 == "A" ] )   ] );
    t.4s <- sum( fours[ names( fours ) %in% unique( codons$Codon[ codons$Three == 4 & codons$L3 == "T" ] )   ] );
    g.4s <- sum( fours[ names( fours ) %in% unique( codons$Codon[ codons$Three == 4 & codons$L3 == "G" ] )   ] );
    c.4s <- sum( fours[ names( fours ) %in% unique( codons$Codon[ codons$Three == 4 & codons$L3 == "C" ] )   ] );
    
    at.4s <- a.4s + t.4s;
    gc.4s <- g.4s + c.4s;
    
    results <- c();
    
    stat1 <- fisher.test( matrix( c( a.4s, t.4s, floor( at.4s / 2 ), floor( at.4s / 2 ) ), nrow = 2 ) );
    stat1$p.value;
    stat1$estimate;
    
    results <- rbind( results, c( stat1$p.value, stat1$estimate ) );
    
    stat2 <- fisher.test( matrix( c( g.4s, c.4s, floor( gc.4s / 2 ), floor( gc.4s / 2 ) ), nrow = 2 ) );
    stat2$p.value;
    stat2$estimate;
    
    results <- rbind( results, c( stat2$p.value, stat2$estimate ) );
    
    ag.4s <- a.4s + g.4s;
    tc.4s <- t.4s + c.4s;
    
    stat3 <- fisher.test( matrix( c( a.4s, g.4s, floor( ag.4s * (100 - gc) / 100 ), floor( ag.4s * gc / 100 ) ), nrow = 2 ) );
    stat3$p.value;
    stat3$estimate;
    
    results <- rbind( results, c( stat3$p.value, stat3$estimate ) );
    
    stat4 <- fisher.test( matrix( c( t.4s, c.4s, floor( tc.4s * (100 - gc) / 100 ), floor( tc.4s * gc / 100 ) ), nrow = 2 ) );
    stat4$p.value;
    stat4$estimate;
    
    results <- rbind( results, c( stat4$p.value, stat4$estimate ) );
    
    atgc.4s <- a.4s + t.4s + g.4s + c.4s;
    stat5 <- fisher.test( matrix( c( at.4s, gc.4s, floor( atgc.4s * (100 - gc) / 100 ), floor( atgc.4s * gc / 100 ) ), nrow = 2 ) );
    stat5$p.value;
    stat5$estimate;
    
    results <- rbind( results, c( stat5$p.value, stat5$estimate ) );
    
    rownames(results) <- c("AT", "GC", "AG", "TC", "ATGC");
    colnames(results) <- c("p.value", "OR");
    
    results[, "OR"] <- 1 / results[, "OR"];
    
    return( results );
}


## -- load the following three functions first --
## -- Models can be : Chen, or Model06 from Martin --
dMuCalc <- function( type, data = dat, iodata = "new", model = "Chen" ){
    results <- list();
    
    ## -- interoperonic --
    if( iodata == "new" ){
        dMu.GC <- ( data[, "ioGle"] - data[, "ioCle"] ) / ( data[, "ioGle"] + data[, "ioCle"] );
        dMu.AT <- ( data[, "ioAle"] - data[, "ioTle"] ) / ( data[, "ioAle"] + data[, "ioTle"] );
        
    } else {
        dMu.GC <- ( data[, "intergenicG"] - data[, "intergenicC"] ) / ( data[, "intergenicG"] + data[, "intergenicC"] );
        dMu.AT <- ( data[, "intergenicA"] - data[, "intergenicT"] ) / ( data[, "intergenicA"] + data[, "intergenicT"] );
        
    }
    
    if( type != "i" ){
        ## -- gc --
        gc.le <- (data[, paste( "G", type, "le" , sep = "")] - data[, paste( "C", type, "le" , sep = "")]) / ( data[, paste( "G", type, "le" , sep = "")] + data[, paste( "C", type, "le" , sep = "")] );
        gc.la <- (data[, paste( "G", type, "la" , sep = "")] - data[, paste( "C", type, "la" , sep = "")]) / ( data[, paste( "G", type, "la" , sep = "")] + data[, paste( "C", type, "la" , sep = "")] );
        # dMu.gc.coding = \Delta \mu' = \tau \Delta \mu
        if (model == "Model06") {
            # -- Model06 --
            paste("Using Model06\n") ;
            dMu.gc.coding <- ( gc.le - gc.la ) * ( 1 + gc.le * gc.la ) / ( 2 - ( gc.le ^ 2 + gc.la ^ 2 ) ) ;
        } else if (model == "Chen") {
            # -- simplest model --
            paste("Using simplest model\n") ;
            dMu.gc.coding <- ( gc.le - gc.la ) / 2;
        }
        
        ## -- at --
        at.le <- (data[, paste( "A", type, "le" , sep = "")] - data[, paste( "T", type, "le" , sep = "")]) / ( data[, paste( "A", type, "le" , sep = "")] + data[, paste( "T", type, "le" , sep = "")] );
        at.la <- (data[, paste( "A", type, "la" , sep = "")] - data[, paste( "T", type, "la" , sep = "")]) / ( data[, paste( "A", type, "la" , sep = "")] + data[, paste( "T", type, "la" , sep = "")] );
        # dMu.at.coding = \Delta \mu' = \tau \Delta \mu
        if (model == "Model06") {
            # -- Model06 --
            paste("Using Model06\n") ;
            dMu.at.coding <- ( at.le - at.la ) * ( 1 + at.le * at.la ) / ( 2 - ( at.le ^ 2 + at.la ^ 2 ) );
        } else if (model == "Chen") {
            # -- simplest model --
            paste("Using simplest model\n") ;
            dMu.at.coding <- ( at.le - at.la ) / 2;
        } 
        
        ## -- factors --
        factor.gc <- lm(  dMu.gc.coding ~ dMu.GC + 0 )$coefficients;  # = \tau
        factor.at <- lm(  dMu.at.coding ~ dMu.AT + 0 )$coefficients;  # = \tau
        
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


## -- created on May 18, 2014 --
## -- last modified : May 18, 2014 --
## -- input dMu.xs = dMu.1s or dMu.2s --
## -- return 
## -- [["GC.le"]], [["GC.la"]] -- zAA for GC at leading and lagging strand respectively --
##    [["AT.le"]], [["AT.la"]] -- zAA for AT at leading and lagging strand respectively --
zAAcalcM06 <- function( dMu.xs, dMu.4s, dMu.io, z.4s ){
    ## ----------------- GC -------------------
    ## -- using Martin's Model 6 --
    gc.le <- dMu.xs[["raw.GC.le"]];
    gc.la <- dMu.xs[["raw.GC.la"]];
    
    tau <- dMu.4s[["GC.fac"]] ;
    dMu.trans <- tau * dMu.io[["GC"]] ; # corresponds to dMu.4s.M6 
    
    ## -- fit 
    fit <- nls( dMu.trans ~ f1 * ( gc.le - gc.la ) * ( 1 + f1^2 * gc.le * gc.la ) / ( 2 - f1^2 * (gc.le^2 + gc.la^2) ) ,
                start = list( f1 = 0.4 ), trace=TRUE, na.action="na.omit" ) ;
    f1.gc <- coefficients( fit ) ; # f1.gc ; 1/f1.gc ; # fraction of free sites f_ns = 1/f1 = 0.410 $$$
    
    ## -- fitted 
    fitted <- fitted.values(fit);
    #plot(dMu.trans[!is.na(dMu.trans)], fitted); abline(0,1)
    #cor.test(dMu.trans[!is.na(dMu.trans)], fitted); # R=0.915
    
    # selection at 1s sites (Eq. 15 & 16):
    z.1s.le.M6.gc <- 2 * (gc.le * f1.gc - dMu.trans) / ( 1 - gc.le^2 * f1.gc^2 ) ;
    z.1s.la.M6.gc <- 2 * (gc.la * f1.gc + dMu.trans) / ( 1 - gc.la^2 * f1.gc^2 ) ;
    
    # selection on AAs (Eq. 17):
    zAA.1s.le.M6.gc <- z.1s.le.M6.gc - z.4s[["GC.le"]]; ## fixed a bug here
    zAA.1s.la.M6.gc <- z.1s.la.M6.gc - z.4s[["GC.la"]]; ## fixed a bug here ; Sep 25, 2014
    
    ## -- prepare output --
    result <- list();
    result[["GC.la"]] <- zAA.1s.la.M6.gc;
    result[["GC.le"]] <- zAA.1s.le.M6.gc;
    
    ## ----------------- AT -------------------
    at.le <- dMu.xs[["raw.AT.le"]];
    at.la <- dMu.xs[["raw.AT.la"]];
    
    tau <- dMu.4s[["AT.fac"]] ;
    dMu.trans <- tau * dMu.io[["AT"]] ; # corresponds to dMu.4s.M6 
    
    ## -- fit 
    fit <- nls( dMu.trans ~ f1 * ( at.le - at.la ) * ( 1 + f1^2 * at.le * at.la ) / ( 2 - f1^2 * (at.le^2 + at.la^2) ) ,
                start = list( f1 = 0.4 ), trace=TRUE, na.action="na.omit" ) ;
    f1.at <- coefficients( fit ) ;  # f1.at ; 1/f1.at ; # fraction of free sites f_ns = 1/f1 = 0.410 $$$
    
    ## -- fitted 
    fitted <- fitted.values(fit);
    #plot(dMu.trans[!is.na(dMu.trans)], fitted); abline(0,1)
    #cor.test(dMu.trans[!is.na(dMu.trans)], fitted); # R=0.915
    
    # selection at 1s sites (Eq. 15 & 16):
    z.1s.le.M6.at <- 2 * (at.le * f1.at - dMu.trans) / ( 1 - at.le^2 * f1.at^2 ) ;
    z.1s.la.M6.at <- 2 * (at.la * f1.at + dMu.trans) / ( 1 - at.la^2 * f1.at^2 ) ;
    
    # selection on AAs (Eq. 17):
    zAA.1s.le.M6.at <- z.1s.le.M6.at - z.4s[["AT.le"]]; ## bug fix ; sep 25, 2014;
    zAA.1s.la.M6.at <- z.1s.la.M6.at - z.4s[["AT.la"]]; ## bug fix ; sep 25, 2014;
    
    ## -- prepare output --
    result[["AT.la"]] <- zAA.1s.la.M6.at;
    result[["AT.le"]] <- zAA.1s.le.M6.at;
    
    ## -- return --
    return( result );
}


## -- last modified : Dec 25, 2013 --
## -- NOTE: simplest model is the same as other models --
## -- also switch leading and lagging strand --
## -- input : z = z scores at ns and 2s sites
##            z.trans = z scores at 4s sites ...
## -- structure of resulting lists --
## -- [["GC.le"]], [["GC.la"]] -- zAA for GC at leading and lagging strand respectively --
##    [["AT.le"]], [["AT.la"]] -- zAA for AT at leading and lagging strand respectively --
##    [["AT.fac"]], [["GC.fac"]] -- estimated factor x of transcription (z4s) for AT and GC respectively --
##    [["fac.raw"]] -- an dataframe shows how fac.AT and fac.GC are determined --
zAAcalc <- function( z, z.trans = z.4s ){
    ## -- first of all, calculate factors for GC and AT respectively --
    fac.raw <- c(); ## a data frame to hold different factor values and corresponding correlations --
    
    for( x in seq(0, 3, by = 0.001) ){
        gc <- cor.test(  z[["GC.le"]] - z.trans[["GC.le"]] * x, z[["GC.la"]] - z.trans[["GC.la"]] * x )$estimate;
        at <- cor.test(  z[["AT.le"]] - z.trans[["AT.le"]] * x, z[["AT.la"]] - z.trans[["AT.la"]] * x )$estimate;
        
        gc.ran <- cor.test(  z[["GC.le"]] - sample( z.trans[["GC.le"]] * x ), z[["GC.la"]] - sample( z.trans[["GC.la"]] * x ) )$estimate;
        at.ran <- cor.test(  z[["AT.le"]] - sample( z.trans[["AT.le"]] * x ), z[["AT.la"]] - sample( z.trans[["AT.la"]] * x ) )$estimate;
        
        fac.raw <- rbind( fac.raw, c( x, gc, at, gc.ran, at.ran) );
    }
    
    fac.raw <- as.data.frame( fac.raw );
    colnames( fac.raw ) <- c( "x", "gc", "at", "gc.ran", "at.ran" );
    
    at.fac <- fac.raw[ order(fac.raw[, "at"], decreasing= T),  ][1, "x"];
    gc.fac <- fac.raw[ order(fac.raw[, "gc"], decreasing= T),  ][1, "x"];
    
    ## -- calculate zAA --
    results <- list();
    results[["GC.le"]] <- z[["GC.le"]] - z.trans[["GC.le"]] * gc.fac;
    results[["GC.la"]] <- z[["GC.la"]] - z.trans[["GC.la"]] * gc.fac;
    
    results[["AT.le"]] <- z[["AT.le"]] - z.trans[["AT.le"]] * at.fac;
    results[["AT.la"]] <- z[["AT.la"]] - z.trans[["AT.la"]] * at.fac;
    
    results[["AT.fac"]] <- at.fac;
    results[["GC.fac"]] <- gc.fac;
    
    results[["fac.raw"]] <- fac.raw;
    
    ## -- switch / swap --
    results[["GC.le.s"]] <- z[["GC.le"]] - z.trans[["GC.la"]] * gc.fac;
    results[["GC.la.s"]] <- z[["GC.la"]] - z.trans[["GC.le"]] * gc.fac;
    
    results[["AT.le.s"]] <- z[["AT.le"]] - z.trans[["AT.la"]] * at.fac;
    results[["AT.la.s"]] <- z[["AT.la"]] - z.trans[["AT.le"]] * at.fac;
    
    
    return( results );
}

## -- last modified: Dec 25, 2013 --
## -- Jan 20, 2014: simplest model --
## -- 
## -- input :
##            dMu  = delta Mu at 1,2 or 4 sites ...
##            dMu.io = delta Mu derived from interoperonic regions, default = dMu.io --
## -- output : list of :
##             [["GC.le"]], [["GC.la"]] = calculated from dMu and raw skews ...
##             [["AT.le"]], [["AT.la"]] = ...
zCalc <- function( dMu, dMuIO = dMu.io, model = "Chen" ){
    results <- list();
    
    if (model == "Model06") {
        # -- Model06 --
        results[[ "GC.le" ]] <- 2 * ( dMu[[ "raw.GC.le" ]] - dMuIO[["GC"]] * dMu[["GC.fac"]] ) / ( 1 - dMu[[ "raw.GC.le" ]] ^2 ) ;
        results[[ "GC.la" ]] <- 2 * ( dMu[[ "raw.GC.la" ]] + dMuIO[["GC"]] * dMu[["GC.fac"]] ) / ( 1 - dMu[[ "raw.GC.la" ]] ^2 ) ;
        
        results[[ "AT.le" ]] <- 2 * ( dMu[[ "raw.AT.le" ]] - dMuIO[["AT"]] * dMu[["AT.fac"]] ) / ( 1 - dMu[[ "raw.AT.le" ]] ^2 ) ;
        results[[ "AT.la" ]] <- 2 * ( dMu[[ "raw.AT.la" ]] + dMuIO[["AT"]] * dMu[["AT.fac"]] ) / ( 1 - dMu[[ "raw.AT.la" ]] ^2 ) ;
        
    } else if (model == "Chen") {
        ## -- simplest model --
        results[[ "GC.le" ]] <- ( dMu[[ "raw.GC.le" ]] - dMuIO[["GC"]] * dMu[["GC.fac"]] ) ;
        results[[ "GC.la" ]] <- ( dMu[[ "raw.GC.la" ]] + dMuIO[["GC"]] * dMu[["GC.fac"]] ) ;
        
        results[[ "AT.le" ]] <- ( dMu[[ "raw.AT.le" ]] - dMuIO[["AT"]] * dMu[["AT.fac"]] ) ;
        results[[ "AT.la" ]] <- ( dMu[[ "raw.AT.la" ]] + dMuIO[["AT"]] * dMu[["AT.fac"]] ) ;
    }
    
    ## -- return -
    return( results );
}


## -- taiStats; Dec 14, 2013  --
## -- last modified : March 6, 2014 --
## -- last modified / verified : May 18, 2014 --
## -- last modified : Nov 3, 2014 --
## -- input : data  = tai10 or tai5  or things similar (table S2) --
## --         dat = counts of letters at 1s, 2s, ns and io regions --
## --         dMuIO = dMuIO --
## --         model = Model06 | Chen, default == Model06
## -- output : a list of lists, each list contains list[[i]] = list( dMu.4s, dMu.2s, dMu.1s, z.4s, z.2s, z.1s, zAA.2s, zAA.2s, zAA.1s ); ...
taiStats <- function( data, dat, dMuIO, z4s, model = "Model06" ){
    tai.stats <- list();
    for( i in 1:max( data[, "quantile"] ) ){
        print(i);
        ## -- prepare data --
        tmp <- data[ data$quantile == i & data$acc %in% rownames(dat), ];
        rownames(tmp) <- tmp$acc;
        tmp <- tmp[ rownames(dat), ];
        
        ## -- calculate dMus ....
        l <- list();
        l[["dMu.4s"]] <- dMuCalc4PrecalcSkews( "4", data = tmp, data2 = dat, model = model );
        l[["dMu.2s"]] <- dMuCalc4PrecalcSkews( "2", data = tmp, data2 = dat, model = model );
        l[["dMu.1s"]] <- dMuCalc4PrecalcSkews( "1", data = tmp, data2 = dat, model = model );
        
        ## -- z --
        l[["z.4s"]] <- zCalc( l[["dMu.4s"]], dMuIO = dMuIO, model= model );
        l[["z.2s"]] <- zCalc( l[["dMu.2s"]], dMuIO = dMuIO, model= model );
        l[["z.1s"]] <- zCalc( l[["dMu.1s"]], dMuIO = dMuIO, model= model );
        
        ## -- zaa --
        if( model == "Model06" ){
            l[["zAA.2s"]] <- zAAcalcM06(  dMu.xs= l[["dMu.2s"]], dMu.4s= l[["dMu.4s"]], dMu.io= dMuIO, z.4s = z4s );
            l[["zAA.1s"]] <- zAAcalcM06(  dMu.xs= l[["dMu.1s"]], dMu.4s= l[["dMu.4s"]], dMu.io= dMuIO, z.4s = z4s );
        } else if( model == "Chen" ){
            l[["zAA.2s"]] <- zAAcalc( l[["z.2s"]], z.trans= l[["z.4s"]] );
            l[["zAA.1s"]] <- zAAcalc( l[["z.1s"]], z.trans= l[["z.4s"]] );
        }
        
        ## -- save --
        tai.stats[[i]] <- l;
    }
    
    return( tai.stats );
}


## -- Dec 13, 2013 --
## -- last modified : March 6, 2014 --
## -- last verified / modified  : May 18, 2014 
## -- input: type = one of "4", "2", "1" --
## --        data = tai10, NOTE: data contains precalculated skews --
## --        data2 = dat  ...
## --        model = Chen | Model06 
## - output : a list : [["dMuAT"]], [["dMuGC"]]
dMuCalc4PrecalcSkews <- function( type, data = tai10, data2 = dat, iodata = "new", model = "Model06" ){
    results <- list();
    
    ## -- interoperonic --
    if( iodata == "new" ){
        dMu.GC <- ( data2[, "ioGle"] - data2[, "ioCle"] ) / ( data2[, "ioGle"] + data2[, "ioCle"] );
        dMu.AT <- ( data2[, "ioAle"] - data2[, "ioTle"] ) / ( data2[, "ioAle"] + data2[, "ioTle"] );
        
    } else {
        dMu.GC <- ( data2[, "intergenicG"] - data2[, "intergenicC"] ) / ( data2[, "intergenicG"] + data2[, "intergenicC"] );
        dMu.AT <- ( data2[, "intergenicA"] - data2[, "intergenicT"] ) / ( data2[, "intergenicA"] + data2[, "intergenicT"] );
    }
    
    
    if( type != "i" ){
        ## -- gc --
        gc.le <- as.numeric(data[, paste( "GCle", type, sep = "")] );
        gc.la <- as.numeric(data[, paste( "GCla", type, sep = "")] );
        
        # dMu.gc.coding = \Delta \mu' = \tau \Delta \mu
        if (model == "Model06") {
            # -- Model06 --
            paste("Using Model06\n") ;
            dMu.gc.coding <- ( gc.le - gc.la ) * ( 1 + gc.le * gc.la ) / ( 2 - ( gc.le ^ 2 + gc.la ^ 2 ) ) ;
        } else if (model == "Chen") {
            # -- simplest model --
            paste("Using simplest model\n") ;
            dMu.gc.coding <- ( gc.le - gc.la ) / 2;
        }
        
        ## -- at --
        at.le <- as.numeric(data[, paste( "ATle", type, sep = "")]);
        at.la <- as.numeric(data[, paste( "ATla", type, sep = "")]);
        # dMu.at.coding = \Delta \mu' = \tau \Delta \mu
        if (model == "Model06") {
            # -- Model06 --
            paste("Using Model06\n") ;
            dMu.at.coding <- ( at.le - at.la ) * ( 1 + at.le * at.la ) / ( 2 - ( at.le ^ 2 + at.la ^ 2 ) );
        } else if (model == "Chen") {
            # -- simplest model --
            paste("Using simplest model\n") ;
            dMu.at.coding <- ( at.le - at.la ) / 2;
        } 
        
        ## -- factors --
        factor.gc <- lm(  dMu.gc.coding ~ dMu.GC + 0 )$coefficients;  # = \tau
        factor.at <- lm(  dMu.at.coding ~ dMu.AT + 0 )$coefficients;  # = \tau
        
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



## -- Dec 14, 2013 --
## -- last modified : may 18, 2014 --
## -- input : ouput of taiStats --
##            strand : le = leading, la = lagging, combined = (le + la) / 2
## -- output: a data frame in which the data can be used for boxplots 
getDataFromTaiStats1 <- function( tai.stats, strand = "le" ){
    
    df <- c();
    for( i in 1:(length(tai.stats)) ){
        
        if( strand == "le" ){
            trans.at <- tai.stats[[i]][["z.4s"]][[ "AT.le" ]];
            trans.gc <- tai.stats[[i]][["z.4s"]][[ "GC.le" ]];
            
            zAA2s.at <- tai.stats[[i]][["zAA.2s"]][[ "AT.le" ]];
            zAA2s.gc <- tai.stats[[i]][["zAA.2s"]][[ "GC.le" ]];
            
            zAA1s.at <- tai.stats[[i]][["zAA.1s"]][[ "AT.le" ]];
            zAA1s.gc <- tai.stats[[i]][["zAA.1s"]][[ "GC.le" ]];
        } else if ( strand == "la" ){
            trans.at <- tai.stats[[i]][["z.4s"]][[ "AT.la" ]];
            trans.gc <- tai.stats[[i]][["z.4s"]][[ "GC.la" ]];
            
            zAA2s.at <- tai.stats[[i]][["zAA.2s"]][[ "AT.la" ]];
            zAA2s.gc <- tai.stats[[i]][["zAA.2s"]][[ "GC.la" ]];
            
            zAA1s.at <- tai.stats[[i]][["zAA.1s"]][[ "AT.la" ]];
            zAA1s.gc <- tai.stats[[i]][["zAA.1s"]][[ "GC.la" ]];
        } else if ( strand == "combined" ){
            trans.at <- ( tai.stats[[i]][["z.4s"]][[ "AT.le" ]] + tai.stats[[i]][["z.4s"]][[ "AT.la" ]] ) / 2 ;
            trans.gc <- (tai.stats[[i]][["z.4s"]][[ "GC.le" ]] + tai.stats[[i]][["z.4s"]][[ "GC.la" ]] ) / 2;
            
            zAA2s.at <- ( tai.stats[[i]][["zAA.2s"]][[ "AT.le" ]] + tai.stats[[i]][["zAA.2s"]][[ "AT.la" ]] ) / 2 ;
            zAA2s.gc <- (tai.stats[[i]][["zAA.2s"]][[ "GC.le" ]] + tai.stats[[i]][["zAA.2s"]][[ "GC.la" ]] ) / 2;
            
            zAA1s.at <- ( tai.stats[[i]][["zAA.1s"]][[ "AT.le" ]] + tai.stats[[i]][["zAA.1s"]][[ "AT.la" ]] ) / 2 ;
            zAA1s.gc <- (tai.stats[[i]][["zAA.1s"]][[ "GC.le" ]] + tai.stats[[i]][["zAA.1s"]][[ "GC.la" ]] ) / 2;
        }
        
        df <- rbind( df, 
                     data.frame( 
                         tai = i, 
                         acc = rownames( dat ),
                         "trans.at" = trans.at, 
                         "trans.gc" = trans.gc,
                         
                         "zAA2.at" = zAA2s.at, 
                         "zAA2.gc" = zAA2s.gc,
                         
                         "zAA1.at" = zAA1s.at, 
                         "zAA1.gc" = zAA1s.gc
                     ) );
    } ## -- end of for( i in 1:10 )
    
    df <- df[ complete.cases(df), ];
    
    return( df );
}
