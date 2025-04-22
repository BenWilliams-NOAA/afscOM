simulate_multiple <- function(nsims, seeds=NA, nyears, nages, nsexes, nregions, nfleets,...){

    land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims))
    disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims))
    caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims))
    faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims))
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions, nsims))
    # naa[1,,,] = init_naa

    f           = array(NA, dim=c(nyears, 1, 1, nregions, nfleets, nsims))
    recruits    = array(NA, dim=c(nyears+1, 1, 1, nregions, nsims))

    survey_preds <- list(
        rpns = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys, nsims)),
        rpws = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys, nsims)),
        acs  = array(NA, dim=c(nyears, nages, nsexes, nregions, nsurveys+nfleets, nsims))
    )

    survey_obs <- list(
        catch = array(NA, dim=c(nyears, 1, 1, nregions, nfleets, nsims)),
        rpns = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys, nsims)),
        rpws = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys, nsims)),
        acs  = array(NA, dim=c(nyears, nages, nsexes, nregions, nsurveys+nfleets, nsims))
    )

    if(all(is.na(seeds))){
        set.seed(1120)
        seeds <- sample(1:(nsims*100), nsims)
    }

    print(paste("Running", nsims, "simulations."))
    for(s in 1:nsims){
        print(s)
        seed <- seeds[s]
        model_options$seed <- seed
        outs <- project(
            ...
        )

        # update state
        land_caa[,,,,,s] <- outs$land_caa
        disc_caa[,,,,,s] <- outs$disc_caa
        caa[,,,,,s] <- outs$caa
        faa[,,,,,s] <- outs$faa
        naa[,,,,s] <- outs$naa
        # output_matrices$recruits[y+1,,1,] <- apply(input$naa_tmp[,1,,,drop=FALSE], 3, sum)

        f[,,,,,s] <- outs$f
        # output_matrices$recruits[y+1,,,] <- rec

        if(model_options$simulate_observations){
            survey_preds$rpns[,,,,,s] <- outs$survey_preds$rpns
            survey_preds$rpws[,,,,,s] <- outs$survey_preds$rpws
            survey_preds$acs[,,,,,s]  <- outs$survey_preds$acs

            survey_obs$catch[,,,,,s] <- outs$survey_obs$catch
            survey_obs$rpns[,,,,,s] <- outs$survey_obs$rpns
            survey_obs$rpws[,,,,,s] <- outs$survey_obs$rpws
            survey_obs$acs[,,,,,s]  <- outs$survey_obs$acs
        }

    }

    outputs <- listN(land_caa, disc_caa, caa, faa, naa, f, recruits, survey_preds, survey_obs, seeds)
    return(outputs)

}
