#' Plot demographic parameters
#'
#' Generate standard plot for all demographic
#' parameter matrices. Output plots will include:
#' - weight-at-age
#' - maturity
#' - natural mortality
#' - selectivity
#' - retention
#' - discard mortality
#' Plots will show values of demographic quantities across
#' time, age, sex, region, and fleet as appropriate.
#'
#' @param dem_params list of demographic parameter matrices
#' @param params vector of names of parameters to plot (not implemented)
#' @param out_dir optional directory to save plots to
#' @param ... additional parameters to be passed to `ggsave`
#'
#' @export
#'
#'
plot_demographic_parameters <- function(dem_params, params=NA, show_plots=TRUE, out_dir=NA, ...){

    waa_plot  <- plot_waa(dem_params$waa)
    mat_plot  <- plot_mat(dem_params$mat)
    mort_plot <- plot_mort(dem_params$mort, is_dmr=FALSE)
    dmr_plot  <- plot_mort(dem_params$dmr,  is_dmr=FALSE)
    sel_plot  <- plot_selret(dem_params$sel, is_selectivity = TRUE)
    ret_plot  <- plot_selret(dem_params$ret, is_selectivity = FALSE)
    survsel_plot <- plot_selret(dem_params$surv_sel, is_selectivity = TRUE)

    if(show_plots){
        show(waa_plot)
        show(mat_plot)
        show(mort_plot)
        show(dmr_plot)
        show(sel_plot)
        show(ret_plot)
        show(survsel_plot)
    }

    if(!is.na(out_dir)){
        ggsave(filename=file.path(out_dir, "waa_plot.png"), plot=waa_plot, ...)
        ggsave(filename=file.path(out_dir, "mat_plot.png"), plot=mat_plot, ...)
        ggsave(filename=file.path(out_dir, "mort_plot.png"), plot=mort_plot, ...)
        ggsave(filename=file.path(out_dir, "dmr_plot.png"), plot=dmr_plot, ...)
        ggsave(filename=file.path(out_dir, "sel_plot.png"), plot=sel_plot, ...)
        ggsave(filename=file.path(out_dir, "ret_plot.png"), plot=ret_plot, ...)
        ggsave(filename=file.path(out_dir, "survsel_plot.png"), plot=survsel_plot, ...)
    }

    return(listN(waa_plot, mat_plot, mort_plot, dmr_plot, sel_plot, ret_plot, survsel_plot))

}

#' Plot Weight-at-Age
#'
#' Generates plot of weight-at-age across time, age, sex,
#' and region as appropriate.
#'
#' @param waa four-dimensional weight-at-age matrix
#'
#' @export
#'
#'
plot_waa <- function(waa){

    dimensions <- get_model_dimensions(waa)

    waa_df <- reshape2::melt(waa) %>%
        group_by(age, sex, region) %>%
        distinct(value, .keep_all=TRUE) %>%
        mutate(time_block = factor(time, labels=c(1:length(unique(time)))))

    ymax <- round(1.2*waa_df %>% pull(value) %>% max, 2)

    plot <- ggplot(waa_df, aes(x=age, y=value, color=sex, linetype=time_block))+
        geom_line(linewidth=1)+
        scale_y_continuous(limits=c(0, ymax))+
        coord_cartesian(expand=0)+
        labs(x="Age", y="Weight", color="Sex", linetype="Time Block", title="Weight-at-Age")+
        theme_bw()

    if(dimensions$nregions > 1){
        plot <- plot + facet_wrap(~region, scales="free_y")
    }

    return(plot)

}

#' Plot Maturity-at-Age
#'
#' Generates plot of maturity-at-age across time, age, sex,
#' and region as appropriate. Only the female maturity ogive
#' is plotted.
#'
#' @param mat four-dimensional maturity-at-age matrix
#'
#' @export
#'
plot_mat <- function(mat){

    dimensions <- get_model_dimensions(mat)

    mat_df <- reshape2::melt(mat) %>%
        group_by(age, sex, region) %>%
        distinct(value, .keep_all=TRUE) %>%
        mutate(time_block = factor(time, labels=c(1:length(unique(time))))) %>%
        filter(sex == "F")

    ymax <- 1.0

    plot <- ggplot(mat_df, aes(x=age, y=value, color=sex, linetype=time_block))+
        geom_line(linewidth=1)+
        scale_y_continuous(limits=c(0, ymax))+
        coord_cartesian(expand=0)+
        labs(x="Age", y="Maturity", color="Sex", linetype="Time Block", title="Female Maturity-at-Age")+
        theme_bw()

    if(dimensions$nregions > 1){
        plot <- plot + facet_wrap(~region, scales="free_y")
    }

    return(plot)

}

#' Plot Mortality
#'
#' Generates plot of mortality-at-age across time, age, sex,
#' and region as appropriate. This function can be used for
#' both natural and discard mortality
#'
#' @param waa four-dimensional mortality-at-age matrix
#'
#' @export
#'
plot_mort <- function(mort, is_dmr=FALSE){

    dimensions <- get_model_dimensions(mort)

    mort_df <- reshape2::melt(mort) %>%
        group_by(age, sex, region) %>%
        distinct(value, .keep_all=TRUE) %>%
        mutate(time_block = factor(time, labels=c(1:length(unique(time)))))

    ymax <- round(1.2*mort_df %>% pull(value) %>% max, 2)

    plot <- ggplot(mort_df, aes(x=age, y=value, color=sex, linetype=time_block))+
        geom_line(linewidth=1)+
        scale_y_continuous(limits=c(0, ymax), expand=c(0.01, 0.01))+
        scale_x_continuous(expand=c(0, 0))+
        # coord_cartesian(expand=0)+
        labs(x="Age", y="Instanteous Mortality", color="Sex", linetype="Time Block", title=ifelse(is_dmr, "Discard Mortality-at-Age", "Natural Mortality-at-Age"))+
        theme_bw()

    if(dimensions$nregions > 1){
        plot <- plot + facet_wrap(~region, scales="free_y")
    }

    return(plot)
}

#' Plot Selectivity-at-Age and Retention-at-Age
#'
#' Generates plot of selectivity-at-age or retention-at-age
#' across time, age, sex, region, and fleet as appropriate.
#'
#' @param selret five-dimensional selectivity-at-age or
#' retention-at-age matrix
#'
#' @export
#'
plot_selret <- function(selret, is_selectivity=TRUE){

    dimensions <- get_model_dimensions(selret)

    selex_df <- reshape2::melt(selret) %>%
        group_by(age, sex, region, fleet) %>%
        distinct(value, .keep_all=TRUE) %>%
        mutate(time_block = factor(time, labels=c(1:length(unique(time)))))

    ymax <- 1

    plot <- ggplot(selex_df, aes(x=age, y=value, color=sex, linetype=time_block))+
        geom_line(linewidth=0.8)+
        scale_y_continuous(limits=c(0, ymax), expand=c(0.01, 0.01))+
        scale_x_continuous(expand=c(0, 0))+
        labs(x="Age", y=ifelse(is_selectivity, "Selectivity", "Retention"), color="Sex", linetype="Time Block", title=ifelse(is_selectivity, "Selectivity", "Retention"))+
        theme_bw()

    if(dimensions$nregions > 1 & dimensions$nfleets <= 1){
        plot <- plot + facet_wrap(~region, scales="free_y")
    }else if(dimensions$nfleets > 1 & dimensions$nregions <= 1){
        plot <- plot + facet_wrap(~fleet, scales="free_y")
    }else if(dimensions$nfleets > 1 & dimensions$nregions > 1){
        plot <- plot + facet_grid(rows=vars(region), cols=vars(fleet))
    }

    return(plot)
}


plot_ssb <- function(ssb, comparison=NA){

    nregions <- ncol(ssb)

    ssb_df <- ssb %>% as_tibble() %>%
        rownames_to_column("time") %>%
        mutate(time=as.numeric(time)) %>%
        pivot_longer(c(2:(nregions+1)), names_to="region", values_to="ssb")

    xmax <- round(ssb_df %>% pull(time) %>% max, -1)
    ymax <- round(1.2*ssb_df %>% pull(ssb) %>% max, -1)

    plot <- ggplot(ssb_df, aes(x=time, y=ssb))+
        geom_line(linewidth=1)+
        scale_y_continuous(limits=c(0, ymax))+
        labs(x="Time", y="Spawning Biomass", title="Spawning Stock Biomass")+
        coord_cartesian(expand=0)+
        theme_bw()+
        theme(
            panel.grid.minor = element_blank()
        )

    if(nregions > 1){
        plot <- plot + facet_wrap(~region)
    }

    return(plot)

 }

 plot_bio <- function(bio){

    nregions <- ncol(bio)

    bio_df <- bio %>% as_tibble() %>%
        rownames_to_column("time") %>%
        mutate(time=as.numeric(time)) %>%
        pivot_longer(c(2:(nregions+1)), names_to="region", values_to="bio")

    xmax <- round(bio_df %>% pull(time) %>% max, -1)
    ymax <- round(1.2*bio_df %>% pull(bio) %>% max, -1)

    plot <- ggplot(bio_df, aes(x=time, y=bio))+
        geom_line(linewidth=1)+
        # scale_x_continuous(breaks=seq(0, xmax, length.out=6))+
        scale_y_continuous(limits=c(0, ymax))+
        labs(x="Time", y="Total Biomass", title="Total Biomass")+
        coord_cartesian(expand=0)+
        theme_bw()+
        theme(
            panel.grid.minor = element_blank()
        )

    if(nregions > 1){
        plot <- plot + facet_wrap(~region)
    }

    return(plot)

 }

 plot_catch <- function(catch){

    nregions <- ncol(catch)

    catch_df <- catch %>% as_tibble() %>%
        rownames_to_column("time") %>%
        mutate(time=as.numeric(time)) %>%
        pivot_longer(c(2:(nregions+1)), names_to="region", values_to="catch")

    xmax <- round(catch_df %>% pull(time) %>% max, -1)
    ymax <- round(1.2*catch_df %>% pull(catch) %>% max, -1)

    plot <- ggplot(catch_df, aes(x=time, y=catch))+
        geom_line(linewidth=1)+
        # scale_x_continuous(breaks=seq(0, xmax, length.out=6))+
        scale_y_continuous(limits=c(0, ymax))+
        labs(x="Time", y="Catch", title="Total Landed Catch")+
        coord_cartesian(expand=0)+
        theme_bw()+
        theme(
            panel.grid.minor = element_blank()
        )

    if(nregions > 1){
        plot <- plot + facet_wrap(~region)
    }

    return(plot)

 }

 plot_f <- function(f){

    nregions <- ncol(f)

    f_df <- f %>% as_tibble() %>%
        rownames_to_column("time") %>%
        mutate(time=as.numeric(time)) %>%
        pivot_longer(c(2:(nregions+1)), names_to="region", values_to="f")

    xmax <- round(f_df %>% pull(time) %>% max, -1)
    ymax <- 1.2*round(f_df %>% pull(f) %>% max, 2)

    plot <- ggplot(f_df, aes(x=time, y=f))+
        geom_line(linewidth=1)+
        # scale_x_continuous(breaks=seq(0, xmax, length.out=6))+
        scale_y_continuous(limits=c(0, ymax))+
        labs(x="Time", y="Fishing Mortality Rate", title="Fishing Mortality")+
        coord_cartesian(expand=0)+
        theme_bw()+
        theme(
            panel.grid.minor = element_blank()
        )

    if(nregions > 1){
        plot <- plot + facet_wrap(~region)
    }

    return(plot)

}


plot_atage <- function(atage){
    atage_df <- reshape2::melt(atage) %>% as_tibble() %>%
        group_by(time, age) %>%
        summarise(naa = sum(value)) %>%
        group_by(time) %>%
        mutate(total_naa = sum(naa)) %>%
        mutate(prop = naa/total_naa) %>%
        mutate(avg_age = weighted.mean(age, prop))


    plot <- ggplot(atage_df)+
        geom_point(aes(x=time, y=age, size=prop, color=prop))+
        geom_line(
            data=atage_df %>% distinct(time, avg_age),
            aes(x=time, y=avg_age),
            linewidth=2,
            color='red',
            show.legend = FALSE
        )+
        coord_cartesian(expand=0.01)+
        theme_bw()

    return(plot)
}






