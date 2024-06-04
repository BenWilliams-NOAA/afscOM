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
#' Plots will show valyes of demographic quantities across
#' time, age, sex, region, and fleet as appropriate.
#'
#' @param dem_params list of demographic parameter matrices
#' @param params vector of names of parameters to plot
#' @param out_dir optional directory to save plots to
#' @param ... additional parameters to be passed to `ggsave`
#'
#' @export plot_demographic_parameters
#'
#' @example
#'
plot_demographic_parameters <- function(dem_params, params=NA, out_dir=NA, ...){
    waa_plot <- plot_waa(dem_params$waa)
    mat_plot <- plot_mat(dem_params$mat)
    mort_plot <- plot_mort(dem_params$mort)
    dmr_plot <- plot_mort(dem_params$dmr)
    sel_plot <- plot_selret(dem_params$sel)
    ret_plot <- plot_selret(dem_params$ret)
    survsel_plot <- plot_selret(dem_params$surv_sel)

    show(waa_plot)
    show(mat_plot)
    show(mort_plot)
    show(dmr_plot)
    show(sel_plot)
    show(ret_plot)
    show(survsel_plot)

    if(!is.na(out_dir)){
        ggsave(filename=file.path(out_dir, "waa_plot.png"), plot=waa_plot, ...)
        ggsave(filename=file.path(out_dir, "mat_plot.png"), plot=mat_plot, ...)
        ggsave(filename=file.path(out_dir, "mort_plot.png"), plot=mort_plot, ...)
        ggsave(filename=file.path(out_dir, "dmr_plot.png"), plot=dmr_plot, ...)
        ggsave(filename=file.path(out_dir, "sel_plot.png"), plot=sel_plot, ...)
        ggsave(filename=file.path(out_dir, "ret_plot.png"), plot=ret_plot, ...)
        ggsave(filename=file.path(out_dir, "survsel_plot.png"), plot=survsel_plot, ...)
    }

}

#' Plot Weight-at-Age
#' 
#' Generates plot of weight-at-age across time, age, sex,
#' and region as appropriate.
#'
#' @param waa four-dimenrionsal weight-at-age matrix
#'
#' @export plot_waa
#'
#' @example
#'
plot_waa <- function(waa){
    waa_df <- reshape2::melt(waa) %>%
        group_by(age, sex, region) %>%
        distinct(value, .keep_all=TRUE) %>%
        mutate(time_block = as.factor(time))

    plot <- ggplot(waa_df, aes(x=age, y=value, color=sex, linetype=time_block))+
        geom_line(linewidth=1)+
        scale_y_continuous(limits=c(0, 7))+
        coord_cartesian(expand=0)+
        labs(x="Age", y="Weight", color="Sex", linetype="Time Block")+
        theme_bw()+
        facet_wrap(~region, scales="free_y")

    return(plot)

}

#' Plot Maturity-at-Age
#' 
#' Generates plot of maturity-at-age across time, age, sex,
#' and region as appropriate. Only the female maturity ogive
#' is plotted.
#'
#' @param mat four-dimenrionsal maturity-at-age matrix
#'
#' @export plot_mat
#'
#' @example
#'
plot_mat <- function(mat){
    mat_df <- reshape2::melt(mat) %>%
        group_by(age, sex, region) %>%
        distinct(value, .keep_all=TRUE) %>%
        mutate(time_block = as.factor(time)) %>%
        filter(sex == "F")

    plot <- ggplot(mat_df, aes(x=age, y=value, color=sex, linetype=time_block))+
        geom_line(linewidth=1)+
        scale_y_continuous(limits=c(0, 1.0))+
        coord_cartesian(expand=0)+
        labs(x="Age", y="Maturity", color="Sex", linetype="Time Block")+
        theme_bw()+
        facet_wrap(~region, scales="free_y")

    return(plot)

}

#' Plot Mortality
#' 
#' Generates plot of mortality-at-age across time, age, sex,
#' and region as appropriate. This function can be used for 
#' both natural and discard mortality
#'
#' @param waa four-dimenionsal mortality-at-age matrix
#'
#' @export plot_mort
#'
#' @example
#'
plot_mort <- function(mort){
    mort_df <- reshape2::melt(mort) %>%
        group_by(age, sex, region) %>%
        distinct(value, .keep_all=TRUE) %>%
        mutate(time_block = as.factor(time))

    plot <- ggplot(mort_df, aes(x=age, y=value, color=sex, linetype=time_block))+
        geom_line(linewidth=1)+
        scale_y_continuous(limits=c(0, 1.0), expand=c(0.01, 0.01))+
        scale_x_continuous(expand=c(0, 0))+
        # coord_cartesian(expand=0)+
        labs(x="Age", y="Instanteous Mortality", color="Sex", linetype="Time Block")+
        theme_bw()+
        facet_wrap(~region, scales="free_y")

    return(plot)
}

#' Plot Selectivity-at-Age and Retention-at-Age
#' 
#' Generates plot of selectivity-at-age or retention-at-age
#' across time, age, sex, region, and fleet as appropriate.
#'
#' @param selret five-dimenrionsal selectivity-at-age or
#' retention-at-age matrix
#'
#' @export plot_selret
#'
#' @example
#'
plot_selret <- function(selret){
    selex_df <- reshape2::melt(selret) %>%
        group_by(age, sex, region, fleet) %>%
        distinct(value, .keep_all=TRUE) %>%
        mutate(time_block = as.factor(time))

    plot <- ggplot(selex_df, aes(x=age, y=value, color=sex, linetype=time_block))+
        geom_line(linewidth=1)+
        scale_y_continuous(limits=c(0, 1.0), expand=c(0.01, 0.01))+
        scale_x_continuous(expand=c(0, 0))+
        labs(x="Age", y="Selectivity", color="Sex", linetype="Time Block")+
        theme_bw()+
        facet_grid(rows=vars(region), cols=vars(fleet))

    return(plot)
}
