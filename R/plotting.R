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
plot_demographic_parameters <- function(dem_params, params=NA, out_dir=NA, ...){

    waa_plot  <- plot_waa(dem_params$waa)
    mat_plot  <- plot_mat(dem_params$mat)
    mort_plot <- plot_mort(dem_params$mort, is_dmr=FALSE)
    dmr_plot  <- plot_mort(dem_params$dmr,  is_dmr=FALSE)
    sel_plot  <- plot_selret(dem_params$sel, is_selectivity = TRUE)
    ret_plot  <- plot_selret(dem_params$ret, is_selectivity = FALSE)
    survsel_plot <- plot_selret(dem_params$surv_sel, is_selectivity = TRUE)

    if(!is.na(out_dir)){
        ggplot2::ggsave(filename=file.path(out_dir, "waa_plot.png"), plot=waa_plot, ...)
        ggplot2::ggsave(filename=file.path(out_dir, "mat_plot.png"), plot=mat_plot, ...)
        ggplot2::ggsave(filename=file.path(out_dir, "mort_plot.png"), plot=mort_plot, ...)
        ggplot2::ggsave(filename=file.path(out_dir, "dmr_plot.png"), plot=dmr_plot, ...)
        ggplot2::ggsave(filename=file.path(out_dir, "sel_plot.png"), plot=sel_plot, ...)
        ggplot2::ggsave(filename=file.path(out_dir, "ret_plot.png"), plot=ret_plot, ...)
        ggplot2::ggsave(filename=file.path(out_dir, "survsel_plot.png"), plot=survsel_plot, ...)
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
        dplyr::group_by(age, sex, region) %>%
        dplyr::distinct(value, .keep_all=TRUE) %>%
        dplyr::mutate(time_block = factor(time, labels=c(1:length(unique(time)))))

    ymax <- round(1.2*waa_df %>% dplyr::pull(value) %>% max, 2)

    plot <- ggplot2::ggplot(waa_df, ggplot2::aes(x=age, y=value, color=sex, linetype=time_block))+
        ggplot2::geom_line(linewidth=1)+
        ggplot2::scale_y_continuous(limits=c(0, ymax))+
        ggplot2::coord_cartesian(expand=0)+
        ggplot2::labs(x="Age", y="Weight", color="Sex", linetype="Time Block", title="Weight-at-Age")+
        ggplot2::theme_bw()

    if(dimensions$nregions > 1){
        plot <- plot + ggplot2::facet_wrap(~region, scales="free_y")
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
        dplyr::group_by(age, sex, region) %>%
        dplyr::distinct(value, .keep_all=TRUE) %>%
        dplyr::mutate(time_block = factor(time, labels=c(1:length(unique(time))))) %>%
        dplyr::filter(sex == "F")

    ymax <- 1.0

    plot <- ggplot2::ggplot(mat_df, ggplot2::aes(x=age, y=value, color=sex, linetype=time_block))+
        ggplot2::geom_line(linewidth=1)+
        ggplot2::scale_y_continuous(limits=c(0, ymax))+
        ggplot2::coord_cartesian(expand=0)+
        ggplot2::labs(x="Age", y="Maturity", color="Sex", linetype="Time Block", title="Female Maturity-at-Age")+
        ggplot2::theme_bw()

    if(dimensions$nregions > 1){
        plot <- plot + ggplot2::facet_wrap(~region, scales="free_y")
    }

    return(plot)

}

#' Plot Mortality
#'
#' Generates plot of mortality-at-age across time, age, sex,
#' and region as appropriate. This function can be used for
#' both natural and discard mortality
#'
#' @param mort four-dimensional mortality-at-age matrix
#' @param is_dmr title to discard mortality, default: FALSE = natural mortality
#'
#' @export
#'
plot_mort <- function(mort, is_dmr=FALSE){

    dimensions <- get_model_dimensions(mort)

    mort_df <- reshape2::melt(mort) %>%
        dplyr::group_by(age, sex, region) %>%
        dplyr::distinct(value, .keep_all=TRUE) %>%
        dplyr::mutate(time_block = factor(time, labels=c(1:length(unique(time)))))

    ymax <- round(1.2*mort_df %>% dplyr::pull(value) %>% max, 2)

    plot <- ggplot2::ggplot(mort_df, ggplot2::aes(x=age, y=value, color=sex, linetype=time_block))+
        ggplot2::geom_line(linewidth=1)+
        ggplot2::scale_y_continuous(limits=c(0, ymax), expand=c(0.01, 0.01))+
        ggplot2::scale_x_continuous(expand=c(0, 0))+
        #ggplot2::coord_cartesian(expand=0)+
        ggplot2::labs(x="Age", y="Instanteous Mortality", color="Sex", linetype="Time Block", title=ifelse(is_dmr, "Discard Mortality-at-Age", "Natural Mortality-at-Age"))+
        ggplot2::theme_bw()

    if(dimensions$nregions > 1){
        plot <- plot + ggplot2::facet_wrap(~region, scales="free_y")
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
#' @param is_selectivity default: TRUE, otherwise retention is output
#' @export
#'
plot_selret <- function(selret, is_selectivity=TRUE){

    dimensions <- get_model_dimensions(selret)

    selex_df <- reshape2::melt(selret) %>%
        dplyr::group_by(age, sex, region, fleet) %>%
        dplyr::distinct(value, .keep_all=TRUE) %>%
        dplyr::mutate(time_block = factor(time, labels=c(1:length(unique(time)))))

    ymax <- 1

    plot <- ggplot2::ggplot(selex_df, ggplot2::aes(x=age, y=value, color=sex, linetype=time_block))+
        ggplot2::geom_line(linewidth=0.8)+
        ggplot2::scale_y_continuous(limits=c(0, ymax), expand=c(0.01, 0.01))+
        ggplot2::scale_x_continuous(expand=c(0, 0))+
        ggplot2::labs(x="Age", y=ifelse(is_selectivity, "Selectivity", "Retention"), color="Sex", linetype="Time Block", title=ifelse(is_selectivity, "Selectivity", "Retention"))+
      ggplot2::theme_bw()

    if(dimensions$nregions > 1 & dimensions$nfleets <= 1){
        plot <- plot + ggplot2::facet_wrap(~region, scales="free_y")
    }else if(dimensions$nfleets > 1 & dimensions$nregions <= 1){
        plot <- plot + ggplot2::facet_wrap(~fleet, scales="free_y")
    }else if(dimensions$nfleets > 1 & dimensions$nregions > 1){
        plot <- plot + ggplot2::facet_grid(rows=ggplot2::vars(region), cols=ggplot2::vars(fleet))
    }

    return(plot)
}

#' Plot Spawning stock biomass
#'
#'
#' @param ssb annual spawning biomass
#' @param compare_ts i don't know what this does - default:NULL
#' @export
#'
plot_ssb <- function(ssb, compare_ts=NULL){

    nregions <- ncol(ssb)

    ssb_df <- ssb %>%
      tibble::as_tibble() %>%
        tibble::rownames_to_column("time") %>%
        dplyr::mutate(time=as.numeric(time)) %>%
        tidyr::pivot_longer(c(2:(nregions+1)), names_to="region", values_to="ssb")

    if(!is.null(compare_ts)){
        ssb_df <- ssb_df %>%
          dplyr::left_join(
            compare_ts %>%
              tibble::as_tibble() %>%
                tibble::rownames_to_column("time") %>%
                dplyr::mutate(time=as.numeric(time)) %>%
                tidyr::pivot_longer(c(2:(nregions+1)), names_to="region", values_to="comp")
        )
    }

    xmax <- round(ssb_df %>% dplyr::pull(time) %>% max, -1)
    ymax <- round(1.2*ssb_df %>% dplyr::pull(ssb) %>% max, -1)

    plot <- ggplot2::ggplot(ssb_df, ggplot2::aes(x=time, y=ssb))+
        ggplot2::geom_line(linewidth=1)+
        ggplot2::scale_y_continuous(limits=c(0, ymax))+
        ggplot2::labs(x="Time", y="Spawning Biomass", title="Spawning Stock Biomass")+
        ggplot2::coord_cartesian(expand=0)+
        ggplot2::theme_bw()+
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank()
        )

    if(!is.null(compare_ts)){
        plot <- plot + ggplot2::geom_line(ggplot2::aes(y=comp), color="red")
    }

    if(nregions > 1){
        plot <- plot + ggplot2::facet_wrap(~region)
    }

    return(plot)

 }

#' Plot Total biomass
#'
#'
#' @param bio annual total biomass
#' @param compare_ts i don't know what this does - default:NULL
#' @export
#'
 plot_bio <- function(bio, compare_ts=NULL){

    nregions <- ncol(bio)

    bio_df <- bio %>%
      tibble::as_tibble() %>%
        tibble::rownames_to_column("time") %>%
        dplyr::mutate(time=as.numeric(time)) %>%
        tidyr::pivot_longer(c(2:(nregions+1)), names_to="region", values_to="bio")

    if(!is.null(compare_ts)){
        bio_df <- bio_df %>%
          dplyr::left_join(
            compare_ts %>%
              tibble::as_tibble() %>%
                tibble::rownames_to_column("time") %>%
                dplyr::mutate(time=as.numeric(time)) %>%
                tidyr::pivot_longer(c(2:(nregions+1)), names_to="region", values_to="comp")
        )
    }

    xmax <- round(bio_df %>% dplyr::pull(time) %>% max, -1)
    ymax <- round(1.2*bio_df %>% dplyr::pull(bio) %>% max, -1)

    plot <- ggplot2::ggplot(bio_df, ggplot2::aes(x=time, y=bio))+
        ggplot2::geom_line(linewidth=1)+
        # ggplot2::scale_x_continuous(breaks=seq(0, xmax, length.out=6))+
        ggplot2::scale_y_continuous(limits=c(0, ymax))+
        ggplot2::labs(x="Time", y="Total Biomass", title="Total Biomass")+
        ggplot2::expand_limits(x=0, y=0) +
        ggplot2::theme_bw()+
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank()
        )

    if(!is.null(compare_ts)){
        plot <- plot + ggplot2::geom_line(ggplot2::aes(y=comp), color="red")
    }

    if(nregions > 1){
        plot <- plot + ggplot2::facet_wrap(~region)
    }

    return(plot)

 }

 #' Plot Catch
 #'
 #'
 #' @param catch annual catch
 #' @param compare_ts i don't know what this does - default:NULL
 #' @export
 #'
 plot_catch <- function(catch, compare_ts=NULL){

    nregions <- ncol(catch)

    catch_df <- catch %>%
      tibble::as_tibble() %>%
        tibble::rownames_to_column("time") %>%
        dplyr::mutate(time=as.numeric(time)) %>%
        tidyr::pivot_longer(c(2:(nregions+1)), names_to="region", values_to="catch")

    if(!is.null(compare_ts)){
        catch_df <- catch_df %>%
          dplyr::left_join(
            compare_ts %>%
              tibble::as_tibble() %>%
                tibble::rownames_to_column("time") %>%
                dplyr::mutate(time=as.numeric(time)) %>%
                tidyr::pivot_longer(c(2:(nregions+1)), names_to="region", values_to="comp")
        )
    }

    xmax <- round(catch_df %>% dplyr::pull(time) %>% max, -1)
    ymax <- round(1.2*catch_df %>% dplyr::pull(catch) %>% max, -1)

    plot <- ggplot2::ggplot(catch_df, ggplot2::aes(x=time, y=catch))+
        ggplot2::geom_line(linewidth=1)+
        # ggplot2::scale_x_continuous(breaks=seq(0, xmax, length.out=6))+
        ggplot2::scale_y_continuous(limits=c(0, ymax))+
        ggplot2::labs(x="Time", y="Catch", title="Total Landed Catch")+
        ggplot2::coord_cartesian(expand=0)+
        ggplot2::theme_bw()+
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank()
        )

    if(!is.null(compare_ts)){
        plot <- plot + ggplot2::geom_line(ggplot2::aes(y=comp), color="red")
    }

    if(nregions > 1){
        plot <- plot + ggplot2::facet_wrap(~region)
    }

    return(plot)

 }

 #' Plot Fishing mortality biomass
 #'
 #'
 #' @param f fishing mortality
 #' @param compare_ts i don't know what this does - default:NULL
 #' @export
 #'
 plot_f <- function(f, compare_ts=NULL){

    nregions <- ncol(f)

    f_df <- f %>%
      tibble::as_tibble() %>%
        tibble::rownames_to_column("time") %>%
        dplyr::mutate(time=as.numeric(time)) %>%
        tidyr::pivot_longer(c(2:(nregions+1)), names_to="region", values_to="f")

    if(!is.null(compare_ts)){
        f_df <- f_df %>%
          dplyr::left_join(
            compare_ts %>%
              tibble::as_tibble() %>%
                tibble::rownames_to_column("time") %>%
                dplyr::mutate(time=as.numeric(time)) %>%
                tidyr::pivot_longer(c(2:(nregions+1)), names_to="region", values_to="comp")
        )
    }

    xmax <- round(f_df %>% dplyr::pull(time) %>% max, -1)
    ymax <- 1.2*round(f_df %>% dplyr::pull(f) %>% max, 2)

    plot <- ggplot2::ggplot(f_df, ggplot2::aes(x=time, y=f))+
        ggplot2::geom_line(linewidth=1)+
        # ggplot2::scale_x_continuous(breaks=seq(0, xmax, length.out=6))+
        ggplot2::scale_y_continuous(limits=c(0, ymax))+
        ggplot2::labs(x="Time", y="Fishing Mortality Rate", title="Fishing Mortality")+
        ggplot2::coord_cartesian(expand=0)+
        ggplot2::theme_bw()+
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank()
        )

    if(!is.null(compare_ts)){
        plot <- plot + ggplot2::geom_line(ggplot2::aes(y=comp), color="red")
    }

    if(nregions > 1){
        plot <- plot + ggplot2::facet_wrap(~region)
    }

    return(plot)

}

 #' Plot Numbers at age
 #'
 #'
 #' @param atage numbers at age
 #' @export
 #'
plot_atage <- function(atage){
    atage_df <- reshape2::melt(atage) %>%
      tibble::as_tibble() %>%
        dplyr::group_by(time, age) %>%
        dplyr::summarise(naa = sum(value)) %>%
        dplyr::group_by(time) %>%
        dplyr::mutate(total_naa = sum(naa),
                      prop = naa/total_naa,
                      avg_age = stats::weighted.mean(age, prop))


    plot <- ggplot2::ggplot(atage_df)+
        ggplot2::geom_point(ggplot2::aes(x=time, y=age, size=prop, color=prop))+
        ggplot2::geom_line(
            data=atage_df %>% dplyr::distinct(time, avg_age),
            ggplot2::aes(x=time, y=avg_age),
            linewidth=2,
            color='red',
            show.legend = FALSE
        )+
        ggplot2::coord_cartesian(expand=0.01)+
      ggplot2::theme_bw()

    return(plot)
}






