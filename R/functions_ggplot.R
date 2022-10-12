UniquePanelCoords <- ggplot2::ggproto(
  "UniquePanelCoords", ggplot2::CoordCartesian,

  num_of_panels = 1,
  panel_counter = 1,
  panel_ranges = NULL,

  setup_layout = function(self, layout, params) {
    self$num_of_panels <- length(unique(layout$PANEL))
    self$panel_counter <- 1
    layout
  },

  setup_panel_params =  function(self, scale_x, scale_y, params = list()) {
    if (!is.null(self$panel_ranges) & length(self$panel_ranges) != self$num_of_panels)
      stop("Number of panel ranges does not equal the number supplied")

    train_cartesian <- function(scale, limits, name, given_range = NULL) {
      if (is.null(given_range)) {
        expansion <- ggplot2:::default_expansion(scale, expand = self$expand)
        range <- ggplot2:::expand_limits_scale(scale, expansion,
                                               coord_limits = self$limits[[name]])
      } else {
        range <- given_range
      }

      out <- list(
        ggplot2:::view_scale_primary(scale, limits, range),
        sec = ggplot2:::view_scale_secondary(scale, limits, range),
        arrange = scale$axis_order(),
        range = range
      )
      names(out) <- c(name, paste0(name, ".", names(out)[-1]))
      out
    }

    cur_panel_ranges <- self$panel_ranges[[self$panel_counter]]
    if (self$panel_counter < self$num_of_panels)
      self$panel_counter <- self$panel_counter + 1
    else
      self$panel_counter <- 1

    c(train_cartesian(scale_x, self$limits$x, "x", cur_panel_ranges$x),
      train_cartesian(scale_y, self$limits$y, "y", cur_panel_ranges$y))
  }
)

coord_panel_ranges <- function(panel_ranges, expand = TRUE, default = FALSE, clip = "on") {
  ggplot2::ggproto(NULL, UniquePanelCoords, panel_ranges = panel_ranges,
                   expand = expand, default = default, clip = clip)
}
# plot
plt.res <- function(bias, variance)
{
  library(latex2exp)
  library(ggplot2)
  bias <- abs(bias)
  mse <- variance + bias^2
  library(tidyr)
  est_bias_m <- reshape2::melt(bias, varnames = c('b', 'param'))%>%
    cbind(name = 'bias')%>%separate(param, c('param', 'dim'), sep = '[.]')

  est_var_m <- reshape2::melt(variance, varnames = c('b', 'param'))%>%
    cbind(name = 'variance')%>%separate(param, c('param', 'dim'), sep = '[.]')

  est_mse_m <- reshape2::melt(mse, varnames = c('b', 'param'))%>%
    cbind(name = 'mse')%>%separate(param, c('param', 'dim'), sep = '[.]')

  est_all_m <- rbind(est_bias_m, est_var_m, est_mse_m)
  xy.labs <- c(c('psi[2]',
                 'psi[3]'),
               c('bias', 'var', 'mse'))
  names(xy.labs) <- c('2', '3',
                      'bias', 'variance', 'mse')

  library(ggplot2)
  # delete the elast 0
  # est_all_m <- subset(est_all_m, subset = dim!=0)
  est_all_m$param <- factor(est_all_m$param, levels = c('AIPW', 'RT', 'EE', 'ELAS'))
  # begin our plots
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  est_all_m <- est_all_m[!est_all_m$dim=='tau',]
  p1 <-
    ggplot(est_all_m, aes(x=b, y=value,
                          color = param,
                          shape = param))+
    geom_line(size = 0.8)+geom_point(stroke = 2)+
    facet_wrap(c('dim','name'), scales = "free",
               labeller =as_labeller(xy.labs,
                                     default = label_parsed))+theme(legend.position = 'top')+
    ylab('')+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text = element_text(size = 10),
          text = element_text(size = 12),
          legend.text=element_text(size=12))+
    coord_panel_ranges(panel_ranges = list(
      list(y=c(-0.01,0.2)),
      list(y=c(0,0.022)),
      list(y=c(0,0.04)),

      list(y=c(-0.01,0.2)),
      list(y=c(0,0.022)),
      list(y=c(0,0.04))))+
    scale_shape_discrete(solid = T,
                         name = '', label = c('aipw',
                                              'rt', 'eff',
                                              'elas'))+
    scale_color_manual(values = c('#67D5B5',
                                  '#EE7785',
                                  '#C89EC4',
                                  '#84B1ED'),
                       name = '', label = c('aipw',
                                            'rt', 'eff',
                                            'elas'))+
    # c('#67D5B5',
    #   '#EE7785',
    #   '#C89EC4',
    #   '#84B1ED')

    theme(plot.margin=unit(c(0.3,0.3,0,0),"cm"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))

  p2 <- ggplot(est_all_m, aes(x=b, y=value,
                              color = param,
                              shape = param))+
    geom_line()+geom_point(stroke = .25)+
    theme(axis.text = element_text(size = 12),
          text = element_text(size = 15),
          legend.text=element_text(size=12))+
    facet_grid(c('dim','name'), scales = "free",
               labeller =as_labeller(xy.labs,
                                     default = label_parsed))+theme(legend.position = 'top')

  library(grid)
  library(gtable)
  gt1 <-  ggplot_gtable(ggplot_build(p1))
  gt2 <-  ggplot_gtable(ggplot_build(p2))
  gt1$grobs[grep('strip-t.+1$', gt1$layout$name)] = gt2$grobs[grep('strip-t', gt2$layout$name)]

  gt.side1 = gtable_filter(gt2, 'strip-r-1')
  gt.side2 = gtable_filter(gt2, 'strip-r-2')
  gt.side3 = gtable_filter(gt2, 'strip-r-3')

  gt1 = gtable_add_cols(gt1, widths=gt.side1$widths[1], pos = -1)
  gt1 = gtable_add_grob(gt1, zeroGrob(), t = 1, l = ncol(gt1), b=nrow(gt1))

  panel_id <- gt1$layout[grep('panel-.+1$', gt1$layout$name),]
  gt1 = gtable_add_grob(gt1, gt.side1, t = panel_id$t[1], l = ncol(gt1))
  gt1 = gtable_add_grob(gt1, gt.side2, t = panel_id$t[2], l = ncol(gt1))
  gt1 = gtable_add_grob(gt1, gt.side3, t = panel_id$t[3], l = ncol(gt1))


  gt1 <- gtable_add_cols(gt1, widths = unit(0.15, 'cm'), pos = -1)
  grid.newpage()
  grid.draw(gt1)

}


#' ggplot for plotting
#' @description
#' `plt.res.combine()` reproduce the plots in Yang et al., (2022), Figure 4
#' @export
plt.res.combine <- function(bias.psi0, bias.psi1,
                            variance.psi0, variance.psi1,
                            AIPW = T)
{
  library(latex2exp)
  library(ggplot2)
  bias.psi0 <- abs(bias.psi0);bias.psi1 <- abs(bias.psi1)
  mse.psi0 <- variance.psi0 + bias.psi0^2
  mse.psi1 <- variance.psi1 + bias.psi1^2
  library(tidyr)

  est_bias_m0 <- reshape2::melt(bias.psi0, varnames = c('b', 'param'))%>%
    cbind(name = 'bias')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.0'))

  est_bias_m1 <- reshape2::melt(bias.psi1, varnames = c('b', 'param'))%>%
    cbind(name = 'bias')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.1'))

  est_var_m0 <- reshape2::melt(variance.psi0, varnames = c('b', 'param'))%>%
    cbind(name = 'variance')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.0'))

  est_var_m1 <- reshape2::melt(variance.psi1, varnames = c('b', 'param'))%>%
    cbind(name = 'variance')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.1'))

  est_mse_m0 <- reshape2::melt(mse.psi0, varnames = c('b', 'param'))%>%
    cbind(name = 'mse')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.0'))

  est_mse_m1 <- reshape2::melt(mse.psi1, varnames = c('b', 'param'))%>%
    cbind(name = 'mse')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.1'))

  est_all_m <- rbind(est_bias_m0, est_bias_m1,
                     est_var_m0, est_var_m1,
                     est_mse_m0, est_mse_m1)
  xy.labs <- c(c('psi[1]==0',
                 'psi[2]==0',
                 'psi[1]==1',
                 'psi[2]== 1'),
               c('bias', 'var', 'MSE'))
  names(xy.labs) <- c('2.0', '3.0',
                      '2.1', '3.1',
                      'bias', 'variance', 'mse')

  library(ggplot2)
  # delete the elast 0
  if(AIPW == F)
  {
    line.color <- c('#EE7785',
                    '#C89EC4',
                    '#84B1ED')
    legend.label <- c('RT', 'Eff', 'Elastic')
    est_all_m <- est_all_m[est_all_m$param!='AIPW',]
    # est_all_m <- subset(est_all_m, subset = dim!=0)
    est_all_m$param <- factor(est_all_m$param, levels = c('RT', 'EE','ELAS'),
                              labels = legend.label)
  }else
  {
    line.color <- c('#EE7785',
                    '#67D5B5',
                    '#C89EC4',
                    '#84B1ED')
    legend.label <- c('RT.AIPW','RT.EE', 'Eff','Elas')
    est_all_m$param <- factor(est_all_m$param, levels = c('AIPW','RT', 'EE','ELAS'),
                              labels = legend.label)
  }

  # begin our plots
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  # est_all_m <- est_all_m[!est_all_m$dim=='tau',]
  est_all_m$name <- factor(est_all_m$name,
                           levels = c('bias', 'variance', 'mse'))
  est_all_m$dim.name <- factor(est_all_m$dim.name, levels = c('2.0','3.0',
                                                              '2.1','3.1'))
  p1 <-
    ggplot(est_all_m, aes(x=b, y=value,
                          color = param,
                          shape = param))+
    geom_line(size = 0.8)+geom_point(stroke = 2)+
    facet_wrap(~name+dim.name, scales = "free",
               labeller =as_labeller(xy.labs,
                                     default = label_parsed))+
    theme(legend.position = 'top')+
    scale_shape_discrete(solid = T,
                         name = '', label = legend.label)+
    scale_color_manual(values = line.color,
                       name = '', label = legend.label)+
    ylab('')+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text = element_text(size = 10),
          text = element_text(size = 12),
          legend.text=element_text(size=12))+
    coord_panel_ranges(panel_ranges = list(
      list(y=c(-0.01,0.15)),
      list(y=c(-0.01,0.15)),
      list(y=c(-0.01,0.15)),
      list(y=c(-0.01,0.15)),

      list(y=c(0,0.022)),
      list(y=c(0,0.022)),
      list(y=c(0,0.022)),
      list(y=c(0,0.022)),

      list(y=c(0,0.03)),
      list(y=c(0,0.03)),
      list(y=c(0,0.03)),
      list(y=c(0,0.03))
    ))+
    # c('#67D5B5',
    #   '#EE7785',
    #   '#C89EC4',
    #   '#84B1ED')
    theme(plot.margin=unit(c(0.3,0.3,0,0),"cm"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))

  p2 <-
    ggplot(est_all_m, aes(x=b, y=value,
                          color = param,
                          shape = param))+
    geom_line()+geom_point(stroke = .25)+
    theme(axis.text = element_text(size = 12),
          text = element_text(size = 15),
          legend.text=element_text(size=12))+
    facet_grid(name~dim.name, scales = "free",
               labeller =as_labeller(xy.labs,
                                     default = label_parsed))+theme(legend.position = 'top')

  library(grid)
  library(gtable)
  gt1 <-  ggplot_gtable(ggplot_build(p1))
  gt2 <-  ggplot_gtable(ggplot_build(p2))
  gt1$grobs[grep('strip-t.+1$', gt1$layout$name)] = gt2$grobs[grep('strip-t', gt2$layout$name)]

  gt.side1 = gtable_filter(gt2, 'strip-r-1')
  gt.side2 = gtable_filter(gt2, 'strip-r-2')
  gt.side3 = gtable_filter(gt2, 'strip-r-3')

  gt1 = gtable_add_cols(gt1, widths=gt.side1$widths[1], pos = -1)
  gt1 = gtable_add_grob(gt1, zeroGrob(), t = 1, l = ncol(gt1), b=nrow(gt1))

  panel_id <- gt1$layout[grep('panel-.+1$', gt1$layout$name),]
  gt1 = gtable_add_grob(gt1, gt.side1, t = panel_id$t[1], l = ncol(gt1))
  gt1 = gtable_add_grob(gt1, gt.side2, t = panel_id$t[2], l = ncol(gt1))
  gt1 = gtable_add_grob(gt1, gt.side3, t = panel_id$t[3], l = ncol(gt1))


  gt1 <- gtable_add_cols(gt1, widths = unit(0.15, 'cm'), pos = -1)
  grid.newpage()
  grid.draw(gt1)


}
