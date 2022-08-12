# This is an interactive script to output conditonal plots - 
# IT DOES NOT RUN ON ITS OW - TIERS TO THE PLOT SPACE USE SCRIPT

(sg_ce_plot <-  plot(ce_sg, plot = F,
                     line_args = list("se" = T,
                                      "color" = pal2[1],
                                      "fill" = pal2[1]))[[1]] + 
  # scale_color_manual(values = palnew[3])+         
  theme_tufte() +
  xlab("Human Mobility") +
  ylab("Space Use")+
  theme(axis.line = element_line(size = .5),
        # axis.text = element_blank(),
        axis.ticks = element_blank(),
        # axis.title = element_blank(),
        # aspect.ratio = 1
  ))

(ghm_ce_plot <-  plot(ce_ghm, plot = F,
                      line_args = list("se" = T,
                                       "color" = pal2[1],
                                       "fill" = pal2[1]))[[1]] +
    # scale_color_manual(values = palnew[3])+         
    theme_tufte() +
    xlab("Human Modification") +
    ylab("Space Use")+
    theme(axis.line = element_line(size = .5),
          # axis.text = element_blank(),
          axis.ticks = element_blank(),
          # axis.title = element_blank(),
          # aspect.ratio = 1
    ))
(int_ce_plot <-  plot(ce_int, plot = FALSE,
                      line_args = list("se"=T ,
                                       "alpha" = 0.2
                                       # size = 4
                                       ))[[1]] +
    theme_tufte() +
    scale_color_manual(values = pal2, name = "Human \n Modification",
                       labels = c("High", "Low")) +
    scale_fill_manual(values = pal2, name = "Human \n Modification",
                      labels = c("High", "Low")) +
    ylab("Space Use")+
    xlab("Human Mobility")+
    theme(
      axis.line = element_line(),
          # axis.text = element_blank(),
          axis.ticks = element_blank(),
          # axis.title = element_blank(),
          aspect.ratio = 1,
          # legend.position = "none"
    ))
