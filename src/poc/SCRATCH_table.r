# library(formattable)
library(tidyverse)
# library(knitr)
# library(kableExtra)
library(ggplot2)
library(ggthemes)
library(patchwork)

df <- data.frame(x = rnorm(10),
                 y = rnorm(10),
                 id = as.factor(1:10))

# https://stackoverflow.com/questions/67441457/change-bars-color-and-orientation-in-a-table
# cb <- function(x) {
#   range <- max(abs(x))
#   width <- round(abs(x / range * 50), 2)
#   ifelse(
#     x > 0,
#     paste0(
#       '<span style="display: inline-block; border-radius: 2px; ', 
#       'padding-right: 2px; background-color: lightgreen; width: ', 
#       width, '%; margin-left: 50%; text-align: left;">', x, '</span>'
#     ),
#     paste0(
#       '<span style="display: inline-block; border-radius: 2px; ', 
#       'padding-right: 2px; background-color: lightpink; width: ', 
#       width, '%; margin-right: 50%; text-align: right; float: right; ">', x, '</span>'
#     )
#   )
# }
# color_bar2 <- function(color = "lightgray", fun = proportion, ...) {
#   fun <- match.fun(fun)
#   formatter("span",
#             style = function(x) style(
#               display = "inline-block",
#               direction = ifelse(x > 0, "rtl", 
#                                  ifelse(x < 0, "ltr", "ltr")),
#               # direction = 
#               "unicode-bidi" = "plaintext",
#               "border-radius" = "4px",
#               "padding-right" = "2px",
#               # "background-color" = csscolor(color),
#               "text-align" = ifelse(x > 0, "right", 
#                                  ifelse(x < 0, "left", "left")),
#               "background-color" = ifelse(x > 0, "lightgreen", 
#                                        ifelse(x < 0, "lightpink", "black")),
#               width = percent(fun(as.numeric(x), ...))
#             ))
# }
# 
# 
# formattable(df)
# formattable(df, list(
#   x = color_bar2("lightgreen")
#   # y = color_bar("lightgreen"))
#   # market_share = color_bar("lightblue"),
#   # revenue = sign_formatter,
#   # profit = sign_formatter)
#   ))


# df %>%
#   mutate(
#     x = cb(round(x, 2)),
#     y = cb(round(y, 2))
#   ) %>%
#   kable(escape = F) %>%
#   kable_styling("hover", full_width = F) %>%
#   column_spec(1:2, width = "3cm") %>%
#   row_spec(0, align = "c")

pal <- c("#E66100", "#5D3A9B")

(p1 <- df %>%
    mutate(sign = case_when(x < 0 ~ "n",
                            x >= 0 ~ "p")) %>%
    ggplot()+
    geom_segment(aes(color = sign, x = 0, xend = x/10, y = id, yend = id),
                 size = 5, alpha = 0.5) +
    scale_color_manual(values = pal) +
    geom_text(aes(y = id, x = 0, label = round(x,2))) +
    ggtitle("Human Mobility") +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)))

(p2 <- df %>%
    mutate(sign = case_when(y < 0 ~ "n",
                            y >= 0 ~ "p")) %>%
    ggplot()+
    geom_segment(aes(color = sign, x = 0, xend = y, y = id, yend = id),
                 size = 5, alpha = 0.5) +
    scale_color_manual(values = pal) +
    geom_text(aes(y = id, x = 0, label = round(y,2))) +
    ggtitle("Human Modification")+
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)))

p1+p2