## set axis tick and title size in all ggplots (global theme)
theme_set(
  theme_minimal()+
    theme(
      axis.text.x  = element_text(size=12),
      axis.text.y  = element_text(size=12),
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14)
    )
)


summary_axis_theme <- theme_minimal()+
  theme(
    axis.text.x  = element_text(size=12),
    axis.title.x = element_text(size=14),
    axis.text.y  = element_text(size=11),
    axis.title.y = element_text(size=10)#,
    #plot.background = element_rect(fill="#FAF7F7", color=NA)
  )


summary_axis_theme2 <- theme_minimal()+
  theme(
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    axis.text.x = element_text(size=18, face="bold"),
    axis.text.y = element_text(size=18, face="bold"),
    plot.background = element_rect(fill="#FAF7F7", color=NA)
  )

comp_axis_theme <- theme_minimal()+
  theme(
    axis.text.x  = element_text(size=12),
    axis.title.x = element_text(size=14),
    axis.text.y  = element_text(size=12),
    axis.title.y = element_text(size=14)
  )


my_ggtheme <- function( font.size=1, font.face="bold", font.x.angle=90, font.y.angle=90,
                        title.x.angle=0, title.y.angle=90,
                        font.x.tit.size=1,
                        font.y.tit.size=1,
                        font.leg.size=1,
                        font.leg.tit.size=1,
                        font.tit.size=1,
                        leg.pos="bottom",
                        hor.x.adj=NULL, ver.x.adj=NULL, base_size=14){
  theme_bw(base_size=base_size) %+replace%
    theme(
      axis.text.x=element_text(angle=font.x.angle, 
                               face=font.face, 
                               size=rel(font.size),
                               hjust=hor.x.adj, 
                               vjust=ver.x.adj),
      axis.text.y=element_text(angle=font.y.angle, 
                               face=font.face, 
                               size=rel(font.size)),
      axis.title.x=element_text(face=font.face, size=rel(font.x.tit.size), angle=title.x.angle),
      axis.title.y=element_text(face=font.face, size=rel(font.y.tit.size), angle=title.y.angle),
      legend.text=element_text(face=font.face, size=rel(font.leg.size)),
      legend.title=element_text(face=font.face, size=rel(font.leg.tit.size)),
      title = element_text(face=font.face, size=rel(font.tit.size)),
      legend.position = leg.pos,
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(1.5, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      strip.background = element_rect(fill = "#17252D", color = "#17252D")#,
      #strip.text = element_text(size = rel(0.85), face = "bold", color = "white", margin = margin(5,0,5,0))
    )
  
}

# heatmap_ggtheme <- function(border.col="grey50"){
#   
#   theme(
#     #axis.title.x = element_blank(),
#     #axis.title.y = element_blank(),
#     axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95, face="bold"),
#     axis.text.y = element_text(face="bold"),
#     panel.grid.major = element_blank(),
#     panel.border = element_rect(colour = border.col, fill=NA, linewidth=1),
#     panel.background = element_blank(),
#     axis.line = element_line(colour = border.col),
#     axis.ticks = element_blank(),
#     legend.justification = c("right","center"),
#     #legend.margin=margin(grid::unit(0,"cm")),
#     legend.position = "right",
#     legend.direction = "vertical")
#   
# }

heatmap_ggtheme <- function(border.col="grey50", base_size=14, font.x.angle=90, 
                            font.x.size=0.7, font.y.size=0.7,
                            font.x.tit.size=1, font.y.tit.size=1, plot_title=1){
theme_bw(base_size=base_size) %+replace%
    theme(
    plot.title = element_text(size = rel(plot_title[1]), face = "bold", margin = margin(0,0,5,0), hjust = 0),
    axis.title = element_text(size = rel(0.85), face = "bold"),
    axis.text.x = element_text(angle = font.x.angle[1], size=rel(font.x.size[1]), vjust = 0.2, hjust = 0.95, face="bold"),
    axis.text.y = element_text(face="bold", size=rel(font.y.size[1])),
    axis.title.x=element_text(face="bold", size=rel(font.x.tit.size[1])),
    axis.title.y=element_text(face="bold", size=rel(font.y.tit.size[1]), angle=90),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = border.col[1], fill=NA, linewidth=1),
    panel.background = element_blank(),
    axis.line = element_line(colour = border.col[1]),
    axis.ticks = element_blank(),
    legend.justification = c("right","center"),
    #legend.margin=margin(grid::unit(0,"cm")),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = rel(0.85), face = "bold"),
    legend.text = element_text(size = rel(0.70), face = "bold"),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.size = unit(1.5, "lines"),
    legend.background = element_rect(fill = "transparent", colour = NA),
    strip.background = element_rect(fill = "#17252D", color = "#17252D"),
    strip.text = element_text(size = rel(0.85), face = "bold", color = "white", margin = margin(5,0,5,0))
    )
  
}

