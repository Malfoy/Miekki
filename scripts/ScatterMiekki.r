library("ggplot2")
library("ggthemes")
library("ggpmisc")
library("ggpubr")

args = commandArgs(trailingOnly=TRUE)




lm_eqn <- function(df){
    m <- lm(df$V2 ~ df$V1, df);
    eq <- substitute(italic("y") == a + b %.% italic("x")*","~~italic(r)^2~"="~r2,
         list(a = format(coef(m)[1], digits = 10),
              b = format(coef(m)[2], digits = 10),
             r2 = format(summary(m)$r.squared, digits = 3))
             )
    as.character(as.expression(eq));
}



ggplot_alternative <- function()
{
df <- read.table(args[1],header=FALSE)
p<-ggplot(df, aes(x=(df$V1), y=df$V2)) + geom_point(alpha=1,color=1) + geom_smooth(method=lm,formula = y ~ x) +
stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) + geom_abline(intercept = 0, slope = 1)+  theme_light()

#~  +stat_density_2d(aes(fill = ..level..), geom = "polygon")

p+labs(x = "Real",y="Estimated") + ggtitle( paste(nrow(df)," points"))

p2<-ggplot(df, aes(x=(df$V3), y=df$V4)) + geom_point(alpha=1,color=1) + geom_smooth(method=lm,formula = y ~ x) +
stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) + geom_abline(intercept = 0, slope = 1)+  theme_light()

#~  +stat_density_2d(aes(fill = ..level..), geom = "polygon")

p2+labs(x = "Real",y="Estimated") + ggtitle( paste(nrow(df)," points"))

figure <- ggarrange(p, p2,
                    labels = c("jaccard", "Intersection"),
                    ncol = 2, nrow = 1)

}



ggsave(
  paste(args[1],".png"),
  ggplot_alternative(),
  device="png",
  width = 20,
  height = 10,
  dpi = 1200
)

#~ ggsave(
#~   paste(args[1],"intersection.pdf"),
#~   ggplot_alternative2(),
#~   device="pdf",
#~   width = 10,
#~   height = 10,
#~   dpi = 1200
#~ )
