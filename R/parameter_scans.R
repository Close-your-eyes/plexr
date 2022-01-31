library(tidyverse)

pl5 <- function(aa,bb,cc,dd,gg, Conc) {
  dd + (aa - dd) / ((1 + (Conc / cc)^bb))^gg
}



#Std. Curve: FI = dd + (aa - dd) / ((1 + (Conc / cc)^-bb))^gg

dd <- 1
aa <- 1e5
cc <- 500
bb <- -1
gg <- 1
conc <- c(seq(0,     0.01,  0.0001),
          seq(0.01,  0.1,   0.001),
          seq(0.1,   1,     0.01),
          seq(1,     10,    0.1),
          seq(10,    100,   1),
          seq(100,   1000,  10),
          seq(1000,  10000, 100),
          seq(10000, 1e5,   1000))

gg_std <-
  ggplot() +
  scale_y_log10(label = fcexpr::sci10, breaks = c(1,10,1e3,1e5)) +
  scale_x_log10(label = fcexpr::sci10, breaks = c(0.001,0.01,0.1,1,10,1e3,1e5)) +
  theme_bw()

## vary dd
out <- do.call(rbind,lapply(c(1:10), function(dd) {
  data.frame(FI = sapply(conc, function(x) pl5(aa,bb,cc,dd,gg,x)), aa=aa, bb=bb, cc=cc, dd=dd, gg=gg, conc=conc)
})) %>%
  dplyr::mutate(dd = as.factor(dd))
dd_plot <- ggplot(out, aes(x = conc, y = FI, color = dd)) +
  geom_line() +
  scale_y_log10(label = fcexpr::sci10, breaks = c(1,10,1e3,1e5)) +
  scale_x_log10(label = fcexpr::sci10, breaks = c(0.001,0.01,0.1,1,10,1e3,1e5)) +
  theme_bw() +
  labs(title = paste0("aa=",aa, ", bb=", bb, ", cc=", cc, ", gg=",gg))

p_arg +
  geom_line(data = out, aes(x = conc, y = FI, color = dd))


## vary aa
out <- do.call(rbind,lapply(c(5e4,1e5,2e5,3e5,5e5,1e6), function(aa) {
  data.frame(FI = sapply(conc, function(x) pl5(aa,bb,cc,dd,gg,x)), aa=aa, bb=bb, cc=cc, dd=dd, gg=gg, conc=conc)
})) %>%
  dplyr::mutate(aa = as.factor(aa))
aa_plot <- ggplot(out, aes(x = conc, y = FI, color = aa)) +
  geom_line() +
  scale_y_log10(label = fcexpr::sci10, breaks = c(1,10,1e3,1e5)) +
  scale_x_log10(label = fcexpr::sci10, breaks = c(0.001,0.01,0.1,1,10,1e3,1e5)) +
  theme_bw() +
  labs(title = paste0("dd=",dd, ", bb=", bb, ", cc=", cc, ", gg=",gg))
aa_plot

## vary bb
#c(-5,-3,-2,-1,-0.5,-0.25,-0.1,1,2,10)
out <- do.call(rbind,lapply(c(-5,-3,-2,-1.5,-1), function(bb) {
  data.frame(FI = sapply(conc, function(x) pl5(aa,bb,cc,dd,gg,x)), aa=aa, bb=bb, cc=cc, dd=dd, gg=gg, conc=conc)
})) %>%
  dplyr::mutate(bb = as.factor(bb))
bb_plot <- ggplot(out, aes(x = conc, y = FI, color = bb)) +
  geom_line() +
  scale_y_log10(label = fcexpr::sci10, breaks = c(1,10,1e3,1e5)) +
  scale_x_log10(label = fcexpr::sci10, breaks = c(0.001,0.01,0.1,1,10,1e3,1e5)) +
  theme_bw() +
  labs(title = paste0("dd=",dd, ", aa=", aa, ", cc=", cc, ", gg=",gg))
bb_plot

## vary cc
out <- do.call(rbind,lapply(, function(cc) {
  data.frame(FI = sapply(conc, function(x) pl5(aa,bb,cc,dd,gg,x)), aa=aa, bb=bb, cc=cc, dd=dd, gg=gg, conc=conc)
})) %>%
  dplyr::mutate(cc = as.factor(cc))
cc_plot <- ggplot(out, aes(x = conc, y = FI, color = cc)) +
  geom_line() +
  scale_y_log10(label = fcexpr::sci10, breaks = c(1,10,1e3,1e5)) +
  scale_x_log10(label = fcexpr::sci10, breaks = c(0.001,0.01,0.1,1,10,1e3,1e5)) +
  theme_bw() +
  labs(title = paste0("dd=",dd, ", aa=", aa, ", bb=", bb, ", gg=",gg))
cc_plot

## vary gg
out <- do.call(rbind,lapply(c(1,1.5,2,2.5,3,5,10), function(gg) {
  data.frame(FI = sapply(conc, function(x) pl5(aa,bb,cc,dd,gg,x)), aa=aa, bb=bb, cc=cc, dd=dd, gg=gg, conc=conc)
})) %>%
  dplyr::mutate(gg = as.factor(gg))
gg_plot <- ggplot(out, aes(x = conc, y = FI, color = gg)) +
  geom_line() +
  scale_y_log10(label = fcexpr::sci10, breaks = c(1,10,1e3,1e5)) +
  scale_x_log10(label = fcexpr::sci10, breaks = c(0.001,0.01,0.1,1,10,1e3,1e5)) +
  theme_bw() +
  labs(title = paste0("dd=",dd, ", aa=", aa, ", bb=", bb, ", cc=",cc))
gg_plot


## vary gg + bb + cc
out <- do.call(rbind,lapply(c(10,50,100,200,500,1000,1e4), function(cc) {
  do.call(rbind,lapply(c(1,1.5,2,2.5,3,5,10), function(gg) {
    do.call(rbind,lapply(c(-5,-3,-2,-1.5,-1), function(bb) {
      data.frame(FI = sapply(conc, function(x) pl5(aa,bb,cc,dd,gg,x)), aa=aa, bb=bb, cc=cc, dd=dd, gg=gg, conc=conc)
    }))
  }))
})) %>%
  dplyr::mutate(bb = as.factor(bb)) %>%
  dplyr::mutate(gg = as.factor(gg)) %>%
  dplyr::mutate(cc = as.factor(cc))

gg_bb_cc_plot1 <-
  ggplot(out, aes(x = conc, y = FI, color = gg)) +
  geom_line() +
  scale_y_log10(label = fcexpr::sci10, breaks = c(1,10,1e3,1e5)) +
  scale_x_log10(label = fcexpr::sci10, breaks = c(0.001,0.1,1,10,1e3,1e5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(title = paste0("dd=",dd, ", aa=", aa)) +
  facet_grid(rows = vars(bb), cols = vars(cc))
gg_bb_cc_plot1

gg_bb_cc_plot2 <-
  ggplot(out, aes(x = conc, y = FI, color = bb)) +
  geom_line() +
  scale_y_log10(label = fcexpr::sci10, breaks = c(1,10,1e3,1e5)) +
  scale_x_log10(label = fcexpr::sci10, breaks = c(0.001,0.1,1,10,1e3,1e5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(title = paste0("dd=",dd, ", aa=", aa)) +
  facet_grid(rows = vars(gg), cols = vars(cc))
gg_bb_cc_plot2


gg_bb_cc_plot3 <-
  ggplot(out, aes(x = conc, y = FI, color = cc)) +
  geom_line() +
  scale_y_log10(label = fcexpr::sci10, breaks = c(1,10,1e3,1e5)) +
  scale_x_log10(label = fcexpr::sci10, breaks = c(0.001,0.1,1,10,1e3,1e5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(title = paste0("dd=",dd, ", aa=", aa)) +
  facet_grid(rows = vars(gg), cols = vars(bb))
gg_bb_cc_plot3









