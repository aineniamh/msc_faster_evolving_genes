res <- read.csv('asymmetry_analysis.txt', header = T, sep=',')
res <- read.csv('fixed_asymmetry_analysis_file.txt', header = T, sep=',')


res = na.omit(res)
library(ggplot2)
res$Change <- res$Omega1/res$Ancestral
res$logChange<- log2(res$Change)
res$Change1<-res$Omega1-res$Omega0
res$ancchange<-log2(res$Ancestral/res$Ancestral)
ggplot(res, aes(mac_gib, Change))+geom_point()+geom_point(data = res, aes(x = mac_gib, y=Omega1, colour = Species))+theme(legend.position = "none",axis.text.x = element_blank())
5678

ggplot(res, aes(Ancestral, Omega1))+geom_point()+geom_line(data=res, aes(Ancestral, Omega1, group=mac_gib))

ggplot(res, aes(mac_gib, logChange))+theme_classic()+
  xlab("Gene tree")+
  geom_hline(yintercept = 0)+
  geom_point(data = res, aes(mac_gib, ancchange), size=2)+
  geom_point(aes(colour=res$Species), size=2)+
  ylab(expression("log"["2"]*"(Fold-change d"["N"]*"d"["S"]*")")) +
  geom_line(data = res, aes(mac_gib, logChange, group=mac_gib_species))+
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        legend.title = element_blank(),
        axis.text.x = element_blank())

ggplot(res, aes(mac_gib, Omega1))+
  #geom_point(data =res, aes(x =mac_gib_species, y= Change, colour = Species))+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_line(data = res, aes(mac_gib, Omega1, group=mac_gib))+
  theme_classic()+
  xlab("Gene tree")+
  ylab(expression("d"["N"]*"d"["S"]))+
  geom_point(data = res, aes(mac_gib, Ancestral), size = 4)+
  geom_jitter(aes(fill = Species), shape=21,height = 0, width = 0.4, size = 3)+
  #scale_y_continuous(limits = c(-5.5, +5.5))+
  #geom_point(data=res, aes(mac_gib, Omega0))+
  theme(legend.position = "bottom",
        text = element_text(size = 20),
        legend.title = element_blank(),
        axis.text.x = element_blank())
        #axis.ticks = element_blank())

wilcox.test(res$Ancestral, res$Omega1)


fisher.test()
res <- subset(res, Omega1 < 6)
res <- subset(res, Omega0 < 10)
