load("~/Documents/10T_K36me2/data/genebody.distribution.mat.gz/cpm/gene_centric/merged/K36me2/tables/PA.K36me2.gene_centric.table.RData")
load("~/Documents/10T_K36me2/data/genebody.distribution.mat.gz/cpm/gene_centric/merged/K36me2/tables/K36M.K36me2.gene_centric.table.RData")
load("~/Documents/10T_K36me2/data/genebody.distribution.mat.gz/cpm/gene_centric/merged/K36me2/tables/SETD2KO.K36me2.gene_centric.table.RData")
load("~/Documents/10T_K36me2/data/genebody.distribution.mat.gz/cpm/gene_centric/merged/K36me2/tables/NSD12DKO.K36me2.gene_centric.table.RData")
#
PA.K36me2.gene_centric.table$condition <- "PA"
SETD2KO.K36me2.gene_centric.table$condition <- "SETD2KO"
K36M.K36me2.gene_centric.table$condition <- "K36M-OE"
NSD12DKO.K36me2.gene_centric.table$condition <- "DKO"
#
ag <- rbind(PA.K36me2.gene_centric.table,SETD2KO.K36me2.gene_centric.table,K36M.K36me2.gene_centric.table,NSD12DKO.K36me2.gene_centric.table)
#
ag$condition <- factor(ag$condition,levels = c("PA","SETD2KO","DKO","K36M-OE"))
#
ag$quartile <- factor(x=ag$quartile,levels=c("4","3","2","1","zero"))

ggplot(data=ag,aes(x = x, y = v, color = quartile)) +
  geom_point(show.legend = FALSE) +
  labs(x="",y="",title = "") +
  guides(color=guide_legend("Expression\nquantile"),size=9,family="Helvetica") +
  scale_color_discrete(labels=c("4","3","2","1","0")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=7,family="Helvetica",colour = "black"),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "blue",size=9,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
        size = 9, color = "black",family = "Helvetica"),
    legend.text=element_text(size=8,family = "Helvetica",color = "black"),
    legend.title=element_text(size=8,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0,'cm'),
    panel.spacing.x = unit(0.1,'cm'),
    legend.position="bottom",
        legend.justification="center",
        legend.box.spacing = unit(-0.001, "cm")) +
  geom_vline(xintercept=c(200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600),linetype="dashed",alpha=0.3) +
  geom_vline(xintercept=c(500,900,1300,1700,2100,2500),linetype="solid",alpha=0.1,size=4.5,color="blue") +
  facet_wrap(~ condition, scales = "free",labeller = labeller(condition = 
    c("PA" = "PA",
      "DKO" = "NSD1/2-DKO",
      "TKO" = "NSD1/2-SETD2-TKO",
      "K36M-OE" = "K36M-OE",
      "SETD2KO" = "SETD2-KO")),ncol = 1,nrow=4)
ggsave(filename = "K36me2.gene_body.plots.png",path = "/Users/padilr1/Documents/10T_K36me_Manuscript/images_not_pushed/K36me2/V2",width = 7,height=17,units="cm",dpi = 600,device = "png")