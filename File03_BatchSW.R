################################################################################
#
# Code to accompany: "The batched stepped wedge design: a design robust to 
#                     delays in cluster recruitment" 
#                     by J Kasza, R Bowden, R Hooper, AB Forbes.
#
# Questions or comments to J Kasza, jessica.kasza@monash.edu
#
#
# BatchedSW File 3: code to generate the nested loop plots 
#
# 2022-04-11
################################################################################

library("looplot")
library("ggpubr")


#Continuous outcomes:
theoreticalpower <- read.csv(file="ContsCombos.csv") 
continuousresults <- read.csv(file="contsresults.csv")

resultssim <- cbind(theoreticalpower, continuousresults[,c(14)])
resultssim <- as.data.frame(resultssim)
names(resultssim)[1] <- "T"
names(resultssim)[2] <- "NBatch"
names(resultssim)[3] <- "Olap"
names(resultssim)[4] <- "K"
names(resultssim)[5] <- "m"
names(resultssim)[6] <- "SharedTime"
names(resultssim)[7] <- "EffectSize"
names(resultssim)[8] <- "icc"
names(resultssim)[9] <- "CAC"
names(resultssim)[10] <- "Power"
names(resultssim)[11] <- "RejPercentage"

resultssim$powerSE <- sqrt((resultssim$RejPercentage/1000)*(1-resultssim$RejPercentage/1000)/1000)
resultssim$powerminusSE <- resultssim$RejPercentage/1000 - 2*sqrt((resultssim$RejPercentage/1000)*(1-resultssim$RejPercentage/1000)/1000)
resultssim$powerplusSE <- resultssim$RejPercentage/1000 + 2*sqrt((resultssim$RejPercentage/1000)*(1-resultssim$RejPercentage/1000)/1000)
resultssim$powerp <- resultssim$RejPercentage/1000


resultssim$icc <- as.factor(resultssim$icc)
resultssim$CAC <- as.factor(resultssim$CAC)
names(resultssim) <- make.unique(names(resultssim))

resultssim015 <- resultssim[resultssim$SharedTime==1 & resultssim$EffectSize == 0.15 ,]
resultssim015 <- subset(resultssim015, select = -c(SharedTime,EffectSize, T, NBatch, K, m, powerSE) )
names(resultssim015) <- make.unique(names(resultssim015))

resultssim0 <- resultssim[resultssim$SharedTime==1 & resultssim$EffectSize == 0 ,]
resultssim0 <- subset(resultssim0, select = -c(SharedTime,EffectSize, T, NBatch, K, m, powerSE) )
names(resultssim0) <- make.unique(names(resultssim0))


NLP_conts015 <- nested_loop_plot(resdf =  resultssim015, 
                                 x = "icc", steps = c("Olap"), 
                                 grid_rows = "CAC",
                                 steps_y_base = 0.4, steps_y_height = 0.05, steps_y_shift = 0.05,
                                 methods = c("powerp", "powerminusSE", "powerplusSE", "Power"),
                                 line_linetypes = c(1, 1, 1, 2),
                                 colors = c("black", "black", "black", "darkred"),
                                 point_shapes = c("powerp" =16, "powerminusSE" = 20, "powerplusSE" = 20, "Power" =17),
                                 point_size = 2, line_size = 0.75,
                                 steps_values_annotate = TRUE, steps_annotation_size = 3,
                                 steps_annotation_nudge = 0.2,
                                 y_name = "Power",
                                 x_name = "ICC",
                                 y_expand_add = c(0.1, NULL), 
                                 legend_name = "Legend", steps_names=c("Overlap"),
                                 legend_labels=c("Simulated power", "Sim. power -2se", "Sim. power +2se", "Theoretical power"),
                                 post_processing = list(
                                   add_custom_theme = list(
                                     legend.position="bottom",
                                     axis.text.x = element_text(angle = -90, 
                                                                vjust = 0.5, 
                                                                size = 8) 
                                   )))

NLP_conts0 <- nested_loop_plot(resdf =  resultssim0, 
                               x = "icc", steps = c("Olap"), 
                               grid_rows = "CAC",
                               steps_y_base = -0.01, steps_y_height = 0.01, steps_y_shift = 0.01,
                               methods = c("powerp", "powerminusSE", "powerplusSE", "Power"),
                               line_linetypes = c(1, 1, 1, 2),
                               colors = c("black", "black", "black", "darkred"),
                               point_shapes = c("powerp" =16, "powerminusSE" = 20, "powerplusSE" = 20, "Power" =17),
                               point_size = 2, line_size = 0.75,
                               steps_values_annotate = TRUE, steps_annotation_size = 3,
                               steps_annotation_nudge = 0.2,
                               y_name = "Type I error rate",
                               x_name = "ICC",
                               hline_intercept = c(0,0.05), 
                               y_expand_add = c(0.02, NULL),
                               y_breaks=c(0, 0.025, 0.05, 0.075),
                               legend_name = "Legend", steps_names=c("Overlap"),
                               legend_labels=c("Simulated error", "Sim. error -2se", "Sim. error +2se", "Theoretical error"),
                               post_processing = list(
                                 add_custom_theme = list(
                                   legend.position="bottom",
                                   axis.text.x = element_text(angle = -90, 
                                                              vjust = 0.5, 
                                                              size = 8) 
                                 )))


contsboth <- ggarrange(NLP_conts0, NLP_conts015, ncol=2)
ggsave(contsboth, file= "NLP_contsALL.pdf", device = "pdf", width = 40, height = 20, units="cm")

ggsave(NLP_conts0, file= "NLP_conts0.pdf", device = "pdf", width = 17, height = 20, units="cm")
ggsave(NLP_conts015, file= "NLP_conts015.pdf", device = "pdf", width = 17, height = 20, units="cm")

####################################################################################
#Binary plots
binaryresults202204 <- read.csv(file="binaryresults_202204.csv")


resultssimbin <- as.data.frame(binaryresults202204)
names(resultssimbin)[2] <- "T"
names(resultssimbin)[3] <- "NBatch"
names(resultssimbin)[4] <- "Olap"
names(resultssimbin)[5] <- "K"
names(resultssimbin)[6] <- "m"
names(resultssimbin)[7] <- "SharedTime"
names(resultssimbin)[8] <- "EffectSize"
names(resultssimbin)[9] <- "ICC"
names(resultssimbin)[10] <- "CAC"
names(resultssimbin)[11] <- "Rej01"
names(resultssimbin)[12] <- "theorpower1"



resultssimbin$powerSE1 <- sqrt((resultssimbin$Rej01/1000)*(1-resultssimbin$Rej01/1000)/1000)
resultssimbin$powerminusSE1 <- resultssimbin$Rej01/1000 - 2*sqrt((resultssimbin$Rej01/1000)*(1-resultssimbin$Rej01/1000)/1000)
resultssimbin$powerplusSE1 <- resultssimbin$Rej01/1000 + 2*sqrt((resultssimbin$Rej01/1000)*(1-resultssimbin$Rej01/1000)/1000)
resultssimbin$powerp1 <- resultssimbin$Rej01/1000


resultssimbin$ICC <- as.factor(resultssimbin$ICC)
names(resultssimbin) <- make.unique(names(resultssimbin))

resultssimbin0025 <- resultssimbin[resultssimbin$EffectSize == 0.025 ,]
resultssimbin0 <- resultssimbin[resultssimbin$EffectSize == 0 ,]


NLP_bin0025_gee <- nested_loop_plot(resdf =  resultssimbin0025, 
                                    x = "ICC", steps = c("Olap"), 
                                    steps_y_base = 0.325, steps_y_height = 0.01, steps_y_shift = 0.01,
                                    methods = c("powerp1", "powerminusSE1", "powerplusSE1", "theorpower1"),
                                    line_linetypes = c(1, 1, 1, 2),
                                    colors = c("black", "black", "black", "darkred"),
                                    point_shapes = c("powerp1" =16, "powerminusSE1" = 20, "powerplusSE1" = 20, "theorpower1" =17),
                                    point_size = 2, line_size = 0.75,
                                    steps_values_annotate = TRUE, steps_annotation_size = 3,
                                    steps_annotation_nudge = 0.2,
                                    y_name = "Power",
                                    x_name = "ICC",
                                    legend_name = "Legend", steps_names=c("Overlap"),
                                    legend_labels=c("Simulated power", "Sim. power -2se", "Sim. power +2se", "Theoretical power"),
                                    post_processing = list(
                                      add_custom_theme = list(
                                        legend.position="bottom",
                                        axis.text.x = element_text(angle = -90, 
                                                                   vjust = 0.5, 
                                                                   size = 8) 
                                      )))



NLP_bin0_gee <- nested_loop_plot(resdf =  resultssimbin0, 
                                 x = "ICC", steps = c("Olap"), 
                                 steps_y_base = 0.01, steps_y_height = 0.01, steps_y_shift = 0.01,
                                 methods = c("powerp1", "powerminusSE1", "powerplusSE1", "theorpower1"),
                                 line_linetypes = c(1, 1, 1, 2),
                                 colors = c("black", "black", "black", "darkred"),
                                 point_shapes = c("powerp1" =16, "powerminusSE1" = 20, "powerplusSE1" = 20, "theorpower1" =17),
                                 point_size = 2, line_size = 0.75,
                                 steps_values_annotate = TRUE, steps_annotation_size = 3,
                                 steps_annotation_nudge = 0.2,
                                 y_name = "Type I error rate",
                                 x_name = "ICC",
                                 hline_intercept = 0.05, 
                                 y_expand_add = c(0.02, NULL), 
                                 legend_name = "Legend", steps_names=c("Overlap"),
                                 legend_labels=c("Simulated error", "Sim. error -2se", "Sim. error +2se", "Theoretical error"),
                                 post_processing = list(
                                   add_custom_theme = list(
                                     legend.position="bottom",
                                     axis.text.x = element_text(angle = -90, 
                                                                vjust = 0.5, 
                                                                size = 8) 
                                   )))


binboth_0025_gee <- ggarrange(NLP_bin0_gee, NLP_bin0025_gee, ncol=2)
ggsave(binboth_0025_gee, file= "NLP_bingee0025ALL.pdf", device = "pdf", width = 40, height = 10, units="cm")

