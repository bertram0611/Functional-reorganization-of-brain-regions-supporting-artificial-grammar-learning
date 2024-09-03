library(tidyverse) # for data manipulation and visualization
library(tidyr) # can use gather to convert wide to long format
library(dplyr) # can use pipe
FC <- read_csv("FC_PreRest_Learn.csv", col_names = FALSE) # This csv file was produced by FC_analysis_Fig4.m
newheaders <- c("PreRest_Ch2_9", "Learn_Ch2_9", "PreRest_Ch9_30", "Learn_Ch9_30","PreRest_Ch9_34","Learn_Ch9_34","PreRest_Ch9_35","Learn_Ch9_35","PreRest_Ch17_38","Learn_Ch17_38","PreRest_Ch39_40","Learn_Ch39_40","PreRest_Ch2_43","Learn_Ch2_43","PreRest_Ch17_43","Learn_Ch17_43","PreRest_Ch19_43","Learn_Ch19_43","PreRest_Ch22_43","Learn_Ch22_43")
colnames(FC) <- newheaders 

FC <- FC %>%
  mutate(ID=c(1:15),
         .before=PreRest_Ch2_9)

FC_long <- FC %>% gather(phase,rvalue,PreRest_Ch2_9:Learn_Ch22_43, factor_key = TRUE)
levels(FC_long$phase)

library(lattice)
# source https://publicifsv.sund.ku.dk/~jufo/courses/rm2017/plotRrepeated.pdf
# https://stackoverflow.com/questions/59693411/r-boxplot-draw-lines-between-each-subject-in-case-of-repeated-measurements

tiff("Individual_Pre_learn_Ch2-9.tiff", units="in", width=2.4, height=2.4, res=600)
xyplot(rvalue ~ phase, group = ID, data = FC_long[1:30,], type = "b", xlab="", ylab="", ylim=c(-1.8,1.8))
dev.off()

tiff("Individual_Pre_learn_Ch9-30.tiff", units="in", width=2.4, height=2.4, res=600)
xyplot(rvalue ~ phase, group = ID, data = FC_long[31:60,], type = "b", xlab="", ylab="", ylim=c(-1.8,1.8))
dev.off()

tiff("Individual_Pre_learn_Ch9-34.tiff", units="in", width=2.4, height=2.4, res=600)
xyplot(rvalue ~ phase, group = ID, data = FC_long[61:90,], type = "b", xlab="", ylab="", ylim=c(-1.8,1.8))
dev.off()

tiff("Individual_Pre_learn_Ch9-35.tiff", units="in", width=2.4, height=2.4, res=600)
xyplot(rvalue ~ phase, group = ID, data = FC_long[91:120,], type = "b", xlab="", ylab="", ylim=c(-1.8,1.8))
dev.off()


xyplot(rvalue ~ phase, group = ID, data = FC_long[121:150,], type = "b")
xyplot(rvalue ~ phase, group = ID, data = FC_long[151:180,], type = "b")
xyplot(rvalue ~ phase, group = ID, data = FC_long[181:210,], type = "b")
xyplot(rvalue ~ phase, group = ID, data = FC_long[211:240,], type = "b")
xyplot(rvalue ~ phase, group = ID, data = FC_long[241:270,], type = "b")
xyplot(rvalue ~ phase, group = ID, data = FC_long[271:300,], type = "b")


# Scatter plot to find potential linear or nonlinear trend
ggplot(data = FC, mapping = aes(x = PreRest_Ch2_9, y = Learn_Ch2_9)) +
  geom_point() +
  geom_smooth(aes(x = PreRest_Ch2_9, y = Learn_Ch2_9), method = "lm",
              data = FC)

ggplot(data = FC, mapping = aes(x = PreRest_Ch9_30, y = Learn_Ch9_30)) +
  geom_point() 

ggplot(data = FC, mapping = aes(x = PreRest_Ch9_34, y = Learn_Ch9_34)) +
  geom_point() 

ggplot(data = FC, mapping = aes(x = PreRest_Ch9_35, y = Learn_Ch9_35)) +
  geom_point() 

ggplot(data = FC, mapping = aes(x = PreRest_Ch17_38, y = Learn_Ch17_38)) +
  geom_point() 

ggplot(data = FC, mapping = aes(x = PreRest_Ch39_40, y = Learn_Ch39_40)) +
  geom_point()

ggplot(data = FC, mapping = aes(x = PreRest_Ch2_43, y = Learn_Ch2_43)) +
  geom_point()

ggplot(data = FC, mapping = aes(x = PreRest_Ch17_43, y = Learn_Ch17_43)) +
  geom_point()

ggplot(data = FC, mapping = aes(x = PreRest_Ch19_43, y = Learn_Ch19_43)) +
  geom_point()

ggplot(data = FC, mapping = aes(x = PreRest_Ch22_43, y = Learn_Ch22_43)) +
  geom_point()


# 
FC$LP_Ch2_9 <- FC$Learn_Ch2_9-FC$PreRest_Ch2_9
FC$LP_Ch9_30 <- FC$Learn_Ch9_30-FC$PreRest_Ch9_30
FC$LP_Ch9_34 <- FC$Learn_Ch9_34-FC$PreRest_Ch9_34
FC$LP_Ch9_35 <- FC$Learn_Ch9_35-FC$PreRest_Ch9_35

FC$LP_Ch17_38 <- FC$Learn_Ch17_38-FC$PreRest_Ch17_38
FC$LP_Ch39_40 <- FC$Learn_Ch39_40-FC$PreRest_Ch39_40


FC$LP_Ch2_43 <- FC$Learn_Ch2_43-FC$PreRest_Ch2_43
FC$LP_Ch17_43 <- FC$Learn_Ch17_43-FC$PreRest_Ch17_43
FC$LP_Ch19_43 <- FC$Learn_Ch19_43-FC$PreRest_Ch19_43
FC$LP_Ch22_43 <- FC$Learn_Ch22_43-FC$PreRest_Ch22_43

# Scatter plot to find potential linear or nonlinear trend
ggplot(data = FC, mapping = aes(x = PreRest_Ch2_9, y = LP_Ch2_9)) +
  geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(aes(x = PreRest_Ch2_9, y = LP_Ch2_9), method = "lm",
              data = FC, se=F)

ggplot(data = FC, mapping = aes(x = PreRest_Ch9_30, y = LP_Ch9_30)) +
  geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(aes(x = PreRest_Ch9_30, y = LP_Ch9_30), method = "lm",
              data = FC, se=F)

ggplot(data = FC, mapping = aes(x = PreRest_Ch9_34, y = LP_Ch9_34)) +
  geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(aes(x = PreRest_Ch9_34, y = LP_Ch9_34), method = "lm",
              data = FC, se=F)

ggplot(data = FC, mapping = aes(x = PreRest_Ch9_35, y = LP_Ch9_35)) +
  geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(aes(x = PreRest_Ch9_35, y = LP_Ch9_35), method = "lm",
              data = FC, se=F)

ggplot(data = FC, mapping = aes(x = PreRest_Ch17_38, y = LP_Ch17_38)) +
  geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(aes(x = PreRest_Ch17_38, y = LP_Ch17_38), method = "lm",
              data = FC, se=F)

ggplot(data = FC, mapping = aes(x = PreRest_Ch39_40, y = LP_Ch39_40)) +
  geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(aes(x = PreRest_Ch39_40, y = LP_Ch39_40), method = "lm",
              data = FC, se=F)

ggplot(data = FC, mapping = aes(x = PreRest_Ch2_43, y = LP_Ch2_43)) +
  geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(aes(x = PreRest_Ch2_43, y = LP_Ch2_43), method = "lm",
              data = FC, se=F)

ggplot(data = FC, mapping = aes(x = PreRest_Ch17_43, y = LP_Ch17_43)) +
  geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(aes(x = PreRest_Ch17_43, y = LP_Ch17_43), method = "lm",
              data = FC, se=F)

ggplot(data = FC, mapping = aes(x = PreRest_Ch19_43, y = LP_Ch19_43)) +
  geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(aes(x = PreRest_Ch19_43, y = LP_Ch19_43), method = "lm",
              data = FC, se=F)

ggplot(data = FC, mapping = aes(x = PreRest_Ch22_43, y = LP_Ch22_43)) +
  geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(aes(x = PreRest_Ch22_43, y = LP_Ch22_43), method = "lm",
              data = FC, se=F)













