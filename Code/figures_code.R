#making dataframes for plotting####

#the 2 data files below are a consolidation of csv files made above



kt_wrasse_theta<-read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/theta_kt_wrasse.csv", header=TRUE)
kt_wrasse_sigma<-read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/rate_kt_wrasse.csv" , header=TRUE)

input_wrasse_theta<-read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/theta_input_wrasse.csv", header=TRUE)
input_wrasee_sigma<- read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/rate_input_wrasse.csv", header=TRUE)

output_wrasse_theta<- read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/theta_output_wrasse.csv", header=TRUE)
output_wrasse_sigma<- read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/rate_output_wrasse.csv", header=TRUE)

kt_stoma_theta<-read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/theta_kt_shrimp.csv", header=TRUE)
kt_stoma_sigma<-read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/rate_kt_shrimp.csv", header=TRUE)

input_stoma_theta<-read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/theta_input_shrimp.csv", header=TRUE)

output_stoma_theta<- read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/theta_output_shrimp.csv", header=TRUE)
output_stoma_sigma<- read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/rate_output_shrimp.csv", header=TRUE)

coupler_stoma_theta<- read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/theta_coupler_shrimp.csv", header=TRUE)
coupler_stoma_sigma<- read.csv("~/Desktop/durophagy/data/ouwie parameters/final parameters/rate_coupler_shrimp.csv", header=TRUE)



#plotting####
#kt rate
kt_wrasse_sigma$regime<-as.factor(kt_wrasse_sigma$regime)
wrasse_kt_rate_plot=ggplot(kt_wrasse_sigma, aes(x=regime, y=sigma.sq)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_fill_manual(values=c('dodgerblue3','seagreen3'))+
  scale_size_manual(values=c(1,1))+
  theme(legend.position="none")+
  ylim(0,.003)
wrasse_kt_rate_plot<- wrasse_kt_rate_plot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                 panel.background = element_blank(), axis.line = element_line(colour = "black")) +stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")

#kt theta
kt_wrasse_theta$regime<-as.factor(kt_wrasse_theta$regime)
wrasse_kt_theta_plot=ggplot(kt_wrasse_theta, aes(x=regime, y=theta)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_fill_manual(values=c('dodgerblue3','seagreen3'))+
  scale_size_manual(values=c(1,1))+
  theme(legend.position="none")+
  ylim(-.4,0)
wrasse_kt_theta_plot<- wrasse_kt_theta_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")

#input rate
input_wrasee_sigma$regime<-as.factor(input_wrasee_sigma$regime)
wrasse_input_rate_plot=ggplot(input_wrasee_sigma, aes(x=regime, y=sigma.sq)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_fill_manual(values=c('dodgerblue3','seagreen3'))+
  scale_size_manual(values=c(1,1))+
  theme(legend.position="none")+
  ylim(0,.0017)
wrasse_input_rate_plot<- wrasse_input_rate_plot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")

#input theta
input_wrasse_theta$regime<-as.factor(input_wrasse_theta$regime)
wrasse_input_theta_plot=ggplot(input_wrasse_theta, aes(x=regime, y=theta)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_fill_manual(values=c('dodgerblue3','seagreen3'))+
  scale_size_manual(values=c(1,1))+
  theme(legend.position="none")+
  ylim(-.55,-.3)
wrasse_input_theta_plot<- wrasse_input_theta_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")


#output rate
output_wrasse_sigma$regime<-as.factor(output_wrasse_sigma$regime)
wrasse_output_rate_plot=ggplot(output_wrasse_sigma, aes(x=regime, y=sigma.sq)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_fill_manual(values=c('dodgerblue3','seagreen3'))+
  scale_size_manual(values=c(1,1))+
  theme(legend.position="none")+
  ylim(0,.0005)
wrasse_output_rate_plot<- wrasse_output_rate_plot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                         panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")


#output theta
output_wrasse_theta$regime<-as.factor(output_wrasse_theta$regime)
wrasse_output_theta_plot=ggplot(output_wrasse_theta, aes(x=regime, y=theta))+ 
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_fill_manual(values=c('dodgerblue3','seagreen3'))+
  scale_size_manual(values=c(1,1))+
  theme(legend.position="none")+
  ylim(-.35,-.1)
wrasse_output_theta_plot<- wrasse_output_theta_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                            panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")


#diff versions of same grid
#grid.arrange(wrasse_kt_rate_plot,wrasse_kt_theta_plot, wrasse_input_rate_plot, wrasse_input_theta_plot, wrasse_output_rate_plot,wrasse_output_theta_plot, ncol=2)

ggarrange(wrasse_kt_theta_plot,wrasse_input_theta_plot,wrasse_output_theta_plot, ncol=3, labels = c("A","B","C"),font.label = list(size = 20))
ggarrange(wrasse_kt_rate_plot, wrasse_input_rate_plot, wrasse_output_rate_plot, ncol=3, labels = c("A","B","C"),font.label = list(size = 20))

#stomatopods

#kt rate
kt_stoma_sigma$regime<-as.factor(kt_stoma_sigma$regime)
stoma_kt_rate_plot=ggplot(kt_stoma_sigma, aes(x=regime, y=sigma.sq)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_fill_manual(values=c('dodgerblue3','seagreen3'))+
  scale_size_manual(values=c(1,1))+
  theme(legend.position="none")+
  ylim(0,.0075)
stoma_kt_rate_plot<- stoma_kt_rate_plot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")


#kt theta
kt_stoma_theta$regime<-as.factor(kt_stoma_theta$regime)
stoma_kt_theta_plot=ggplot(kt_stoma_theta, aes(x=regime, y=theta)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_fill_manual(values=c('dodgerblue3','seagreen3'))+
  scale_size_manual(values=c(1,1))+
  theme(legend.position="none")+
  ylim(.72,.9)
stoma_kt_theta_plot<- stoma_kt_theta_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                  panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")
#input theta
input_stoma_theta$regime<-as.factor(input_stoma_theta$regime)
stoma_input_theta_plot=ggplot(input_stoma_theta, aes(x=regime, y=theta)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_fill_manual(values=c('dodgerblue3','seagreen3'))+
  scale_size_manual(values=c(1,1))+
  theme(legend.position="none")+
  ylim(-.2,-.07)
stoma_input_theta_plot<- stoma_input_theta_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")

#output rate
output_stoma_sigma$regime<-as.factor(output_stoma_sigma$regime)
output_rate_plot=ggplot(output_stoma_sigma, aes(x=regime, y=sigma.sq)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_fill_manual(values=c('dodgerblue3','seagreen3'))+
  scale_size_manual(values=c(1,1))+
  theme(legend.position="none")+
  ylim(0,.0045)
stoma_output_rate_plot<- output_rate_plot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")

#output theta
output_stoma_theta$regime<-as.factor(output_stoma_theta$regime)
stoma_output_theta_plot=ggplot(output_stoma_theta, aes(x=regime, y=theta)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_color_manual(values=c('dodgerblue3','seagreen3'))+
  scale_fill_manual(values=c(1,1))+
  theme(legend.position="non")+
  ylim(-1.,-.9)
stoma_output_theta_plot<- stoma_output_theta_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")


#coupler rate
coupler_stoma_sigma$regime<-as.factor(coupler_stoma_sigma$regime)
stoma_coupler_rate_plot=ggplot(coupler_stoma_sigma, aes(x=regime, y=sigma.sq)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_fill_manual(values=c('dodgerblue3','seagreen3'))+
  scale_size_manual(values=c(1,1))+
  theme(legend.position="none")+
  ylim(.00005,.0002)
stoma_coupler_rate_plot<- stoma_coupler_rate_plot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                         panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")

#coupler theta
coupler_stoma_theta$regime<-as.factor(coupler_stoma_theta$regime)
stoma_coupler_theta_plot=ggplot(coupler_stoma_theta, aes(x=regime, y=theta)) +
  geom_violin(aes(fill=regime)) + 
  geom_smooth(method=lm)+
  scale_shape_manual(values=c(15, 16))+ 
  scale_color_manual(values=c('dodgerblue3','seagreen3'))+
  scale_fill_manual(values=c(1,1))+
  theme(legend.position="non")+
  ylim(0.0325,.039)
stoma_coupler_theta_plot <- stoma_coupler_theta_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=22), axis.text.y = element_text(size=22),
                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_summary(fun=median, geom="point", size=5, color="red")+stat_summary(fun=mean, geom="point", size=5, color="gold")



ggarrange(stoma_kt_theta_plot, stoma_output_theta_plot, stoma_coupler_theta_plot, stoma_input_theta_plot, ncol=4, labels = c("A","B","C","D"),font.label = list(size = 20))
          
ggarrange(stoma_kt_rate_plot, stoma_output_rate_plot, stoma_coupler_rate_plot, ncol=3,labels = c("A","B","C"),font.label = list(size = 20))
