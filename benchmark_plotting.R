library(ggplot2)

data <- read.csv("plate_operations_rollout_benchmark.csv")

opts <- data$operations[data$nsims==100]
seq <- order(order(opts,decreasing=F))
data$opt=rep(opts,each=5)
data$order <- rep(seq,each=5)
data$nsims_fact=as.factor(data$nsims)


ggplot(data=data)+
  geom_point(aes(x=order,y=operations,color=nsims_fact))+
  theme_classic()+
  xlab("Run")+
  ylab("Number of Required Pipette Operations")

ggplot(data=data)+
  geom_point(aes(x=nsims,y=time))+
  geom_smooth(aes(x=nsims,y=time),method=lm)

model <- lm(time~nsims,data=data)

summary(model)
