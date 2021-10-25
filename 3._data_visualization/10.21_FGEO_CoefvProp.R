# PLOTS OF COEF V PROP GLM AND GAM

# output study A

## CHOOSE ONE
#GLM coef net A
dat<-read.csv("data/output.study.net.A.CM_9.29.21.csv", sep=",",header=TRUE) 
#GLM coef mod 2
dat<-read.csv("data/output.study.net.A_9.29.21.csv", sep=",",header=TRUE) 

# plot prop change A against coef A: different
ggplot(data = dat, mapping = aes(x=prop.change.A, y = coef_A))+
  geom_point()+
  ylab("coef_A")+
  xlab("prop.change.A")+
  geom_abline(intercept = 0, slope = 1, linetype="dashed")+
  xlim(-5,1)+
  ylim(-5,1)

# plot net prop change A against prop change A: very linear
ggplot(data = dat, mapping = aes(x=net.prop.change.A, y = prop.change.A))+
  geom_point()+
  ylab("prop.change.A")+
  xlab("net.prop.change.A")+
  geom_abline(intercept = 0, slope = 1, linetype="dashed")+
  xlim(-5,1)+
  ylim(-5,1)

# output study CM
ggplot(data = dat, mapping = aes(x=prop.change.CMH, y = coef_CMH))+
  geom_point()+
  ylab("coef_CMH")+
  xlab("prop.change.CMH")+
  geom_abline(intercept = 0, slope = 1, linetype="dashed")+
  xlim(-1,1)+
  ylim(-1,1)

# plot net prop change CM against prop change CM:
ggplot(data = dat, mapping = aes(x=net.prop.change.CMH, y = prop.change.CMH))+
  geom_point()+
  ylab("prop.change.CMH")+
  xlab("net.prop.change.CMH")+
  geom_abline(intercept = 0, slope = 1, linetype="dashed")+
  xlim(-1,1)+
  ylim(-1,1)

