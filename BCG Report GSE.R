# BGC Project

# Modelling an influenza vaccination intervention using 
# different modelling strategies

# General stoschastic epidemic

source("BGC Report GSE functions.R")

# Paramaters
seed<-1917+0:2
n.all<-c(200,1000,5000)


m<-1 # Initially infected
lambda <- 0.5 # Lambda of Poisson process (no of contacts)
gamma <- 4.8 # Gamma of Exp model of infection time
# Number of vaccinated each time point
vacc.ant <- matrix(c(rep(0,3),rep(round(n.all*0.05),2)),nrow=3) 
vacc.start.T <- c(0,15,25)
vacc.int.T <- rep(1,3)

sim.res<-list();i<-1;j<-1
for (j in 1:3){
  for (i in 1:3){
    print(c(i,j))
    sim.res<-c(sim.res,
               list(
                 GSE.sim(n=n.all[i],
                         m=m,lambda=lambda,gamma=gamma,
                         vacc.ant=vacc.ant[i,j],
                         vacc.start.T=vacc.start.T[j],
                         vacc.int.T=vacc.int.T[j],
                         seed=seed[i])[[2]]))
  }
}

sim.no<-1000
sim.peak.all<-list();i<-1;j.length<-3;all.data<-list()
for (j in 1:j.length){
  sim.peak<-list()
  for (i in 1:sim.no){
    print(c(j,i))
    temp.data<-GSE.sim(n=n.all[2],
                         m=m,lambda=lambda,gamma=gamma,
                         vacc.ant=vacc.ant[2,j],
                         vacc.start.T=vacc.start.T[j],
                         vacc.int.T=vacc.int.T[j],
                         seed=seed[2]+i)[[2]]
    max.no<-which(max(temp.data$I)==temp.data$I)[1]
    temp.data$cum.I<-temp.data$R[dim(temp.data)[1]]
    temp.data$max.T<-temp.data$time[dim(temp.data)[1]]
    
    if(i==1){sim.peak<-temp.data[max.no,]}else{sim.peak<-rbind(sim.peak,temp.data[max.no,])}
    if(i==1){all.data<-list(temp.data)}else{all.data<-c(all.data,list(temp.data))}
  }
  if(i==1){sim.peak.all<-list(sim.peak)}else{sim.peak.all<-c(sim.peak.all,list(sim.peak))}
}

# Plot the results
interv<-c("No vaccine intervention",
          paste("Vaccination started at t =",vacc.start.T[2]),
          paste("Vaccination started at t =",vacc.start.T[3]))

png(file="PeakT.01.png",width=623,height=600)
par(mfrow=c(3,1))
for (k in 1:3){
  hist(sim.peak.all[[k]]$time,breaks=100,
       xlim=c(0,max(sim.peak.all[[1]]$time)),
       main=interv[k],
       ylab="Simulations",
       xlab="Peak of epidemic")
}
dev.off()

png(file="PeakI.01.png",width=623,height=600)
par(mfrow=c(3,1))
for (k in 1:3){
  hist(sim.peak.all[[k]]$I,breaks=100,
       xlim=c(0,max(sim.peak.all[[1]]$I)),
       main=interv[k],
       ylab="Simulations",
       xlab="Infected at peak")
}
dev.off()



png(file="MaxT.01.png",width=623,height=600)
par(mfrow=c(3,1))
for (k in 1:3){
  hist(sim.peak.all[[k]]$max.T,breaks=100,
       xlim=c(0,max(sim.peak.all[[1]]$max.T)),
       main=interv[k],
       ylab="Simulations",
       xlab="Length of epidemic")  
}
dev.off()

png(file="CumI.01.png",width=623,height=600)
par(mfrow=c(3,1))
for (k in 1:3){
  hist(sim.peak.all[[k]]$cum.I,breaks=100,
       xlim=c(0,max(sim.peak.all[[1]]$cum.I)),
       main=interv[k],
       ylab="Simulations",
       xlab="Total number of infected")
  
}
dev.off()

png(file="Scatter.01.png",width=623,height=600)
par(mfrow=c(3,1))
for (k in 1:3){
  with(sim.peak.all[[k]],
       plot(max.T,cum.I,
            ylim=c(0,max(sim.peak.all[[1]]$cum.I)),
            xlim=c(0,max(sim.peak.all[[1]]$max.T)),
            main=interv[k],
            ylab="Infected",
            xlab="Length of epidemic")
  )
}
dev.off()


# Plot the results
png(file="DiffPop01.png",width=623,height=600)
par(mfrow=c(3,1))
for(k in 1:3){
  n<-sim.res[[k]]$S[1]+sim.res[[k]]$I[1]+sim.res[[k]]$R[1]
  plot(sim.res[[k]]$time,sim.res[[k]]$S,type="l",col="Green",ylim=c(0,n),ylab="Population size",xlab="Time")
  lines(x=sim.res[[k]]$time,sim.res[[k]]$I, col="Red")
  lines(x=sim.res[[k]]$time,sim.res[[k]]$R, col="Blue")
#  lines(x=sim.res[[k]]$time,sim.res[[k]]$V, col="Black")
  legend("topright", c("Susceptibles", "Infecteds", "Recovereds"), 
         pch = 1, col = c("Green","Red","Blue"))
}
dev.off()

png(file="DiffPop02.png",width=623,height=600)
par(mfrow=c(3,1))
for(k in c(2,5,8)){
  n<-sim.res[[k]]$S[1]+sim.res[[k]]$I[1]+sim.res[[k]]$R[1]
  plot(sim.res[[k]]$time,sim.res[[k]]$S,type="l",col="Green",ylim=c(0,n),ylab="Population size",xlab="Time")
  lines(x=sim.res[[k]]$time,sim.res[[k]]$I, col="Red")
  lines(x=sim.res[[k]]$time,sim.res[[k]]$R, col="Blue")
  lines(x=sim.res[[k]]$time,sim.res[[k]]$V, col="Black")
  legend("topright", c("Susceptibles", "Infecteds", "Recovereds", "Vaccinated"), 
         pch = 1, col = c("Green","Red","Blue","Black"))
}
dev.off()

png(file="DiffPop12.png",width=623,height=600)
par(mfrow=c(3,3))
for(k in 1:3){
  n<-sim.res[[k]]$S[1]+sim.res[[k]]$I[1]+sim.res[[k]]$R[1]
  plot(sim.res[[k]]$time,sim.res[[k]]$S,type="l",col="Green",ylim=c(0,n),ylab="Population size",xlab="Time")
  lines(x=sim.res[[k]]$time,sim.res[[k]]$I, col="Red")
  lines(x=sim.res[[k]]$time,sim.res[[k]]$R, col="Blue")
  
  if(k==1){title(interv[1])}
  
  l<-k+3
  n<-sim.res[[l]]$S[1]+sim.res[[l]]$I[1]+sim.res[[l]]$R[1]
  plot(sim.res[[l]]$time,sim.res[[l]]$S,type="l",col="Green",ylim=c(0,n),ylab="Population size",xlab="Time")
  lines(x=sim.res[[l]]$time,sim.res[[l]]$I, col="Red")
  lines(x=sim.res[[l]]$time,sim.res[[l]]$R, col="Blue")
  lines(x=sim.res[[l]]$time,sim.res[[l]]$V, col="Black")

  if(k==1){title(interv[2])}
  
  l<-k+6
  n<-sim.res[[l]]$S[1]+sim.res[[l]]$I[1]+sim.res[[l]]$R[1]
  plot(sim.res[[l]]$time,sim.res[[l]]$S,type="l",col="Green",ylim=c(0,n),ylab="Population size",xlab="Time")
  lines(x=sim.res[[l]]$time,sim.res[[l]]$I, col="Red")
  lines(x=sim.res[[l]]$time,sim.res[[l]]$R, col="Blue")
  lines(x=sim.res[[l]]$time,sim.res[[l]]$V, col="Black")
  
  if(k==1){title(interv[3])}
  
  if(k==3){legend("right", c("Susceptibles", "Infecteds", "Recovereds", "Vaccinated"), 
         pch = 1, col = c("Green","Red","Blue","Black"))}
}
dev.off()

png(file="Scatter.02.png",width=623,height=600)
par(mfrow=c(3,1))
for (k in 1:3){
  with(sim.peak.all[[k]],
       plot(time,I,
            ylim=c(0,max(sim.peak.all[[1]]$I)),
            xlim=c(0,max(sim.peak.all[[1]]$time)),
            main=interv[k],
            ylab="Infected at peak",
            xlab="Time of peak")
  )
}
dev.off()

head(sim.peak.all[[1]])
k<-3
def.outb<-20
# Total infected
mean(sim.peak.all[[k]]$cum.I[sim.peak.all[[k]]$cum.I>def.outb]/1000)
sd(sim.peak.all[[k]]$cum.I[sim.peak.all[[k]]$cum.I>def.outb]/1000)
# Peak infected
mean(sim.peak.all[[k]]$I[sim.peak.all[[k]]$cum.I>def.outb]/1000)
sd(sim.peak.all[[k]]$I[sim.peak.all[[k]]$cum.I>def.outb]/1000)
# Peak time
mean(sim.peak.all[[k]]$time[sim.peak.all[[k]]$cum.I>def.outb])
sd(sim.peak.all[[k]]$time[sim.peak.all[[k]]$cum.I>def.outb])
# P(Outbreak)
mean(sim.peak.all[[k]]$cum.I>50)
