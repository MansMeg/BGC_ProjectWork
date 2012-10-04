# BGC Project

# Modelling an influenza vaccination intervention using 
# different modelling strategies

# General stoschastic epidemic
# Paramaters
n<-5000 # Population
m<-5 # Initially infected
lambda <- 1.5 # Lambda of Poisson process (no of contacts)
gamma <- 2 # Gamma of Exp model of infection time
vacc.ant <- 1000 # Number of vaccinated each time point 
vacc.start.T <- 5
vacc.int.T <- .1

# Create data
sim.data<-data.frame(inf=rep(0,n),
                     rec=rep(0,n),
                     inf.t=rep(NA,n),
                     inf.length=rep(NA,n),
                     contacts=rep(NA,n),
                     rec.t=rep(NA,n))

# Adding index case(s)
index.no<-sample(x=1:n,size=m)
sim.data[index.no,]<-matrix(data=rep(c(1,0,0,NA,NA,NA),m),nrow=m,byrow=T)

# Counting variables
t.int<-vacc.start.T
vacc.nr<-1

# Simulation
while(sum(sim.data$inf>sim.data$rec)>0){
  # Identifies the next infectious person
  inf.no<-which(sim.data$inf.t==min(sim.data$inf.t[sim.data$inf==1 & sim.data$rec==0],na.rm=T))
  if(length(inf.no)==0){break}
    
  # Vaccination
  if(min(sim.data[inf.no,3])>t.int && vacc.ant>0){
                       to.vacc<-((vacc.nr-1)*vacc.ant+1):(vacc.nr*vacc.ant)
  	
  	sim.data[to.vacc,2]<-1

 	 	# Infected after vaccination
  		after.vac<-sim.data[to.vacc,1]==1 & sim.data[to.vacc,3]>t.int
  		sim.data[to.vacc[after.vac],1]<-0
  		sim.data[to.vacc[after.vac],3]<-NA
 
  	print(paste("Vaccinated",vacc.ant,"persons at t =",t.int,sep=" "))
  	t.int<-t.int+vacc.int.T
  	vacc.nr<-vacc.nr+1
 	  	
  }

  # Identifies the next infectious person
  inf.no<-which(sim.data$inf.t==min(sim.data$inf.t[sim.data$inf==1 & sim.data$rec==0],na.rm=T))
  if(length(inf.no)==0){break}
  
  # Simulates infection time and no of contacts
  sim.data[inf.no,4]<-rexp(length(inf.no),rate=1/gamma)
  sim.data[inf.no,5]<-rpois(length(inf.no),lambda*sim.data[inf.no,4])
  
  # Categorizing the infectious persons as "recovered"
  sim.data[inf.no,2]<-1
  
  # Randomizing new infectious persons
  for(j in 1:length(inf.no)){
    if(sim.data[inf.no[j],5]>0){
      new.inf<-ceiling(n*runif(sim.data[inf.no[j],5]))
      new.inf<-new.inf[sim.data[new.inf,2]==0] # Only selecting "new" infections
      
      # Simulating infection times for newly infected
      if(length(new.inf)>0){
        sim.data[new.inf,1]<-1
        
        inf.t<-runif(length(new.inf),
                   min=sim.data[inf.no[j],3],
                   max=sim.data[inf.no[j],3]+sim.data[inf.no[j],4])
        sim.data[new.inf,3]<-pmin(sim.data[new.inf,3],inf.t,na.rm=T)
      }
    }  
  }
}

# Calculates the recovery time
inf<-!is.na(sim.data$inf.t)
#sim.data$rec.t[inf]<-ifelse(sim.data$inf.t[inf]+sim.data$inf.length[inf]>sim.data$rec.t[inf],sim.data$rec.t[inf],sim.data$inf.t[inf]+sim.data$inf.length[inf])
sim.data$rec.t<-sim.data$inf.t+sim.data$inf.length


#sim.data$inf.t[sim.data$inf.length==0]<-NA
#sim.data$inf.length[sim.data$inf.length==0]<-NA

# Summarizing the data based on time at rec and inf
inf.df<-as.data.frame(aggregate(sim.data$inf.t,by=list(sim.data$inf.t),FUN=length))
names(inf.df)<-c("time","new.inf")

rec.df<-as.data.frame(aggregate(sim.data$rec.t,by=list(sim.data$rec.t),FUN=length))
names(rec.df)<-c("time","new.rec")

sim.data.aggr<-merge(x=inf.df,y=rec.df,by="time",all=T)
sim.data.aggr[is.na(sim.data.aggr)]<-0

# Calculating S(t), I(t) and R(t)
sim.data.aggr$S<-n-cumsum(sim.data.aggr$new.inf)
sim.data.aggr$R<-cumsum(sim.data.aggr$new.rec)
sim.data.aggr$I<-cumsum(sim.data.aggr$new.inf)-sim.data.aggr$R
sim.data.aggr$V<-0
for(i in 1:(vacc.nr-1)){
	sim.data.aggr$V	[sim.data.aggr$time>vacc.start.T+(i-1)*vacc.int.T]<-i*vacc.ant}

# Plot the results
plot(sim.data.aggr$time,sim.data.aggr$S,type="l",col="Blue",ylim=c(0,n),ylab="Population size",xlab="Time")
lines(x=sim.data.aggr$time,sim.data.aggr$I,col="Red")
lines(x=sim.data.aggr$time,sim.data.aggr$R, col="Green")
lines(x=sim.data.aggr$time,sim.data.aggr$V, col="Yellow")


