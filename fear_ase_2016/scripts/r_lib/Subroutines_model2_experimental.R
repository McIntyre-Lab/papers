#Library of functions
#this runs on Lauren's data and calculates intervals
#It loops though the catagorizations
#Model
#This model is a modification of the previous model
#because the structure of the data changed.
#theta G "average" rate in the sample. 
#theta is centered in 1-p to take into accoount the bias
#in the estimation of theta
#See specification of the model below

#x_i RNA hits on A
#l_i RNA no hits on A, hits on G
#y_i DNA hits on A
#k_i DNA no hits on A, hits on G
#i is biorep

#model for p 
#y_i|gamma_i~poisson(gamma_i)
#gamma_i|k_i,p~gamma(k_i,p/(1-p))
#k_i~Poisson(delta), j=1,...J_i
#delta~gamma(a_delta,b_delta)
#p~beta(v,v)

#Model for theta
#Once a sample p_1,p_2,... if p is given
#x_i|eta_i~poisson(eta_i)
#eta_i~gamma(l_i,theta/(1+theta))
#l_i~Poisson(lamda)
#lambda~gamma(a_lambda,b_lambda)
#theta_i~beta((1-p_i)t,p_i t)

#That is the estimation of p is not afected by the RNA data
#The only things that change in the Gibbs sample
#(compared to model2)
#is that when updating p all the terms 
#depending on theta are ignored
#is that it when updating p


#The main function in this library is:
#gibbs_sampler
###########
r_eta_i<-function(n_i,th){
	return(  rgamma(1,n_i,rate=(1-th)^(-1))  )
		}
r_eta<-function(nx,theta){
	dimen<-length(nx)
	output<-rep(0,dimen)
	for(i in 1:dimen)
	output[i]<-r_eta_i(as.numeric(nx[i]),theta)
	return(output)
	}
#r_eta(nx=n_i,th=sum(x)/sum(nx))
#########################################################################
r_lambda<-function(l_d,Ix,a_lambda,b_lambda){
	return(rgamma(1,l_d+a_lambda,rate=b_lambda+Ix ))
	}
#r_lambda(l_d=lv_d,Ix,a_lambda,b_lambda)
#########################################################################
r_theta_MH<-function(th,p,l_d,eta_d,tt,sigma_theta=0.005){
	u	<- th/(1-th)
 	#We have two different metropolis-hasting schemes depending on the value of t
	if(tt<0){  #for relatively small t
	u_candid<-rgamma(1,l_d+(1-p)*tt,eta_d)
	alpha<-((1+u)/(1+u_candid))^tt
			}else{		#for large t
	
	u_candid<- rnorm(1,mean=u,sd=sigma_theta)  #u_cand is a sample of f(u)propto 1/(1+u)^t
	#print(paste("u candidate=", u_candid))
	alpha<-ifelse(0<u_candid,
	exp(tt*log((1+u)/(1+u_candid))+dgamma(u_candid,shape=l_d+(1-p)*tt,rate = eta_d, log = TRUE)-dgamma(u,shape=l_d+(1-p)*tt,rate = eta_d, log = TRUE)),
	0)
    }
    #print(alpha)
	accept<-ifelse(alpha>runif(1),1,0)
	output<-ifelse(accept == 1,u_candid/(1+u_candid),th)
	return(list(theta=output,accept=accept ))
	}

#r_theta_MH(th=theta,p=.5,l_d=lv_d,eta_d=sum(eta_i),t=10,sigma_theta=0.005)
###########################################
log_complete_cond_p<-function(p,theta,kd,gammad,tt,v){
return((kd+v-1)*log(p/(1-p))+2*(v-1)*log(1-p)-gammad*p/(1-p))
}
#log_complete_dond_p(p,theta,kd,sum(gamma_i),t,v)


r_p_MH<-function(p,theta,kd,gamma_i,tt,v,sigma_MH=.001){
	#if(is.na(p)){print("p NA")}
	p_cand<-rnorm(1,p,sd=sigma_MH)
    #if(is.na(p_cand)){print("p_cand_na")}
    if(p_cand>0 && p_cand<1){
    sum_gamma<-	sum(gamma_i)
    log_alpha<-log_complete_cond_p(p_cand,theta,kd,sum_gamma,tt,v)-log_complete_cond_p(p,theta,kd,sum_gamma,tt,v)
    accept<-ifelse(log_alpha>log(runif(1)),1,0)
	output<-ifelse(accept == 1,p_cand,p)
}else{output<-p;accept<-0}
	return(list(p=output,accept=accept ))
	}

#r_p_MH(p=0.5,theta,kd=kv_d,gamma_i,t,v,sigma_MH=.05)    
    

##n steps of Gibbs sampler
#theta=theta_i
gibbs_n_steps<-function(
n_steps,p,theta,
eta_i,lambda,gamma_i,delta,
n_i,lv_d,
m_i,kv_d,kd,
a_lambda,b_lambda,a_delta,b_delta,tt,v,sigma_MH,sigma_theta){
accept_counter_p 	<-	0
accept_counter_theta<-	0
Ix<-length(n_i)
Iy<-length(m_i)
	#print(paste("p in gibbs_nsteps",p))
	#I generate an observation of p from the second model 
	for(ii in 1:n_steps){
	gamma_i	<- r_eta(nx=m_i,    rep(p,length(y))    )
	#delta	<- r_lambda(kv_d,Iy,a_delta,b_delta)
	tem	<-  r_p_MH(p,theta,kd=kv_d,gamma_i,tt,v,sigma_MH=sigma_MH)
	p	<- tem$p
	accept_counter_p<-accept_counter_p+tem$accept
	#print(paste("p step",ii))
	}
	
	#With this p as hyperparameter in the first
	#model a value of theta is generated 
	for(ii in 1:n_steps){
	eta_i	<-r_eta(nx=n_i,th=theta)
	#lambda	<-r_lambda(l_d=lv_d,Ix,a_lambda,b_lambda)
	tem<-r_theta_MH(th=theta,p,l_d=lv_d,eta_d=sum(eta_i),tt=tt,sigma_theta=sigma_theta)
	theta <-tem$theta
	accept_counter_theta<- accept_counter_theta+tem$accept
    #print(paste("theta step",ii))
    #print(tem$accept)
    #print(theta)
	}
	
return(list(p=p,theta=theta,eta=eta_i, lambda=lambda,gamma=gamma_i,delta=delta,accepted_p=accept_counter_p,accepted_theta=accept_counter_theta))
	}
	
	
#gibbs_n_steps(n_steps=100,p,theta=theta_i,eta_ij,lambda_i,gamma_ij,delta_i,n_ij,nx,lv_id,Jx,m_ij,ny,kv_id,Jy,kdd,a_lambda,b_lambda,a_delta,b_delta,t,v)#############################################

#################################################################################
#################################################################################
#################################################################################
#x: vector of RNA hits on A
#nx:tecrep index. eg, if nx[1]=1,nx[1]=1, then x[1] and x[2] are techrep counts of the same biological sample. 
#n_ij: Total number of RNA counts
#Analogously y,ny,m_ij are DNA counts  
gibbs_sampler<-function(
##data 
x,n_i,y,m_i,
##Hyperparameters
tt=		1000,   #prior variance of theta_i is p(1-p)/(t+1)
v=		100,    #prior variance of p is 1/(4(v+1))
sigma_MH=.04,  	#Standard deviacion in the 
				#proposed distribution in the
                #Metropolis Hasting algorithm when simulating p
sigma_theta=0.04,#Standard deviations of the sigma_thetas in the MH. I uses these only when t is large
#Parameers of the gamma distribution of lambda and lambda delta
a_lambda=1000,
b_lambda=1,
a_delta=1000,
b_delta=1,
##Gibbs Sampler parameters
burn_in=1000,storing_every=10,m_simulations=100,figures=1){

x<-as.numeric(x)
n_i<-as.numeric(n_i)
Ix<-length(x)
xv_d<-	sum(x)
lv_i<-n_i-x
lv_d<-sum(lv_i)
nv_d<-sum(n_i)

y<-as.numeric(y)
m_i<-as.numeric(m_i)

Iy<-	length(y)
yv_d<-	sum(y)
m_d		<-sum(m_i)
kv_i	<-m_i-y
kv_d	<-sum(kv_i)

#Row estimation of p
#1-ydd/mdd

##Setting initial values
eta_i	<-x

theta_star<-ifelse(xv_d>0,1-xv_d/nv_d,1-1/(nv_d+1))#Row 
theta_star<-ifelse(theta_star == 1,0.9,theta_star)
theta_star<-ifelse(theta_star == 0,0.1,theta_star)
#estimations of theta
lambda  <-	mean(as.numeric(lv_i))
pstar	<-	ifelse(sum(y)>0,1-sum(y)/sum(m_i),1-1/(sum(m_i)+1))
pstar	<-	ifelse(pstar== 1,0.9,pstar)
pstar	<-	ifelse(pstar== 0,0.1,pstar)
p		<-	pstar	
gamma_i	<-	as.numeric(y)   
delta	<-	mean(as.numeric(kv_i))
##Starting with the function

temp		<-list(p=p,theta=theta_star,eta=eta_i,lambda=lambda,gamma=gamma_i,delta=delta,accepted_p=0,accepted_theta=0)
	
	accepted_MH_p<-	0
	accepted_MH_theta<-	0
	
	theta_vector<- rep(0,nrow=m_simulations)
	eta_matrix	<- matrix(0,nrow=m_simulations,ncol=length(eta_i))
	lambda_vector<- rep(0,nrow=m_simulations)

	gamma_matrix	<- matrix(0,nrow=m_simulations,ncol=length(gamma_i))
	delta_vector<- rep(0,nrow=m_simulations)
	
	p_vector	<- rep(0,m_simulations)
  	
#burn in
accepted_MH_p<-0
min_accepted_MH_theta<-0
max_accepted_MH_theta<-1
counter<-1
#sigma_MH=.04
#sigma_theta=rep(0.04,10)

while((accepted_MH_p<0.25 | accepted_MH_p>0.75|  min_accepted_MH_theta<0.25  | max_accepted_MH_theta>0.75) & counter<100){

#burn_in=100;set.seed(0)
temp<-gibbs_n_steps(n_steps=burn_in,p=temp$p,theta=temp$theta,eta_i=temp$eta,lambda=temp$lambda,gamma_i=temp$gamma,delta=temp$delta,n_i,lv_d,m_i,kv_d,kd,a_lambda,b_lambda,a_delta,b_delta,tt,v,sigma_MH,sigma_theta)
#############################################
accepted_MH_p		<-temp$accepted_p/burn_in
accepted_MH_theta	<-temp$accepted_theta/burn_in
if(accepted_MH_p<0.25){sigma_MH<-max(sigma_MH-0.01,0.001)}
if(accepted_MH_p>0.75){sigma_MH<-sigma_MH+0.01}
min_accepted_MH_theta<-min(accepted_MH_theta)
max_accepted_MH_theta<-max(accepted_MH_theta)
if(min_accepted_MH_theta<0.25){sigma_theta<-max(.001,sigma_theta-0.01)}
if(max_accepted_MH_theta>0.75){sigma_theta<-sigma_theta+0.01}
counter<-counter+1
print("Adjusting acceptance rates in burn in")
print("Proportion of Accepted in burn in")
print(paste("p:",accepted_MH_p))
print("theta")
print(accepted_MH_theta)
print(paste("sigma_MH_p",sigma_MH))
print(paste("sigma_MH_theta",sigma_theta[1]))
}
accepted_MH_p		<-temp$accepted_p
accepted_MH_theta	<-temp$accepted_theta
#print(paste("accepted_MH_p",accepted_MH_p))
#print("end of burn in")
######################

for(m in 1:m_simulations){
temp<-gibbs_n_steps(n_steps=storing_every,p=temp$p,theta=temp$theta,eta_i=temp$eta,lambda=temp$lambda,gamma_i=temp$gamma,delta=temp$delta,n_i,lv_d,m_i,kv_d,kd,a_lambda,b_lambda,a_delta,b_delta,tt,v,sigma_MH,sigma_theta)

accepted_MH_p		<-accepted_MH_p+temp$accepted_p
accepted_MH_theta	<-accepted_MH_theta+temp$accepted_theta
	
	theta_vector[m]		<- temp$theta
	eta_matrix[m,]		<- temp$eta
	lambda_vector[m]	<- temp$lambda
	gamma_matrix[m,]	<- temp$gamma
	delta_vector[m]		<- temp$delta
	p_vector[m]			<- temp$p
	}
print("end of simulation")
if(figures == 1){
#par(mfrow=c(3,2))

lims<-pstar+c(-1,1)*1.965*sqrt(pstar*(1-pstar)/m_d)
plot(c(1:m_simulations,1,1),c(p_vector,lims),type="n")
lines(1:m_simulations,p_vector,type="l")
abline(h=quantile(p_vector,c(.025,.975)),col="red",lwd=3)
abline(h=mean(p_vector),col="red",lwd=3)
abline(h=pstar,col="blue",lwd=3,lty=2)
abline(h=lims,col="blue",lwd=3,lty=2)

lims<-theta_star+c(-1,1)*1.965*sqrt(theta_star*(1-theta_star)/nv_d)
plot(c(1:m_simulations,1,1),c(theta_vector,lims),type="n")
lines(1:m_simulations,theta_vector,type="l")
abline(h=quantile(theta_vector,c(.025,.975)),col="red",lwd=3)
abline(h=mean(theta_vector),col="red",lwd=3)
abline(h=theta_star,col="blue",lwd=3,lty=2)
abline(h=lims,col="blue",lwd=3,lty=2)

plot(1:m_simulations,lambda_vector,type="l")
plot(1:m_simulations,delta_vector,type="l")

print("Total Proportion of Accepted (burn in+ skipped and sample chain)")
print(paste("p:",accepted_MH_p/(burn_in+m_simulations*storing_every) ))
print("theta")
print(accepted_MH_theta/(burn_in+m_simulations*storing_every))
print("bias(p)=mean(p_i) -1/2")
print(mean(p_vector)-.5)
print("sample mean - estimated theta")
print(dif_theta<-theta_star-mean(theta_vector))
print("% of RNA   theta est    difference")
print(cbind(theta_star,mean(theta_vector),dif_theta))
print("% of counts , p_model    %-pmode")
print(paste(pstar,mean(p_vector),pstar-mean(p_vector)))

plot(c(0,0),c(1,1),type="n",xlim=c(0,1.5),ylim=c(0,1))
legend("bottom",legend=c("Standard 95% Conf I","95% Cred Int"),col=c("blue","red"),lty=c(3,2),lwd=3)
text(c(0),c(.8),paste("v=",v,", t=",tt	, ", sigma_MH=",sigma_MH),pos=4)
text(c(0),c(.6),paste("a_lambda=",a_lambda,", b_lambda=",b_lambda),pos=4)
text(c(0),c(.4),paste("a_delta=",a_delta,", b_delta=",b_delta),pos=4)
}
return(list(p=p_vector,theta=theta_vector))
}

