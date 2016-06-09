
#This is the Negative binomial version (and then extension) of
#the Poisson Gamma model in
#PGmodelextended_withboth_library4 and described in 
#section 3 of model_extension_withboth5.PDF
#Here we consider that we have six bioreps

#The model has the following structure 
#x, y and z are negative binomial with means
#      tester(x)			line(y)								both(z)
#Line 1	beta_1 q_1			alpha_y	beta_1  q_2					(alpha+1)beta_1 q_3
#Line 2	gamma 	beta_2 q_1 	alpha_y	beta_2  q_2	 gamma delta	(alpha+1)beta_2 q_3
#and overdispersion parameter phi
#To induce the negative binomial distribution 
#we introduce a set latent rvs nu
#See NBmodel_Withboth_testing.R where data according to the model above is generated
#and the functions in this library are tested

library("MCMCpack")
ralpha=function(alpha,y..,z..,beta1,beta2,gamma_,delta,nu,q1,q2,a_alpha=1/2,b_alpha=1/2,sd_MH=0.05){
	out=list(alpha=alpha,accept=0)
	candidate=rnorm(1,alpha,sd=sd_MH)
	if(candidate>0){
		logaccept=	
		z..*log((candidate+1)/(alpha+1))+
		y..*log(candidate/alpha)+
		(alpha-candidate)*(sum(beta1*nu$y1)*q1[2]+gamma_*delta*sum(beta2*nu$y2)*q2[2]+(sum(beta1*nu$z1)*q1[3]+sum(beta2*nu$z2)*q2[3]))+#likelihood
		#-(a_alpha+1)*log(candidate/alpha)-b_alpha*(1/candidate-1/alpha)#inverse gamma prior for alpha
		(alpha-candidate)*b_alpha+(a_alpha-1)*log(candidate/alpha)  #gamma prior alpha
			
		#print(paste("logaccept",logaccept))
	if(log(runif(1))<logaccept) out=list(alpha=candidate,accept=1);
	}
	return(out)
	}





rphi_MH=function(phi,nu,I1,I2,a_phi,b_phi,sd_MH=0.05){
	out=list(phi=phi,accept=0)
	candidate=rnorm(1,phi,sd=sd_MH)
	if(candidate>0){
		
		logaccept=	
		(1/candidate-1/phi)*sum(log(unlist(nu)))-
		(1/candidate-1/phi)*sum(   (unlist(nu)))+
		3*(I1+I2)*(	(1/candidate)*log(1/candidate)-(1/phi)*log(1/phi)-
					lgamma(1/candidate)+lgamma(1/phi))+
		#(a_phi-1)*log(candidate/phi)-b_phi*(candidate-phi)	#gamma prior for phi
		-(a_phi+1)*log(candidate/phi)-b_phi*(1/candidate-1/phi)#Inverse gamma prior
	if(log(runif(1))<logaccept) out=list(phi=candidate,accept=1);
	}
	return(out)
	}

#############################################################
gibbs_NB_niters_withboth=function(
n_iter=10,x1,y1,z1,x2,y2,z2,
all,q1,q2=q1,a,b,I1,I2,alpha_sd_MH=0.5,phi_sd_MH=0.5){
acceptaceratealpha=acceptaceratephi=0;
y..=sum(c(y1,y2));z..=sum(c(z1,z2));

for(ii in 1:n_iter){	
(all$gamma<-rgamma(1,a$gamma+sum(x2)+sum(y2),
rate=with(all,b$gamma+
				q2[1]*		sum(beta2*nu$x2)+
				q2[2]*alpha*sum(beta2*nu$y2)*delta)))

###updating betas and nus
for(i in 1:I1){
(all$beta1[i]<-rgamma(1,a$beta+x1[i]+y1[i]+z1[i],
rate=with(all,b$beta+q1[1]*nu$x1[i] +q1[2]*nu$y1[i]*alpha +q1[3]*nu$z1[i]*(alpha+1))))

all$nu$x1[i]=rgamma(1,1/all$phi+x1[i],with(all,1/phi+q1[1]*			beta1[i]))
all$nu$y1[i]=rgamma(1,1/all$phi+y1[i],with(all,1/phi+q1[2]*alpha*	beta1[i]))
all$nu$z1[i]=rgamma(1,1/all$phi+z1[i],with(all,1/phi+q1[3]*(alpha+1)*beta1[i]))
}

for(i in 1:I2){
(all$beta2[i]<-rgamma(1,a$beta+x2[i]+y2[i]+z2[i],
rate=with(all,b$beta+q2[1]*nu$x2[i]*gamma    +q2[2]*nu$y2[i]*alpha*gamma*delta+q2[3]*nu$z2[i]*(alpha+1) )))

(all$nu$x2[i]<-rgamma(1,1/all$phi+x2[i],with(all,1/phi+q2[1]*			beta2[i]*gamma)))
(all$nu$y2[i]<-rgamma(1,1/all$phi+y2[i],with(all,1/phi+q2[2]*alpha*		beta2[i]*gamma*delta)))
(all$nu$z2[i]<-rgamma(1,1/all$phi+z2[i],with(all,1/phi+q2[3]*(alpha+1)*	beta2[i])))
}

all$delta=rgamma(1,a$delta+sum(y2),rate=with(all,b$delta+q2[2]*alpha*gamma*sum(beta2*nu$y2)))

(tem=with(all,ralpha(alpha,y..,z..,beta1,beta2,gamma,delta,nu,q1,q2,a$alpha,b$alpha,sd_MH=alpha_sd_MH)))
all$alpha=tem$alpha
acceptaceratealpha=acceptaceratealpha+tem$accept

tem=rphi_MH(all$phi,all$nu,I1,I2,a$phi,b$phi,sd_MH=phi_sd_MH)
all$phi=tem$phi
acceptaceratephi=acceptaceratephi+tem$accept
}
return(list(all=all,acceptaceratephi=acceptaceratephi/n_iter,acceptaceratealpha=acceptaceratealpha/n_iter));
}









##############################################
#
#xs=c(1,4,0,2,0,3)
#ys=c(0,0,0,3,5,2)
#zs=c(0,12,34,21,32,12)

## fusion_id        qsim_both        qsim_line      qsim_tester       M_num_reps 
#      "r101"      "F10001_SI
#xs=c( 2,  	7, 10, 		15,  	6, 53)
#ys=c(5, 	12, 3, 		16, 	36, 55)
#zs=c(118, 	159,205, 	278, 	508, 640)
#q_=c(0.08333333, 0.08333333, 0.83333333)
#############################################
gibbs_NegBin_withboth=function(nsim=100,nburnin=100,nlag=10,
xs,		#Vector of reads assigned to "tester" allele
ys,		#Vector of reads assigned to "line" allele
zs,		#Vector of reads assigned to both alleles
is=rep(1:3,2),		#index for biorep
js=c(rep(1,3),rep(2,3)),	#index for Line1=1, vigin=Line2
a_alpha=0.01,
b_alpha=0.01,
a_beta=0.01,
b_beta=0.01, 
a_gamma=0.01,
b_gamma=0.01, 
a_delta=0.01,			#The interaction term
b_delta=0.01,			
a_phi=1,
b_phi=100,				#phi is apriori small	
alpha_sd_MH=0.1,
phi_sd_MH=0.1,			#Initial SD for the normal proposal to sample from phi via MH
q1=c(0.2,0.2,0.6),		#Bias correction hyperparameter in Line 1
q2=q1,					#Bias correction hyperparameter in Line 2
plots=FALSE					#If true plots are depicted
){
I1=length(unique(is[which(js==1)]))
I2=length(unique(is[which(js==2)]))
a=	list(alpha=a_alpha,beta=a_beta,gamma=a_gamma,delta=a_delta,phi=a_phi)
b=	list(alpha=b_alpha,beta=b_beta,gamma=b_gamma,delta=b_delta,phi=b_phi)

#Initiating parameters		
x1=xs[which(js==1)]; if(sum(x1)==0)	x1[sample(1:I1,1)]=1;
y1=ys[which(js==1)]; if(sum(y1)==0)	y1[sample(1:I1,1)]=1;
z1=zs[which(js==1)]; if(sum(z1)==0)	z1[sample(1:I1,1)]=1;
x2=xs[which(js==2)]; if(sum(x2)==0)	x2[sample(1:I2,1)]=1;
y2=ys[which(js==2)]; if(sum(y2)==0)	y2[sample(1:I2,1)]=1;
z2=zs[which(js==2)]; if(sum(z2)==0)	z2[sample(1:I2,1)]=1;
y..=sum(ys)
z..=sum(zs)
#Checking that I have not a vector of all zeros (this messes with the convergence, if so I add a 1 at random)
if(sum(x1)==0)	x1[sample(1:I1)]=1


beta1=pmax(apply(rbind(x1,y1,z1)/q1,MARGIN=2,FUN=mean),0.5)
beta2=pmax(apply(rbind(x2,y2,z2)/q2,MARGIN=2,FUN=mean),0.5)
#beta2/beta1

gamma_=	mean(((x2+y2)/beta2)/(pmax(x1+y1,1)/beta1))
		#((x2.+y2.)/z2.)/((x1.+y1.)/z1.)
delta=1
qys=c(q1[2],q2[2])
qxs=c(q1[1],q2[1])
alpha=mean((ys/qys[js])/(pmax(xs/qxs[js],0.5)))

#(2/3)*((y../q_[2]+(z1.+z2.*gamma_)/q_[3])/(x../q_[1])-1/2)
#mean(c((y../q_[2])/(x../q_[1]),2*((z../q_[3])/(x../q_[1]))-1))
#max(0.1,mean(c(mean(ys/(mu*q_[2])),2*mean(zs)/(mu*q_[3])-1)));

beta1=pmax(apply(rbind(x1,y1/alpha,z1)/q1,MARGIN=2,FUN=mean),0.5)
beta2=pmax(apply(rbind(x2/(gamma_),y2/(alpha*gamma_*delta),z2)/q2,MARGIN=2,FUN=mean),0.5)

nu=list(x1=rep(1,I1),y1=rep(1,I1),z1=rep(1,I1),
		x2=rep(1,I2),y2=rep(1,I2),z2=rep(1,I2))


		#mean((x2./I2)/max(.5,(x1./I1)),(y2./I2)/max(.1,(y1./I1)),(z1./I1)/
		#max(.1,(z2./I2))        )
		#mean(x2./(q_[1]*I2),y2./(alpha*q_[2]*I2),z1./(((alpha+1)/2)*q_[3]*I1))/mu

#print(paste("alphax=",alphax," alphay=",alphay," gamma=",gamma_))
#print("beta1");print(beta1);
#print("beta2");print(beta2);

#print("Expected counts assuming initial values")
#print(
#rbind(	c(q_[1]*sum(beta1)*alphax,		 q_[2]*sum(beta1)*		alphay,			q_[3]*sum(beta1)),
#		c(q_[1]*sum(beta2)*gamma_*alphax,q_[2]*sum(beta2)*gamma_*alphay*delta,	q_[3]*sum(beta2)))
#)

#Printing counts
totals=rbind(c(sum(x1),sum(y1)),c(sum(x2),sum(y2)))
rownames(totals)=c("row1","row2")
colnames(totals)=c("x","y")
#totals/q_[c(1,2)]
phi=b_phi/(a_phi+1)#Mode of an inverse gamma



all=list(beta1=beta1,beta2=beta2,alpha=alpha,gamma=gamma_,delta=delta,nu=nu,phi=phi)




#print("Sample Proportion")
sample.props<-addmargins(prop.table(totals))
totalswboth=cbind(totals,c(sum(z1),sum(z2)))
#colnames(totalswboth)=c("x","y","z")
#print("Total counts")
#print(totalswboth)
#print("Total counts divided by q")
#print(cbind(totalswboth[,1]/q_[1],totalswboth[,2]/q_[2],totalswboth[,3]/q_[3]))


####Adjusting the variance of the normal proposal for  the Metropolis Hasting to sample both alpha and phi
#algorith to sample alpha in the Gibbs sampler
acceptaceratealpha=	acceptaceratephi=	1;
maxiter=100
counter=0
while((	acceptaceratealpha<0.25 | acceptaceratealpha>0.75 |
		acceptaceratephi<0.25 	| acceptaceratephi	>0.75) & counter<maxiter){

tem=gibbs_NB_niters_withboth(n_iter=100,x1,y1,z1,x2,y2,z2,all,q1,q2,a,b,I1,I2,alpha_sd_MH=alpha_sd_MH,phi_sd_MH=phi_sd_MH)
all=tem$all
acceptaceratephi=		tem$acceptaceratephi
acceptaceratealpha=		tem$acceptaceratealpha

if(acceptaceratephi<0.25){phi_sd_MH=phi_sd_MH/2}
if(acceptaceratephi>0.75){phi_sd_MH=2*phi_sd_MH}

if(acceptaceratealpha<0.25){alpha_sd_MH=alpha_sd_MH/2}
if(acceptaceratealpha>0.75){alpha_sd_MH=2*alpha_sd_MH}
counter=counter+1

print(paste("Iteration for MH alpha and phi", counter))
print(paste("acceptace rate alpha", acceptaceratealpha))
print(paste("next SD for MH for alpha ", alpha_sd_MH))
print(paste("acceptace rate phi", acceptaceratephi))
print(paste("next SD for MH for phi ", phi_sd_MH))
}



tem=gibbs_NB_niters_withboth(n_iter=nburnin,x1,y1,z1,x2,y2,z2,all,q1,q2,a,b,I1,I2,alpha_sd_MH=alpha_sd_MH,phi_sd_MH=phi_sd_MH)
all=tem$all

tem$acceptaceratealpha
tem$acceptaceratephi



phi_chain=nux11_chain=nuy11_chain=delta_chain=beta11_chain=beta21_chain=alpha_chain=gamma_chain=rep(NA,nsim)
theta_chain=		matrix(NA,nrow=nsim,ncol=4);colnames(theta_chain)=c("11","12","21","22");


for(ii in 1:nsim){
	all=gibbs_NB_niters_withboth(n_iter=nlag,x1,y1,z1,x2,y2,z2,all,q1,q2,a,b,I1,I2,phi_sd_MH=0.5)$all
	
	delta_chain[ii]=	all$delta	
	beta11_chain[ii]=	all$beta1[1]
	beta21_chain[ii]=	all$beta2[1]
	alpha_chain[ii]=	all$alpha
	gamma_chain[ii]=	all$gamma
	nux11_chain[ii]=	all$nu$x1[1]
	nuy11_chain[ii]=	all$nu$y1[1]
	phi_chain[ii]=		all$phi
	total=	(1+all$gamma)+all$alpha*(1+all$gamma*all$delta)
	theta_chain[ii,"11"]=1	 		         			/total;
	theta_chain[ii,"12"]=all$alpha         	 			/total;
	theta_chain[ii,"21"]=    		all$gamma			/total;
	theta_chain[ii,"22"]=all$alpha*	all$gamma*all$delta	/total;
	}
	
#print(paste(" gamma=",mean.gamma<-mean(gamma_chain)))	
#print(paste(" alphax=",mean.alphax<-mean(alphax_chain)," alphay=",mean.alphay<-mean(alphay_chain)))
#print(paste("interaction =",mean.delta<-mean(delta_chain)))	


#print(paste("posterior means beta1=",beta1.mean," beta2=",beta2.mean))	

#print(paste("correction bias q=",q_))
#print("Expected values")
#expected=
#rbind(beta1.mean*c(q_[1]*mean.alphax		   ,q_[2]*mean.alphay,			q_[3]),
#	  beta2.mean*c(q_[1]*mean.gamma*mean.alphax,q_[2]*mean.gamma*mean.alphay*mean.delta,q_[3]))

#print(expected)
#expected.countswocorrection=
#rbind(beta1.mean*c(mean.alphax,				mean.alphay,					1),
#	  beta2.mean*c(mean.alphax*mean.gamma,	mean.gamma*mean.alphay*mean.delta,1))
	  
##checking convergence	
if(plots==TRUE){
	par(ask=TRUE)
	plot(mcmc(cbind(gamma_chain)))
	plot(mcmc(cbind(alpha_chain,delta_chain)))
	plot(mcmc(cbind(beta11_chain,beta21_chain)))
	plot(mcmc(theta_chain))
	plot(mcmc(cbind(nux11_chain,nuy11_chain,phi_chain)))	
	par(ask=FALSE)
}
me=			apply(theta_chain,MARGIN=2,mean)
q025=		apply(theta_chain,MARGIN=2,function(xx){quantile(xx,0.025)})
q975=		apply(theta_chain,MARGIN=2,function(xx){quantile(xx,0.975)})


#print("estimated props")

tem=rbind(me[c("11","12")],me[c("21","22")]);

rownames(tem)=c("row1","row2")
colnames(tem)=c("x","y")
#print(tem);
#print(round(
estimated.props<-addmargins(tem)
#,3))



theta_q025=rbind(q025[c("11","12")],q025[c("21","22")])
theta_q975=rbind(q975[c("11","12")],q975[c("21","22")])


summaryt=matrix(NA,nrow=6,ncol=6)
colnames(summaryt)=c("sampleprop","mean","q_025","q_975","p_value","flagdifonehalf")
rownames(summaryt)=c("d1","1d","11o1d","21o2d","11od1","12od2")

td1=theta_chain[,"11"]+theta_chain[,"21"]
t1d=theta_chain[,"11"]+theta_chain[,"12"]
tem=cbind(	td1,
			t1d,
			theta_chain[,"11"]/t1d,
			theta_chain[,"21"]/(theta_chain[,"21"]+theta_chain[,"22"]),
			theta_chain[,"11"]/td1,
			theta_chain[,"12"]/(theta_chain[,"12"]+theta_chain[,"22"])
			)
colnames(tem)=rownames(summaryt)
head(tem)

cts=cbind(c(sum(x1),sum(x2)),c(sum(y1),sum(y2)))			
sp=c(
sum(cts[,1])/sum(cts),
sum(cts[1,])/sum(cts),
cts[1,1]/sum(cts[1,]),			
cts[2,1]/sum(cts[2,]),
cts[1,1]/sum(cts[,1]),
cts[1,2]/sum(cts[,2]))
			
for(i in 1:6){
	tem2=tem[,i]	
	summaryt[i,"sampleprop"]=		sp[i]
	summaryt[i,"mean"]=				mean(tem2)
	summaryt[i,c("q_025","q_975")]=	quantile(tem2,c(0.025,0.975))
	pval=	2*min(c(mean(tem2<0.5),mean(tem2>0.5)))
	summaryt[i,"p_value"]=pval
	summaryt[i,"flagdifonehalf"]=ifelse(pval<0.05,1,0)	
}

Bayes.pval_independence=2*min(c(mean(delta_chain<1),mean(delta_chain>1)))

op <- options(warn = (-1)) # suppress warnings
independenceChisqpvalue=chisq.test(x=totals)$p.value
options(op) # reset the default value, if you want

out=list(
tests=summaryt,
counts=totalswboth,
independencetest=list(Chisqpvalue=independenceChisqpvalue,bayesianpvalue=Bayes.pval_independence),
proportions=list(sample=sample.props,estimated=estimated.props,estimated_q025=theta_q025,estimated_q975=theta_q975)
)
return(out)
	}
