#Functions for model 
#This model is to test AI
#The model is the following
#y_i: Counts from paternal reference for rep i
#x_i: counts from maternal reference for rep i
#alpha: Treatment effect
#mu: overall mean. this is a nuisance parameter
#beta_i: replication i specific effect
#s_i:  Abundance: not random
#q:  not random, bias computed with simulation
#Model
#y_i~Poisson(mu alpha beta_i s_i q)
#x_i~Poisson(mu beta_i s_i (1-q))
#mu  ~gamma(a_mu,b_mu)
#alpha ~gamma(a_mu=1,b_mu=1)
#beta_i ~gamma(a_beta=1,b_beta=1) for i=1,...,I


x=c(21,43,32)
y=c(12,54,23)
both=c(32,65,43)

#It computes the abundances
#both: is the number of counts that align to both references
abundance=function(x,y,both){
    (x+y+both)/mean(x+y+both)
}

#It samples mu given the rest of the parameters
#sumcounts:
#q_: controls for bias
#when we have DNA information
#q_ is the proportion of ''DNA'' reads 
#from the paternal reads
#sumcounts: total number of RNA reads \sum_i (x_i+y_i)
#s: abundances

r_mu=function(sumcounts,alpha,betas,s,q_,a_mu,b_mu){
    a_post= a_mu+sumcounts 
    b_post= b_mu+(alpha*q_+1-q_)*sum(s*betas) 
    rgamma(1,shape=a_post,rate=b_post) 
}

#It samples alpha
#sumy: sum over all paternal reference reads=sum_i y_i
ralpha=function(sumy,mu,betas,s,q_,a_alpha,b_alpha){
    a_post= a_alpha+sumy
    b_post= b_alpha+mu*q_*sum(s*betas) 
    rgamma(1,shape=a_post,rate=b_post)  
}

#It samples the vector of betas
rbetas=function(x,y,alpha,mu,s,q_,a_beta,b_beta){
    betas=rep(NA,length(x))
    for(i in 1:length(betas)){
        a_post=  a_beta+x[i]+y[i]
        b_post=  b_beta+mu*s[i]*(alpha*q_+1-q_) 
        betas[i]= rgamma(1,shape=a_post,rate=b_post) 
    } 
    return(betas) 
}

#It performs n steps in the object all
#It updates all the parameters of the model
#the object all is a list that contains 
#the different parameters of the model and observed values
gibbs_nsteps=function(nsteps,all){
    for(i in 1:nsteps){
        all$alpha= ralpha(all$sumy,all$mu,all$betas,all$s,all$q,all$a_alpha,all$b_alpha)
        all$mu=  r_mu(all$sumcounts,all$alpha,all$betas,all$s,all$q,all$a_mu,all$b_mu)
        all$betas= rbetas(all$x,all$y,all$alpha,all$mu,all$s,all$q,all$a_beta,all$b_beta)
    } 
    return(all)
}


##Simulation from the Poisson-Gamma model
#with q_ fixed
#nsim: number of MCMC samples to be generated
#nburnin: Gibbs iterations to be burn and discard
#lag: the number of gibbs steps between sampled 
#imputed values
#it returns a list of MCMC posterior samples of
#the imputed values of the model
gibbs_poissongamma=function(nsim=10,nburnin=10,lag=10,
                            x,y,both,
                            a_mu=1/2,
                            b_mu=1/2,
                            a_alpha=1/2,
                            b_alpha=1/2,
                            a_beta=1/2,
                            b_beta=1/2,
                            q_=0.5,
                            abundance=TRUE
                            ){
    if(abundance==TRUE){
        s=   abundance(x,y,both)
        print("considering abundances")
        print(s)
    }else{
        s=  rep(1,length(x))
        #print("considering abundances all equal to 1")
    }
    sumcounts= sum(c(x,y))
    sumy=  sum(y)
    I=   length(x)
    mu=   sumcounts/(2*I)   
    alpha=  1
    betas=  (y/(mu*alpha*s*q_)+x/(mu*s*(1-q_)))/2

    all=list(x=x,y=y,sumcounts=sumcounts,sumy=sumy,mu=mu,alpha=alpha,betas=betas,
    q=q_,s=s,a_mu=a_mu,b_mu=b_mu,a_alpha=a_alpha,b_alpha=b_alpha,a_beta=a_beta,b_beta=b_beta)

    all=gibbs_nsteps(nburnin,all)

    mus_chain=alphas_chain= rep(NA,nsim)
    betas_chain=   matrix(NA,ncol=length(x),nrow=nsim)
    for(i in 1:nsim){
        all=   gibbs_nsteps(lag,all)
        mus_chain[i]= all$mu
        betas_chain[i,]=all$betas
        alphas_chain[i]=all$alpha
    }
    list(alphas=alphas_chain,mus=mus_chain,betas=betas_chain)
}

#gibbs(nsim=10,nburnin=10,lag=10,
#x,y,both,
#a_mu=1/2,
#b_mu=1/2,
#a_alpha=1/2,
#b_alpha=1/2,
#a_beta=1/2,
#b_beta=1/2,
#q_=0.5)

#Here we assume a random distribution for the now hyper-parameter 
#q. We estimate q from the DNA and insert it into the 
#RNA model. The RNA model is the same as in Rita's paper
#x_DNA: vector of DNA reads from the maternal reference
#y_DNA: vector of DNA reads from the paternal reference
#q:  starting value of q for the Gibbs sampler
gibbs_n_steps_for_q=function(n_steps,q,x_DNA,y_DNA,sigma_MH=0.1,v=1){
    p=q
    accept_counter_p=0
    m_i=x_DNA+y_DNA
    kd=sum(x_DNA)
    for(ii in 1:n_steps){
        gamma_i <- r_eta(nx=m_i,rep(p,length(y)))
        tem=r_p_MH(p=p,theta=NA,kd=kd,gamma_i=gamma_i,tt=NA,v,sigma_MH=sigma_MH)
        p=  tem$p
        accept_counter_p<-accept_counter_p+tem$accept 
    }
    return(list(acceptrate=accept_counter_p/n_steps,q=p))
}
#gibbs_n_steps_for_q(n_steps=100,q=0.1,x_DNA=c(20,20),y_DNA=c(100,100),sigma_MH=0.1,v=1)

#x=sample(200:300,3)
#y=sample(200:300,3)
#x_DNA=sample(200:300,3)
#y_DNA=sample(200:300,3)
#Here we assume the Poisson Gamma model but q is random and sampled 
#from the posterior distribution of q
#in the model
#y_DNA|x_DNA,q~Negative binomial(x_DNA,q)
#q~beta(v,v) 

#Simulation from Poisson Gamma model with q random
#It can either simulate q
#from the model y_XNA|x_DNA NegBinomial(y_DNA,q), p~Beta(v,v)
#q: plays the role of paternal alleles
#v: is the parameter of the beta distribution 
#when v is 1, q has a uniform prior distribution
#When qs=vector ther qs are not samples but instead
#the sample i of the imputed parameter in the gibbs
#step i is conditional on the value of qs[i]
#It returns a list of vectors of the imputed parameters
gibbs_poissongamma_q_random=function(nsim=10,nburnin=10,lag=10,
                                     x,y,both,x_DNA,y_DNA,
                                     a_mu=1/2,
                                     b_mu=1/2,
                                     a_alpha=1/2,
                                     b_alpha=1/2,
                                     a_beta=1/2,
                                     b_beta=1/2,
                                     abundance=TRUE, 
                                     sigma_MH=0.1,  #SD of the normal dist to simulate from p in the MH step
                                     v=1,    #parameter of the beta prior distribution for p 
                                     qs=NULL    #vector of qs to use in the model     #If they are not given they are sampled from the NB DNA model
                                    ){

    if(length(qs)<nsim & length(qs)!=0){print("error: length(qs) smaller than n_sim")}

    if(abundance==TRUE){
        s=   abundance(x,y,both)
        print("considering abundances")
        print(s)
    }else{
        s=  rep(1,length(x))
        #print("considering abundances all equal to 1")
    }
    sumcounts= sum(c(x,y))
    sumy=  sum(y)
    I=   length(x)
    mu=   sumcounts/(2*I)   
    alpha=  1

    if(length(qs)==0){
        q_=sum(x_DNA)/sum(x_DNA+y_DNA)
        accepted_MH_p<-0

        while(accepted_MH_p<0.25 | accepted_MH_p>0.75){
            tem= gibbs_n_steps_for_q(n_steps=100,q_,x_DNA,y_DNA,sigma_MH=sigma_MH,v)
            accepted_MH_p=tem$acceptrate
            if(accepted_MH_p<0.25){sigma_MH<-max(sigma_MH-0.01,0.001)}
            if(accepted_MH_p>0.75){sigma_MH<-sigma_MH+0.01}
            q_=tem$q
        }

        q_ =gibbs_n_steps_for_q(n_steps=nburnin,q_,x_DNA,y_DNA,sigma_MH=sigma_MH,v)$q

    }else{q_=mean(qs);}

    betas=  (y/(mu*alpha*s*q_)+x/(mu*s*(1-q_)))/2

    all=list(x=x,y=y,sumcounts=sumcounts,sumy=sumy,mu=mu,alpha=alpha,betas=betas,
    q=q_,s=s,a_mu=a_mu,b_mu=b_mu,a_alpha=a_alpha,b_alpha=b_alpha,a_beta=a_beta,b_beta=b_beta)

    all=gibbs_nsteps(nburnin,all)

    mus_chain=alphas_chain=qs_chain= rep(NA,nsim)
    betas_chain=   matrix(NA,ncol=length(x),nrow=nsim)

    for(i in 1:nsim){
        all$q=   ifelse(length(qs)==0,gibbs_n_steps_for_q(lag,all$q,x_DNA,y_DNA,sigma_MH,v)$q,qs[i])
        all=   gibbs_nsteps(lag,all)
        mus_chain[i]= all$mu
        betas_chain[i,]=all$betas
        alphas_chain[i]=all$alpha
        qs_chain[i]= all$q
    }
    list(alphas=alphas_chain,mus=mus_chain,betas=betas_chain,qs=qs_chain)
}
