library(stabledist)
library("libstable4u")
library(truncnorm)
library("TruncatedNormal")
library(invgamma)
library(MASS)
library(stats)
library(Rcpp)
library(RcppArmadillo)
library(popbio)
library(mvtnorm)
library("splines2")
library(ggplot2)
library(viridis)
library(patchwork)
library(MatrixExtra)
library("Rfast")
library(MatrixExtra)
#library(parallel)
library(mgcv)

# set seed -----------------------------------------------------------------

# for some reason, I set seed number twice.

set.seed(1)

# load data ---------------------------------------------------------------
load("/wind_speed_subset.RData")
#load("/Users/myungsooyoo/Desktop/My_stuff/research/extreme/windstorm/data/wind_speed_subset.RData")
time.length<-ncol(data.mat)
N<-nrow(data.mat)

H_mat<-diag(N)
n_data<-N

# data transform --------------------------------------------------------
data.mat<-notLog( notExp( (data.mat) ) -1 )
locationwise_mean<-mean(data.mat,na.rm=TRUE)

## mean centerering
data.mat<-(data.mat-locationwise_mean)

sigmoid_trans<-function(x){
  result<-1/(1+exp(-x))
  return(result)
}

# b spline ----------------------------------------------------------------
## first make P using spline and extract the coefficient
x<-seq(0,1,length.out=time.length)
aaa<-bSpline(x, degree=2,df=10,intercept = TRUE)
phi_spline_mat <- aaa
q<-dim(phi_spline_mat)[2]
rm(aaa)


# load basis --------------------------------------------------------------
load("basis_100.RData")
#load("/Users/myungsooyoo/Desktop/My_stuff/research/extreme/simul/simul_feb2025/windstorm_basis/basis_100.RData")
J<-ncol(Phi_mat)
M_mat_base<-matrix(0,J,J)

## below is to make transition matrix SPARSE
temp<-250
idx<-which(dist_btw_knot<=temp,arr.ind =TRUE)

M_mat2<-matrix(0,nrow=J,ncol=J)
for(i in 1:nrow(idx)){
  M_mat2[idx[i,1],idx[i,2]]<-1
}

for(j in 1:nrow(M_mat2)){
  print(j)
  print(which(M_mat2[j,]==1))
}
M_mat3<-M_mat2
positions<- which(M_mat3==1,arr.ind = TRUE)
M_mat3<- cbind(positions,NA)
M_mat3<-as.data.frame(M_mat3)
for(j in 1: nrow(M_mat3)){
  if(M_mat3[j,1]==M_mat3[j,2]){
    M_mat3[j,3]<-"a"
  }
}

for(j in 1:nrow(M_mat3)){
  if(is.na(M_mat3[j,3]) ){
    row_temp<-M_mat3[j,1]
    col_temp<-M_mat3[j,2]
    if(row_temp==col_temp -1){
      M_mat3[j,3]<-"b"
    } else if (row_temp==col_temp+1){
      M_mat3[j,3]<-"d"
    } else if (row_temp> col_temp){
      M_mat3[j,3]<-"e"
    } else{
      M_mat3[j,3]<-"c"
    }
  }
}

M_mat2_2<-M_mat2
for(j in 1: nrow(M_mat3)){
  M_mat2_2[M_mat3[j,1],M_mat3[j,2]]<- M_mat3[j,3]
}

idx_a<-which(M_mat2_2=="a")
idx_b<-which(M_mat2_2=="b")
idx_c<-which(M_mat2_2=="c")
idx_d<-which(M_mat2_2=="d")
idx_e<-which(M_mat2_2=="e")
idx_zero<-which(M_mat2_2==0)

rm(idx,M_mat2,M_mat2_2,M_mat3)
rm(positions,col_temp,i)
rm(row_temp,temp)
rm(j)
rm(x)
# functions ----------------------------------------------------------
log_mixture_density2<-function(residual,
                               lambda1,
                               kappa1,
                               nu1,
                               lambda2,
                               kappa2,
                               nu2,
                               prob){
  pars_1<-c(lambda1,kappa1,nu1,0)## alpha, beta, sigma, mu
  temp1<-stable_pdf(residual,pars_1) 
  temp0<-dnorm(residual,mean=0,sd=nu2*sqrt(2))
  result<-prob*temp1 +(1-prob)*temp0
  result<-log(result)
  return(result)
}
log_mixture_density_all3<-function(residual,
                                   lambda1,
                                   kappa1,
                                   nu1,
                                   nu2,
                                   prob){
  pars_1<-c(lambda1,kappa1,nu1,0)## alpha, beta, sigma, mu
  temp1<-stable_pdf(residual,pars_1) 
  temp0<-dnorm(residual,mean=0,sd=nu2*sqrt(2))
  result<-prob * temp1 +(1-prob)*temp0
  result<-log(result)
  return(result)
}

cppFunction('SEXP residual_vec(arma::vec alpha_vec, arma::mat M_mat, arma::vec alpha_vec_before) {
            arma::vec result = alpha_vec-M_mat * alpha_vec_before;
            return Rcpp::wrap( result );}',
            depends='RcppArmadillo')

cppFunction('SEXP residual_mat_rcpp(arma::mat alpha_mat, arma::mat M_mat, arma::mat alpha_mat_before) {
            arma::mat result = alpha_mat- M_mat * alpha_mat_before;
            return Rcpp::wrap( result );}',
            depends='RcppArmadillo')

cppFunction('SEXP alpha_update(double J, double C_m, arma::mat y_mat, double length_t,
            double sigma_square_d,
            double epsilon_nuisance, 
            arma::mat Phi_mat, 
            arma::mat Phi_t_H_t_H_Phi,
            arma::mat Phi_t_H_t,
            arma::mat alpha_i_minus,
            arma::mat alpha_i,
            arma::mat delta_mat,
            arma::mat M_mat,
            arma::mat H_mat) {
            
            for(int t=1; t<=(length_t+1); ++t){
               if(t==1){
                arma::mat D=(M_mat.t() * M_mat )/ epsilon_nuisance +(1/C_m)*arma::eye( J,J );
                arma::vec b= ( M_mat.t()  *(alpha_i_minus.col(t)-delta_mat.col(t-1)) )/ epsilon_nuisance;
                arma::mat D_inverse =arma::inv_sympd(D);
                alpha_i.col(t-1)=arma::mvnrnd( D_inverse * b,D_inverse);
             
              } else if(t==(length_t+1)){
                arma::mat D= ( 
                  (1/epsilon_nuisance)*arma::eye( J,J ) +
                   (Phi_t_H_t_H_Phi * 1/sigma_square_d )
                  );
                arma::vec b= (
                  (M_mat * alpha_i.col(t-2))/epsilon_nuisance +
                  delta_mat.col(t-2)/epsilon_nuisance+ 
                  (Phi_t_H_t* (y_mat.col(t-2)  ) *1/sigma_square_d )
                );
                arma::mat D_inverse =arma::inv_sympd	(D);
                alpha_i.col(t-1)=arma::mvnrnd( D_inverse * b,D_inverse);

              } else  {
                arma::mat D= (
                  (1/epsilon_nuisance)*arma::eye( J,J )+
                   (M_mat.t() * M_mat)/epsilon_nuisance+
                  (Phi_t_H_t_H_Phi * 1/sigma_square_d )
                  );
                arma::vec b= ( 
                  (M_mat * alpha_i.col(t-2))/epsilon_nuisance +
                  delta_mat.col(t-2)/epsilon_nuisance+ 
                  ( M_mat.t()  * alpha_i_minus.col(t) ) /epsilon_nuisance -                 
                  ( M_mat.t()  * delta_mat.col(t-1) ) /epsilon_nuisance +                 
                  (Phi_t_H_t * ( y_mat.col(t-2) ) *1/sigma_square_d)
                  );
                arma::mat D_inverse =arma::inv_sympd	(D);
                alpha_i.col(t-1)=arma::mvnrnd( D_inverse * b,D_inverse);

              }
            }
            return Rcpp::wrap(alpha_i );}',
            depends='RcppArmadillo')


beta_function_new<-function(x,
                            current,
                            proposal,
                            proposal_beta_sd,
                            q,
                            current_log_likelihood_all,
                            candidate_log_likelihood_all){
  current_log_likelihood<- sum( dnorm(current[,x],mean=0,sd=sqrt(10^6),log=TRUE) )+
    current_log_likelihood_all[x]+
    sum(log(dtruncnorm(proposal[,x],a=-6,b=6,mean=current[,x],sd=proposal_beta_sd)))
  candidate_log_likelihood<-sum( dnorm(proposal[,x],mean=0,sd=sqrt(10^6),log=TRUE) )+
    candidate_log_likelihood_all[x]+
    sum(log(dtruncnorm(current[,x],a=-6,b=6,mean=proposal[,x],sd=proposal_beta_sd)))
  log_ratio<-candidate_log_likelihood-current_log_likelihood
  return(log_ratio)
}


stable_function<-function(x,
                          current_log_likelihood_all,
                          candidate_log_likelihood_all,
                          current_log_nu1,
                          current_log_nu2,
                          current_lambda_1,
                          current_kappa_1,
                          proposal_log_nu1,
                          proposal_log_nu2,
                          proposal_lambda_1,
                          proposal_kappa_1,
                          hyper_alpha_nu1,
                          hyper_beta_nu1,
                          hyper_alpha_nu2,
                          hyper_beta_nu2,
                          epsilon_epsilon,
                          proposal_sd_lambda,
                          proposal_sd_kappa){
  current_log_likelihood<-current_log_likelihood_all[x]+
    dinvgamma(exp(current_log_nu1[x]),hyper_alpha_nu1,rate=hyper_beta_nu1,log=TRUE)+current_log_nu1[x]+
    dinvgamma(exp(current_log_nu2[x]),hyper_alpha_nu2,rate=hyper_beta_nu2,log=TRUE)+current_log_nu2[x]+
    log( dtruncnorm(proposal_lambda_1[x],a=0+epsilon_epsilon,b=2-epsilon_epsilon,mean=current_lambda_1[x],sd=proposal_sd_lambda) )+
    log( dtruncnorm(proposal_kappa_1[x],a=-1+epsilon_epsilon,b=1-epsilon_epsilon,mean=current_kappa_1[x],sd=proposal_sd_kappa) )+
    log(dtruncnorm(current_lambda_1[x],a=0,b=2,mean=1,sd=0.14))
  candidate_log_likelihood<-candidate_log_likelihood_all[x]+
    dinvgamma(exp(proposal_log_nu1[x]),hyper_alpha_nu1,rate=hyper_beta_nu1,log=TRUE)+proposal_log_nu1[x]+
    dinvgamma(exp(proposal_log_nu2[x]),hyper_alpha_nu2,rate=hyper_beta_nu2,log=TRUE)+proposal_log_nu2[x]+
    log( dtruncnorm(current_lambda_1[x],a=0+epsilon_epsilon,b=2-epsilon_epsilon,mean=proposal_lambda_1[x],sd=proposal_sd_lambda) )+
    log( dtruncnorm(current_kappa_1[x],a=-1+epsilon_epsilon,b=1-epsilon_epsilon,mean=proposal_kappa_1[x],sd=proposal_sd_kappa) )+
    log(dtruncnorm(proposal_lambda_1[x],a=0,b=2,mean=1,sd=0.14))
  
  log_ratio<-candidate_log_likelihood-current_log_likelihood
  return(log_ratio)
}

# mcmc parameter mat and initial ------------------------------------------
iterations<-20000
burnin<-10000

# storage -----------------------------------------------------------------

parameter.mat.sigma.square.d <-numeric(iterations)

## beta
parameter.mat.beta.list<-list()
for(i in 1: J){
  parameter.mat.beta.list[[i]]<-matrix(NA,q,iterations)
}

## alpha
post.alpha.mat<-list()
for(i in 1:length((burnin+1) : iterations) ){
  post.alpha.mat[[i]]<-matrix(NA,nrow=J,ncol=time.length+1)
}

## delta
post.delta.mat<-list()
for(i in 1:length((burnin+1) : iterations) ){
  post.delta.mat[[i]]<-matrix(NA,nrow=J,ncol=time.length)
}

## M_mat
post.a.mat<-matrix(NA,nrow=J,ncol=iterations)
post.b.mat<-numeric(iterations)
post.c.mat<-numeric(iterations)
post.d.mat<-numeric(iterations)
post.e.mat<-numeric(iterations)


### stable erros
lambda1_mat<-matrix(NA,nrow=J,ncol=iterations)
nu1_mat<-matrix(NA,nrow=J,ncol=iterations)
kappa1_mat<-matrix(NA,nrow=J,ncol=iterations)
nu2_mat<-matrix(NA,nrow=J,ncol=iterations)

count<-0
meanmean_pred<-matrix(0,nrow=N,ncol=time.length)
M2_mat<-matrix(0,nrow=N,ncol=time.length)




# initial and hyperparameter-----------------------------------------------------------------
## sigma square
hyper_alpha_d<-0.1
hyper_beta_d<-0.1

hyper_alpha_nu1<-2
hyper_beta_nu1<-100*(hyper_alpha_nu1-1)

hyper_alpha_nu2<-2
hyper_beta_nu2<-20*(hyper_alpha_nu2-1)

post.alpha.mat[[1]][,1]<-rnorm(J)
post.alpha.mat[[1]][,2:(time.length+1)]<-solve(t(H_mat %*% Phi_mat)%*%(H_mat %*% Phi_mat)) %*% t(H_mat %*%Phi_mat) %*% data.mat
###
post.alpha.mat.previous<-post.alpha.mat[[1]]

post.delta.mat[[1]]<-matrix(rnorm(J* (time.length)),nrow=J,ncol=time.length)
###
post.delta.mat.previous<-post.delta.mat[[1]]

lambda1_mat[,1]<-runif(J,1,1.5)
nu1_mat[,1]<-rep(100,J)
kappa1_mat[,1]<-runif(J,-1,1)
nu2_mat[,1]<-rep(0.5,J)

for(j in 1:J){
  parameter.mat.beta.list[[j]][,1]<-rnorm(q)
}


post.a.mat[,1]<-1
post.b.mat[1]<-1
post.c.mat[1]<-1
post.d.mat[1]<-1
post.e.mat[1]<-1


# MH jump -----------------------------------------------------------------


proposal_sd_nu<-0.23
proposal_sd_lambda<-0.23
proposal_sd_kappa<-0.23
proposal_beta_sd<-1.8

proposal_delta_sd<-4.1


accept.nu_lambda_ratio<-matrix(0,nrow=J,ncol=iterations)
accept.beta.ratio<-matrix(0,nrow=J,ncol=iterations)

# mcmc start --------------------------------------------------------------

epsilon_epsilon<-0.0000000001
epsilon_nuisance<- 0.5 
shape<-hyper_alpha_d + ( (time.length) *n_data)/2

Diagonal_J<-Diagonal(J)

#### storage for the predict.time
Phi_t_H_t_H_Phi<- t(Phi_mat)%*% t(H_mat)%*%H_mat %*% Phi_mat
Phi_t_H_t<- t(Phi_mat)%*% t(H_mat)

print("mcmc start")


for(j in 2:iterations){

  parameter.mat.sigma.square.d[j]<-rinvgamma(1, shape,rate=(sum( ( (data.mat-H_mat%*%Phi_mat%*%post.alpha.mat.previous[,-1]) )^2 ) +2*hyper_beta_d )/2 )

  ## a vec
  M_rest<-M_mat_base
  M_rest[idx_b]<-post.b.mat[j-1]
  M_rest[idx_c]<-post.c.mat[j-1]
  M_rest[idx_d]<-post.d.mat[j-1]
  M_rest[idx_e]<-post.e.mat[j-1]
  b.temp<-numeric(J)
  D.temp<-matrix(0,J,J)
  for(t in 2: (time.length+1)){
    b.temptemp<-diag(post.alpha.mat.previous[,t-1]) %*% (post.alpha.mat.previous[,t]-M_rest %*% post.alpha.mat.previous[,t-1] - post.delta.mat.previous[,t-1])
    b.temptemp <- b.temptemp /epsilon_nuisance
    b.temp<-b.temp+b.temptemp
    D.temptemp<- diag( (post.alpha.mat.previous[,t-1] )^2 ) 
    D.temptemp<-D.temptemp/epsilon_nuisance
    D.temp<-D.temp+D.temptemp
  }
  D.temp<-D.temp+ 1/1000* diag(J)
  D_inv<-solve(D.temp)
  post.a.mat[,j]<- c(rmvnorm(1, mu= D_inv %*% b.temp, sigma=D_inv))
  ## b
  M_rest<-M_mat_base
  M_rest[idx_a]<-post.a.mat[,j]
  M_rest[idx_c]<-post.c.mat[j-1]
  M_rest[idx_d]<-post.d.mat[j-1]
  M_rest[idx_e]<-post.e.mat[j-1]
  design_mat<-M_mat_base
  design_mat[idx_b]<-1
  b.temp<-0
  D.temp<-0
  design_mat_cross_prod<-t(design_mat) %*% design_mat
  for(t in 2: (time.length+1)){
    b.temptemp<-t( post.alpha.mat.previous[,t-1] ) %*% t(design_mat) %*%( 
      post.alpha.mat.previous[,t]-M_rest %*% post.alpha.mat.previous[,t-1] -post.delta.mat.previous[,t-1]
    )
    b.temptemp <- b.temptemp /epsilon_nuisance
    b.temp<-b.temp+b.temptemp
    #
    d.temptemp<-t( post.alpha.mat.previous[,t-1] ) %*% design_mat_cross_prod%*%post.alpha.mat.previous[,t-1] 
    d.temptemp<-d.temptemp/epsilon_nuisance
    D.temp<-D.temp+d.temptemp
  }
  D.temp<- D.temp+ 1/1000
  post.b.mat[j]<-rnorm(1, mean=b.temp/D.temp,sd= sqrt(1/D.temp))
  ## c
  M_rest<-M_mat_base
  M_rest[idx_a]<-post.a.mat[,j]
  M_rest[idx_b]<-post.b.mat[j]
  M_rest[idx_d]<-post.d.mat[j-1]
  M_rest[idx_e]<-post.e.mat[j-1]
  design_mat<-M_mat_base
  design_mat[idx_c]<-1
  b.temp<-0
  D.temp<-0
  design_mat_cross_prod<-t(design_mat) %*% design_mat
  for(t in 2: (time.length+1)){
    b.temptemp<-t( post.alpha.mat.previous[,t-1] ) %*% t(design_mat) %*%( 
      post.alpha.mat.previous[,t]-M_rest %*% post.alpha.mat.previous[,t-1] -post.delta.mat.previous[,t-1]
    )
    b.temptemp <- b.temptemp /epsilon_nuisance
    b.temp<-b.temp+b.temptemp
    #
    d.temptemp<-t( post.alpha.mat.previous[,t-1] ) %*% design_mat_cross_prod%*%post.alpha.mat.previous[,t-1] 
    d.temptemp<-d.temptemp/epsilon_nuisance
    D.temp<-D.temp+d.temptemp
  }
  D.temp<- D.temp+ 1/1000
  post.c.mat[j]<-rnorm(1, mean=b.temp/D.temp,sd= sqrt(1/D.temp))
  ## d
  M_rest<-M_mat_base
  M_rest[idx_a]<-post.a.mat[,j]
  M_rest[idx_b]<-post.b.mat[j]
  M_rest[idx_c]<-post.c.mat[j]
  M_rest[idx_e]<-post.e.mat[j-1]
  design_mat<-M_mat_base
  design_mat[idx_d]<-1
  b.temp<-0
  D.temp<-0
  design_mat_cross_prod<-t(design_mat) %*% design_mat
  
  for(t in 2: (time.length+1)){
    b.temptemp<-t( post.alpha.mat.previous[,t-1] ) %*% t(design_mat) %*%( 
      post.alpha.mat.previous[,t]-M_rest %*% post.alpha.mat.previous[,t-1] -post.delta.mat.previous[,t-1]
    )
    b.temptemp <- b.temptemp /epsilon_nuisance
    b.temp<-b.temp+b.temptemp
    #
    d.temptemp<-t( post.alpha.mat.previous[,t-1] ) %*% design_mat_cross_prod%*%post.alpha.mat.previous[,t-1] 
    d.temptemp<-d.temptemp/epsilon_nuisance
    D.temp<-D.temp+d.temptemp
  }
  D.temp<- D.temp+ 1/1000
  post.d.mat[j]<-rnorm(1, mean=b.temp/D.temp,sd= sqrt(1/D.temp))
  ## e
  M_rest<-M_mat_base
  M_rest[idx_a]<-post.a.mat[,j]
  M_rest[idx_b]<-post.b.mat[j]
  M_rest[idx_c]<-post.c.mat[j]
  M_rest[idx_d]<-post.d.mat[j]
  design_mat<-M_mat_base
  design_mat[idx_e]<-1
  b.temp<-0
  D.temp<-0
  design_mat_cross_prod<-t(design_mat) %*% design_mat
  
  for(t in 2: (time.length+1)){
    b.temptemp<-t( post.alpha.mat.previous[,t-1] ) %*% t(design_mat) %*%( 
      post.alpha.mat.previous[,t]-M_rest %*% post.alpha.mat.previous[,t-1] -post.delta.mat.previous[,t-1]
    )
    b.temptemp <- b.temptemp /epsilon_nuisance
    b.temp<-b.temp+b.temptemp
    #
    d.temptemp<-t( post.alpha.mat.previous[,t-1] ) %*% design_mat_cross_prod%*%post.alpha.mat.previous[,t-1] 
    d.temptemp<-d.temptemp/epsilon_nuisance
    D.temp<-D.temp+d.temptemp
  }
  D.temp<- D.temp+ 1/1000
  post.e.mat[j]<-rnorm(1, mean=b.temp/D.temp,sd= sqrt(1/D.temp))
  ## current_M_mat
  current_M_mat<-M_mat_base
  current_M_mat[idx_a]<-post.a.mat[,j]
  current_M_mat[idx_b]<-post.b.mat[j]
  current_M_mat[idx_c]<-post.c.mat[j]
  current_M_mat[idx_d]<-post.d.mat[j]
  current_M_mat[idx_e]<-post.e.mat[j]
  # alpha:gibbs
  post.alpha.mat.current<-matrix(0,J,time.length+1)
  post.alpha.mat.current<-alpha_update(J=J,
                                       C_m=1000,
                                       y_mat=data.mat,
                                       length_t=time.length,
                                       sigma_square_d=parameter.mat.sigma.square.d[j],
                                       epsilon_nuisance = epsilon_nuisance,
                                       Phi_mat=Phi_mat,
                                       Phi_t_H_t_H_Phi=Phi_t_H_t_H_Phi,
                                       Phi_t_H_t=Phi_t_H_t,
                                       alpha_i_minus=post.alpha.mat.previous,
                                       alpha_i=post.alpha.mat.current,
                                       delta_mat=post.delta.mat.previous,
                                       M_mat=current_M_mat,
                                       H_mat=H_mat)
  if(j >= (burnin+1)){
    post.alpha.mat[[j-burnin]]<-post.alpha.mat.current
  }
  
  ### delta_t
  ### delta_t
  ### delta_t
  ### delta_t
  current_p_mat<-t( simplify2array(lapply(1:J, function(x) c(sigmoid_trans(phi_spline_mat %*%parameter.mat.beta.list[[x]][,j-1]))) ) )
  
  proposal<-matrix(rnorm(J*time.length,mean=c(post.delta.mat.previous),sd=proposal_delta_sd),J,time.length)
  
  current_log_likelihood<-matrix( dnorm(c(post.alpha.mat.current[,-1]),
                                        mean=c(current_M_mat%*% post.alpha.mat.current[,-c(time.length+1)]+post.delta.mat.previous),
                                        sd=sqrt(epsilon_nuisance),
                                        log=TRUE),
                                  J,
                                  time.length
  )
  candidate_log_likelihood<-matrix( dnorm(c(post.alpha.mat.current[,-1]),
                                          mean=c(current_M_mat%*% post.alpha.mat.current[,-c(time.length+1)]+proposal),
                                          sd=sqrt(epsilon_nuisance),
                                          log=TRUE),
                                    J,
                                    time.length
  )
  current_prior<-t(mapply(log_mixture_density_all3,residual=asplit(post.delta.mat.previous,1),
                          lambda1=as.list(lambda1_mat[,j-1]),
                          kappa1=as.list(kappa1_mat[,j-1]),
                          nu1=as.list(nu1_mat[,j-1]),
                          nu2=as.list(nu2_mat[,j-1]),
                          prob=asplit(current_p_mat,1)))

  candidate_prior<-t(mapply(log_mixture_density_all3,residual=asplit(proposal,1),
                            lambda1=as.list(lambda1_mat[,j-1]),
                            kappa1=as.list(kappa1_mat[,j-1]),
                            nu1=as.list(nu1_mat[,j-1]),
                            nu2=as.list(nu2_mat[,j-1]),
                            prob=asplit(current_p_mat,1)))
  current_log_likelihood<-current_log_likelihood+current_prior
  candidate_log_likelihood<-candidate_log_likelihood+candidate_prior
  log_ratio<-c(candidate_log_likelihood-current_log_likelihood)
  update_idx<-which(log(runif(J*time.length)) <=log_ratio)
  post.delta.mat.current<-post.delta.mat.previous
  post.delta.mat.current[update_idx]<-c(proposal)[update_idx]
  if(j >=(burnin+1)){
    post.delta.mat[[j-burnin]]<-post.delta.mat.current
  }
  ### beta
  ### beta
  ### beta
  ### beta
  proposal<-  simplify2array(lapply(1:J, function(x) rtruncnorm(1, a=-6,b=6,
                                                                mean=parameter.mat.beta.list[[x]][,j-1],
                                                                sd=proposal_beta_sd))
  )
  
  current<- simplify2array(lapply(1:J, function(x) parameter.mat.beta.list[[x]][,j-1]))
  
  proposal_p <-sigmoid_trans(phi_spline_mat%*%proposal)
  current_p <- sigmoid_trans(phi_spline_mat%*%current) 
  ## current
  current_log_likelihood_all<-apply(mapply(log_mixture_density_all3,residual=asplit(post.delta.mat.current,1),
                                           lambda1=as.list(lambda1_mat[,j-1]),
                                           kappa1=as.list(kappa1_mat[,j-1]),
                                           nu1=as.list(nu1_mat[,j-1]),
                                           nu2=as.list(nu2_mat[,j-1]),
                                           prob=asplit(current_p,2))
                                    ,2,sum)
  ## candidate

  candidate_log_likelihood_all<-apply(mapply(log_mixture_density_all3,residual=asplit(post.delta.mat.current,1),
                                             lambda1=as.list(lambda1_mat[,j-1]),
                                             kappa1=as.list(kappa1_mat[,j-1]),
                                             nu1=as.list(nu1_mat[,j-1]),
                                             nu2=as.list(nu2_mat[,j-1]),
                                             prob=asplit(proposal_p,2))
                                      ,2,sum)
  
  log_ratio<-unlist(lapply(X=1:J, beta_function_new,current=current,proposal=proposal,
                           proposal_beta_sd=proposal_beta_sd,
                           q=q,
                           current_log_likelihood_all=current_log_likelihood_all,
                           candidate_log_likelihood_all=candidate_log_likelihood_all))
  update_idx<-which(log(runif(J)) <=log_ratio)
  for( i in update_idx){
    parameter.mat.beta.list[[i]][,j]<-proposal[,i]
    accept.beta.ratio[i,j]<-1
  }
  for( i in which( (1:J) %in% update_idx==FALSE ) ){
    parameter.mat.beta.list[[i]][,j]<-current[,i]
  }
  ## stable errors
  ## stable errors
  ## stable errors
  ## stable errors
  ## stable errors
  ## lambda_1, nu_1, lambda_2, nu_2
  current_p_mat<-t( simplify2array(lapply(1:J, function(x) c(sigmoid_trans(phi_spline_mat %*%parameter.mat.beta.list[[x]][,j]))) ) )
  
  ## proposal
  
  proposal_lambda_1<-unlist(lapply(1:J, function(x) rtruncnorm(1, a=0+epsilon_epsilon,b=2-epsilon_epsilon,mean=lambda1_mat[x,j-1],sd=proposal_sd_lambda)))
  proposal_kappa_1<-unlist(lapply(1:J, function(x) rtruncnorm(1, a=-1+epsilon_epsilon,b=1-epsilon_epsilon,mean=kappa1_mat[x,j-1],sd=proposal_sd_kappa)))
  proposal_log_nu1<-unlist(lapply(1:J, function(x) rnorm(1,mean=log(nu1_mat[x,j-1]), sd=proposal_sd_nu)))
  proposal_log_nu2<-unlist(lapply(1:J, function(x) rnorm(1,mean=log(nu2_mat[x,j-1]), sd=proposal_sd_nu)))
  
  current_lambda_1<-unlist( lapply(1:J, function(x) lambda1_mat[x,j-1]))
  current_kappa_1<-unlist( lapply(1:J, function(x) kappa1_mat[x,j-1]))
  current_log_nu1<-unlist( lapply(1:J, function(x) log(nu1_mat[x,j-1])))
  current_log_nu2<-unlist( lapply(1:J, function(x) log(nu2_mat[x,j-1])))
  # current
  current_log_likelihood_all <-apply(mapply(log_mixture_density_all3,residual=asplit(post.delta.mat.current,1),
                                            lambda1=as.list(current_lambda_1),
                                            kappa1=as.list(current_kappa_1),
                                            nu1=as.list(exp(current_log_nu1)),
                                            nu2=as.list(exp(current_log_nu2)),
                                            prob=asplit(current_p_mat,1)),2,sum)
  candidate_log_likelihood_all <-apply(mapply(log_mixture_density_all3,residual=asplit(post.delta.mat.current,1),
                                              lambda1=as.list(proposal_lambda_1),
                                              kappa1=as.list(proposal_kappa_1),
                                              nu1=as.list(exp(proposal_log_nu1)),
                                              nu2=as.list(exp(proposal_log_nu2)),
                                              prob=asplit(current_p_mat,1)),2,sum)
  
  ###

  log_ratio<-unlist(lapply(1:J,stable_function, current_log_likelihood_all=current_log_likelihood_all,
                           candidate_log_likelihood_all=candidate_log_likelihood_all,
                           current_log_nu1=current_log_nu1,
                           current_log_nu2=current_log_nu2,
                           current_lambda_1=current_lambda_1,
                           current_kappa_1=current_kappa_1,
                           proposal_log_nu1=proposal_log_nu1,
                           proposal_log_nu2=proposal_log_nu2,
                           proposal_lambda_1=proposal_lambda_1,
                           proposal_kappa_1=proposal_kappa_1,
                           hyper_alpha_nu1=hyper_alpha_nu1,
                           hyper_beta_nu1=hyper_beta_nu1,
                           hyper_alpha_nu2=hyper_alpha_nu2,
                           hyper_beta_nu2=hyper_beta_nu2,
                           epsilon_epsilon=epsilon_epsilon,
                           proposal_sd_lambda=proposal_sd_lambda,
                           proposal_sd_kappa)
  )
  
  update_idx<-which(log(runif(J)) <=log_ratio)
  
  nu1_mat[,j]<-exp(current_log_nu1)
  nu2_mat[,j]<-exp(current_log_nu2)
  lambda1_mat[,j]<-lambda1_mat[,j-1]
  kappa1_mat[,j]<-kappa1_mat[,j-1]
  nu1_mat[update_idx,j]<-exp(proposal_log_nu1[update_idx])
  nu2_mat[update_idx,j]<-exp(proposal_log_nu2[update_idx])
  lambda1_mat[update_idx,j]<-proposal_lambda_1[update_idx]
  kappa1_mat[update_idx,j]<-proposal_kappa_1[update_idx]
  accept.nu_lambda_ratio[update_idx,j]<-1
  ###
  post.alpha.mat.previous<-post.alpha.mat.current
  post.delta.mat.previous<-post.delta.mat.current
  ###
  print(j)
  print(parameter.mat.sigma.square.d[j])

  
  if(j %% 30==0){
    print(paste("lambda1=",lambda1_mat[1,j]))
    print(paste("nu1=",nu1_mat[1,j]))
    print(paste("nu2=",nu2_mat[1,j]))
    print(paste("hit ratio for nu1=", mean(accept.nu_lambda_ratio[1,j:(j-29)]) ) )
    print(paste("hit ratio for nu2=", mean(accept.nu_lambda_ratio[2,j:(j-29)]) ) )
    print(paste("hit ratio for beta1=", mean(accept.beta.ratio[1,j:(j-29)]) ) )
    print(paste("hit ratio for beta2=", mean(accept.beta.ratio[2,j:(j-29)]) ) )

  }
  #### prediction
  if(j>=burnin){
    temptemp<-matrix(NA,N,time.length)
    for(t in 1: time.length){
      gaussian_error<-rnorm(N,mean=0,sd=sqrt(parameter.mat.sigma.square.d[j]))
      temptemp[,t]<-c(
        Phi_mat%*%post.alpha.mat.current[,t+1] +gaussian_error
      )
    }
    temp<-temptemp
    count<-count+1
    delta<- temp-meanmean_pred
    meanmean_pred<-meanmean_pred+delta/count
    delta2<-temp-meanmean_pred
    M2_mat <-M2_mat + delta *delta2
  }
}



# remove unnecessary objects------------------------------------------------------------------
rm(gaussian_error,current_log_nu2)
rm(current_log_likelihood,current_log_likelihood_all)
rm(candidate_log_likelihood,candidate_log_likelihood_all)
rm(post.alpha.mat.current,post.alpha.mat.previous)
rm(H_mat,delta,delta2)
rm(current_M_mat,current,current_p,current_p_mat)
rm(current_prior,candidate_prior)
rm(Diagonal_J)
rm(Phi_t_H_t,Phi_t_H_t_H_Phi)
rm(b.temp,b.temptemp,D_inv,D.temp,d.temptemp,D.temptemp)
rm(design_mat,design_mat_cross_prod)
rm(M_rest)
rm(post.delta.mat.current,post.delta.mat.previous)
rm(temptemp)
rm(temp)
rm(i,j)
rm(update_idx)
rm(proposal_log_nu2)
rm(log_ratio)
rm(M2_mat)
rm(proposal,proposal_p)
rm(current_kappa_1,current_lambda_1,current_log_nu1,epsilon_epsilon)
rm(proposal_sd_nu,proposal_sd_lambda,proposal_sd_kappa,proposal_log_nu1)
rm(proposal_lambda_1,proposal_kappa_1,proposal_delta_sd,proposal_beta_sd)

# save objects only after burnin -----------------------------------------------------------

### sigma square d
parameter.mat.sigma.square.d<-parameter.mat.sigma.square.d[-c(1:burnin)]
### lambda 1, stable parameters..

kappa1_mat<-kappa1_mat[,-c(1:burnin)]
lambda1_mat<-lambda1_mat[,-c(1:burnin)]
nu1_mat<-nu1_mat[,-c(1:burnin)]
nu2_mat<-nu2_mat[,-c(1:burnin)]


### beta
dim(parameter.mat.beta.list[[1]])
for(j in 1:J){
  parameter.mat.beta.list[[j]]<-parameter.mat.beta.list[[j]][,-c(1:burnin)]
}


## post M mat

post.a.mat<-post.a.mat[,-c(1:burnin)]
post.b.mat<-post.b.mat[-c(1:burnin)]
post.c.mat<-post.c.mat[-c(1:burnin)]
post.d.mat<-post.d.mat[-c(1:burnin)]
post.e.mat<-post.e.mat[-c(1:burnin)]
gc()

### delta

## alpha
lower_alpha<-apply(simplify2array(post.alpha.mat), 1:2, quantile, prob = c(0.025))
upper_alpha<-apply(simplify2array(post.alpha.mat), 1:2, quantile, prob = c(0.975))
meanmean_alpha<- mean.list(post.alpha.mat)
#rm(post.alpha.mat)

## accept ratio
accept.nu_lambda_ratio<-accept.nu_lambda_ratio[,-c(1:burnin)]
accept.beta.ratio<-accept.beta.ratio[,-c(1:burnin)]
### delta accept ratio
accept.delta.ratio<-matrix(1,nrow=J*time.length,ncol=length(post.delta.mat)-1)

k=2
for(k in 2: length(post.delta.mat)){
  idx<-which( c( post.delta.mat[[k]]) ==c( post.delta.mat[[k-1]]) ) ## same
  accept.delta.ratio[idx,k-1]<-0
}

apply(accept.delta.ratio,1,mean) 

apply(accept.nu_lambda_ratio,1,mean)
which(apply(accept.nu_lambda_ratio,1,mean)==0)
which(apply(accept.nu_lambda_ratio,1,mean)==0)


apply(accept.beta.ratio,1,mean)

# save --------------------------------------------------------------------
#### save .RData file
#save.image("windstorm_stable.RData")



