#'@title Calibrated Ratio Estimator under Double Sampling Design
#' @description Population ratio estimator under two-phase random sampling design has gained enormous popularity in present era. This package provides functions for estimation calibrated population ratio under two phase sampling design, including the approximate variance of the ratio estimator. The improved ratio estimator can be applicable for both the case, when auxiliary data is available at unit level or aggregate level for first phase sampled. Calibration weight of each unit of the second phase sample was calculated. Single and combined inclusion probabilities were also estimated for both phases under two phase random sampling. The improved ratio estimator's percentage coefficient of variation was also determined as a measure of accuracy.
#' @param N Population size
#' @param FSU First stage sampling units
#' @param SSU Second stage sampling units
#' @import MASS stats
#'
#' @return
#'  \itemize{
#'   \item CalEstimate: Estimate value of calibration estimator
#'   \item CalVariance: Variance of calibration estimator
#'   \item CV: Coefficient of variance
#'   \item SampleSize: Sample Size of FSU and SSU
#'   \item DesignWeight: Design weight vector
#'   \item InclusionProb: Inclusion probability vector
#'   \item Correlation: Correlation value
#' }
#'
#' @export
#'
#' @examples
#' f1<-rnorm(100,20,5)
#' f2<-rnorm(100,20,5)
#' fsu<-cbind(f1,f2)
#' s1<-rnorm(50,20,5)
#' s2<-rnorm(50,20,5)
#' s3<-rnorm(50,20,5)
#' s4<-rnorm(50,20,5)
#' ssu<-cbind(s1,s2,s3,s4)
#' RCRatio(N=1000, FSU=fsu, SSU=ssu)
#'
#' @references
#' \itemize{
#'\item Islam, S., Chandra, H., Sud, U.C., Basak, P., Ghosh, N. and Ojasvi, P.R. (2021). A Revised Calibration Weight based Ratio Estimator in Two-phase Sampling: A Case when Unit Level Auxiliary Information is Available for the First-phase Sample, Journal of Indian Society of Agricultural Statistics, 75(2), 147–156.
#'\item Ozgul, N. (2021). New improved calibration estimator based on two auxiliary variables in stratified two-phase sampling. Journal of Statistical Computation and Simulation, 91(6), 1243-1256.
#' }

RCRatio<-function(N,FSU,SSU ){
  fsu<-as.data.frame(FSU)
  ssu<-as.data.frame(SSU)
  n1=dim(fsu)[1] # first phase sample size
  n2=dim(ssu)[1] # second phase sample size
  d2i<-n1/n2; di<-N/n2
  d1i=N/n1;#di=n1/n2;d1i;di

  pi_1i=1/d1i  #First phase inclusion probability

  pi_1ij=(n1*(n1-1))/(N*(N-1)) #joined inclusion probability
  pi_2ij=(n2*(n2-1))/(n1*(n1-1))

  pi=1/di #overall inclusion probability


  x1_fsu<-fsu[,1]
  x2_fsu<-fsu[,2]

  y_ssu<-ssu[,1]
  z_ssu<-ssu[,3]

  x1_ssu<-ssu[,2]
  x2_ssu<-ssu[,4]

  #Pearson's correlation coefficient

  Rho_yz<-cor(y_ssu,z_ssu)
  Rho_yx1<-cor(y_ssu,x1_ssu)
  Rho_yx2<-cor(y_ssu,x2_ssu)
  Rho_zx1<-cor(z_ssu,x1_ssu)
  Rho_zx2<-cor(z_ssu,x2_ssu)
  cor<-cbind(yz=Rho_yz,yx1=Rho_yx1,yx2=Rho_yx2,zx1=Rho_zx1,zx2=Rho_zx2)

  t.xs1<-as.matrix(apply(fsu,2,sum))

  t.xs2<-as.matrix(apply(ssu[,c(1,4)],2,sum))

  # Direct ratio estimator
  Est_R<- sum(y_ssu)/sum(z_ssu) #Simple est_R

  A_prime=d1i*t.xs1[1]
  B_prime=d1i*t.xs1[2]

  sum_da=di*t.xs2[1]
  sum_db=di*t.xs2[2]

  a=x1_ssu
  b=x2_ssu

  sum_daa=di*(t(a)%*%a)
  sum_dbb=di*(t(b)%*%b)
  sum_dab=di*(t(a)%*%b)
  #calibration weight

  w_cap1<- rep(0,n2)
  for (i in 1:n2){
    c1=((di*a[i])*(sum_dbb*(A_prime-sum_da)-sum_dab*(B_prime-sum_db)))
    c2=(sum_daa*sum_dbb-sum_dab*sum_dab)
    c3=((di*b[i])*(sum_daa*(B_prime-sum_db)-sum_dab*(A_prime-sum_da)))
    w_cap1[i]=di+c1/c2+c3/c2

  }
  #proposed calibration estimator
  R_cap2<- (t(w_cap1)%*%y_ssu)/(t(w_cap1)%*%z_ssu)

  #Variance of proposed calibration estimator
  sum_daa=di*(t(a)%*%a)
  sum_dbb=di*(t(b)%*%b)
  sum_dab=di*(t(a)%*%b)

  sum_dyy=di*(t(y_ssu)%*%y_ssu)
  sum_dzz=di*(t(z_ssu)%*%z_ssu)
  sum_dya=di*(t(y_ssu)%*%a)
  sum_dyb=di*(t(y_ssu)%*%b)
  sum_dza=di*(t(z_ssu)%*%a)
  sum_dzb=di*(t(z_ssu)%*%b)

  t_cap_z<-di*(t(z_ssu)%*%z_ssu)

  c11<-N*(N-n1)/(n1*t_cap_z^2)
  c22<-N^2*(n1-n2)/(n1*n2*t_cap_z^2)

  R_cap<-Est_R
  u_cap<-y_ssu-R_cap*z_ssu
  q<-sum_daa*sum_dbb-sum_dab^2

  l1<-(sum_dyb*sum_dyy-sum_dya*sum_dab)/q
  l2<-(sum_dya*sum_dbb-sum_dyb*sum_dab)/q
  l3<-(sum_dzb*sum_daa-sum_dza*sum_dab)/q
  l4<-(sum_dza*sum_dbb-sum_dzb*sum_dab)/q

  v_cap=u_cap+b%*%(l3*R_cap-l1)+a%*%(l4*R_cap-l2)

  v_cap_mean=mean(v_cap)
  vec_1=c(rep(1,n2))


  SS_v=(1/(n2-1))*(vec_1%*%(v_cap-vec_1*v_cap_mean)^2)

  y_cap_mean=mean(y_ssu)
  z_cap_mean=mean(z_ssu)

  SS_y=(1/(n2-1))*(vec_1%*%(y_ssu-vec_1*y_cap_mean)^2)
  SS_z=(1/(n2-1))*(vec_1%*%(z_ssu-vec_1*z_cap_mean)^2)
  Rho_yz<-cor(y_ssu,z_ssu)

  SS_u=SS_y+R_cap^2*SS_z-2*R_cap*Rho_yz*sqrt(SS_y)*sqrt(SS_z)

  Var_RC2_SRS<-c11*SS_u+c22*SS_v

  CV=(sqrt(Var_RC2_SRS)/R_cap2)*100
  samplesize<-cbind(FS=n1,SS=n2)
  designweight<-cbind(FS=d1i, SS=d2i, Overall=di)
  incprob<-cbind(pi_1i,pi_1ij,pi_2ij,pi_2ij)
  cor<-cbind(yz=Rho_yz,yx1=Rho_yx1,yx2=Rho_yx2,zx1=Rho_zx1,zx2=Rho_zx2)
  out<-list(CalEstimate=R_cap2,CalVariance=Var_RC2_SRS,CV=CV,SampleSize=samplesize, DesignWeight=designweight, InclusionProb=incprob,Correlation=cor )

  return(out)
}
