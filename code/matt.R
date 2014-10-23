#Nested logit coding assignment
source('header.R')

#Load data
val<-readRDS(paste0(varSave, 'choice_values.rds'))
lambda<-readRDS(paste0(varSave, 'lambda_values.rds'))

# Simulate 10,000 choices corresponding to V_ik and lambda_k values: 
N <- 10^5
epsilon <- rgev(n = nrow(lambda)*N, xi = 0, mu = 0, beta = 1)
eta <- rgev(n = length(unique(val$choice))*N, xi = 0, mu = 0, beta = 1)

# 1 under the traditional specification where epsilon_k resolves before eta_i
IV <- left_join(val, lambda, by = "bucket") %>%
      group_by(bucket) %>%
      summarise(lambda = unique(lambda), iv = log(sum(exp(val/lambda)))) %>%
      mutate(iv = lambda * iv)
nest_choice <- data.frame(
  "bucket" = apply(IV$iv + matrix(epsilon,nrow = 3,byrow = TRUE),2,which.max)
)

# nest_choice <- ldply(seq(N), function(x){
#   util = IV$iv + rgev(n = nrow(IV), xi = 0, mu = 0, beta = 1)
#   which.max(util)
# }) 

item_choice1 <- data.frame(
  "bucket" = nest_choice$bucket,
  "choice" = LETTERS[
    apply(
      matrix(
        left_join(nest_choice,val,by="bucket")$val,
        nrow = length(unique(val$choice)),
        ncol = N) + 
        matrix(eta,nrow = 3,byrow = TRUE),
      2,which.max)
    ]
)

# 2 under the single-shot specification where epsilon_k and eta_i both resolve up front  
a <- matrix(eta,nrow = 3,byrow = TRUE)
U2 <- apply(
  temp<-
  val$val + 
  matrix(rep(epsilon,3),nrow = 9, byrow = TRUE) + 
  rep(lambda$lambda,3) * 
  a[rep(seq(nrow(a)),each = 3),],
  2,
  which.max
)
item_choice2 <- data.frame(
  "bucket" = val$bucket[U2],
  "choice" = val$choice[U2]
)

# Estimate choice probabilities and calculate covariance matrix using bootstrap
# 1 under the traditional specification where epsilon_k resolves before eta_i
# 2 under the single-shot specification where epsilon_k and eta_i both resolve up front

k_i<- array(val[,1:2])
Prob <- function(k_i, data){
  P <- mdply(k_i, 
        function(k,i) with(data,sum(bucket == k & choice == i)/nrow(data))
    )
  colnames(P)<- c("bucket", "choice", "prob")
  return(P)
}
Prob1 <- Prob(k_i,item_choice1)
Prob2 <- Prob(k_i,item_choice2)

# Bootstrap covariance matrix
b.cov <- function(data, num, size){
  resamples <- lapply(1:num, function(i) data[sample(1:nrow(data),size,replace = TRUE),])
  r.prob <- sapply(resamples, function(x) Prob(k_i,x)[,3])
  covProb <- data.frame(cov(t(r.prob)),row.names = c("A1","B1","C1","A2","B2","C2","A3","B3","C3"))
  names(covProb) <- row.names(covProb)
  return(covProb)
} 

sigma1 <- b.cov(item_choice1,100,100) %>%
      as.matrix
sigma2 <- b.cov(item_choice2,100,100) %>%
      as.matrix



# True probability
Prob_bucket <- function (bucket){
  Pbucket = exp(IV$iv[bucket])/sum(exp(IV$iv))
  return (data.frame("bucket" = 1:3,Pbucket))
}


Prob_choice_bucket <- function (choice,bucket){
  Pchoice <- exp(val$val[val$choice == choice & val$bucket == bucket])/sum(exp(val$val[val$bucket==bucket]))
  return (Pchoice)
}

trueProb <- left_join(
  ddply(val,
        .(bucket,choice),
        summarise,
        PchoiceBucket = Prob_choice_bucket(choice,bucket)),
  Prob_bucket(1:3),
  by = "bucket"
) %>% 
  mutate(prob = PchoiceBucket * Pbucket) %>%
  select(bucket,choice,prob)

wald <- function (theta, sigma, H0,rem,n = N){
  W<- n * t(theta[-rem] - H0[-rem])%*%
    (sigma[-rem,-rem]) %*%
    (theta[-rem] - H0[-rem])
  return(W)
}


(Truth <- arrange(trueProb,choice,bucket))
W1 <- wald(Prob1$prob, sigma1, Truth$prob,rem=-1,n=N)
W1 > qchisq(0.95,8) # Fail to reject the null hypothesis of equal probabilities

W2 <- wald(Prob2$prob,sigma2,Truth$prob,rem=c(1,2,4),n=N)
W2 > qchisq(0.95,6) # Fail to reject the null hypothesis of equal probabilities

