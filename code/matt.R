#Nested logit coding assignment
source('header.R')

#Load data
val<-readRDS(paste0(varSave, 'choice_values.rds'))
lambda<-readRDS(paste0(varSave, 'lambda_values.rds'))

# Simulate 10,000 choices corresponding to epsilon and eta values for each choice/bucket: 
N <- 10^5
nBuckets <- nrow(lambda)
nChoices <- length(unique(val$choice))
set.seed(1023)
epsilon <- rgev(n = nBuckets*N, xi = 0, mu = 0, beta = 1)
eta <- rgev(n = nChoices*N, xi = 0, mu = 0, beta = 1)
k_i<- array(val[,1:2])

#################################################################################################
# 1 under the traditional specification where epsilon_k resolves before eta_i
# We estimate the inclusive value of each bucket to choose a bucket for each simulation.
# Then we estimate the value of each choice given the bucket selected for each simulation,
# and select the best choice given the bucket.
#################################################################################################
IV <- left_join(val, lambda, by = "bucket") %>%
      group_by(bucket) %>%
      summarise(lambda = unique(lambda), iv = log(sum(exp(val/lambda)))) %>%
      mutate(iv = lambda * iv)

nest_choice <- data.frame(
  "bucket" = apply(IV$iv + 
                    matrix(epsilon,nrow = nBuckets,byrow = TRUE),
                   2,which.max)
)

item_choice1 <- data.frame(
  "bucket" = nest_choice$bucket,
  "choice" = LETTERS[
    apply(
      matrix(
        left_join(nest_choice,val,by="bucket")$val,
        nrow = nChoices,
        ncol = N) + 
        matrix(eta,nrow = nChoices,byrow = TRUE),
      2,which.max)
  ]
)

#################################################################################################
# 2 under the single-shot specification where epsilon_k and eta_i both resolve up front 
# We directly estimate the value of each bucket-choice combination for each simulation, 
# and select the maximum one.
#################################################################################################
a <- matrix(eta,nrow = nChoices,byrow = TRUE)
U2 <- apply(val$val + 
              matrix(rep(epsilon,nChoices),nrow = nBuckets * nChoices, byrow = TRUE) + 
              rep(lambda$lambda,nChoices) * 
              a[rep(seq(nChoices),each = nBuckets),],
            2,which.max
            )

#Alternative U2 calc. I think it's a little easier to read? What do you think? The results are identical.
#eps <- matrix(epsilon,nrow = nChoices,byrow = TRUE)
#U3<-apply(mdply(k_i,
#                function(k, i) 
#                  val$val[val$bucket==k & val$choice==i] + 
#                  eps[k,] + 
#                  lambda$lambda[k]*a[i,]
#                )%>% 
#                select(c(-1,-2)),
#          2,which.max)
      


item_choice2 <- data.frame(
  "bucket" = val$bucket[U2],
  "choice" = val$choice[U2]
)

#################################################################################################
# Estimate the probability of selecting each bucket-choice combination under each of the methods above.
# Also estimate the covariance matrix of estimates with bootstrap.
# Use 100 resamplings of size 100 each for the bootstrap.
#################################################################################################
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
  covProb <- cov(t(r.prob))
  rownames(covProb)<- c("A1","A2","A3","B1","B2","B3","C1","C2","C3")
  colnames(covProb)<- c("A1","A2","A3","B1","B2","B3","C1","C2","C3")
  return(covProb)
} 

sigma1 <- b.cov(item_choice1,100,100)
sigma2 <- b.cov(item_choice2,100,100)

#################################################################################################
# Compute the true, theoretical probability. 
# We first calculate the probability of selecting each bucket.
# We then calculate the probabilty of selecting each choice within the bucket.
# The final probability of selecting each bucket-choice combination is the above two probabilities
# multiplied together.
#################################################################################################
Prob_bucket <- function (bucket){
  Pbucket = exp(IV$iv[bucket])/sum(exp(IV$iv))
  return (data.frame("bucket" = 1:3,Pbucket))
}


Prob_choice_bucket <- function (choice, bucket){
  Pchoice <- exp(val$val[val$choice == choice & val$bucket == bucket])/sum(exp(val$val[val$bucket==bucket]))
  return (Pchoice)
}

trueProb <- left_join(
  ddply(val,.(bucket,choice),summarise,PchoiceBucket = Prob_choice_bucket(choice,bucket)),
  Prob_bucket(1:3),
  by = "bucket"
) %>% 
  mutate(prob = PchoiceBucket * Pbucket) %>%
  select(bucket, choice, prob) %>%
  arrange(choice, bucket)

#################################################################################################
# We calculate the Wald statistic for our estimates using each of the two methods above.
# The hypothesis we are testing is that our estimates are equal to the true, theoretical probs.
# We expect to not reject the hypothesis in the first method, and to reject the hypothesis
# in the second case. 
# Our expectations are confirmed.
#################################################################################################
wald <- function (theta, sigma, H0){
  W<- t(theta - H0)%*%
    ginv(sigma) %*%
    (theta - H0)
  return(W)
}

W1 <- wald(Prob1$prob, sigma1, trueProb$prob)
1 - pchisq(W1,9) # Fail to reject the null hypothesis of equal probabilities

W2 <- wald(Prob2$prob,sigma2,trueProb$prob)
1 - pchisq(W2,9) # Reject the null hypothesis of equal probabilities

