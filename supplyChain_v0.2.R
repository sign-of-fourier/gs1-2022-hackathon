library(scales)


# Markov Chain Simulation of Supply Chain
# Ordering to Transportation to Consumption back to Ordering (next time step)


# input parameters
n_products<-25
n_time_steps<-100
n_orders <- 10
n_sources <- 4
n_transports <- 7
n_observations <- 100
n_stores<- 10


# time step
tm<-.0001


# variance and covariance of ordering factors

order_distribution<-runif(n_orders)
order_distribution<-order_distribution/sum(order_distribution)

src_distribution<-runif(n_sources)
src_distribution<-src_distribution/sum(src_distribution)

ORD_TO_SRC<-matrix(runif(n_sources*n_orders), nrow=n_sources, ncol=n_orders)
SRC_TO_TNS<-matrix(runif(n_sources*n_transports), nrow=n_transports, ncol=n_sources)
TNS_TO_CNS<-matrix(runif(n_orders*n_transports), nrow=n_stores, ncol=n_transports)
CNS_TO_ORD<-matrix(runif(n_orders*n_transports), nrow=n_stores, ncol=n_stores)

solve(order_distribution, src_distribution)

SRC_VCV <-matrix(runif(n_sources^2), ncol=n_sources, nrow=n_sources)
SRC_VCV<-SRC_VCV %*% t(SRC_VCV)*sqrt(tm)
diag(SRC_VCV)<-diag(SRC_VCV)+.1
eigen(SRC_VCV)$values

TNS_VCV <-matrix(runif(n_transports^2), nrow=n_transports, ncol=n_transports)
TNS_VCV<-TNS_VCV %*% t(TNS_VCV)*sqrt(tm)
eigen(TNS_VCV)$values
TRN_MIN<-c(.2, -.4)
TRN_MAX<-c(.1, .7)
TRN_MEAN<-c(.1, .4)


CNS_VCV <-matrix(runif(n_orders^2), nrow=n_stores, ncol=n_stores)
CNS_VCV<-CNS_VCV %*% t(CNS_VCV)*sqrt(tm)
diag(CNS_VCV)<-diag(CNS_VCV)+.11
eigen(CNS_VCV)$values

Inventory_Target_Base<-matrix(runif(n_stores*n_products), nrow=n_products, ncol=n_stores) 

Inventory_Target<-matrix(runif(n_stores*n_products), nrow=n_products, ncol=n_stores) 
n_observations<-100


# How ordering is distributed accross the producers

ORD_TO_SRC<-matrix(runif(n_stores*n_sources), nrow=n_sources, ncol=n_stores)
for(a in 1:n_stores){
    ORD_TO_SRC[,a]<-ORD_TO_SRC[,a]/sum(ORD_TO_SRC[,a])
}


# This matrix defines how the production is expected to be 
# distributed accross the transportation channels

SRC_TO_TNS<-matrix(runif(n_sources*n_transports), nrow=n_transports, ncol=n_sources)
for(a in 1:n_sources){
    SRC_TO_TNS[,a]<-SRC_TO_TNS[,a]/sum(SRC_TO_TNS[,a])
}


# estimated demand would be inserted here
# for this demonstration use random numbers

expected_consumption<-matrix(runif(n_stores*n_products)*10, nrow=n_products, ncol=n_stores)


SupplyChain<-function(Inventory_Target, lab, n_observations){

    # initialize Monte Carlo matricies

    MC<-array(0, dim=c(n_observations, n_time_steps, n_products, n_stores))
    MC_Understock<-array(0, dim=c(n_observations, n_time_steps, n_products, n_stores))
    MC_SLS<-array(0, dim=c(n_observations, n_time_steps, n_products, n_stores))
    MC_TSP<-array(0, dim=c(n_observations, n_time_steps, n_products, n_transports))
    MC_SRC<-array(0, dim=c(n_observations, n_time_steps, n_products, n_sources))

    for(a in 1:n_observations){


        Inventory_after_sales<-array(0,dim=c(n_time_steps, nrow=n_products, ncol=n_stores))
        Understock<-array(0,dim=c(n_time_steps, nrow=n_products, ncol=n_stores))
        Sales<-array(0,dim=c(n_time_steps, nrow=n_products, ncol=n_stores))
        Trnspt_Cost<-array(0,dim=c(n_time_steps, nrow=n_products, ncol=n_transports))
        Mfctr_Cost<-array(0,dim=c(n_time_steps, nrow=n_products, ncol=n_sources))

        # Initilize Current Inventory
        # Then calculate orders (deterministic)

        Inventory_after_sales[1,,]<-matrix(rnorm(n_orders), nrow=n_products, ncol=n_stores)+Inventory_Target
        ORD<-matrix(sapply(Inventory_Target-Inventory_after_sales[1,,], max, 0)+expected_consumption, nrow=n_products, ncol=n_stores)


        for(b in 2:n_time_steps){

            # Sources; How much is produced based on how much is ordered
            # with some unknown variation due to possible production problems

            SRC<- t(ORD_TO_SRC%*%t(ORD)) +  matrix(rnorm(n_sources), nrow=n_products, ncol=n_sources) %*% chol(SRC_VCV)

            Mfctr_Cost[b-1,,]<- SRC

            # Transportation (conditional on source; markov chain or conditional distribution)
            # How much is transported also depends on various unknown delays/failures

            TNS<- t(SRC_TO_TNS%*%t(SRC)) +   matrix(rnorm(n_transports*n_products), nrow=n_products, ncol=n_transports) %*% chol(TNS_VCV)
            TNS<-matrix(sapply(TNS, max, 0), nrow=n_products, ncol=n_transports)

            Trnspt_Cost[b-1,,]<- TNS

            RSTK<-ORD
            for(d in 1:n_products){
                srd<-sum(ORD[d,])
                if(srd > 0){
                    RSTK[d,]<-ORD[d,]/srd * sum(TNS[d,])
                }
            }


            # Invetory after New product arrives        
            # Existing inventory + distribution: t(TNS_TO_CNS%*%t(TNS)) + random consumption

            Inventory_aft_ship_bf_sales<-Inventory_after_sales[b-1,,] + RSTK
            CNS<- Inventory_aft_ship_bf_sales - matrix(runif(n_products*n_stores), nrow=n_products, ncol=n_stores) %*% chol(CNS_VCV) - expected_consumption
            Understock[b,,]<-CNS * (CNS < 0)

            CNS<-matrix(sapply(CNS, max, 0), nrow=n_products, ncol=n_stores)

            Sales[b,,]<-Inventory_aft_ship_bf_sales-CNS

            Inventory_after_sales[b,,]<-CNS

            # Reordering (Invetory Target is a paramter to optimize or set given a risk level)

            ORD<-matrix(sapply(Inventory_Target-CNS+expected_consumption, max, 0), nrow=n_products, ncol=n_stores)
        }
        MC[a,,,]<-Inventory_after_sales
        MC_SRC[a,,,]<-Mfctr_Cost
        MC_TSP[a,,,]<-Trnspt_Cost
        MC_SLS[a,,,]<-Sales
        MC_Understock[a,,,]<-Understock
    }    
    

    return(c(
                # Total Profit Sales - Transported goods * .1 - Inventory Produced *.2 - Invetory Stored *.1
                # Real equation could easily be more complicated with loss of generality

                (sum(MC_SLS)-sum(MC_TSP)*.1-sum(MC_SRC)*.2-sum(MC)*.1)/n_observations,

                # unmet demand
                sum(MC_Understock)/n_observations
            )
    )
}

# look at profit/understock tradeoff for various inventory targets

X<-seq(.25, 3, .125)
earnings<-1:length(X)

E<-length(X)
F<-length(X)
for(a in 1:length(X)){
     print(paste(X[a], 'x Base Inventory Level', sep=''))

     S<-SupplyChain(Inventory_Target_Base*X[a], paste('Inventory Factor', X[a]), 12)
     E[a]<-S[1]
     F[a]<-S[2]
     print(paste("Profit:", dollar(S[1])))
     print(paste("Unmet Demand",dollar(S[2])))
}


par(mfrow=c(3,1))
par(mar = c(5, 4, 4, 4) + 0.3)            
plot(X, E, main="Total Profit vs. Shortfall", type='l', col='blue', 
ylab='Profit', xlab='Inventory Target Factor')
par(new = TRUE)

plot(X, -F, col='red', xlab='', ylab='', axes=FALSE, lty=3, type='l')
axis(side = 4, at = pretty(range(-F)))

mtext("Understock", side = 4, line = 2)    



legend('topright', legend=c("Profit", "Understock"),
       col=c("blue", "red"), lty=c(1,3), cex=0.8)



plot(Inventory_after_sales[,1,1], ylim=c(0, 1), xlab='Time', ylab='Inventory', type='l', main='Inventory Product 1 Over Time')

for(a in 2:5){
    lines(MC[a,,1,1])
}

lines(    apply(MC[,,1,1], 2, mean), col='red', lwd=2.0, lty=1 )

legend("top", legend=c("Expected Level", "Random Path"),
       col=c("red", "black"), lty=c(1,1), cex=0.8)



# Average Inventory Level
sum(MC[,n_time_steps,1,])/n_observations

# Amount of Understock; Product 1
sum(MC_Understock[,,1,])/n_observations

# Total Amount of Understock
sum(MC_Understock)/n_observations


hist( c(MC_Understock)[c(MC_Understock < -0.01)], main='Distribution of Understock',
xlab='Amount of Understock at Single Point in Time' )

# Total Profit for Product 1
sum(MC_SLS[,,1,])-sum(MC_TSP[,,1,])*.1-sum(MC_SRC[,,1,])*.2-sum(MC[,,1,])*.1



