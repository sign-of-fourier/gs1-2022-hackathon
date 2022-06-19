library(lubridate)

stores<-read.csv('stores.csv')
transactions<-read.csv('transactions.csv')

data<-merge(transactions, stores, by=c('store_nbr'))

data$day_of_week<-wday(data$date, week_start=1)

data$a_week_ago<-as.Date(data$date)-7
last_week<- data[,c('a_week_ago', 'transactions', 'store_nbr')]
names(last_week)<-c('date', 'last_week_transactions', 'store_nbr')
last_week$date<-as.character(last_week$date)
week_over_week<-merge(data, last_week, by=c('store_nbr', 'date'))

min(unique(data$date))
max(unique(data$date))

train<-week_over_week[week_over_week$date<'2016-01-01',]
test<-week_over_week[week_over_week$date>='2016-01-01',]


dim(test)
dim(train)
dim(data)

model1<-glm(transactions~type+cluster, train[train$day_of_week==1, ], family='gaussian')
p<-predict(model1, test[test$day_of_week==1,])

1-mean((p-test$transactions[test$day_of_week==1])^2)/
  mean((test$transactions[test$day_of_week==1]-mean(test$transactions[test$day_of_week==1]))^2)

# results are very autocorrelated
# look at difference based on using "last week's transactions"

model1<-glm(transactions~type+cluster+last_week_transactions, train[train$day_of_week==1, ], family='gaussian')
p<-predict(model1, test[test$day_of_week==1,])

1-mean((p-test$transactions[test$day_of_week==1])^2)/
  mean((test$transactions[test$day_of_week==1]-mean(test$transactions[test$day_of_week==1]))^2)

for(a in 1:7){
    model<-glm(transactions~type+cluster+last_week_transactions, train[train$day_of_week==a, ], family='gaussian')
    p<-predict(model, test[test$day_of_week==a,])
    print(paste('r squared for model based on day of week', a, 
            1-mean((p-test$transactions[test$day_of_week==a])^2)/
                mean((test$transactions[test$day_of_week==a]-mean(test$transactions[test$day_of_week==a]))^2)
        )
    )
}


