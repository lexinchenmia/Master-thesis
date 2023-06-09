#Research Master Thesis
#title:Social Contagion in Homophilic Value-Based Social Networks: 
#Two-Step Simulations of Network Generation and Optimal Seeding Applied in the Threshold Model
#Lexin Chen

library(igraph)
library(plotrix)

#weight distribution under different homophily strength
fun0 <- function(x) exp(0*x)
fun1 <- function(x) exp(-x)
fun2 <- function(x) exp(-2*x)
fun3 <- function(x) exp(-3*x)

plot(fun1, 0, 10, col="red", lwd = 2, xlab="value difference",ylab="weight/raw probability",main="Weight Distribution under Different Homophily Strength")
plot(fun2, 0, 10, add=TRUE, col="green",lwd = 2)
plot(fun3, 0, 10, add=TRUE, col="gray",lwd = 2)
plot(fun0, 0, 10, add=TRUE, col="blue",lwd = 2)
legend(6, 0.8, legend=c("ρ = 0", "ρ = 1","ρ = 2","ρ = 3"),
       fill=c("blue","red","green","gray"), cex=0.8)


#beta distribution plots

par(mfrow=c(3,1))

x= rbeta(10000,5,1)
hist(x, main="A Conservative Society", freq=FALSE)
lines(density(x), col='red', lwd=3)
abline(v = c(mean(x),median(x)),  col=c("green", "red"), lty=c(2,2), lwd=c(3, 3))


x= rbeta(10000,5,5)
hist(x, main="Symmetrical", freq=FALSE)
lines(density(x), col='red', lwd=3)
abline(v = c(mean(x),median(x)),  col=c("green", "red"), lty=c(2,2), lwd=c(3, 3))


x= rbeta(10000,0.1,5)
hist(x, main="Positive or Right Skewness", freq=FALSE)
lines(density(x), col='red', lwd=3)
abline(v = c(mean(x),median(x)),  col=c("green", "red"), lty=c(2,2), lwd=c(3, 3))


#functions
samps <- function(i, weights, items) {
  if (sum(weights[, i] > 0) > 4) {
    return(sample(items, prob = weights[, i] / sum(weights[, i]), size = 5, replace = FALSE))
  } else {
    return(c(i - 5, i - 4, i - 3, i - 2, i - 1))
  }
}


homophilic_linkage <- function(Distri, Rho, population) {
  weights <- matrix(0, nrow = population, ncol = population)
  for (i in 1:length(Distri)) {
    for (j in 1:length(Distri)) {
      weights[i, j] <- round(exp(-Rho * abs(Distri[i] - Distri[j])),2)
    }
  }
  diag(weights) <- 0
  items <- 1:population
  linklist <- list()
  for (i in 1:length(Distri)) {
    sampsnow <- samps(i, weights, items)
    for (j in 1:length(sampsnow)) {
      pairnow <- sort(c(i, sampsnow[j]))
      linklist <- c(linklist, list(pairnow))
    }
  }
  linklist <- unique(do.call(rbind, linklist))
  linksPerNodeList <- list()
  for (i in 1:population) {
    mylinks <- c(i)
    for (j in 1:nrow(linklist)) {
      if (i %in% linklist[j, ]) {
        if (linklist[j, 1] == i) {
          mylinks <- c(mylinks, linklist[j, 2])
        } else {
          mylinks <- c(mylinks, linklist[j, 1])
        }
      }
    }
    linksPerNodeList <- c(linksPerNodeList, list(mylinks))
  }
  return(linksPerNodeList)
}

normalize <- function(original){
  normalized_value = ((original - 0) / (10 - 0)) * 0.8
  return(normalized_value)
}

trans_edgelist <- function(test){
  test_df<- data.frame(from = 1, to = 1)
  for (i in 1:length(test)){
    for (j in test[[i]]){
      edge <- data.frame(from = i, to = j)
      test_df <- rbind(test_df, edge)
    }
  }
  test_df <- test_df[which(test_df$from != test_df$to),]
  
  for (i in 1:nrow(test_df)){
    if (test_df[i,1] > test_df[i,2]){
      test_df[i,] <- c(test_df[i,2],test_df[i,1])
      
    }
  }
  
  test_df <- unique(test_df)
  
  return(test_df)
}

#task 1
#fix d, alter the homophily strength from 0 to 3 with break as 0.1
#see mean density and transitivity changes


homo_c <- data.frame(homophily = NA, density = NA, transitivity = NA)
for (i in seq(0, 3, 0.3)){
  
  temp_den <- c()
  temp_tran <- c()
  for (j in 1:100){
    conservative_value = round(rbeta(1000,5,5)*10,2)
    
    test <- homophilic_linkage(Distri = conservative_value, Rho = i, population = 1000)
    
    test_df <- trans_edgelist(test)
    
    g <- graph(t(test_df),directed=F)
    g<-simplify(g)
    temp_den <- c(temp_den, graph.density(g))
    temp_tran <- c(temp_tran, transitivity(g))
  }
  temp_df <- data.frame(homophily = i, density = mean(temp_den), transitivity = mean(temp_tran))
  homo_c <- rbind(homo_c,temp_df)
}
homo_c <- homo_c[-1,]

plot(homo_c$homophily,homo_c$density,type="l",col="red",ylim=c(0.0089,0.0355),
     xlab="homophily strength",ylab="",main="Conservative")
lines(homo_c$homophily,homo_c$transitivity,col="green")
legend("topleft", legend=c("density", "transitivity"),
       col=c("red", "green"), lty=1:2, cex=0.8)


homo_c <- data.frame(homophily = NA, density = NA, transitivity = NA)
for (i in seq(0, 3, 0.3)){
  
  temp_den <- c()
  temp_tran <- c()
  for (j in 1:100){
    conservative_value = round(rbeta(1000,5,2)*10,2)
    
    test <- homophilic_linkage(Distri = conservative_value, Rho = i, population = 1000)
    
    test_df <- trans_edgelist(test)
    
    g <- graph(t(test_df),directed=F)
    g<-simplify(g)
    temp_den <- c(temp_den, graph.density(g))
    temp_tran <- c(temp_tran, transitivity(g))
  }
  temp_df <- data.frame(homophily = i, density = mean(temp_den), transitivity = mean(temp_tran))
  homo_c <- rbind(homo_c,temp_df)
}
homo_c <- homo_c[-1,]

homo_o <- data.frame(homophily = NA, density = NA, transitivity = NA)
for (i in seq(0, 3, 0.3)){
  
  temp_den <- c()
  temp_tran <- c()
  for (j in 1:100){
    open_value = round(rbeta(1000,2,5)*10,2)
    
    test <- homophilic_linkage(Distri = open_value, Rho = i, population = 1000)
    
    test_df <- trans_edgelist(test)
    
    g <- graph(t(test_df),directed=F)
    g<-simplify(g)
    temp_den <- c(temp_den, graph.density(g))
    temp_tran <- c(temp_tran, transitivity(g))
  }
  temp_df <- data.frame(homophily = i, density = mean(temp_den), transitivity = mean(temp_tran))
  homo_o <- rbind(homo_o,temp_df)
}
homo_o <- homo_o[-1,]

par(mfrow = c(1,3))

plot(hc$homophily,hc$density,type="l",col="red",ylim=c(0.0089,0.0355),
     xlab="homophily strength",ylab="",main="Conservative")
lines(hc$homophily,hc$transitivity,col="green")
legend("topleft", legend=c("density", "transitivity"),
       col=c("red", "green"), lty=1:2, cex=0.8)

plot(hm$homophily,hm$density,type="l",col="red",ylim=c(0.0089,0.0355),
     xlab="homophily strength",ylab="",main="Moderate")
lines(hm$homophily,hm$transitivity,col="green")
legend("topleft", legend=c("density", "transitivity"),
       col=c("red", "green"), lty=1:2, cex=0.8)

plot(ho$homophily,ho$density,type="l",col="red",ylim=c(0.0089,0.0355),
     xlab="homophily strength",ylab="",main="Progressive")
lines(ho$homophily,ho$transitivity,col="green")
legend("topleft", legend=c("density", "transitivity"),
       col=c("red", "green"), lty=1:2, cex=0.8)

summary(lm(homo_c$density~homo_c$homophily))

summary(lm(hc$density~hc$homophily))


#task 2
#fix d, homophily strength, change social openness
#choose one network each from three types of lists
#compare the mean social openness and network density and transitivity

#open
so_o <- data.frame(so = NA, density = NA, transitivity = NA)
for (i in seq(1, 5, 0.2)){
  
  temp_den <- c()
  temp_tran <- c()
  for (j in 1:100){
    open_value = round(rbeta(1000,i,5)*10,2)
    
    test <- homophilic_linkage(Distri = open_value, Rho = 1.5, population = 1000)
    
    test_df <- trans_edgelist(test)
    
    g <- graph(t(test_df),directed=F)
    g<-simplify(g)
    temp_den <- c(temp_den, graph.density(g))
    temp_tran <- c(temp_tran, transitivity(g))
  }
  temp_df <- data.frame(so = i, density = mean(temp_den), transitivity = mean(temp_tran))
  so_o <- rbind(so_o,temp_df)
  print(i)
}
so_o <- so_o[-1,]

#ylim is free to change
plot(so_o$so,so_o$density,type="l",col="red",ylim=c(0.0099,0.02),
     xlab="mean social openness(reverse)",ylab="",main="Progressive")
lines(so_o$so,so_o$transitivity,col="green")
legend("left", legend=c("density", "transitivity"),
       col=c("red", "green"), lty=1:2, cex=0.8)


#conservative
so_c <- data.frame(so = NA, density = NA, transitivity = NA)
for (i in seq(1, 5, 0.2)){
  
  temp_den <- c()
  temp_tran <- c()
  for (j in 1:100){
    open_value = round(rbeta(1000,5,i)*10,2)
    
    test <- homophilic_linkage(Distri = open_value, Rho = 1.5, population = 1000)
    
    test_df <- trans_edgelist(test)
    
    g <- graph(t(test_df),directed=F)
    g<-simplify(g)
    temp_den <- c(temp_den, graph.density(g))
    temp_tran <- c(temp_tran, transitivity(g))
  }
  temp_df <- data.frame(so = i, density = mean(temp_den), transitivity = mean(temp_tran))
  so_c <- rbind(so_c,temp_df)
  print(i)
}
so_c <- so_c[-1,]

#ylim is free to change
plot(so_c$so,so_c$density,type="l",col="red",
     xlab="mean social openness",ylab="",main="Conservative",ylim=c(0.009938,0.01843))
lines(so_c$so,so_c$transitivity,col="green")
legend("left", legend=c("density", "transitivity"),
       col=c("red", "green"), lty=1:2, cex=0.8)

lm1 <- lm(so_c$transitivity~so_c$so)
summary(lm1)

lm2 <- lm(so_c$density~so_c$so)
summary(lm2)

lm3 <- lm(so_o$transitivity~log(so_o$so))
summary(lm3)

lm4 <- lm(so_o$density~so_o$so)
summary(lm4)

#zoom-in plot
twoord.plot(so_c$so,so_c$density,so_c$so,so_c$transitivity,
            xlab="Mean Social Openness",ylab = "Density",rylab="Transitivity",
            main="Conservative")




#task 3
#seeding and diffusion
#choose 3 networks from three social types
#remember to give process plots

ThModel<-function(node_seed,network,threshold){ 
  #prepare input for the 'calculate_value' function#
  adj_matrix <- igraph::as_adjacency_matrix(network, type = "both")
  en<-data.frame(row = NA, col = NA)
  for (i in 1:vcount(network)){
    for (j in 1:vcount(network)){
      if (adj_matrix[i,j]>0){
        temp <- data.frame(row = i, col = j)
        en <- rbind(en, temp)
      }
    }
  }
  en <- en[-1,]
  each_neighbors <- split(en[, 2], en[, 1]) #get the neigbhour list of each node
  
  nNode<-vcount(network)
  node_status <- rep.int(0, nNode) 
  neighbour_status<-rep.int(0, nNode)  ##percentage of adopted neighbours
  new_infected <- list()
  day_total_infected <- c(0) ### Total number of active people by end of each day
  
  
  ### Day 1 ####
  day <- 1
  node_status[as.numeric((node_seed))] <- 1 
  new_infected[[day]] <-node_seed
  day_total_infected[day]=sum(node_status == 1)
  
  day <- day+1
  
  NotAdopted <- which(node_status == 0)
  Adopted <- which(node_status == 1)
  
  for (i in 1:nNode){
    if (i %in% NotAdopted){
      neighbor <- each_neighbors[[i]]
      neighbour_status[i] <- sum(node_status[neighbor])/length(neighbor)
    }
  }
  
  new_infected[[day]] <- setdiff(which(neighbour_status > threshold), Adopted)
  node_status[new_infected[[day]]] <- 1  #update the staus to 1 for those newly adopted
  day_total_infected[day] <- sum(node_status)
  
  
  
  ####
  
  while (day_total_infected[day] != day_total_infected[day-1]){
    day <- day+1
    
    NotAdopted <- which(node_status == 0)
    Adopted <- which(node_status == 1)
    
    for (i in 1:nNode){
      if (i %in% NotAdopted){
        neighbor <- each_neighbors[[i]]
        neighbour_status[i] <- sum(node_status[neighbor])/length(neighbor)
      }
    }
    
    new_infected[[day]] <- setdiff(which(neighbour_status > threshold), Adopted)
    node_status[new_infected[[day]]] <- 1  #update the staus to 1 for those newly adopted
    day_total_infected[day] <- sum(node_status)
    
    
  }
  #return(day_total_infected)
  return(list(day_total_infected,new_infected))
}


#3.0 seeding strategy
#running diffusion
record <- data.frame(run = NA,
                     cov = NA, 
                     time = NA,
                     choice = NA)
net_list <- list()
for (i in 1:100){
  m_value = round(rbeta(1000, 4.5 , 5)*10,2)#the second and third parameters are altered to collect data of different social types
  
  test <- homophilic_linkage(Distri = m_value, Rho = 1, population = 1000)
  
  test_df <- trans_edgelist(test)    
  g <- graph(t(test_df),directed=F)
  g<-simplify(g)
  
  m_threshold = normalize(m_value)
  
  V(g)$threshold <- m_threshold
  
  deg <- as.vector(degree(g))
  
  # define the new range
  new_min <- 1
  new_max <- 1000
  
  # normalize the range (1, 17) to (1, 1000)
  normalized_deg <- ((deg - min(deg)) / (max(deg) - min(deg))) * (new_max - new_min) + new_min
  
  
  eigen <- unlist(eigen_centrality(g)[[1]])
  cen_m <- data.frame(deg = deg,
                      deg_r = normalized_deg,
                      clo = closeness(g),
                      clo_r = rank(closeness(g)),
                      bet = betweenness(g),
                      bet_r = rank(betweenness(g)),
                      eigen = eigen,
                      eig_r = rank(eigen))
  cen_m$sum_r <- cen_m$deg_r + cen_m$clo_r + cen_m$bet_r + cen_m$eig_r
  cen_m$dbm <- cen_m$deg_r + cen_m$bet_r + cen_m$eig_r
  cen_m$db <- cen_m$deg_r + cen_m$bet_r
  cen_m$de <- cen_m$deg_r + cen_m$eig_r
  
  cen_m$value <- m_value
  cen_m$name <- 1:1000
  
  
  #seeds
  cen_sorted <- cen_m[order(-cen_m$sum_r),]
  seeds_sum <- cen_sorted$name[1:20]
  
  cen_sorted <- cen_m[order(-cen_m$deg_r),]
  seeds_deg <- cen_sorted$name[1:20]
  
  cen_sorted <- cen_m[order(-cen_m$clo_r),]
  seeds_clo <- cen_sorted$name[1:20]
  
  cen_sorted <- cen_m[order(-cen_m$bet_r),]
  seeds_bet <- cen_sorted$name[1:20]
  
  cen_sorted <- cen_m[order(-cen_m$eig_r),]
  seeds_eig <- cen_sorted$name[1:20]
  
  cen_sorted <- cen_m[order(-cen_m$dbm),]
  seeds_dbm <- cen_sorted$name[1:20]
  
  cen_sorted <- cen_m[order(-cen_m$db),]
  seeds_db <- cen_sorted$name[1:20]
  
  cen_sorted <- cen_m[order(-cen_m$de),]
  seeds_de <- cen_sorted$name[1:20]
  
  #model
  #1 = sum, 2=deg, 3=clo, 4=bet, 5=eig, 6=dbm, 7=db, 8=de
  m_sum <- ThModel(seeds_sum, g, m_threshold)
  temp <- data.frame(run = i,
                     cov = max(m_sum[[1]]),
                     time = length(m_sum[[2]])-1,
                     choice = 1)
  record <- rbind(record, temp)
  
  m_deg <- ThModel(seeds_deg, g, m_threshold)
  temp <- data.frame(run = i,
                     cov = max(m_deg[[1]]),
                     time = length(m_deg[[2]])-1,
                     choice = 2)
  record <- rbind(record, temp)
  
  m_clo <- ThModel(seeds_clo, g, m_threshold)
  temp <- data.frame(run = i,
                     cov = max(m_clo[[1]]),
                     time = length(m_clo[[2]])-1,
                     choice = 3)
  record <- rbind(record, temp)
  
  m_bet <- ThModel(seeds_bet, g, m_threshold)
  temp <- data.frame(run = i,
                     cov = max(m_bet[[1]]),
                     time = length(m_bet[[2]])-1,
                     choice = 4)
  record <- rbind(record, temp)
  
  m_eig <- ThModel(seeds_eig, g, m_threshold)
  temp <- data.frame(run = i,
                     cov = max(m_eig[[1]]),
                     time = length(m_eig[[2]])-1,
                     choice = 5)
  record <- rbind(record, temp)
  
  m_dbm <- ThModel(seeds_dbm, g, m_threshold)
  temp <- data.frame(run = i,
                     cov = max(m_dbm[[1]]),
                     time = length(m_dbm[[2]])-1,
                     choice = 6)
  record <- rbind(record, temp)
  
  m_db <- ThModel(seeds_db, g, m_threshold)
  temp <- data.frame(run = i,
                     cov = max(m_db[[1]]),
                     time = length(m_db[[2]])-1,
                     choice = 7)
  record <- rbind(record, temp)
  
  m_de <- ThModel(seeds_de, g, m_threshold)
  temp <- data.frame(run = i,
                     cov = max(m_de[[1]]),
                     time = length(m_de[[2]])-1,
                     choice = 8)
  record <- rbind(record, temp)
  
  
  net_list[[i]] <- g
  
  print(i)
}
record <- record[-1,]



##
cen_sub <- cen_m[,c(1,3,5,7)]
corm <- cor(cen_sub)

p_values <- matrix(nrow = ncol(corm), ncol = nrow(corm))
for (i in 1:nrow(corm)) {
  for (j in 1:ncol(corm)) {
    p_values[i,j] <- cor.test(cen_sub[,i], cen_sub[,j])$p.value
  }
}

cor.test(cen_sub[,3], cen_sub[,4])





#task 4
#analyze diffusion process by method/social type
#these are datasets generated above. For replication, you could generate and name the data yourselves
c1 <- read.csv("record_c1.csv",header=T)
c1.5 <- read.csv("record_c1.5.csv",header=T)
c2 <- read.csv("record_c2.csv",header=T)
c2.5 <- read.csv("record_c2.5.csv",header=T)
c3 <- read.csv("record_c3.csv",header=T)
c3.5 <- read.csv("record_c3.5.csv",header=T)
c4 <- read.csv("record_c4.csv",header=T)
c4.5 <- read.csv("record_c4.5.csv",header=T)
m5 <- read.csv("record_m5.csv",header=T)
p1 <- read.csv("record_p1.csv",header=T)
p1.5 <- read.csv("record_p1.5.csv",header=T)
p2 <- read.csv("record_p2.csv",header=T)
p2.5 <- read.csv("record_p2.5.csv",header=T)
p3 <- read.csv("record_p3.csv",header=T)
p3.5 <- read.csv("record_p3.5.csv",header=T)
p4 <- read.csv("record_p4.csv",header=T)
p4.5 <- read.csv("record_p4.5.csv",header=T)

all <- list(p1,p1.5,p2,p2.5,p3,p3.5,p4,p4.5,m5,c4.5,c4,c3.5,c3,c2.5,c2,c1.5,c1)

#4.1 by method
#1 = sum, 2=deg, 3=clo, 4=bet, 5=eig, 6=dbm, 7=db, 8=de
#4.1.1 sum
mean_sum <- c()
sd_sum <- c()
for (i in 1:17){
  data <- all[[i]]
  data <- data[which(data$choice==1),]
  mean_sum <-c(mean_sum,mean(data$cov))
  sd_sum <- c(sd_sum,sd(data$cov))
  
}


#4.1.2 deg
mean_deg <- c()
sd_deg <- c()
for (i in 1:17){
  data <- all[[i]]
  data <- data[which(data$choice==2),]
  mean_deg <-c(mean_deg,mean(data$cov))
  sd_deg <- c(sd_deg,sd(data$cov))
  
}

#4.1.3 clo
mean_clo <- c()
sd_clo <- c()
for (i in 1:17){
  data <- all[[i]]
  data <- data[which(data$choice==3),]
  mean_clo <-c(mean_clo,mean(data$cov))
  sd_clo <- c(sd_clo,sd(data$cov))
  
}

#4.1.4 bet
mean_bet <- c()
sd_bet <- c()
for (i in 1:17){
  data <- all[[i]]
  data <- data[which(data$choice==4),]
  mean_bet <-c(mean_bet,mean(data$cov))
  sd_bet <- c(sd_bet,sd(data$cov))
  
}

#4.1.5 eig
mean_eig <- c()
sd_eig <- c()
for (i in 1:17){
  data <- all[[i]]
  data <- data[which(data$choice==5),]
  mean_eig <-c(mean_eig,mean(data$cov))
  sd_eig <- c(sd_eig,sd(data$cov))
  
}

#4.1.6dbm
mean_dbm <- c()
sd_dbm <- c()
for (i in 1:17){
  data <- all[[i]]
  data <- data[which(data$choice==6),]
  mean_dbm <-c(mean_dbm,mean(data$cov))
  sd_dbm <- c(sd_dbm,sd(data$cov))
  
}

#4.1.7 db
mean_db <- c()
sd_db <- c()
for (i in 1:17){
  data <- all[[i]]
  data <- data[which(data$choice==7),]
  mean_db <-c(mean_db,mean(data$cov))
  sd_db <- c(sd_db,sd(data$cov))
  
}

#4.1.8 de
mean_de <- c()
sd_de <- c()
for (i in 1:17){
  data <- all[[i]]
  data <- data[which(data$choice==8),]
  mean_de <-c(mean_de,mean(data$cov))
  sd_de <- c(sd_de,sd(data$cov))
  
}

t4_method <- data.frame(m_sum = mean_sum,
                        sd_sum = sd_sum,
                        m_deg = mean_deg,
                        sd_deg = sd_deg,
                        m_clo = mean_clo,
                        sd_clo = sd_clo,
                        m_bet = mean_bet,
                        sd_bet = sd_bet,
                        m_eig = mean_eig,
                        sd_eig = sd_eig,
                        m_dbm = mean_dbm,
                        sd_dbm = sd_dbm,
                        m_db = mean_db,
                        sd_db = sd_db,
                        m_de = mean_de,
                        sd_de = sd_de,
                        so = so)

#mean coverage plot
dat<-t4_method[,c(1,3,5,7,9,11,13,15)]
matplot(so,dat,type="b",pch=1,col=1:8,ylab="mean coverage",xlab="social openness",main="Mean Coverage Plot of Various Strategies")
legend("topright", legend = c("sum","deg","clo","bet","eig","deg+bet+eig","deg+bet","deg+eig"), col=1:8, pch=1)
axis(side = 1, at=1:9)

#standard deviation plot
dat2<-t4_method[,c(2,4,6,8,10,12,14,16)]
matplot(so,dat2,type="b",pch=1,col=1:8,ylab="coverage standard deviation",xlab="social openness",main="Coverage Standard Deviation Plot of Various Strategies")
legend("topright", legend = c("sum","deg","clo","bet","eig","deg+bet+eig","deg+bet","deg+eig"), col=1:8, pch=1)
axis(side = 1, at=1:9)

#maximal coverage plot
library(readxl)
data <- read_excel("record t3.xlsx")#it is another record in which I only record part information,but it is also done using the same codes

plot(data$mean,data$coverage, xlab = "social openness",xlim = c(1,9),
     ylab = "maximal coverage", main = "maximal coverage plot",pch = 19)



#full-rate plot
full_rate <- matrix(1,nrow=17,ncol=8)

for (i in 1:17){
  data <- all[[i]]
  for (j in 1:8){
    temp<-data[which(data$choice==j),]
    full_rate[i,j]<- sum(temp$cov==1000)
  }
  
}
matplot(so,full_rate,type="b",pch=1,col=1:8,ylab="full coverage rate",xlab="social openness",main="Full Coverage Rate Plot of Various Strategies")
legend("topright", legend = c("sum","deg","clo","bet","eig","deg+bet+eig","deg+bet","deg+eig"), col=1:8, pch=1)
axis(side = 1, at=1:9)

#efficiency plot
eff <- matrix(1,nrow=4,ncol=8)
for (i in 1:4){
  data <- all[[i]]
  for (j in 1:8){
    temp<-data[which(data$choice==j),3]
    eff[i,j]<- mean(1000/temp)
  }
  
}


matplot(c(1,1.5,2,2.5),eff,type="b",pch=1,col=1:8,ylab="mean efficiency",xlab="social openness",main="Average Efficiency Plot of Various Strategies")
legend("topright", legend = c("sum","deg","clo","bet","eig","deg+bet+eig","deg+bet","deg+eig"), col=1:8, pch=1,inset=-0.05, xpd=TRUE)







