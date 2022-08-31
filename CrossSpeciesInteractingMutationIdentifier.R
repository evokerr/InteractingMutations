c.D.df <- read.csv("AllCocultures_DexceptUA3 (2).csv",header=TRUE)
c.M.df <- read.csv("AllCocultures_M.csv",header=TRUE)
m.D.df <- read.csv("AllMonocultures_D (1).csv",header=TRUE)
m.M.df <- read.csv("AllMonocultures_M.csv",header=TRUE)

#This function gets mutations that are found in multiple co-cultures for a given species
get.shared.mutations <- function(df, critical.sorting.column.number) {
  culture.names <- unique(df$culture)
  shared.mutations <- character(0)
  for(i in 1:(length(culture.names)-1)) {
    for(j in (i+1):length(culture.names)) {
      intersecting.mutations <- intersect(df[df$culture==culture.names[i],critical.sorting.column.number],df[df$culture==culture.names[j],critical.sorting.column.number])
      if(length(intersecting.mutations)>0) {
        shared.mutations<-c(shared.mutations,intersecting.mutations)
      }
    }
  }
  return(unique(shared.mutations))
}

#This function returns a list of pairs of mutations in an "initiating" species 
#and "responding" species, where the same pair exists in multiple co-cultures and
#the mutation from the responding species is not in any relevant monoculture.
get.pairs.list <- function(c.init.df, c.res.df, m.res.df, c.critical.sorting.column.number, m.critical.sorting.column.number) {
  init.shared.mutations <- get.shared.mutations(c.init.df, c.critical.sorting.column.number)
  res.shared.mutations <- get.shared.mutations(c.res.df, c.critical.sorting.column.number)
  only.coculture.res.shared.mutations <- character(0)
  for(m in res.shared.mutations) {
    if(!(m %in% unique(m.res.df[,m.critical.sorting.column.number]))) {
      only.coculture.res.shared.mutations <- c(only.coculture.res.shared.mutations,m)
    }
  }
  pairs.list <- vector(mode="list",length=0)
  for(m.i in init.shared.mutations) {
    for(m.r in only.coculture.res.shared.mutations) {
      init.cocultures <- unique(c.init.df[c.init.df[,c.critical.sorting.column.number]==m.i,]$culture)
      res.cocultures <- unique(c.res.df[c.res.df[,c.critical.sorting.column.number]==m.r,]$culture)
      cocultures.with.focal.pair <- intersect(init.cocultures, res.cocultures)
      if(length(cocultures.with.focal.pair)>1) {
        pairs.list<-append(pairs.list,"")
        pairs.list[[length(pairs.list)]] <- c(m.i,m.r)
      }
    }
  }
  return(pairs.list)
}

#This function returns a list of pairs of times where the "initiator" mutation
#and "responder" mutation in all relevant co-cultures first pass a critical
#frequency.
find.times <- function(pair.of.mutations, c.init.df, c.res.df, critical.sorting.column.number, critical.frequency) {
  times<-c(0,100,300,500,780,1000)
  init.cocultures <- unique(c.init.df[c.init.df[,critical.sorting.column.number]==pair.of.mutations[1],]$culture)
  res.cocultures <- unique(c.res.df[c.res.df[,critical.sorting.column.number]==pair.of.mutations[2],]$culture)
  cocultures.with.focal.pair <- intersect(init.cocultures, res.cocultures)
  times.list <- vector(mode="list",length=0)
  for(cc in cocultures.with.focal.pair) {
    new.c.init.df <- c.init.df[c.init.df$culture==cc,]
    new.c.res.df <- c.res.df[c.res.df$culture==cc,]
    filtered.new.c.init.df<-new.c.init.df[new.c.init.df[,critical.sorting.column.number]==pair.of.mutations[1],]
    min.t.init<-1000
    for(r in 1:length(filtered.new.c.init.df$culture)) {
      for(c in 11:16){
        if(filtered.new.c.init.df[r,c]>critical.frequency) {
          t.init <- times[c-10]
          break
        }
      }
      if(t.init < min.t.init) {
        min.t.init <- t.init
      }
    }
    filtered.new.c.res.df<-new.c.res.df[new.c.res.df[,critical.sorting.column.number]==pair.of.mutations[2],]
    max.t.res <- 0
    for(r in 1:length(filtered.new.c.res.df$culture)) {
      for(c in 11:16){
        if(filtered.new.c.res.df[r,c]>critical.frequency) {
          t.res <- times[c-10]
          break
        }
      }
      if(t.res > max.t.res) {
        max.t.res <- t.res
      }
    }
    times.list<-append(times.list,0)
    times.list[[length(times.list)]] <- c(min.t.init,max.t.res)
  }
  return(times.list)
}

#This function takes a list of pairs of mutations and performs a permutation test
#on the list of times of appearance of the mutations across replicate co-cultures
#and assembles a vector of p-values with the same indices as the pairs.list argument.
get.p.values.list <- function(pairs.list, c.init.df, c.res.df, critical.sorting.column.number, critical.frequency, n.perms, strictly.greater) {
  p.vec <- numeric(0)
  for(p in 1:length(pairs.list)) {
    t.list <- find.times(pairs.list[[p]], c.init.df, c.res.df, critical.sorting.column.number, critical.frequency)
    p.val <- permutation.test(t.list, n.perms, strictly.greater)
    p.vec <- c(p.vec,p.val)
  }
  return(p.vec)
}

#Returns a p-value for the fraction of permutations on the responder times for
#which the metric of the bootstrapped data is more extreme than the actual
#metric.
permutation.test <- function(list.of.pairs, n.perms, strictly.greater) {
  actual.metric <- compute.metric(list.of.pairs)
  bootstrapped.metrics <- integer(n.perms)
  for(r in 1:n.perms) {
    permuted.second.indices <- sample(1:length(list.of.pairs),length(list.of.pairs))
    bootstrapped.list.of.pairs <- vector(mode = "list", length = length(list.of.pairs))
    for(i in 1:length(list.of.pairs)) {
      bootstrapped.list.of.pairs[[i]] <- c(list.of.pairs[[i]][1],list.of.pairs[[permuted.second.indices[[i]]]][2])
    }
    bootstrapped.metrics[r] <- compute.metric(bootstrapped.list.of.pairs)
  }
  counter <- 0
  for(r in 1:n.perms) {
    if(strictly.greater) {
      if(actual.metric < bootstrapped.metrics[r]) {
        counter <- counter + 1
      }
    } else {
      if(actual.metric <= bootstrapped.metrics[r]) {
        counter <- counter + 1
      }
    }
  }
  return(counter/n.perms)
}

#Computes a score of +1 for an ordering of times of first initiator then responder
#mutation, -1 for the reverse, and 0 for a tie, where each co-culture's score is
#summed to get the metric.
compute.metric <- function(list.of.pairs) {
  metric <- 0
  for(e in 1:length(list.of.pairs)) {
    if(list.of.pairs[[e]][1]<list.of.pairs[[e]][2]){
      metric <- (metric + 1)
    } else {
      if(list.of.pairs[[e]][1]==list.of.pairs[[e]][2]) {
        metric <- metric
      } else {
        metric <- (metric - 1)
      }
    }
  }
  return(metric)
}

#This function takes a list of pairs of mutations and checks the fraction of co-cultures
#on the list where the initiator mutation occurs at or before the responder mutation
#and assembles a vector of prop-values with the same indices as the pairs.list argument.
get.prop.ordering.list <- function(pairs.list, c.init.df, c.res.df, critical.sorting.column.number, critical.frequency) {
  p.vec <- numeric(0)
  for(p in 1:length(pairs.list)) {
    t.list <- find.times(pairs.list[[p]], c.init.df, c.res.df, critical.sorting.column.number, critical.frequency)
    counter <- 0
    for(t.ind in 1:length(t.list)) {
      if(t.list[[t.ind]][1] <= t.list[[t.ind]][2]) {
        counter <- counter + 1
      }
    }
    p.vec <- c(p.vec,counter/length(t.list))
  }
  return(p.vec)
}


mut.pairs.D.M <- get.pairs.list(c.init.df=c.D.df, 
                            c.res.df=c.M.df, 
                            m.res.df=m.M.df, 
                            c.critical.sorting.column.number = 43, 
                            m.critical.sorting.column.number = 29) 

p.values.D.M <- get.p.values.list(pairs.list=mut.pairs.D.M, 
                              c.init.df=c.D.df,
                              c.res.df=c.M.df,
                              critical.sorting.column.number = 43, 
                              critical.frequency = 0.0, 
                              n.perms = 1000, 
                              strictly.greater = FALSE)

prop.values.D.M <- get.prop.ordering.list(pairs.list=mut.pairs.D.M, 
                              c.init.df=c.D.df,
                              c.res.df=c.M.df,
                              critical.sorting.column.number = 43, 
                              critical.frequency = 0.0)

mut.pairs.M.D <- get.pairs.list(c.init.df=c.M.df, 
                            c.res.df=c.D.df, 
                            m.res.df=m.D.df, 
                            c.critical.sorting.column.number = 43, 
                            m.critical.sorting.column.number = 29) 

p.values.MD <- get.p.values.list(pairs.list=mut.pairs.M.D, 
                              c.init.df=c.M.df,
                              c.res.df=c.D.df,
                              critical.sorting.column.number = 43, 
                              critical.frequency = 0.0, 
                              n.perms = 1000, 
                              strictly.greater = TRUE)

prop.values.M.D <- get.prop.ordering.list(pairs.list=mut.pairs.M.D, 
                                      c.init.df=c.M.df,
                                      c.res.df=c.D.df,
                                      critical.sorting.column.number = 43, 
                                      critical.frequency = 0.0)


mut.pairs.D.M
p.values.D.M
prop.values.D.M

mut.pairs.M.D
p.values.M.D
prop.values.M.D


