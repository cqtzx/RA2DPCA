PCA2D <- function(x,k){#2DPCA
  x_dim <- dim(x)
  width <- x_dim[2]; p <- x_dim[3]
  cov <- matrix(0,width,width)
  for(i in 1:p){
    cov <- cov+t(x[,,i])%*%(x[,,i])
  }
  ei <- eigen(cov)
  V <- ei$vectors 
  D <- ei$values 
  W <- V[,1:k]
  W
}

RE_2DPCA <- function(x_train,ix_clean,k,ss){
  di <- dim(x_train)
  rowl <- di[1]; coll <- di[2]; nSub <- di[3] 
  E.mean <- matrix(0,rowl,coll)
  for(i in 1:nSub){
    E.mean <- E.mean + x_train[,,i]
  }
  x_mean <- E.mean/nSub
  x_centered <- array(0,c(rowl,coll,nSub))
  for(i in 1:nSub){
    x_centered[,,i]=x_train[,,i]-x_mean
  }
  W <- PCA2D(x_centered,k)
  x_reco <- array(0,dim(x_train))
  err <- rep(0,length(ss))
  j <- 1
  sss <- length(ix_clean)
  for(iK in ss){
    w <- W[,1:iK]
    for(iSub in 1:nSub){
      x_reco[,,iSub] <- x_centered[,,iSub]%*%w%*%t(w)+x_mean
   
    }
    temp <- x_train-x_reco
    sum <- 0
    for(i in 1:sss){
      sum <- sum+F2norm(temp[,,ix_clean[i]])
    }
    err[j] <- sum/sss
    j <- j+1
  }

  err
}

RA_2DPCA <- function(x,gamma1,lambda1,k){ 
  di <- dim(x)
  r <- di[1]; c <- di[2]; n <- di[3]
  W <- diag(1,c,k)
  E <- FF <- matrix(0,n,r)
  C <- matrix(0,c,c)
  err <- err1 <- h <- ff <- 1
  Lambda1 <- matrix(lambda1,n,r)

  dd <- Omega1 <- Lambda1

  EE <- matrix(0,c,c)
  Y <- array(0,c(r,c,n))
  while(err1 > 1e-4){

    ff.old <- ff
    for(i in 1:n){
      for(j in 1:r){
        Omega1[i,j] <- max(0,Lambda1[i,j]-F2norm(x[j,,i]-t(x[j,,i])%*%W%*%t(W))/gamma1)#F2norm(x[j,,i]%*%W)/gamma
        #Omega2[i,j] <- max(0,Lambda2[i,j]-F2norm((x[j,,i])%*%W)/gamma2)
        dd[i,j] <- Omega1[i,j]/(F2norm((x[j,,i]))*F2norm(x[j,,i]-t(x[j,,i])%*%W%*%t(W))+0.0001)
      }
    }
    EE <- matrix(0,c,c)
    for(i in 1:n){
      for(j in 1:r){
        EE <- EE+(dd[i,j]*(x[j,,i]))%*%t(x[j,,i])
      }
    }
    W1 <- eigen(EE)
    W <- W1$vectors[,1:k] 
    ff <- 0
    for(d in 1:n){
      ff <- ff+F2norm(x[,,d]-x[,,d]%*%W%*%t(W))
    }
    err <- rbind(err,F2norm(ff-ff.old)/F2norm(ff))
    err1 <- F2norm(ff-ff.old)/F2norm(ff)
    h <- h+1
    if(h>100){
      break;
    }
  }
  W
}


RE_RA_2DPCA <- function(x_train,ix_clean,gamma1,lambda1,k,ss){
  di <- dim(x_train)
  rowl <- di[1]; coll <- di[2]; nSub <- di[3] 
  E.mean <- matrix(0,rowl,coll)
  for(i in 1:nSub){
    E.mean <- E.mean + x_train[,,i]
  }
  x_mean <- E.mean/nSub
  x_centered <- array(0,c(rowl,coll,nSub))
  for(i in 1:nSub){
    x_centered[,,i]=x_train[,,i]-x_mean
  }
  W <- RA_2DPCA(x_centered,gamma1,lambda1,k)
  x_reco <- array(0,dim(x_train))
  err <- rep(0,length(ss))
  j <- 1
  sss <- length(ix_clean)
  for(iK in ss){
    w <- W[,1:iK]
    for(iSub in 1:nSub){
      x_reco[,,iSub] <- x_centered[,,iSub]%*%w%*%t(w)+x_mean
    
    }
    temp <- x_train-x_reco
    sum <- 0
    for(i in 1:sss){
      sum <- sum+F2norm(temp[,,ix_clean[i]])
    }
    err[j] <- sum/sss
    j <- j+1
  }
  #plot(kSet,err)
  err
}


OMRA_2DPCA <- function(x,gamma1,lambda1,k){ #GC-2DPCA-l_2,p
  di <- dim(x)
  r <- di[1]; c <- di[2]; n <- di[3]
  W <- diag(1,c,k)
  E <- FF <- matrix(0,n,r)
  C <- matrix(0,c,c)
  err <- err1 <- h <- ff <- 1
  Lambda1  <- matrix(lambda1,n,r)
  
  dd <- Omega1 <- Lambda1

  EE <- matrix(0,c,c)
  Y <- array(0,c(r,c,n))
  M <- matrix(0,r,c)
  while(err1 > 1e-4){

    ff.old <- ff
    for(i in 1:n){
      for(j in 1:r){
        Omega1[i,j] <- max(0,Lambda1[i,j]-F2norm(x[j,,i]-M[j,]-t(x[j,,i]-M[j,])%*%W%*%t(W))/gamma1)#F2norm(x[j,,i]%*%W)/gamma
       
        dd[i,j] <- Omega1[i,j]/(F2norm((x[j,,i])-M[j,])*F2norm(x[j,,i]-M[j,]-t(x[j,,i]-M[j,])%*%W%*%t(W))+0.0001)
      }
    }
    #M
    for(j in 1:r){
      sumM <- matrix(0,1,c)
      sumd <- 0
      for(i in 1:n){
        sumM <- sumM+dd[i,j]*x[j,,i]
        sumd <- sumd+dd[i,j]
      }
      M[j,] <- sumM/sumd
    }
    EE <- matrix(0,c,c)
    for(i in 1:n){
      for(j in 1:r){
        EE <- EE+(dd[i,j]*(x[j,,i]-M[j,]))%*%t(x[j,,i]-M[j,])
      }
    }
    W1 <- eigen(EE)
    W <- W1$vectors[,1:k] ;
    ff <- 0
    for(d in 1:n){
      ff <- ff+F2norm(x[,,d]-M-(x[,,d]-M)%*%W%*%t(W))
    }
    err <- rbind(err,F2norm(ff-ff.old)/F2norm(ff))
    err1 <- F2norm(ff-ff.old)/F2norm(ff)
    h <- h+1
    if(h>100){
      break;
    }
  }
  list(W,M)
}

RE_OMRA_2DPCA <- function(x_train,ix_clean,gamma1,lambda1,k,ss){
  di <- dim(x_train)
  rowl <- di[1]; coll <- di[2]; nSub <- di[3] #很多张输入进去一起重构的
  E.mean <- matrix(0,rowl,coll)
  for(i in 1:nSub){
    E.mean <- E.mean + x_train[,,i]
  }
  x_mean <- E.mean/nSub
  x_centered <- array(0,c(rowl,coll,nSub))
  WW <- OMRA_2DPCA(x_train,gamma1,lambda1,k);M<-WW[[2]]
  for(i in 1:nSub){
    x_centered[,,i]=x_train[,,i]-M
  }
  W<-WW[[1]];
  x_reco <- array(0,dim(x_train))
  err <- rep(0,length(ss))
  j <- 1
  sss <- length(ix_clean)
  for(iK in ss){
    w <- W[,1:iK]
    for(iSub in 1:nSub){
      x_reco[,,iSub] <- x_centered[,,iSub]%*%w%*%t(w)+M
    }
    temp <- x_train-x_reco
    sum <- 0
    for(i in 1:sss){
      sum <- sum+F2norm(temp[,,ix_clean[i]])
    }
    err[j] <- sum/sss
    j <- j+1
  }
  err
}
