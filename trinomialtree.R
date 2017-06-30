trinomialtree<-function(TypeFlag=c("EC","EP","AC","AP"),S0,K,T,r,sigma,N,div){
  #constants for trinomial calculation
  dt = T/N
  nu = r-div-0.5*sigma^2
  dx = sigma*sqrt(3*dt)
  pu = 0.5*((sigma^2*dt+nu^2*dt^2)/dx^2+nu*dt/dx)
  pm = 1-(sigma^2*dt+nu^2*dt^2)/dx^2
  pd = 0.5*((sigma^2*dt+nu^2*dt^2)/dx^2-nu*dt/dx)
  disc = exp(-r*dt)
  firstRow = 1
  r = nRows = lastRow = 2*N+1
  firstCol = 1
  middleRow = s = nCols = lastCol = N+1
  St = ValueofOption = matrix(0,nrow = nRows,ncol=nCols,dimnames = list(paste("NumUps",N:-N,sep="="),paste("T",0:N,sep="=")))
  
  St[middleRow,firstCol] = S0           ##initalize asset prices at maturity
  for(j in 1:(nCols-1)){
    for(i in (middleRow-j+1):(middleRow+j-1)){
      St[i-1,j+1] = St[i,j]*exp(dx)
      St[i,j+1] = St[i,j]
      St[i+1,j+1]=St[i,j]*exp(-dx)
    }
  }
  if(TypeFlag=="EC"){                     ##initialize option value at maturity
    for(i in 1:nRows){
      ValueofOption[i,lastCol] = max(0,St[i,lastCol]-K)
    }
    for (j in (nCols-1):1){              ##Step back through the tree
      for(i in (nCols-j+1):(nCols+j-1)){
        ValueofOption[i,j] = disc*(pu*ValueofOption[i-1,j+1]+pm*ValueofOption[i,j+1]+pd*ValueofOption[i+1,j+1])
      }
    }
  } else if(TypeFlag=="EP"){                     ##initialize option value at maturity
    for(i in 1:nRows){
      ValueofOption[i,lastCol] = max(0,K-St[i,lastCol])
    }
    for (j in (nCols-1):1){              ##Step back through the tree
      for(i in (nCols-j+1):(nCols+j-1)){
        ValueofOption[i,j] = disc*(pu*ValueofOption[i-1,j+1]+pm*ValueofOption[i,j+1]+pd*ValueofOption[i+1,j+1])
      }
    }
  } else if(TypeFlag=="AC"){
    for(i in 1:nRows){
      ValueofOption[i,lastCol] = max(0,St[i,lastCol]-K)
    }
    for(j in (nCols-1):1){
      for(i in (nCols-j+1):(nCols+j-1)){
        ValueofOption[i,j] = disc*(pu*ValueofOption[i-1,j+1]+pm*ValueofOption[i,j+1]+pd*ValueofOption[i+1,j+1])
        ValueofOption[i,j] = max(ValueofOption[i,j],St[i,j]-K)
      }
    }
  } else if(TypeFlag=="AP"){
    for(i in 1:nRows){
      ValueofOption[i,lastCol] = max(0,K-St[i,lastCol])
    }
    for(j in (nCols-1):1){
      for(i in (nCols-j+1):(nCols+j-1)){
        ValueofOption[i,j] = max(disc*(pu*ValueofOption[i-1,j+1]+pm*ValueofOption[i,j+1]+pd*ValueofOption[i+1,j+1]),K-St[i,j])
        ValueofOption[i,j] = max(ValueofOption[i,j],K-St[i,j])
      }
    }
  }

  return(ValueofOption[middleRow,firstCol])
}
