NR      <- 10 #Number of richness gradients
NAm     <-  9 #Number of amplitudes of fluctuation
Am      <- seq(0.1, 0.9, length.out = NAm) #Values of fluctuation amplitudes
NFreq   <- 9
Freq    <- c(1, 2, 4, 8, 16, 32, 64, 128, 256 )
Nrep    <- 10  #Number of replications at each gradient
nyear   <- 5
d_per_y <- 365
newx    <- 1:d_per_y

Vol2ESD <- function(V) (6*V/pi)**.3333333333333333333333333

#Functions to get data
get_data <- function(v,i,k=1,l, type = 'poly'){
  #v: amplitude
  #i: richness
  #k: number of repetition
  #l: frequency
  if (type == 'poly'){
    file <- paste0('P',v,'_',i,'_',k,'_',l,'.out')
    dat  <- read.table(file)
    
    #Obtain sizes
    z <- dat[1,]
    
    #Get size info
    N <- ncol(z)
    
    #The number of size classes
    N <- (N-2)/2

    s <- z[3:(2+N)]
      
    #Obtain data
    dat <- dat[-1,]
    colnames(dat)[1] <- as.character(z[1,1])
    colnames(dat)[2] <- as.character(z[1,2])
    colnames(dat)[3:ncol(dat)] <- as.character(z[1,3:ncol(z)])
    
    return(list(data=dat, sizes=as.numeric(s)))
  }else if (type == 'mono'){
    #Extract all the information of monocultures
    file <- paste0('M',v,'_',i,'_',k,'_',l,'_1.out')
    z    <- read.table(file, header=T)
    
    #Vector storing size information
    LNVs <- numeric(i)
    dat <- matrix(NA, nr=nrow(z), nc=1+(ncol(z)-1) * i)
    
    for (n in 1:i){
      file <- paste0('M',v,'_',i,'_',k,'_',l,'_',n,'.out')
      d <- read.table(file)
      
      z <- d[1,]

      #Obtain size of this monoculture
      LNVs[n]   <- as.numeric(z[3])

      #Remove the first row and obtain the real data
      z <- d[-1,]
      names(z) <- c('Day', 'NO3', 'PHY', 'NPP')
      for (m in 1:ncol(z)) z[,m] <- as.numeric(z[,m])

      if (n == 1){
        dat[,1] <- z$Day
        for (k in 2:ncol(z)) dat[,k] <- z[,k]
         
      }else{
        #Count the starting column
        ic <- 1 + (ncol(z)-1)*(n-1) + 1
        for (k in 2:ncol(z)) dat[,k+ic-2] <- z[,k]
      }
    }     
    return(list(sizes=LNVs, data=dat))
  }else{
    stop("Incorrect type!")
  }
}

#Get data of the final year 
get_LY_PP <- function(v,i,k=1,l, type = 'poly'){
   cff    <- get_data(v,i,k,l, type)
   
   if (type == 'poly') {
     dat <- cff$data

     #Force  dat to be numerical values
     for (u in 1:ncol(dat)) dat[,u] <- as.numeric(as.character(dat[,u]))

     Nout   <- nrow(dat)
     NY     <- (Nout-1)/nyear
     ff     <- (Nout-NY+1):Nout #Final year
     
     NO3s   <- as.numeric(dat$NO3[ff])

     #Mean NO3 annual mean concentration
     NO3_mean <- mean(NO3s)
 
     #Obtain phyto. biomass
     wx <- which(names(dat) == 'NO3')

     #Obtain number of size classes
     NSize <- length(cff$sizes)

     #Obtain biomass
     B <- dat[ff, (wx+1):(wx+NSize)]

     #Obtain production
     P <- dat[ff, (wx+1+NSize):(wx+2*NSize)]

     #Calculate total phytoplankton biomass
     PHYt <- apply(B, 1, sum)

     #Calculate total primary production
     PPt <- apply(P, 1, sum)

     #Obtain average biomass and production for each species
     Bavg <- apply(B, 2, mean)
     PPavg<- apply(P, 2, mean)

     PHYm   <- mean(PHYt) #Mean total phyto. biomass
     PPm   <- mean(PPt) #Mean total primary production
     
     #Mean of mean trait weighted by biomass
     LNVs <- dat$LNV[ff]
     
     #Proportion of phyto. biomass in each day
     p <- PHYt/sum(PHYt)
     avgL <- sum(p*LNVs)
     
     #Variance of mean trait weighted by biomass
     VarL <- sum(p*(LNVs - avgL)^2)

     return(list(sizes = cff$sizes,
                       NO3 = NO3_mean,
                      PHY = PHYm, 
                      PPt = PPm,
                      Biomass= Bavg,
                      NPP=PPavg))    
   }else if (type == 'mono'){
     sizes  <- as.numeric(cff$sizes)
     cff    <- cff$data
     Nout   <- nrow(cff)
     NY     <- (Nout-1)/nyear
     ff     <- (Nout-NY+1):Nout #Obtain the time indexes for the final year
     
     #Calculate the annual mean biomass and productivity of each species in monocultures 
     Nsp <- length(sizes)
     out <- matrix(NA, nr = Nsp, nc = 3)
     for (kk in 1:Nsp){
       #Column index for cff
       w <- 1+3*(kk-1)
       out[kk,1] <- mean(cff[ff, w+1]) #NO3
       out[kk,2] <- mean(cff[ff, w+2]) #Phyto. biomass
       out[kk,3] <- mean(cff[ff, w+3]) #PP
     }
     out <- as.data.frame(out)
     colnames(out) <- c('NO3', 'PHY', 'PP')
     return(list(sizes=sizes, data=out))
   }
}

#Calculate PHY ~ richness  at two different amplitudes
#Calculate average biomass and productivity of the final year for each gradient
#Phytoplankton biomass plot
PHYv <- array(NA, dim = c(NAm, NR-1, Nrep, NFreq))
PPv    <- PHYv
mu     <- PHYv

#Array storing the size information of each treatment/replicate
SIZES <- array(list(), dim= c(NAm, NR-1, Nrep, NFreq))

#Average NO3 concentration
NO3 <- PHYv

#Selection effect
Sel_B <- PHYv  #Based on biomass
Sel_P <- PHYv  #Based on production

#Complementarity effect
Comp_B <- PHYv  #Based on biomass
Comp_P <- PHYv  #Based on production

for (v in 1:NAm){
  for (i in 2:NR){
    for (k in 1:Nrep){
      for (l in 1:NFreq){
        z          <- get_LY_PP(v, i, k-1, l)
       PHYv[v,i-1,k,l] <- z$PHY   #average phyto. biomass
        PPv[v,i-1,k,l] <- z$PPt    #Average primary production

        #Obtain size info from polycultures
        p_size <- z$sizes

      #Calculate effects of selection and complementarity 
      #following Loreau and Hector (2001)
      
      #Obtain monocultures
      z1 <- get_LY_PP(v, i, k-1, l, type = 'mono')
      
      #Obtain sizes of monocultures
      m_size <- z1$sizes

      #Monoculture biomass
      M <- z1$data[,2]
      
      #Monoculture production
      MP <- z1$data$PP

      #sort monoculture biomass and production
      M <- M[order(m_size)]
      MP <- MP[order(m_size)]

      #Store size information
      SIZES[v,i-1,k,l][[1]] <- m_size[order(m_size)]

      #Polyculture biomass extracted from z
      Y <- as.numeric(z$Biomass)

      #sort polyculture biomass
      Y <- Y[order(p_size)]
      
      #Compute relative yield based on biomass
      dRY <- Y/M - 1/i
      
      Comp_B[v,i-1,k,l] <- i*mean(dRY)*mean(M)

      Sel_B[v,i-1,k,l] <- i*cov(dRY, M)

      #Polyculture production extracted from z
      Y <- as.numeric(z$NPP)

      #sort polyculture production
      Y <- Y[order(p_size)]
      
      #Compute relative yield based on biomass
      dRY <- Y/MP - 1/i
      
      Comp_P[v,i-1,k,l] <- i*mean(dRY)*mean(MP)

      Sel_P[v,i-1,k,l] <- i*cov(dRY, MP)

    }
   }
  }
}

save(PHYv, PPv, Comp_B, Sel_B,Comp_P, Sel_P, SIZES, file = 'SelCOMP.Rdata')