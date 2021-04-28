if (!require(plotly)) install.packages('plotly')
if (!require(ggplot2)) install.packages('ggplot2')
library(plotly)
library(ggplot2)

options(scipen=999)

# formatka do testowania
TEST=Firefly_Algorithm(Iterations = 1000, 
                       Bounds = 2, 
                       RandomCoef = 0.3, 
                       MaxGeneration = 30, 
                       Function = OPT_Sphere)


# Generowanie wartoÅci funkcji celu (Light intensity) // zostawiÄ parametr plot="NO", byÅ dla celÃ³w testownaia 
Get_Light_Intensity = function(A,Function=OPT_Sphere, plot="NO"){
  
  if (is.null(nrow(A))){
    return(list(as.numeric(Function(A[1],A[2]))))
  }else{
    Light_Intensity=matrix(ncol=1,nrow=nrow(A),NA)
    for (i in 1:nrow(A)){
      Light_Intensity=Function(A[,1],A[,2])
    }
    if (plot!="NO"){
      temp <- Opt_Value
      Graph=plot_ly(x=A[,1], y=A[,2], z=Opt_Value, type="scatter3d", mode="markers" , color=temp, alpha=0.6) #scatter3d, mesh3d
      return(list(as.matrix(cbind(A,Light_Intensity)),Graph))
    }else{
      return(list(as.matrix(cbind(A,Light_Intensity))))
    }
  }
}

Poziomica=function(from=-5,to=5,by=0.1,Function){
  
  nr=(to-from)/by+1
  
  Warstwice=matrix(nrow=nr^2,ncol=5,NA)
  sequence=seq(to=to,from=from,by=by)
  
  Warstwice[,1]=rep(sequence,nr)
  Warstwice[,2]=rep(sequence,each=nr)
  
  for(i in 1:(nr^2)){
    Warstwice[i,3]=Function(Warstwice[i,1],Warstwice[i,2])
  }
  return(Warstwice)
}

# Core algorithm
Firefly_Algorithm=
  function(
    Iterations,
    MaxAttraction=1,
    MaxGeneration=20,
    LighAbsorbtion=1,
    RandomCoef=0.2,
    Fireflies_Count=20,
    Bounds=2,
    Function=OPT_Sphere){
    
    Fireflies_Coord=matrix(nrow=Fireflies_Count,ncol=2,NA)
    colnames(Fireflies_Coord)=c("X","Y")
    Fireflies_Coord=apply(Fireflies_Coord, c(1,2), function(x){-Bounds+2*Bounds*runif(1,0,1)})
    
    # Wylicza jasnoÅÄ poczÄtkowych robakÃ³w 
    Fireflies_Coord=Get_Light_Intensity(Fireflies_Coord,Function)[[1]]
    
    RandNorm=rnorm(Fireflies_Count)
    REPO=matrix(ncol=3,nrow=Iterations, NA)
    
    for (G in 1:Iterations){
      for (i in 1:Fireflies_Count){
        for (j in 1:i){
          
          # PorÃ³wnywanie jasnoÅci ÅwietlikÃ³w i przybliÅ¼anie do najjaÅniejszego
          if (Fireflies_Coord[j,"Light_Intensity"]>Fireflies_Coord[i,"Light_Intensity"]){
            
            Dist=dist(x = Fireflies_Coord[c(i,j),])
            
            Fireflies_Coord[i,"X"]=Fireflies_Coord[i,"X"]+MaxAttraction*(exp(-1*LighAbsorbtion*(Dist^2)))*
              (Fireflies_Coord[j,"X"]-Fireflies_Coord[i,"X"])+RandomCoef*rnorm(1)
            
            Fireflies_Coord[i,"Y"]=Fireflies_Coord[i,"Y"]+MaxAttraction*(exp(-1*LighAbsorbtion*(Dist^2)))*
              (Fireflies_Coord[j,"Y"]-Fireflies_Coord[i,"Y"])+RandomCoef*rnorm(1)
            
          }
          
          Fireflies_Coord[i,3]=Get_Light_Intensity(Fireflies_Coord[i,],Function)[[1]]
          Fireflies_Coord[j,3]=Get_Light_Intensity(Fireflies_Coord[j,],Function)[[1]]
        }
      }
      # Zapisanie koordynatÃ³w najjaÅniejszych ÅwietlikÃ³w
      REPO[G,]=as.vector(Fireflies_Coord[which.max(Fireflies_Coord[,3]),])
    }
    
    # Wykres i jego wymiary (zamiast stosowaÄ zakres i zakres2 moÅ¼na przekazaÄ customowe argumenty do funkcji Poziomica, zamiast defaultowych)
    Warstwice=Poziomica(Function=Function)
    zakres=max(Warstwice[,1],Warstwice[,2])
    zakres2=min(Warstwice[,1],Warstwice[,2])
    
    # poczÄtkowy odsetek ÅwietlikÃ³w, ktÃ³re znajdujÄ siÄ w obrÄbie wykresu (czyli domyÅlnie tych, ktÃ³re zbiegajÄ do optimum)
    Success=1
    
    # Zabieg techniczny. Å»eby narysowaÄ wartstwice i punkty na jednym wykresie, ggplot potrzebuje Å¼eby pochodziÅy z tego samego df
    for (i in 1:Fireflies_Count){
      if (Fireflies_Coord[i,1]<zakres && Fireflies_Coord[i,1]>zakres2 && Fireflies_Coord[i,2]<zakres && Fireflies_Coord[i,2]>zakres2){
        Warstwice[i,4]=Fireflies_Coord[i,1]
        Warstwice[i,5]=Fireflies_Coord[i,2]
      }else{
        Success=Success-(1/Fireflies_Count)
      }
    }
    
    # Z jakiegoÅ powodu nie dziaÅa nadanie gradientu kolorÃ³w, ktÃ³ry ma wskazywaÄ na wartoÅc funkcji celu. 
    Warstwice=as.data.frame(Warstwice)
    colnames(Warstwice)=c("X","Y","Z","F_X","F_Y")
    
    v = ggplot(Warstwice, aes(X, Y, z = Z, fill=Z)) + 
      stat_contour(aes()) +
      scale_fill_gradient(low = "blue",
                          high = "red",
                          guide = "colourbar") +
      geom_point(aes(x=F_X, y=F_Y)) +
      ggtitle(paste("Warstwice, (",Success*100,"% Of ÅwietlikÃ³w w zakresie)",sep = "")) +
      theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank())
    
    
    print(v)
    print(Fireflies_Coord[which.max(Fireflies_Coord[,3]),])
    return(list(REPO,v))
  }


