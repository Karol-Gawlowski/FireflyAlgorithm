library(plotly)
library(ggplot2)

options(scipen=999)
set.seed(1)

# Test the algorithm
TEST=Firefly_Algorithm(Iterations = 10000, 
                       Bounds = 5, # limits for initialization of fireflies' positions
                       RandomCoef = 0.3, 
                       LightAbsorbtion=2, # >1 the more attracted the fireflies' are to each other 
                       Fireflies_Count=30,
                       Function = OPT_Xin_She_Yang)

# Generate the light intensity (test function value)
# plot parameter is used for debugging
Get_Light_Intensity = function(A,
                               Function=OPT_Sphere, 
                               plot="NO"){
  
  if (is.null(nrow(A))){
    return(list(as.numeric(Function(A[1],A[2]))))
  }else{
    Light_Intensity=Function(A[,1],A[,2])
    if (plot!="NO"){
      temp <- Opt_Value
      Graph=plot_ly(x=A[,1], y=A[,2], z=Opt_Value, type="scatter3d", mode="markers" , color=temp, alpha=0.6) #scatter3d, mesh3d
      return(list(as.matrix(cbind(A,Light_Intensity)),Graph))
    }else{
      return(list(as.matrix(cbind(A,Light_Intensity))))
    }
  }
}

# prepare level curves
LevelCurve_Function=function(from=-5,
                             to=5,
                             by=0.1,
                             Function=OPT_Sphere){
  
  nr=(to-from)/by+1
  
  LevelCurves=matrix(nrow=nr^2,ncol=5,NA)
  sequence=seq(to=to,from=from,by=by)
  
  LevelCurves[,1]=rep(sequence,nr)
  LevelCurves[,2]=rep(sequence,each=nr)
  
  for(i in 1:(nr^2)){
    LevelCurves[i,3]=Function(LevelCurves[i,1],LevelCurves[i,2])
  }
  return(LevelCurves)
}

# Core algorithm
Firefly_Algorithm=
  function(
    Iterations,
    MaxAttraction=1,
    LightAbsorbtion=1,
    RandomCoef=0.2,
    Fireflies_Count=20,
    Bounds=5,
    Function=OPT_Sphere){

    # Initiatial Fireflies' coordinates
    Fireflies_Coord=matrix(runif(Fireflies_Count*2,-Bounds,Bounds), ncol=2)
    colnames(Fireflies_Coord)=c("X","Y")
    
    # Initial Fireflies' brightness intensity
    Fireflies_Coord=Get_Light_Intensity(Fireflies_Coord,Function)[[1]]
    
    #
    REPO=matrix(ncol=3,nrow=Iterations, NA)
    
    for (G in 1:Iterations){
      for (i in 1:Fireflies_Count){
        for (j in 1:i){
          
          # Compare light intensities and move dimmer Fireflies to the brighter ones
          if (Fireflies_Coord[j,"Light_Intensity"]>Fireflies_Coord[i,"Light_Intensity"]){
            
            Dist=dist(x = Fireflies_Coord[c(i,j),])
            
            Fireflies_Coord[i,"X"]=Fireflies_Coord[i,"X"]+MaxAttraction*(exp(-1*LightAbsorbtion*(Dist^2)))*
              (Fireflies_Coord[j,"X"]-Fireflies_Coord[i,"X"])+RandomCoef*rnorm(1)
            
            Fireflies_Coord[i,"Y"]=Fireflies_Coord[i,"Y"]+MaxAttraction*(exp(-1*LightAbsorbtion*(Dist^2)))*
              (Fireflies_Coord[j,"Y"]-Fireflies_Coord[i,"Y"])+RandomCoef*rnorm(1)
            
            # New brightness intensity
            Fireflies_Coord[i,3]=Get_Light_Intensity(Fireflies_Coord[i,],Function)[[1]]
            Fireflies_Coord[j,3]=Get_Light_Intensity(Fireflies_Coord[j,],Function)[[1]]
          }
        }
      }
      
      # Save coordinates of the lightest Fireflies
      REPO[G,]=as.vector(Fireflies_Coord[which.max(Fireflies_Coord[,3]),])
    }
    
    # get bounds
    LevelCurves=LevelCurve_Function(from = -Bounds,
                                    to = Bounds,
                                    by=0.2,
                                    Function=Function)

    # for plot centering, these bounds depend on the number of level curves generated
    U_Bound=max(LevelCurves[,c(1,2)])
    L_Bound=min(LevelCurves[,c(1,2)])
    
    # initial prop of fireflies converging to optimum
    Success=1
    
    # Prep data for ggplot
    for (i in 1:Fireflies_Count){
      # overlay fireflies' positions on the level curves
      # note that the way these coordinates are stored is such, just for the way ggplot takes in the objects
      if (Fireflies_Coord[i,1]<U_Bound && Fireflies_Coord[i,1]>L_Bound && Fireflies_Coord[i,2]<U_Bound && Fireflies_Coord[i,2]>L_Bound){
        LevelCurves[i,4]=Fireflies_Coord[i,1]
        LevelCurves[i,5]=Fireflies_Coord[i,2]
      }else{
        # count the proportion of fireflies that are within the bounds
        Success=Success-(1/Fireflies_Count)
      }
    }
    
    # prepare and print the plot
    LevelCurves=as.data.frame(LevelCurves)
    colnames(LevelCurves)=c("X","Y","Z","F_X","F_Y")
    
    v = ggplot(LevelCurves, aes(X, Y, z = Z)) + 
      stat_contour_filled(aes()) +
      # scale_fill_gradient(low = "blue",
      #                     high = "red",
      #                     guide = "colourbar",
      #                     aesthetics = "fill") +
      geom_point(aes(x=F_X, y=F_Y)) +
      ggtitle("Optimization results") +
      theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank())+
      xlim(-2.5,2.5)+ #hardcoded, because the functions used to test, have optima near the origin
      ylim(-2.5,2.5)
    
    print(v)
    print(Fireflies_Coord[which.max(Fireflies_Coord[,3]),])
    return(list(REPO,v))
  }


