# Functions to test optimization
# Note - where necessary: we replace f(x) with -f(x) 

OPT_Sphere = function(x,y){         # (0,0)
  -(x^2+y^2)
}

OPT_Rosenbrock = function(x,y){    # (1,1)
  -(100*(y-x^2)^2+(1-x)^2)
}

OPT_McCormick = function(x,y){       # (-0.54719,-1.54719)
  -(sin(x+y)+(x-y)^2 -1.5*x+2.5*y+1)
}

OPT_Xin_She_Yang = function(x,y){    # (0,0)
  ((abs(x)+abs(y))*exp(-(x^2+y^2)))
}