M=500
n=100
xweibull<-seq(1,M)
for (i in 1:M){
  xweibull[i]<-max(rweibull(n,shape=40,scale=3))
}
xweibull
logverxweibull<-function(x){
  return(-sum(log_(fdensidad(xweibull,x[1],x[2],x[3]))))
}
logverxweibull2<-function(a,b,c){
  return(sum(log_(fdensidad(xweibull,a,b,c))))
}
a0<-median(xweibull)
b0<-sd(xweibull)

logverperfilnina<-function(c){
  lv<-function(x){
    return(-sum(log_(fdensidad(xweibull,x[1],x[2],c))))
  } ###lv es la funcion de logverosimilitud evaluada en x=(a,b) y c,
  ####queremos maximizar para a,b
  a0<-median(xweibull) ##estimador inicial de a
  b0<-sd(xweibull)##estimador inicial de b
  emv<-optim(c(a0,b0),lv)$par
  return(logverxweibull2(emv[1],emv[2],c))
}
c0<-c(-0.8,-0.1,0.1,0.8)
cinicial=0
for (i in 1:4){
  if((sum((1+c0[i]*(xweibull-a0)/b0)<=0)==0)){
    cinicial=c0[i]
  }
}
if(cinicial==0){
  print("No hay valor inicial disponible")
  return()
}
sum(log_(fdensidad(xweibull,a0,b0,cinicial)))
emv<-optim(c(a0,b0,cinicial),logverxweibull)$par
print(emv)
Rpnina <- function(c){
  return(exp(logverperfilnina(c)-logverperfilnina(emv[3])))
}
h<-Vectorize(Rpnina)
cmin<-emv[2]/(emv[1]-max(xweibull)) 
cmax<-emv[2]/(emv[1]-min(xweibull))
plot.function(h,
              from = cmin+0.1,
              to = cmax-0.1, lwd = 2.5,xlim=c(-1, 2),ylim=c(0, 1),
              col = "springgreen3",
              main = "Verosimilitud relativa perfil de c",
              ylab = "perfil",
              xlab = "Valores")
niv=0.1465
cmin<-emv[3]
lfemv=logverperfilnina(emv[3])
tol=0.01
lfmin=logverperfilnina(cmin)
lv<-function(x){
  return(-sum(log_(fdensidad(xweibull,x[1],x[2],cmin))))
}
print("Aqui")
while(exp(lfmin-lfemv)>tol){
  #print(c(cmin,exp(lfmin-lfemv)))
  if(abs(exp(lfmin-lfemv)-niv)<0.01){
    c1<-cmin
  }
  points(cmin,exp(lfmin-lfemv),col="springgreen3",lwd=0.5)
  cmin=cmin-0.001
  lv<-function(x){
    return(-sum(log_(fdensidad(xweibull,x[1],x[2],cmin))))
  } 
  if(sum((1+cmin*(xweibull-emv[1])/emv[2])<=0)>=1){
    break
  }
  emv<-optim(c(emv[1],emv[2]),lv)$par
  lfmin<-logverxweibull2(emv[1],emv[2],cmin)
}
print("Aqui")
emv<-optim(c(a0,b0,cinicial),logverxweibull)$par
#cmax<-emv[2]/(emv[1]-min(xweibull))
cmax=emv[3]
lfmax=logverperfilnina(cmax)
lv<-function(x){
  return(-sum(log_(fdensidad(xweibull,x[1],x[2],cmax))))
}
while(exp(lfmax-lfemv)>tol){
  points(cmax,exp(lfmax-lfemv),col="springgreen3",lwd=0.5)
  cmax=cmax+0.001
  #print(cmax)
  if(abs(exp(lfmax-lfemv)-niv)<0.01){
    c2<-cmax
  }
  lv<-function(x){
    return(-sum(log_(fdensidad(xweibull,x[1],x[2],cmax))))
  } 
  if(sum((1+cmax*(xweibull-emv[1])/emv[2])<=0)>=1){
    break
  }
  emv<-optim(c(emv[1],emv[2]),lv)$par
  lfmax<-logverxweibull2(emv[1],emv[2],cmax)
}
#uno = uniroot(ints, c(cmin+0.1,0.5))$root
#dos = uniroot(ints, c(1,2))$root
segments(x0 = c1, y0 = niv, x1 = c2, y1 = niv, col = "blue",
         lwd = 2.5)
emv<-optim(c(a0,b0,cinicial),logverxweibull)$par
abline(v=emv[3],lty=2,col="red")
pp_plot(xweibull,emv[1],emv[2],emv[3],0.95)
