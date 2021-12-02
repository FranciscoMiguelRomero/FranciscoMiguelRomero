datosnina<-c(66,88,58.5,170.5,112,65,34,27.5,169,148,35.78,74.33)
datosnino<-c(98,47.5,77,37.83,74,55,75.5,165,64)
datosnorm<-c(41,110.5,130.4,50,60,188.5,57,73,114,199,89.22)
datos<-c(66,88,41,98,47.5,58.5,110.5,130.4,77,50,60,37.83,170.5,188.5,74,
         112,65,34,57,55,73,75.5,27.5,165,169,148,64,35.78,74.33,114,199,89.22)
#####funcion que aparece en el exponente de la distribucion
fexp<-function(x,a,b,c){
  if (c==0){
    return(exp(-(x-a)/b))
  }
  return((1+c*(x-a)/b)^(-1/c))
}
############funcion de densidad
fdensidad<-function(x,a,b,c){
  return((1/b)*(fexp(x,a,b,c))^(c+1)*exp(-fexp(x,a,b,c)))
}
##############funcion de distribucion
fdistribucion<-function(x,a,b,c){
  return(exp(-fexp(x,a,b,c)))
}
##############Grafica PP
pp_plot = function(X, ag, bg,cg, confidence){
  n=length(X)
  X = fdistribucion(sort(X),ag,bg,cg)
  Y = seq(1/(n+1), 1, length.out = n)
  print(sum(abs(X-Y)))
  # puntos de la muestra
  plot(X, Y,xlim=c(0, 1),ylim=c(0, 1),
       main = "Gráfica PP para DGVE",
       xlab = "Probabilidades teóricas",
       ylab = "Probabilidades empíricas",
       pch = 19,
       cex = 0.5)
  # identidad
  abline(a = 0, b = 1, col = "red", lwd = 2)
  # bandas de confianza
  points(qbeta((1 - confidence)/2, 1:n, n + 1 - 1:n), Y,
         type = "l",
         lty = 2)
  points(qbeta((1 + confidence)/2, 1:n, n + 1 - 1:n), Y,
         type = "l",
         lty = 2)
}
##########################################################
##############Ajuste DGVE es la función que ajusta un modelo DGVE al parámetro 
#############datos, grafíca la verosimilitud perfil. No recominedo usarla para
#############datos generales  pues en cuanto falla algo la funcion se detiene
AjusteDGVE<-function(datosnina){
  #####Negativo de Logverosimilitud evaluada en el vector x, usada para la función
  ######optim. Se considera toda la información de la funcion de densidad
  logverdatosnina<-function(x){
    return(-sum(log(fdensidad(datosnina,x[1],x[2],x[3]))))
  }
  #####Logverosimilitud evaluada en a,b,c
  #####Se considera toda la información de la funcion de densidad
  logverdatosnina2<-function(a,b,c){
    return(sum(log(fdensidad(datosnina,a,b,c))))
  }
  ######a0 y b0 son parametros iniciales para optim
  a0<-median(datosnina)
  b0<-sd(datosnina)
  
  ######logverperfilnina es la logverosimilitud perfil evaluada en c.
  ###### la funcion lv dentro de ella es la funcion logverdatosnina con parametro 
  ####x=(x1,x2) y c fijo se define para poder minimizarla con optim, usando como
  #####parametros iniciales a0 y b0.
  logverperfilnina<-function(c){
    lv<-function(x){
      return(-sum(log(fdensidad(datosnina,x[1],x[2],c))))
    } ###lv es la funcion de logverosimilitud evaluada en x=(a,b) y c,
    ####queremos maximizar para a,b
    a0<-median(datosnina) ##estimador inicial de a
    b0<-sd(datosnina)##estimador inicial de b
    emv<-optim(c(a0,b0),lv)$par
    return(logverdatosnina2(emv[1],emv[2],c))
  }
  
  #####c0 son candidatos a considerar como valores iniciales del emv de c en la 
  ##### funcion optim, se toma el mayor que satisface la condicion
  ####(1+c0[i]*(y-a0)/b0)>0 para y todos los valores en datos. cinicial inicia en 0
  #### y se actualiza por el candidato, si ninguno de los propuestos satisface la condicion
  #### entonces regresa que no hay valor inicial disponible
  c0<-c(-0.8,-0.1,0.1,0.8)
  cinicial=0
  for (i in 1:4){
    if((sum((1+c0[i]*(datosnina-a0)/b0)<=0)==0)){
      cinicial=c0[i]
    }
  }
  if(cinicial==0){
    print("No hay valor inicial disponible")
    return()
  }
  ##### EMV encontrado con optim con los valores iniciales a0, b0 y cinicial
  emv<-optim(c(a0,b0,cinicial),logverdatosnina)$par
  print(emv)
  
  #####Verosimilitud Perfil Relativa evaluada en c
  Rpnina <- function(c){
    return(exp(logverperfilnina(c)-logverperfilnina(emv[3])))
  }
  ####vectorize porque aveces plot.function da problemas cuando la función tiene 
  ####un vector como exponente.
  h<-Vectorize(Rpnina)
  #####Teoricamente son los valores de c que satisfacen x_(n)=a-b/c
  #####y x_(1)=a-b/c para a y b los EMV, pero no sirven de mucho pues al cambiar c
  #### se deben cambiar los valores de a y b, solo nos interesa el cuadro para gráficar
  #### h de forma correcta mas adelante
  cmin<-emv[2]/(emv[1]-max(datosnina)) 
  cmax<-emv[2]/(emv[1]-min(datosnina))
  plot.function(h,
                from = emv[3],
                to = emv[3], lwd = 2.5,xlim=c(-1, 2),ylim=c(0, 1),
                col = "springgreen3",
                main = "Verosimilitud relativa perfil de c",
                ylab = "perfil",
                xlab = "Valores")
  #####nivel de verosimilitud
  niv=0.1465
  #######se toma un cmin inicila (EMV de c), lfemv es la logverosimilitud perfil de c
  #evaluada en el EMV, se toma una tolerancia de 0.01. El siguiente codigo lo que hace
  #es grafícar la verosimilitud perfil de c de forma puntual, empieza en cmin, grafíca
  ###(cmin,exp(lfmin-lfemv)), actualiza cmin=cmin-0.001 (Puede ser mas fina si quieren)
  ####, si (1+c0[i]*(y-ev[1])/emv[2])>0 se deja de cumplir entonces sale del while
  ####si no optimiza lv y guarda los EMV de a y b que serviran como valores iniciales
  #para el optim en el siguiente ciclo, actualiza lfmin=logverperfilnina(cmin)
  ### y repite el proceso hasta que la condicion exp(lfmin-lfemv)>tol deje de cumplirse
  # es decir la verosimilitud perfil relativa deje de ser mayor a tol. Durante los ciclos
  #tambien busca el extremo izquierdo del intervalo de verosimilitud de nivel niv, al cual
  #le llama c1
  cmin<-emv[3]
  lfemv=logverperfilnina(emv[3])
  tol=0.01
  lfmin=logverperfilnina(cmin)
  while(exp(lfmin-lfemv)>tol){
    #print(c(cmin,exp(lfmin-lfemv)))
    if(abs(exp(lfmin-lfemv)-niv)<0.01){
      c1<-cmin
    }
    points(cmin,exp(lfmin-lfemv),col="springgreen3",lwd=0.5)
    cmin=cmin-0.001
    lv<-function(x){
      return(-sum(log(fdensidad(datosnina,x[1],x[2],cmin))))
    } 
    if(sum((1+cmin*(datosnina-emv[1])/emv[2])<=0)>=1){
      break
    }
    emv<-optim(c(emv[1],emv[2]),lv)$par
    lfmin<-logverdatosnina2(emv[1],emv[2],cmin)
  }
  
  
  
  #######Lo mismo que todo lo anterior pero ahora el movimiento es "hacia la derecha"
  emv<-optim(c(a0,b0,cinicial),logverdatosnina)$par
  #cmax<-emv[2]/(emv[1]-min(datosnina))
  cmax=emv[3]
  lfmax=logverperfilnina(cmax)
  
  lv<-function(x){
  return(-sum(log(fdensidad(datosnina,x[1],x[2],cmax))))
  }
  while(exp(lfmax-lfemv)>tol){
    points(cmax,exp(lfmax-lfemv),col="springgreen3",lwd=0.5)
    cmax=cmax+0.001
    #print(cmax)
    if(abs(exp(lfmax-lfemv)-niv)<0.01){
      c2<-cmax
    }
    lv<-function(x){
      return(-sum(log(fdensidad(datosnina,x[1],x[2],cmax))))
    } 
    if(sum((1+cmax*(datosnina-emv[1])/emv[2])<=0)>=1){
      break
    }
    emv<-optim(c(emv[1],emv[2]),lv)$par
    lfmax<-logverdatosnina2(emv[1],emv[2],cmax)
  }
  #####Finalmente Gráfica el intervalo de verosimilitud sobre la gráfica
  #y la linea vertical del EMV de c
  segments(x0 = c1, y0 = niv, x1 = c2, y1 = niv, col = "blue",
           lwd = 2.5)
  emv<-optim(c(a0,b0,cinicial),logverdatosnina)$par
  abline(v=emv[3],lty=2,col="red")
  ####Gráfica PP
  pp_plot(datosnina,emv[1],emv[2],emv[3],0.95)}
####Ajuste DGVE aplicada a datos
AjusteDGVE(datos)

