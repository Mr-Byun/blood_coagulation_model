using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

using("deSolve", "MESS", "gtools", "randomcoloR")

# library(deSolve)
# library(MESS) #AUC
# library(gtools) #fold change

Blood <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dTF <- 1/blood*(- (k2*TF*VII-k1*TF_VII) - (k4*TF*VIIa-k3*TF_VIIa))
    dTF_VII <- 1/blood*( (k2*TF*VII-k1*TF_VII) )
    dVII <- 1/blood*(- (k2*TF*VII-k1*TF_VII) - (k5*TF_VIIa*VII) - (k6*Xa*VII) - (k7*IIa*VII))
    dTF_VIIa <- 1/blood*((k4*TF*VIIa-k3*TF_VIIa) - (k9*TF_VIIa*X-k8*TF_VIIa_X) - (k12*TF_VIIa*Xa-k11*TF_VIIa_Xa) - (k14*TF_VIIa*IX-k13*TF_VIIa_IX) + (k15*TF_VIIa_IX) - (k37*TF_VIIa*Xa_TFPI) - (k42*TF_VIIa*ATIII))
    dVIIa <- 1/blood*(-(k4*TF*VIIa-k3*TF_VIIa) + (k5*TF_VIIa*VII) + (k6*Xa*VII) + (k7*IIa*VII) - (k61*VIIa*X-k62*VIIa_X) + (k63*VIIa_X) + (n1*VIIa_X*XIa))
    dXa <- 1/blood*(-(k12*TF_VIIa*Xa-k11*TF_VIIa_Xa) + (k22*IXa_VIIIa_X) - (k28*Xa*Va-k27*Xa_Va) - (k34*Xa*TFPI-k33*Xa_TFPI) - (k38*Xa*ATIII) + (k60*IXa_X) + (k63*VIIa_X) - (k83*Xa*VIII-k84*Xa_VIII) + (k85*Xa_VIII) + (n1*VIIa_X*XIa))
    dIIa <- 1/blood*( (k16*Xa*II) + (k32*mIIa*Xa_Va) - (k41*IIa*ATIII) - (k64*IIa*Tmod-k65*IIa_Tmod) )
    dTF_VIIa_X <- 1/blood*( (k9*TF_VIIa*X-k8*TF_VIIa_X) - (k10*TF_VIIa_X) )
    dX <- 1/blood*( - (k9*TF_VIIa*X-k8*TF_VIIa_X) - (k21*IXa_VIIIa*X-k20*IXa_VIIIa_X) + (k25*IXa_VIIIa_X) - (k58*IXa*X-k59*IXa_X) - (k61*VIIa*X-k62*VIIa_X) + (k82*APC_IXa_VIIIa_X) )
    dTF_VIIa_Xa <- 1/blood*( (k10*TF_VIIa_X) + (k12*TF_VIIa*Xa-k11*TF_VIIa_Xa) - (k36*TF_VIIa_Xa*TFPI-k35*TF_VIIa_Xa_TFPI) )
    dIX <- 1/blood*( - (k14*TF_VIIa*IX-k13*TF_VIIa_IX) - (k55*XIa*IX-k56*XIa_IX) )
    dTF_VIIa_IX <- 1/blood*( (k14*TF_VIIa*IX-k13*TF_VIIa_IX) - (k15*TF_VIIa_IX) )
    dIXa <- 1/blood*( (k15*TF_VIIa_IX) - (k19*IXa*VIIIa-k18*IXa_VIIIa) + (k25*IXa_VIIIa_X) + (k25_1*IXa_VIIIa) - (k40*IXa*ATIII) + (k57*XIa_IX) - (k58*IXa*X-k59*IXa_X) + (k60*IXa_X) + (k79*APC_IXa_VIIIa) + (k82*APC_IXa_VIIIa_X) )
    dII <- 1/blood*( - (k16*Xa*II) - (k30*Xa_Va*II-k29*Xa_Va_II) )
    dVIII <- 1/blood*( - (k17*IIa*VIII) - (k83*Xa*VIII-k84*Xa_VIII) )
    dVIIIa <- 1/blood*( (k17*IIa*VIII) - (k19*IXa*VIIIa-k18*IXa_VIIIa) - (k24*VIIIa-k23*VIIIa1_L*VIIIa2) - (k73*APC*VIIIa-k74*APC_VIIIa) + (k85*Xa_VIII) )
    dIXa_VIIIa <- 1/blood*( (k19*IXa*VIIIa-k18*IXa_VIIIa) - (k21*IXa_VIIIa*X-k20*IXa_VIIIa_X) + (k22*IXa_VIIIa_X) - (k25_1*IXa_VIIIa) - (k77*APC*IXa_VIIIa-k78*APC_IXa_VIIIa) )
    dIXa_VIIIa_X <- 1/blood*( (k21*IXa_VIIIa*X-k20*IXa_VIIIa_X) - (k22*IXa_VIIIa_X) - (k25*IXa_VIIIa_X) - (k80*APC*IXa_VIIIa_X-k81*APC_IXa_VIIIa_X) )
    dVIIIa1_L <- 1/blood*( (k24*VIIIa-k23*VIIIa1_L*VIIIa2) + (k25*IXa_VIIIa_X) + (k25_1*IXa_VIIIa) )
    dVIIIa2 <- 1/blood*( (k24*VIIIa-k23*VIIIa1_L*VIIIa2) + (k25*IXa_VIIIa_X) + (k25_1*IXa_VIIIa) )
    dV <- 1/blood*( - (k26*IIa*V) )
    dVa <- 1/blood*( (k26*IIa*V) - (k28*Xa*Va-k27*Xa_Va) - (k71*APC*Va-k72*APC_Va) )
    dXa_Va <- 1/blood*( (k28*Xa*Va-k27*Xa_Va) - (k30*Xa_Va*II-k29*Xa_Va_II) + (k31*Xa_Va_II) )
    dXa_Va_II <- 1/blood*( (k30*Xa_Va*II-k29*Xa_Va_II) - (k31*Xa_Va_II) )
    dmIIa <- 1/blood*( (k31*Xa_Va_II) - (k32*mIIa*Xa_Va) - (k39*mIIa*ATIII) )
    dTFPI <- 1/blood*( - (k34*Xa*TFPI-k33*Xa_TFPI) - (k36*TF_VIIa_Xa*TFPI-k35*TF_VIIa_Xa_TFPI) )
    dXa_TFPI <- 1/blood*( (k34*Xa*TFPI-k33*Xa_TFPI) - (k37*TF_VIIa*Xa_TFPI) )
    dTF_VIIa_Xa_TFPI <- 1/blood*( (k36*TF_VIIa_Xa*TFPI-k35*TF_VIIa_Xa_TFPI) + (k37*TF_VIIa*Xa_TFPI) )
    dATIII <- 1/blood*( - (k38*Xa*ATIII) - (k39*mIIa*ATIII) - (k40*IXa*ATIII) - (k41*IIa*ATIII) - (k42*TF_VIIa*ATIII) - (k50*XIIa*ATIII) - (k54*XIa*ATIII) )
    dXa_ATIII <- 1/blood*( (k38*Xa*ATIII) )
    dmIIa_ATIII <- 1/blood*( (k39*mIIa*ATIII) )
    dIXa_ATIII <- 1/blood*( (k40*IXa*ATIII) )
    dIIa_ATIII <- 1/blood*( (k41*IIa*ATIII) )
    dTF_VIIa_ATIII <- 1/blood*( (k42*TF_VIIa*ATIII) )
    dXII <- 1/blood*( - (k43*CA*XII) - (k47*XII*K-k48*XII_K) )
    dXIIa <- 1/blood*( (k43*CA*XII) - (k44*XIIa*PK-k45*XIIa_PK) + (k46*XIIa_PK) + (k49*XII_K) - (k50*XIIa*ATIII) - (k51*XIIa*XI-k52*XIIa_XI) + (k53*XIIa_XI) )
    dPK <- 1/blood*( - (k44*XIIa*PK-k45*XIIa_PK) )
    dXIIa_PK <- 1/blood*( (k44*XIIa*PK-k45*XIIa_PK) - (k46*XIIa_PK) )
    dK <- 1/blood*( (k46*XIIa_PK) - (k47*XII*K-k48*XII_K) + (k49*XII_K) )
    dXII_K <- 1/blood*( (k47*XII*K-k48*XII_K) - (k49*XII_K) )
    dXIIa_ATIII <- 1/blood*( (k50*XIIa*ATIII) )
    dXI <- 1/blood*( - (k51*XIIa*XI-k52*XIIa_XI) )
    dXIa <- 1/blood*( (k53*XIIa_XI) - (k54*XIa*ATIII) - (k55*XIa*IX-k56*XIa_IX) + (k57*XIa_IX) )
    dXIIa_XI <- 1/blood*( (k51*XIIa*XI-k52*XIIa_XI) - (k53*XIIa_XI) )
    dXIa_ATIII <- 1/blood*( (k54*XIa*ATIII) )
    dXIa_IX <- 1/blood*( (k55*XIa*IX-k56*XIa_IX) - (k57*XIa_IX) )
    dIXa_X <- 1/blood*( (k58*IXa*X-k59*IXa_X) - (k60*IXa_X) )
    dVIIa_X <- 1/blood*( (k61*VIIa*X-k62*VIIa_X) - (k63*VIIa_X) - (n1*VIIa_X*XIa) )
    dTmod <- 1/blood*( - (k64*IIa*Tmod-k65*IIa_Tmod) )
    dIIa_Tmod <- 1/blood*( (k64*IIa*Tmod-k65*IIa_Tmod) - (k66*IIa_Tmod*PC-k67*IIa_Tmod_PC) + (k68*IIa_Tmod_PC) )
    dPC <- 1/blood*( - (k66*IIa_Tmod*PC-k67*IIa_Tmod_PC) )
    dIIa_Tmod_PC <- 1/blood*( (k66*IIa_Tmod*PC-k67*IIa_Tmod_PC) - (k68*IIa_Tmod_PC) )
    dAPC <- 1/blood*( (k68*IIa_Tmod_PC) - (k71*APC*Va-k72*APC_Va) - (k73*APC*VIIIa-k74*APC_VIIIa) + (k75*APC_Va) + (k76*APC_VIIIa) - (k77*APC*IXa_VIIIa-k78*APC_IXa_VIIIa) + (k79*APC_IXa_VIIIa) - (k80*APC*IXa_VIIIa_X-k81*APC_IXa_VIIIa_X) + (k82*APC_IXa_VIIIa_X) )
    dFg <- 1/blood*( - (k69*IIa*Fg/(k70+Fg)) )
    dF <- 1/blood*( (k69*IIa*Fg/(k70+Fg)) )
    dAPC_Va <- 1/blood*( (k71*APC*Va-k72*APC_Va) - (k75*APC_Va) )
    dAPC_VIIIa <- 1/blood*( (k73*APC*VIIIa-k74*APC_VIIIa) - (k76*APC_VIIIa) )
    dVa_deg <- 1/blood*( (k75*APC_Va) )
    dVIIIa_deg <- 1/blood*( (k76*APC_VIIIa) + (k79*APC_IXa_VIIIa) + (k82*APC_IXa_VIIIa_X) )
    dAPC_IXa_VIIIa <- 1/blood*( (k77*APC*IXa_VIIIa-k78*APC_IXa_VIIIa) - (k79*APC_IXa_VIIIa) )
    dAPC_IXa_VIIIa_X <- 1/blood*( (k80*APC*IXa_VIIIa_X-k81*APC_IXa_VIIIa_X) - (k82*APC_IXa_VIIIa_X) )
    dXa_VIII <- 1/blood*( (k83*Xa*VIII-k84*Xa_VIII) - (k85*Xa_VIII) )
    dF12 <- 1/blood*( (k16*Xa*II) + (k31*Xa_Va_II) - (k86*F12) )
    dF12_deg <- 1/blood*( (k86*F12) )
    
    
    list(c(dTF,dTF_VII,dVII,dTF_VIIa,dVIIa,dXa,dIIa,dTF_VIIa_X,dX,dTF_VIIa_Xa,dIX,dTF_VIIa_IX,dIXa,dII,dVIII,dVIIIa,dIXa_VIIIa,dIXa_VIIIa_X,dVIIIa1_L,dVIIIa2,dV,dVa,dXa_Va,dXa_Va_II,dmIIa,dTFPI,dXa_TFPI,dTF_VIIa_Xa_TFPI,dATIII,dXa_ATIII,dmIIa_ATIII,dIXa_ATIII,dIIa_ATIII,dTF_VIIa_ATIII,dXII,dXIIa,dPK,dXIIa_PK,dK,dXII_K,dXIIa_ATIII,dXI,dXIa,dXIIa_XI,dXIa_ATIII,dXIa_IX,dIXa_X,dVIIa_X,dTmod,dIIa_Tmod,dPC,dIIa_Tmod_PC,dAPC,dFg,dF,dAPC_Va,dAPC_VIIIa,dVa_deg,dVIIIa_deg,dAPC_IXa_VIIIa,dAPC_IXa_VIIIa_X,dXa_VIII,dF12,dF12_deg))
  })
}

parameters <- c(CA=0,k50=0.0000000108,k40=0.0000002817,k42=0.00000032905,k54=0.00000048,k38=0.0000010557,k39=0.00000355,k41=0.0000039177,k43=0.0000085153,k7=0.000011527,k63=0.000012944,k23=0.000033,k16=0.000037641,k35=0.00010062,k86=0.00010959,k33=0.00018016,k5=0.00033895,k53=0.00035,k26=0.00040234,k34=0.00045,k60=0.0011427,k25=0.0013358,k20=0.0013766,k25_1=0.0013946,TF_ss=0.001743,k3=0.0019496,k2=0.0041242,k1=0.004334,k18=0.0050725,k27=0.006,k24=0.009,k6=0.0099754,k14=0.010569,k82=0.013215,k37=0.025387,k79=0.029888,k76=0.03236,k12=0.033,k85=0.0345,k9=0.036246,k21=0.048795,n1=0.05,k65=0.050085,k61=0.059664,k4=0.07568,k77=0.077518,k73=0.079343,k44=0.091853,VIIa_ss=0.1,k30=0.10002,k71=0.10436,k66=0.10524,k64=0.11573,k19=0.1175,k80=0.12893,k28=0.12914,k55=0.13082,k58=0.13304,k51=0.13632,k47=0.15,k83=0.15,k68=0.18579,k32=0.21872,k75=0.24028,k36=0.25639,Tmod_ss=0.5,VIII_ss=0.7,blood=1,k8=1.38,k81=1.4002,k17=1.4489,k74=1.9895,k72=2.0649,k78=2.087,k15=2.3887,TFPI_ss=2.5,k84=3.15,k49=4.1832,k10=8.9988,k11=9.5,VII_ss=10,k57=10.566,k52=18.705,k67=19.338,V_ss=20,k13=20.671,k31=29.479,k56=30.668,XI_ss=31,k22=42.714,k46=47.928,k48=57.614,PC_ss=60,k59=83.207,k62=84.66,IX_ss=90,k69=90.212,k29=149.92,X_ss=160,XII_ss=340,PK_ss=450,II_ss=1400,ATIII_ss=3400,k45=3600,Fg_ss=9000,k70=10740)

#Initial value of VIIa and Xa
#VIIa = 0.058333
#VIII = 0.40833
#Xa = 0

state <- c(TF=0.0010167,TF_VII=0,VII=5.8333,TF_VIIa=0,VIIa=0.058333,Xa=0,IIa=0,TF_VIIa_X=0,X=93.333,TF_VIIa_Xa=0,IX=52.5,TF_VIIa_IX=0,IXa=0,II=816.66,VIII=0.40833,VIIIa=0,IXa_VIIIa=0,IXa_VIIIa_X=0,VIIIa1_L=0,VIIIa2=0,V=11.667,Va=0,Xa_Va=0,Xa_Va_II=0,mIIa=0,TFPI=1.4583,Xa_TFPI=0,TF_VIIa_Xa_TFPI=0,ATIII=1983.3,Xa_ATIII=0,mIIa_ATIII=0,IXa_ATIII=0,IIa_ATIII=0,TF_VIIa_ATIII=0,XII=198.33,XIIa=0,PK=262.5,XIIa_PK=0,K=0,XII_K=0,XIIa_ATIII=0,XI=18.083,XIa=0,XIIa_XI=0,XIa_ATIII=0,XIa_IX=0,IXa_X=0,VIIa_X=0,Tmod=0.29166,IIa_Tmod=0,PC=35,IIa_Tmod_PC=0,APC=0,Fg=5250,F=0,APC_Va=0,APC_VIIIa=0,Va_deg=0,VIIIa_deg=0,APC_IXa_VIIIa=0,APC_IXa_VIIIa_X=0,Xa_VIII=0,F12=0,F12_deg=0)

run = 3600
times <- seq(1, run, by = 1)



draw <- function(var, conc)
{
  f_state <- state
  if(missing(var)) var <- "VIIa"

  col_vector <- distinctColorPalette(length(conc), runTsne = T)
  
  plot(NULL, lwd = 2, ylab = "Thrombin (nM)", xlab = "time (sec)", xlim = c(1, run), ylim = c(0,400), type = "l")
  for(i in conc)
  {
    f_state[var] = i
    out <- as.data.frame(ode(y = f_state, times = times, func = Blood, parms = parameters))
    
    colour <- col_vector[which(conc == i)]
    IIa_12mIIa <- out$IIa + 1.2 * out$mIIa
    lines(times, IIa_12mIIa, lwd = 2, col = colour)
  }
  
  legend("topright", legend = conc, col = col_vector, lty = 1, cex = 0.7)
  title(main = paste0("Thrombin Activation per concentration of ", var, " in normal person"))
  
  f_state["VIII"] = 0
  plot(NULL, lwd = 2, ylab = "Thrombin (nM)", xlab = "time (sec)", xlim = c(1, run), ylim = c(0,400), type = "l")
  for(i in conc)
  {
    f_state[var] = i
    out <- as.data.frame(ode(y = f_state, times = times, func = Blood, parms = parameters))
    
    colour <- col_vector[which(conc == i)]
    IIa_12mIIa <- out$IIa + 1.2 * out$mIIa
    lines(times, IIa_12mIIa, lwd = 2, col = colour)
  }
  
  legend("topright", legend = conc, col = col_vector, lty = 1, cex = 0.7)
  title(main = paste0("Thrombin Activation per concentration of ", var, " in VIII-deficient patient"))
}