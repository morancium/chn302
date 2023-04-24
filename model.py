# INPUT Style:
# Feed -> F
# Feed components -> Zc
# Feed Temperature -> Tf
import argparse
from sympy import *
from utils import *
import json
import numpy as np
import decimal
decimal.getcontext().prec = 100
with open('data.json', 'r') as f:
    data = json.load(f)


global Vf,Lf,v_bar,l_bar,Vj,Lj
v_bar=[]
l_bar=[]
Vj=[]
# Lj=np.zeros(No_trays)

ap = argparse.ArgumentParser()
ap.add_argument("-Z", "--Feed_conc",help="concentration of feed Component",nargs="+",type=float)
ap.add_argument("-f", "--feed_location",help="Where feed is insterted",type=int)
ap.add_argument("-F", "--feed_flowrate",help="feed flowrate",type=float)
ap.add_argument("-Tf", "--Feed_temperature",help="Temperature of the feed",type=float)
ap.add_argument("-P", "--Pressure",help="Pressure of the column",type=float)
ap.add_argument("-D", "--Distilate",help="Distilate flowrate",type=float)
ap.add_argument("-N", "--No_trays",help="Number of trays",type=int)
ap.add_argument("-R", "--Reflux",help="Reflux Ratio",type=float)
ap.add_argument("-Type", "--Type",help="Type of condenser. 0 for partial and 1 for total",type=int)

args = ap.parse_args()
Feed_conc=args.Feed_conc
Feed_temperature=args.Feed_temperature
Pressure=args.Pressure
feed_flowrate=args.feed_flowrate
feed_location=args.feed_location
Distilate=args.Distilate
No_trays=args.No_trays
Type=args.Type
Reflux=args.Reflux
nc=len(Feed_conc)
component=["benzene","toulene","o-xylene"]
Lj=np.zeros(No_trays+2)
# Vj=np.random.uniform(low=feed_flowrate/2, high=feed_flowrate, size=(No_trays+2))
Vj=np.ones((No_trays+2))*(Distilate+10)
Tj=np.ones((No_trays+2))*Feed_temperature

# print(Tj)
# STEP 1
def Flash_calc():
    def f(x):
        # fun=0
        fun= decimal.Decimal(0)
        for i in range(0,nc):
            fun+=Feed_conc[i]/(1+x*(Ki(data["components"][component[i]]["antoine"],Feed_temperature,Pressure)-1))
            # print(Ki(data["components"][component[i]]["antoine"],Feed_temperature,Pressure))
            # fun+=x**i-1
        # print("funf",fun)
        return fun-1
    xans=newton_raphson(f,1)
    # print(xans)
    Vf=xans*feed_flowrate
    Lf=feed_flowrate-Vf
    yi=np.zeros(nc)
    for i in range(0,nc):
        xi=Feed_conc[i]/(1+(Ki(data["components"][component[i]]["antoine"],Feed_temperature,Pressure)-1)*xans)
        yi[i]=Ki(data["components"][component[i]]["antoine"],Feed_temperature,Pressure)*xi
        v_bar.append(Vf*yi[i])
        l_bar.append(Lf*xi)
        pass
    # print(yi)
    return Vf,Lf,v_bar,l_bar,yi
# print(Flash_calc())


# STEP 2
def step2():
    Vf,_,_,_,_=Flash_calc()
    f=feed_location

    for i in range(0,f-1):
        Lj[i]=Vj[i+1]-Distilate
    Lj[f-1]=Vj[f]+Vf-Distilate
    for i in range(f,No_trays+1):
        Lj[i]=Vj[i+1]+feed_flowrate-Distilate
    Lj[-1]=feed_flowrate-Distilate
    Vj[0]=Distilate
    return Lj


# STEP 3
def step3():
    # print(Tj)
    Aji=np.zeros((No_trays+2,nc))
    Lj=step2()
    for j in range (0,No_trays+2):
        for i in range (0,nc):
            Aji[j][i]=Lj[j]/(Vj[j]*Ki(data["components"][component[i]]["antoine"],Tj[j],Pressure))
            # print(Aji[j][i])
    # print(Aji)
    return Aji


# STEP 4
def step4():
    Aji_bar=[]
    gammaji=[]
    lji=np.zeros((No_trays+2,nc))
    Aji=step3()
    Li=np.zeros((nc,No_trays+2))
    ro=np.zeros((nc,No_trays+2))
    lower=np.zeros((nc,No_trays+1))
    for i in range(nc):
        for j in range(No_trays+2):
            ro[i][j]=-1-Aji[j][i]
            if(j!=No_trays+1):
                lower[i][j]=Aji[j][i]
        Li[i][feed_location-1]=v_bar[i]
        Li[i][feed_location]=l_bar[i]
        Aji_bar.append(np.diag(ro[i],0)+np.diag(lower[i],-1)+np.diag(np.ones_like(lower[i]),1))
        gammaji.append(np.linalg.solve(Aji_bar[i],-Li[i]))
    Aji_bar=np.array(Aji_bar)
    gammaji=np.array(gammaji).transpose()
    # print(gammaji)
    for j in range(No_trays+2):
        for i in range(nc):
            lji[j][i]=Aji[j][i]*gammaji[j][i]
    # print(lji)
    # print(gammaji.shape)
    # print(Li)
    return lji, gammaji
# step4()

# STEP 5
def step5():
    Aji=step3()
    bi_calc=np.zeros((nc))
    di_calc=np.zeros((nc))
    di_corr=np.zeros((nc))
    xji=np.zeros((No_trays+2,nc))
    yji=np.zeros((No_trays+2,nc))
    lji,gammaji=step4()
    for i in range(nc):
        bi_calc[i]=lji[No_trays+1][i]
        if Type:
            di_calc[i]=lji[0][i]/(Reflux)
        else:
            di_calc[i]=lji[0][i]/Aji[0][i]
    # print(bi_calc)
    # print(di_calc)
    def g(x):
        g_phi=0
        for i in range(nc):
            g_phi+=feed_flowrate*Feed_conc[i]/(1+x*bi_calc[i]/di_calc[i])
        g_phi=g_phi-Distilate
        # print(g_phi)
        return g_phi
    # print(bi_calc)
    phi=newton_raphson(g,0)
    for i in range(nc):
        di_corr[i]=feed_flowrate*Feed_conc[i]/(1+phi*bi_calc[i]/di_calc[i])
    
    for j in range (No_trays+2):
        denomx=0
        denomy=0
        for i in range(nc): 
            # denomx+=lji[j][i]*di_corr[i]/gammaji[0][i]
            # denomy+=gammaji[j][i]*di_corr[i]/gammaji[0][i]
            denomx+=lji[j][i]*di_corr[i]/di_calc[i]
            denomy+=gammaji[j][i]*di_corr[i]/di_calc[i]
            
        for i in range(nc):
            xji[j][i]=(lji[j][i]*di_corr[i]/di_calc[i])/denomx
            yji[j][i]=(gammaji[j][i]*di_corr[i]/di_calc[i])/denomy
    # print(xji)
    return xji,yji,di_corr

# STEP 6
def step6():
    # Kb method
    xji,_,_=step5()
    kji=np.zeros((No_trays+2,nc))
    kbj=np.zeros((No_trays+2))
    alphaji=np.zeros((No_trays+2,nc))
    Tj_new=np.zeros((No_trays+2))
    for j in range(No_trays+2):
        for i in range(nc):
            kji[j][i]=Ki(data["components"][component[i]]["antoine"],Tj[j],Pressure)
    for j in range(No_trays+2):
        kbj[j]=kji[j][0]
    for j in range(No_trays+2):
        for i in range(nc):
            alphaji[j][i]=kji[j][i]/kbj[j]
    # print(alphaji)
    # print(kbj)
    a,b,c=data["components"][component[0]]["antoine"]
    for j in range(No_trays+2):
        denom=0
        for i in range(nc):
            denom+=alphaji[j][i]*xji[j][i]
        Tj_new[j]=b/(a-log(Pressure/denom,10))-c
    return Tj_new

# STEP 7
def step7():
    xji,yji,di_corr=step5()
    Vj_new=np.zeros(No_trays+2)
    Vj_new[0]=Distilate
    Vj_new[1]=Reflux*Distilate+Distilate
    Vf,_,_,_,yi=Flash_calc()
    summ1=0
    summ2=0
    h0=0
    HD=0
    hB=0
    HF=0
    for i in range (nc):
        summ1+=((Hji(data["components"][component[i]]["Cp_liquid"],Tj[1])/1000+(data["components"][component[i]]["Hv"])))*xji[0][i]
        h0+=(Hji(data["components"][component[i]]["Cp_liquid"],Tj[0])/1000)*xji[0][i]
        summ2+=((Hji(data["components"][component[i]]["Cp_liquid"],Tj[1])/1000+(data["components"][component[i]]["Hv"])))*di_corr[i]/Distilate
        if Type:
            HD+=(Hji(data["components"][component[i]]["Cp_liquid"],Tj[0])/1000)*di_corr[i]/Distilate
        else:
            HD+=(Hji(data["components"][component[i]]["Cp_liquid"],Tj[0])/1000+(data["components"][component[i]]["Hv"]))*di_corr[i]/Distilate
        hB+=Hji(data["components"][component[i]]["Cp_liquid"],Tj[-1])*xji[-1][i]/1000
        HF+=yi[i]*(Hji(data["components"][component[i]]["Cp_liquid"],Tj[0])/1000+data["components"][component[i]]["Hv"])
        
    Qc=Reflux*Distilate*(summ1-h0)+Distilate*(summ2-HD)
    QR=(feed_flowrate-Distilate)*hB+Distilate*HD+Qc-feed_flowrate*HF
    
    for j in range (1,feed_location-1):
        summ_Num=0
        summ_denom=0
        hj=0
        for i in range (nc):
            summ_Num+=(Hji(data["components"][component[i]]["Cp_liquid"],Tj[j+1])/1000+(data["components"][component[i]]["Hv"]))*di_corr[i]/Distilate
            summ_denom+=(Hji(data["components"][component[i]]["Cp_liquid"],Tj[j+1])/1000+(data["components"][component[i]]["Hv"]))*xji[j][i]
            hj+=Hji(data["components"][component[i]]["Cp_liquid"],Tj[j])*xji[j][i]/1000
        Vj_new[j+1]=(Distilate*(HD-summ_Num)+Qc)/summ_denom-hj[j]+Distilate
        
    summ_feedtray_d=0
    summ_feedtray_f=0
    summ_feedtray_f_liq=0
    hf_1=0
    for i in range (nc):
        summ_feedtray_d+=(Hji(data["components"][component[i]]["Cp_liquid"],Tj[feed_location])/1000+(data["components"][component[i]]["Hv"]))*di_corr[i]/Distilate
        summ_feedtray_f+=(Hji(data["components"][component[i]]["Cp_liquid"],Tj[feed_location])/1000+(data["components"][component[i]]["Hv"]))*yi[i]
        summ_feedtray_f_liq+=(Hji(data["components"][component[i]]["Cp_liquid"],Tj[feed_location])/1000+(data["components"][component[i]]["Hv"]))*xji[feed_location-1][i]
        hf_1+=Hji(data["components"][component[i]]["Cp_liquid"],Tj[feed_location-1])*xji[feed_location-1][i]/1000
    Vj_new[feed_location]=(Distilate*(HD-summ_feedtray_d)+Vf*(summ_feedtray_f-HF)+Qc)/(summ_feedtray_f_liq-hf_1)+Distilate-Vf
        
    for j in range(feed_location,No_trays+1):
        summ_Num1=0
        summ_denom1=0
        Hj=0
        for i in range(nc):
            summ_Num1+=Hji(data["components"][component[i]]["Cp_liquid"],Tj[j])*xji[-1][i]/1000
            summ_denom1+=Hji(data["components"][component[i]]["Cp_liquid"],Tj[j])*yji[j+1][i]/1000
            Hj+=(Hji(data["components"][component[i]]["Cp_liquid"],Tj[j+1])/1000+(data["components"][component[i]]["Hv"]))*yji[j+1][i]
        Vj_new[j+1]=((feed_flowrate-Distilate)*(summ_Num1-hB)+QR)/(Hj-summ_denom1)
        
    return Vj_new
# Flash_calc()7
# print(step7())
# print(Feed_conc)

Vj_new=step7()
Tj_new=step6()
tolT = np.linalg.norm(Tj-Tj_new)
tolV = np.linalg.norm(Vj-Vj_new)
while (tolT>1e-6 and tolV>1e-6 ):
    # tolT = np.linalg.norm(Tj-Tj_new)
    tolV = np.linalg.norm(Vj-Vj_new)
    Tj=Tj_new
    Vj=Vj_new
    Vj_new=step7()
    # Tj_new=step6()
    print(Vj_new)
    # print(Tj_new)
# main()