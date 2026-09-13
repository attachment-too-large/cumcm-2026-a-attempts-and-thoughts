# -*- coding: utf-8 -*-
import numpy as np, openpyxl, io, sys
from scipy.optimize import curve_fit

wb = openpyxl.load_workbook(r"C:\Users\qing1\Desktop\数A\附件\附件1.xlsx", data_only=True)
ws = wb.active
rows = list(ws.iter_rows(min_row=2, values_only=True))
t  = np.array([float(r[0]) for r in rows])
Ta = np.array([float(r[1]) for r in rows])
Ca = np.array([float(r[2]) for r in rows])
out=[]
out.append("N=%d  t: %.1f..%.1f  Ta0=%.4f TaEnd=%.4f  Ca0=%.6f CaEnd=%.6f"%(len(t),t[0],t[-1],Ta[0],Ta[-1],Ca[0],Ca[-1]))

def fT(t,Tset,tau):  return Tset-(Tset-28.0)*np.exp(-t/tau)
def fC(t,Cset,tau):  return Cset-(Cset-0.01963)*np.exp(-t/tau)

pT,_ = curve_fit(fT,t,Ta,p0=[50.2,1800.0])
pC,_ = curve_fit(fC,t,Ca,p0=[0.0509,2800.0])
rT = Ta-fT(t,*pT); rC = Ca-fC(t,*pC)
out.append("Tset=%.4f tauT=%.2f  RMSE=%.4f  R2=%.6f"%(pT[0],pT[1],np.sqrt(np.mean(rT**2)),1-np.sum(rT**2)/np.sum((Ta-Ta.mean())**2)))
out.append("Cset=%.6f tauC=%.2f  RMSE=%.3e R2=%.6f"%(pC[0],pC[1],np.sqrt(np.mean(rC**2)),1-np.sum(rC**2)/np.sum((Ca-Ca.mean())**2)))

# per-parameter standard errors
def se(f,p,t,y):
    J=np.zeros((len(t),len(p)))
    for i in range(len(p)):
        dp=p.copy(); h=abs(p[i])*1e-6+1e-12; dp[i]+=h
        J[:,i]=(f(t,*dp)-f(t,*p))/h
    s2=np.sum((y-f(t,*p))**2)/(len(t)-len(p))
    cov=s2*np.linalg.inv(J.T@J)
    return np.sqrt(np.diag(cov))
out.append("se(T): %s"%np.array2string(se(fT,pT,t,Ta),precision=4))
out.append("se(C): %s"%np.array2string(se(fC,pC,t,Ca),precision=6))
open(r"C:\Users\qing1\Desktop\数A\new_approach\step1_fit.txt","w",encoding="utf-8").write("\n".join(out))
print("done")
