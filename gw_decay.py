import numpy as np
from constants import *
import math

#kepler third law to obtain a in Rsun as function of P (in days) and m1, m2 in Msun
def kepler3_a(P,m1,m2):
    return float(((P*24.*3600.)**2*cgrav*(m1+m2)*Msun/(4.*math.pi**2))**(1./3.)/Rsun);

#kepler third law to obtain P in days as function of a (in Rsun) and m1, m2 in Msun
def kepler3_P(a,m1,m2):
    return float((4.*math.pi**2*(a*Rsun)**3/(cgrav*(m1+m2)*Msun))**(1./2.)/(24.*3600.));

#merger time due to GWs in a circular orbit
def beta_gw(m1,m2):
    return 64./5.*cgrav**3/clight**5*m1*m2*(m1+m2)*Msun**3

#merger time (in Gyr) for initial orbital separation a (Rsun) and masses m1,m2 (Msun)
#Peters (1964)
def T_merger_a(a,m1,m2):
    return (a*Rsun)**4/(4.*beta_gw(m1,m2))/(secyear*1e9)

#merger time (in Gyr) for initial orbital period P (days) and masses m1,m2 (Msun)
def T_merger_P(P,m1,m2):
    return T_merger_a(kepler3_a(P,m1,m2),m1,m2)

prefact_dadt = -64./5.*cgrav**3/clight**5
v73div24 = 73./24.
v37div96 = 37./96.
#da_dt from Peters (1964), eq. (5.6)
def da_dt(a,e,m1,m2):
    #return -64./5.*(cgrav**3*m1*m2*(m1+m2))/(clight**5*a**3*(1.-e**2)**(3.5)) * \
    #        (1.+73./24.*e**2+37./96.*e**4)
    return prefact_dadt*(m1*m2*(m1+m2))/(a**3*(1.-e**2)**(3.5)) * \
            (1.+v73div24*e**2+v37div96*e**4)

prefact_dedt = -304./15.*cgrav**3/clight**5
v121div304 = 121./304.
#de_dt from Peters (1964), eq. (5.7)
def de_dt(a,e,m1,m2):
    #return -304./15.*e*(cgrav**3*m1*m2*(m1+m2))/(clight**5*a**4*(1.-e**2)**(2.5)) * \
    #        (1.+121./304.*e**2)
    return prefact_dedt*e*(m1*m2*(m1+m2))/(a**4*(1.-e**2)**(2.5)) * \
            (1.+v121div304*e**2)

#Merger time (in Gyr) for initial orbital separation a (Rsun) and eccentricity e.
#This integrates de/dt and da/dt
#Timestep is defined as dt_control*min(a/(da/dt),e/(de/dt)).
#Evolution is terminated once the separation is equal to end_condition*max(R_1sch,R_2sch)
#Timestep limit for eccentricity is only used while e>emin
def T_merger_a_e(a,e,m1,m2,dt_control=0.01,end_condition=100.0,emin=1e-3):
    a = a*Rsun
    m1 = m1*Msun
    m2 = m2*Msun
    final_sep = end_condition*(2*cgrav/clight**2*max(m1,m2))
    t = 0
    k = 0
    v1div6 = 1.0/6.0
    while a>final_sep:
        k += 1
        #####Euler is too slow
        #dadt = da_dt(a,e,m1,m2)
        ##if e>0:
        #dedt = de_dt(a,e,m1,m2)
        #if e > 1e-3:
        #    dt = abs(dt_control*min(a/dadt,e/dedt))
        #else:
        #    dt = abs(dt_control*a/dadt)
        #a += dadt*dt
        #if e>0:
        #    e = min(e+dedt*dt,0)
        #t += dt

        #second order RK
        k1 = da_dt(a,e,m1,m2)
        l1 = de_dt(a,e,m1,m2)
        if e > 1e-3:
            dt = abs(dt_control*min(a/k1,e/l1))
        else:
            dt = abs(dt_control*a/k1)
        k2 = da_dt(a+0.5*k1*dt,max(e+0.5*l1*dt,0),m1,m2)
        l2 = de_dt(a+0.5*k1*dt,max(e+0.5*l1*dt,0),m1,m2)
        k3 = da_dt(a+0.5*k2*dt,max(e+0.5*l2*dt,0),m1,m2)
        l3 = de_dt(a+0.5*k2*dt,max(e+0.5*l2*dt,0),m1,m2)
        k4 = da_dt(a+k3*dt,max(e+l3*dt,0),m1,m2)
        l4 = de_dt(a+k3*dt,max(e+l3*dt,0),m1,m2)
        kf = (k1+2*k2+2*k3+k4)*v1div6
        lf = (l1+2*l2+2*l3+l4)*v1div6
        a += kf*dt
        if e>0:
            e = max(e+lf*dt,0)
        t += dt

    return t/secyear/1e9

#Merger time (in Gyr) for initial orbital period P (days) and eccentricity e.
def T_merger_P_e(P,e,m1,m2,dt_control=0.01,end_condition=100.0,emin=1e-3):
    return T_merger_a_e(kepler3_a(P,m1,m2),e,m1,m2,dt_control,end_condition,emin)

