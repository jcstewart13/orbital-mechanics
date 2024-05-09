
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import streamlit as st

m=100*1000
mu=398601
lunarmu=4905
solarmu=132.71*(10**9)
cd=1
s=10**(-5)

def T(y, t):
    mpayload=100
    mtotal=y[6]
    mfuel=mtotal-mpayload
    direction=np.array([y[3],y[4],y[5]])/(np.sqrt(y[3]**2+y[4]**2+y[5]**2))
    if mfuel > 0:
        if t<50:
            thr=INIT_THRUST*direction*y[6]/100000
        elif t<1200:
            thr=FIN_THRUST*direction*y[6]/100000
            
        else:
            thr = [0,0,0]
    else:
        thr=[0,0,0]
    return np.array(thr)

def lunar(y, t):
    lunar0 = 0-lunarmu*y[0]/(np.sqrt((y[0]-384400)**2+y[1]**2+y[2]**2)**3)
    lunar1 = 0-lunarmu*y[1]/(np.sqrt((y[0]-384400)**2+y[1]**2+y[2]**2)**3)
    lunar2 = 0-lunarmu*y[2]/(np.sqrt((y[0]-384400)**2+y[1]**2+y[2]**2)**3)
    return np.array([lunar0, lunar1, lunar2])

def solar(y, t):
    solar0 = 0-solarmu*y[0]/(np.sqrt((y[0]-149*10**6)**2+y[1]**2+y[2]**2)**3)
    solar1 = 0-solarmu*y[1]/(np.sqrt((y[0]-149*10**6)**2+y[1]**2+y[2]**2)**3)
    solar2 = 0-solarmu*y[2]/(np.sqrt((y[0]-149*10**6)**2+y[1]**2+y[2]**2)**3)

    return np.array([solar0, solar1, solar2])

def J2(y, t):
    magni=np.sqrt(y[0]**2+y[1]**2+y[2]**2)
    J2=[0,0,0]
    J2[0]=-mu*0.0010826*6371**2*(3/(2*magni**4)-(15/2)*(y[0]**2)/(magni**6))*y[0]/magni
    J2[1]=-mu*0.0010826*6371**2*(3/(2*magni**4)-(15/2)*(y[1]**2)/(magni**6))*y[1]/magni
    J2[2]=-mu*0.0010826*6371**2*(3*y[2]/(magni**5)+(3/(2*magni**4)-(15/2)*(y[2]**2)/(magni**6))*y[2]/magni)
    return J2

def Isp(y, t):
    mpayload=100
    mtotal=y[6]
    mfuel=mtotal-mpayload
    if mfuel>5000:
        isp=3
    elif mfuel>1300:
        isp=3.5
    elif mfuel>5:
        isp=4
    else:
        isp=2.6
    return isp

def drag(y, t):
    alt=np.sqrt(y[0]**2+y[1]**2+y[2]**2)-6378
    if alt < 0:
        alt = 0
    d=np.exp(-alt/8.5)*1.225*(1000)**3
    q=1/2*(y[3]**2+y[4]**2+y[5]**2)*d
    sd=cd*s*q

    direction=np.array([y[3],y[4],y[5]])/(np.sqrt(y[3]**2+y[4]**2+y[5]**2))
    vd=-direction*sd

    return vd


def df(t, y):
    isp = Isp(y,t)
    thr = T(y,t)
    vd = drag(y,t)
    J2out = J2(y,t)
    lunarout = lunar(y,t)
    solarout = solar(y,t)
    dy1dt = y[3]
    dy2dt = y[4]
    dy3dt = y[5]
    dy4dt = 0-mu*y[0]/(np.sqrt(y[0]**2+y[1]**2+y[2]**2)**3)+thr[0]/y[6]+vd[0]/y[6]+J2out[0]+lunarout[0]+solarout[0]
    dy5dt = 0-mu*y[1]/(np.sqrt(y[0]**2+y[1]**2+y[2]**2)**3)+thr[1]/y[6]+vd[1]/y[6]+J2out[1]+lunarout[0]+solarout[0]
    dy6dt = 0-mu*y[2]/(np.sqrt(y[0]**2+y[1]**2+y[2]**2)**3)+thr[2]/y[6]+vd[2]/y[6]+J2out[2]+lunarout[0]+solarout[0]

    dmdt = 0-(np.sqrt(thr[0]**2+thr[1]**2+thr[2]**2)/isp)
    
    return np.array([dy1dt, dy2dt, dy3dt, dy4dt, dy5dt, dy6dt, dmdt])

alt = 1
normvelocity0 = 0.5

with st.sidebar:
    theta = st.number_input(label="Launch Angle", 
                      min_value=0.0,
                      max_value=0.01,
                      value=0.001,
                      step=.0001,
                      format="%0.4f"
                      )
    INIT_THRUST = st.number_input(label = "Intial Thrust for the First 50 Seconds",
                        min_value=0,
                        max_value=2500,
                        value=1000,
                        step=100)
    FIN_THRUST = st.number_input(label = "Intial Thrust for the Last 1150 Seconds",
                        min_value=0,
                        max_value=2500,
                        value=1000,
                        step=100)
    input_rt =  st.number_input(label = "How many seconds do you want your code to run? (Recommended value is around 75000)",
                        min_value=25000,
                        max_value=150000,
                        value=75000,
                        step=100)

st.write("""This code applies differential calculus to simulate a rocket launched on Earth to orbit around Earth. The code implements thrust, lunar and solar gravitational pulls, J2 (Earth is not a perfect sphere), specific impulse, and drag. 

Try playing around with the thrust and launch angle to successfully simulate a rocket into orbit. Click “run” to see your 3D and 2D (altitude vs. time) graphs. If your 2D graph is a line going up, your thrust is too large or your launch angle is too small. If your rocket crashes back into the Earth (lithobrakes), at some point in time the altitude in your 2D graph will be negative. The 2D graph of a successful orbit will be in the shape of a sine curve. 
""")

y0 = np.array([6378+alt,0,0,normvelocity0,0,-theta*normvelocity0,m])
sol=scipy.integrate.solve_ivp(df, [0,input_rt], y0, max_step = 500)

st.subheader('altitude (km) vs time (s)')
fig, ax = plt.subplots()
ax.plot(sol.t, np.sqrt(sol.y[0]**2+sol.y[1]**2+sol.y[2]**2)-6378)
st.write(fig)

st.subheader('3d plot')
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')
ax.plot(sol.y[0], sol.y[1], sol.y[2])
st.write(fig)


