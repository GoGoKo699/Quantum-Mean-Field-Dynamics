import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
rng = np.random.default_rng(241111)
T = 1.0
sigma = 0.8
samples = 24000

def v(t):
    return 0.2*np.sin(np.pi*t) + 0.1*np.abs(t-0.4)

def drift(x, t):
    return -0.7*np.sin(x) + 0.4*v(t)*np.cos(x)

def generator(x, t):
    return drift(x,t)*np.cos(x)-0.5*sigma*sigma*np.sin(x)

records=[]
for H in [1.0,0.25,0.0625,0.015625]:
    left=T-H
    for level in [3,4,5,6,7,8,9]:
        h=2.0**(-level)
        tau=left+H*rng.random(samples)
        sign=np.where(tau<left+H/2,1.0,-1.0)
        xf=np.full(samples,0.3)
        xc=xf.copy()
        coarse_noise=np.zeros(samples)
        for i in range(2**level):
            t=i*h
            dt=np.clip(tau-t,0.0,h)
            dw=np.sqrt(dt)*rng.standard_normal(samples)
            xf += drift(xf,t)*dt + sigma*dw
            coarse_noise += dw
            if i%2:
                tc=(i-1)*h
                dtc=np.clip(tau-tc,0.0,2*h)
                xc += drift(xc,tc)*dtc+sigma*coarse_noise
                coarse_noise.fill(0.0)
        tf=np.floor(tau/h)*h
        tc=np.floor(tau/(2*h))*(2*h)
        diff=(H/2)*sign*(generator(xf,tf)-generator(xc,tc))
        std=float(np.std(diff,ddof=1))
        rms=float(np.sqrt(np.mean(diff*diff)))
        records.append({'H':H,'level':level,'h':h,'n_samples':samples,'std':std,'rms':rms,'std_over_Hh':std/(H*h),'rms_over_Hh':rms/(H*h),'mean':float(np.mean(diff))})
out={'model':'dX=(-0.7 sin X+0.4 v(t) cos X)dt+0.8dW, g=sin X, v=0.2sin(pi t)+0.1|t-0.4|','purpose':'Test L2 fine-coarse random-time Dynkin-detail correction proportional to H*h, with partial terminal increments shared exactly in floating-point simulation. Finite numerical check, not proof or quantum execution.','records':records}
(ROOT/'random_time_coupling_check.json').write_text(json.dumps(out,indent=2)+'\n')
for H in [1.0,0.25,0.0625,0.015625]:
    rows=[r for r in records if r['H']==H]
    slope=float(np.polyfit(np.log([r['h'] for r in rows]),np.log([r['std'] for r in rows]),1)[0])
    ratios=[r['std_over_Hh'] for r in rows]
    print(json.dumps({'H':H,'std_loglog_slope':slope,'ratio_min':min(ratios),'ratio_max':max(ratios)}))
