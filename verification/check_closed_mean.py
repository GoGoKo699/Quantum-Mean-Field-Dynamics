"""Check feedback reconstruction on a classically solvable mean-field diffusion.

dX = lam * E[sin X] dt + sigma dW, X_0=x0.
The exact mean is known. This is a validation case, not an advantage benchmark.
The reconstruction uses exact expectations plus deterministic bounded coefficient
perturbations. It checks global feedback, interpolation and scalar readout;
it does not emulate amplitude estimation or certify a quantum gate count.
"""
from pathlib import Path
import json
import math
import numpy as np
from scipy.integrate import quad

OUT = Path(__file__).resolve().parent
LAM, SIGMA, X0, T = 0.8, 0.7, 0.4, 1.0
G0 = math.sin(X0)
B0 = LAM + SIGMA**2 / 2

def exact_mu(t):
    return 2*np.arctan(np.tan(X0/2)*np.exp(
        2*LAM/SIGMA**2*(1-np.exp(-SIGMA**2*t/2))))

def exact_mean(t):
    return np.exp(-SIGMA**2*t/2)*np.sin(exact_mu(t))

def exact_terminal_cos():
    return float(np.exp(-SIGMA**2*T/2)*np.cos(exact_mu(T)))

def reconstruct(target_nodes, delta, iteration):
    """Midpoint coefficients, each with a deterministic error <= delta."""
    n = len(target_nodes)-1
    ans = np.empty_like(target_nodes)
    ans[0] = target_nodes[0]
    ans[n] = target_nodes[n] + delta*math.sin(1.7+iteration)
    width = n
    while width >= 2:
        left = np.arange(0, n, width)
        right = left + width
        mid = left + width//2
        detail = target_nodes[mid] - (target_nodes[left]+target_nodes[right])/2
        error = delta*np.sin(0.37*(mid+1)+0.91*iteration)
        ans[mid] = (ans[left]+ans[right])/2 + detail + error
        width //= 2
    return np.clip(ans, -1.0, 1.0)

feedback = []
for eps in [2.0**-4, 2.0**-6, 2.0**-8]:
    eta = eps/(8*math.exp(LAM*T))
    j = max(1, math.ceil(math.log2(8*B0*T/eta)))
    n = 2**j
    h = T/n
    delta = eta/(8*(j+1))
    grid = np.linspace(0, T, n+1)
    r = 0
    while 2*(LAM*T)**r/math.factorial(r) > eps/8:
        r += 1
    v = np.full(n+1, G0)
    rows = []
    for it in range(1, r+1):
        integ = np.concatenate(([0.0], np.cumsum((v[:-1]+v[1:])*h/2)))
        phi = np.exp(-SIGMA**2*grid/2)*np.sin(X0+LAM*integ)
        updated = reconstruct(phi, delta, it)
        detail_error = float(np.max(np.abs(updated-phi)))
        max_slope = float(np.max(np.abs(np.diff(updated)))/h)
        curve_error = float(np.max(np.abs(updated-exact_mean(grid))))
        assert detail_error <= (j+1)*delta*(1+1e-10)
        assert max_slope <= B0 + 4*delta/h + delta/T + 1e-10
        rows.append({"iteration":it,"node_reconstruction_error":detail_error,
                     "max_slope":max_slope,"node_feedback_error":curve_error})
        v = updated
    mu_approx = X0 + LAM*float(np.sum((v[:-1]+v[1:])*h/2))
    terminal = math.exp(-SIGMA**2*T/2)*math.cos(mu_approx)
    error = abs(terminal-exact_terminal_cos())
    assert error < eps
    feedback.append({"epsilon":eps,"eta":eta,"temporal_levels":j,
        "nodes":n+1,"picard_iterations":r,"coefficient_error_budget":delta,
        "terminal_cos_estimate":terminal,"exact_terminal_cos":exact_terminal_cos(),
        "terminal_abs_error":error,"iterations":rows})

# A separate deterministic check of the random-time Dynkin identity.
knots = np.linspace(0, 1, 9)
values = 0.2*np.sin(np.pi*knots) + 0.1*np.abs(knots-0.4)
prefix = np.concatenate(([0.0], np.cumsum(np.diff(knots)*(values[:-1]+values[1:])/2)))
def curve(t):
    return float(np.interp(t, knots, values))
def curve_integral(t):
    idx = min(int(t*8), 7)
    dt = t-knots[idx]
    slope = (values[idx+1]-values[idx])*8
    return float(prefix[idx]+values[idx]*dt+slope*dt*dt/2)
def phi(t):
    return math.exp(-SIGMA**2*t/2)*math.sin(X0+LAM*curve_integral(t))
def generator_expectation(t):
    mu = X0+LAM*curve_integral(t)
    return math.exp(-SIGMA**2*t/2)*(LAM*curve(t)*math.cos(mu)
                                  -SIGMA**2/2*math.sin(mu))
identities=[]
for width in [1.0, 0.25, 0.015625]:
    left=T-width; mid=left+width/2
    val=phi(mid)-(phi(left)+phi(T))/2
    def integrate(a,b):
        pts=[float(x) for x in knots if a<x<b]
        return quad(generator_expectation,a,b,points=pts,epsabs=1e-13,epsrel=1e-13)[0]
    dynkin=(integrate(left,mid)-integrate(mid,T))/2
    assert abs(dynkin-val)<1e-11
    identities.append({"H":width,"detail_from_values":val,
        "detail_from_generator":dynkin,"absolute_discrepancy":abs(dynkin-val)})

result={"scope":"Exact-mean feedback reconstruction with bounded deterministic coefficient perturbations; no quantum simulation or advantage claim.",
        "model":{"lambda":LAM,"sigma":SIGMA,"x0":X0,"T":T},
        "feedback":feedback,"dynkin_identity":identities}
(OUT/"closed_mean_results.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"terminal_errors":[x["terminal_abs_error"] for x in feedback],
                  "max_dynkin_discrepancy":max(x["absolute_discrepancy"] for x in identities)}))

