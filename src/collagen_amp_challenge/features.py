
import math
import numpy as np
from collections import Counter

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA)
BOMAN_SCALE = {'L':4.92,'I':4.92,'V':4.04,'F':2.98,'M':2.35,'W':2.33,'A':1.81,'C':1.28,'G':0.94,'Y':-0.14,'P':0.0,'T':-2.57,'S':-3.40,'H':-4.66,'Q':-5.54,'K':-5.55,'N':-6.64,'E':-6.81,'D':-8.72,'R':-14.92}
EISENBERG = {'A':0.62,'R':-2.53,'N':-0.78,'D':-0.90,'C':0.29,'Q':-0.85,'E':-0.74,'G':0.48,'H':-0.40,'I':1.38,'L':1.06,'K':-1.50,'M':0.64,'F':1.19,'P':0.12,'S':-0.18,'T':-0.05,'W':0.81,'Y':0.26,'V':1.08}
APV = {'A':0.307,'R':0.106,'N':0.240,'D':0.479,'C':0.165,'Q':0.248,'E':0.449,'G':0.265,'H':0.202,'I':0.198,'L':0.246,'K':0.111,'M':0.265,'F':0.246,'P':0.327,'S':0.281,'T':0.242,'W':0.172,'Y':0.185,'V':0.200}
ALPHA = {'A':1.41,'C':0.70,'D':1.01,'E':1.51,'F':1.13,'G':0.57,'H':1.00,'I':1.08,'K':1.16,'L':1.34,'M':1.30,'N':0.90,'P':0.59,'Q':1.17,'R':0.98,'S':0.77,'T':0.83,'V':1.06,'W':1.08,'Y':0.69}
PKA={'Nterm':8.6,'Cterm':3.6,'C':8.5,'D':3.9,'E':4.1,'H':6.5,'K':10.8,'R':12.5,'Y':10.1}

def canonical(seq):
    s = str(seq).strip().upper()
    return s if s and set(s) <= AA_SET else None

def charge(seq,pH=7.0):
    s=seq
    pos=1/(1+10**(pH-PKA['Nterm']))
    for a in 'HKR': pos += s.count(a)/(1+10**(pH-PKA[a]))
    neg=1/(1+10**(PKA['Cterm']-pH))
    for a in 'CDEY': neg += s.count(a)/(1+10**(PKA[a]-pH))
    return pos-neg

def boman(seq): return -sum(BOMAN_SCALE[a] for a in seq)/len(seq)
def apv(seq): return sum(APV[a] for a in seq)/len(seq)
def alpha(seq): return sum(ALPHA[a] for a in seq)/len(seq)

def hmoment(seq, angle=100.0, window=11):
    w=min(len(seq),window); th=math.radians(angle); best=0.0
    for st in range(len(seq)-w+1):
        vc=vs=0.0
        for k,a in enumerate(seq[st:st+w],1):
            h=EISENBERG[a]; ang=th*k
            vc += h*math.cos(ang); vs += h*math.sin(ang)
        best=max(best, math.hypot(vc,vs)/w)
    return best

def feature_names():
    names=["length","charge","boman","hmoment","apv","alpha","basic_frac","acidic_frac","aromatic_frac","hydrophobic_frac","proline_frac","glycine_frac"]
    names += [f"aa_{a}" for a in AA]
    names += [f"dipep_{a}{b}" for a in AA for b in AA]
    return names

def featurize(seq):
    s=canonical(seq)
    if s is None: raise ValueError(f"Noncanonical sequence: {seq}")
    n=len(s); c=Counter(s)
    vals=[
        n,charge(s),boman(s),hmoment(s),apv(s),alpha(s),
        sum(c[a] for a in "KRH")/n,
        sum(c[a] for a in "DE")/n,
        sum(c[a] for a in "FWY")/n,
        sum(c[a] for a in "AILMFWVY")/n,
        c["P"]/n,c["G"]/n
    ]
    vals += [c[a]/n for a in AA]
    dc=Counter(s[i:i+2] for i in range(len(s)-1))
    denom=max(1,n-1)
    vals += [dc[a+b]/denom for a in AA for b in AA]
    return np.asarray(vals,dtype=np.float32)
