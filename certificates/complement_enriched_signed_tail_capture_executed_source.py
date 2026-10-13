import numpy as np
from flint import arb,arb_mat,ctx
ctx.prec=192
p='li-positivity/certificates/';a=np.load(p+'complement_residual_enrichment_action_matrices.npz');ph=np.load(p+'complement_physical_frame_locator.npz')
def exact(m):return arb_mat([[arb(float(v)) for v in row] for row in m])
def fro(m):return sum((m[i,j]**2 for i in range(m.nrows()) for j in range(m.ncols())),arb(0)).sqrt()
G=exact(a['GY']);G=(G+G.transpose())/2;C=exact(a['CY']);H=exact(ph['Hcenter']);H=(H+H.transpose())/2;g=arb(14115739)/10**9;W=(G+g*arb_mat(np.eye(56).astype(int).tolist())).inv();R=H-C*W*C.transpose()
rv=np.array([[float(R[i,j].mid()) for j in range(56)] for i in range(56)]);v=exact(np.linalg.eigh(rv)[1][:,-1,None]);r=(v.transpose()*R*v)[0,0]/(v.transpose()*v)[0,0]
e=arb(32)/10**10+2*(arb(1)/20000)*fro(W*C.transpose())+(arb(1)/20000)**2*fro(W)
b=arb(979236349)/2000000000000+arb(11)*10**7/99999999*(r-e)
assert b>arb(14)/100
print('Verified robust capture barrier',b,'with Frobenius W error',e)
