# New norm-certificate question: fine actual r8 cell masses.
# Coarse zeros are established once in this run; only nonzero coarse
# cells are refined. No previously successful trace loop is replayed.
from flint import arb, acb, ctx
from pathlib import Path
import json, hashlib, time
ctx.prec=160
root=Path(__file__).resolve().parents[1]
rows=[(arb(n).log(),2*arb(p).log()/arb(n).sqrt()) for n,p in [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2)]]
pi=arb.pi();logpi=pi.log()
def upperpart(k,den):
 lo=arb(k)/den; hi=arb(k+1)/den
 c=(lo+hi)/2;h=(hi-lo)/2;X=c+arb(0,h)
 L=arb(10) if k==0 else min(arb(10),2/(lo*lo)+pi/(2*lo))
 m=acb(arb(1)/4,c/2).digamma().real-logpi+arb(0,L*h)
 delta=sum((coef*(a*X).cos() for a,coef in rows),arb(0))-m
 U=delta.upper()
 return U if U>0 else arb(0)
start=time.monotonic(); coarsezeros=[];active=[]
for k in range(80000):
 U=upperpart(k,8)
 if U>0:active.append(k)
 else:coarsezeros.append(k)
# Fine cells width1/256 are integrated; four form each kernel cell1/64.
masses={};finez=0;fineactive=0;masssum=arb(0)
for k in active:
 for j in range(32):
  f=32*k+j;U=upperpart(f,256)
  if U>0:
   key=f//4;mass=U/(256*pi)
   masses[key]=masses.get(key,arb(0))+mass
   masssum+=mass;fineactive+=1
  else:finez+=1
cells=[{'k':k,'denominator':64,'measure_upper_ball':str(v),'measure_upper_float':float(v.upper())} for k,v in sorted(masses.items())]
assert masssum<14
record={'status':'PASS new actual norm-cell execution','precision_bits':160,'coarse_width':'1/8','fine_mass_width':'1/256','kernel_cell_width':'1/64','coarse_zero_indices':coarsezeros,'coarse_active_indices':active,'coarse_zero_count':len(coarsezeros),'coarse_active_count':len(active),'fine_active_count':fineactive,'fine_zero_count':finez,'kernel_cells':cells,'total_measure_upper_ball':str(masssum),'Fourier_tail':'r8=0 above10000 byGC1-GC2','gamma_method':'fullacbcenterpluscompleteST2radius','elapsed_seconds':time.monotonic()-start,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(root/'certificates/complement_signed_symbol_norm_cells.json').write_text(json.dumps(record,separators=(',',':'))+'\n')
print('PASS newnormrun coarseactive',len(active),'kernelcells',len(cells),'fineactive',fineactive,'measureupper',masssum,'elapsed',record['elapsed_seconds'])
