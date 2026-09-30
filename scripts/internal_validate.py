from pathlib import Path
from rapidfuzz import fuzz
AA=set('ACDEFGHIKLMNPQRSTVWY');R=Path(__file__).resolve().parents[1]
def read(p):
 s=[];cur=''
 for ln in Path(p).read_text().splitlines():
  if ln.startswith('>'):
   if cur:s.append(cur);cur=''
  elif ln.strip():cur+=ln.strip().upper()
 if cur:s.append(cur)
 return s
lib=read(R/'generate/library.fasta');top=read(R/'generate/top.fasta')
assert len(lib)==50000 and len(set(lib))==50000
assert len(top)==100 and len(set(top))==100 and set(top)<=set(lib)
assert all(8<=len(s)<=50 and set(s)<=AA for s in lib)
mx=max(fuzz.ratio(a,b)/100 for i,a in enumerate(top) for b in top[i+1:])
assert mx<0.8
print({'library':len(lib),'top':len(top),'max_internal_similarity':mx})
