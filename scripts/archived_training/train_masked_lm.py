from pathlib import Path
import random, json, joblib, sys
import numpy as np, pandas as pd, torch
from torch import nn
from torch.utils.data import Dataset, DataLoader, random_split

ROOT=Path('/mnt/data/collagen_amp_ml'); MODEL=ROOT/'models'; OUT=ROOT/'outputs'; PUB=Path('/mnt/data/ML-AMP')
SEED=42; random.seed(SEED); np.random.seed(SEED); torch.manual_seed(SEED)
AA='ACDEFGHIKLMNPQRSTVWY'; PAD=0; MASK=21; vocab={a:i+1 for i,a in enumerate(AA)}; inv={i+1:a for i,a in enumerate(AA)}
MAXLEN=30

# AMP sequences, challenge-like length range
amp=pd.read_csv(PUB/'unlabelled_positive.csv').Sequence.astype(str).str.upper()
amp=[s for s in amp if 18<=len(s)<=30 and set(s)<=set(AA)]
# collagen candidates excluding repeat/cuticle; sample to balance sources
coll=pd.read_csv('/mnt/data/amp_challenge_provenance_refined/master_all_candidates_with_refined_provenance.csv')
coll=coll[~coll.Provenance_Refined.isin(['collagen-like/repeat','nematode/cuticle collagen'])].Lowest_Peptide.astype(str).tolist()
random.shuffle(coll); coll=coll[:len(amp)]
seqs=amp+coll; random.shuffle(seqs)

class MaskDS(Dataset):
 def __init__(self,seqs):self.seqs=seqs
 def __len__(self):return len(self.seqs)
 def __getitem__(self,idx):
  s=self.seqs[idx]; ids=[vocab[a] for a in s]+[PAD]*(MAXLEN-len(s))
  pos=random.randrange(len(s)); target=ids[pos]; ids[pos]=MASK
  return torch.tensor(ids,dtype=torch.long),torch.tensor(pos),torch.tensor(target-1)

class MaskedBiLSTM(nn.Module):
 def __init__(self):
  super().__init__(); self.emb=nn.Embedding(22,48,padding_idx=PAD); self.lstm=nn.LSTM(48,64,batch_first=True,bidirectional=True); self.head=nn.Sequential(nn.Linear(128,96),nn.ReLU(),nn.Linear(96,20))
 def forward(self,x,pos):
  h,_=self.lstm(self.emb(x)); b=torch.arange(x.size(0),device=x.device); return self.head(h[b,pos])

ds=MaskDS(seqs); nval=max(1000,int(.1*len(ds))); ntr=len(ds)-nval
tr,val=random_split(ds,[ntr,nval],generator=torch.Generator().manual_seed(SEED))
tl=DataLoader(tr,batch_size=512,shuffle=True,num_workers=0); vl=DataLoader(val,batch_size=512,num_workers=0)
model=MaskedBiLSTM(); opt=torch.optim.AdamW(model.parameters(),lr=2e-3,weight_decay=1e-4); lossfn=nn.CrossEntropyLoss()
best=1e9; hist=[]
for ep in range(1,7):
 model.train(); tot=correct=n=0; ls=0
 for x,p,y in tl:
  opt.zero_grad(); z=model(x,p); loss=lossfn(z,y); loss.backward(); opt.step(); ls+=loss.item()*len(y); correct+=(z.argmax(1)==y).sum().item();n+=len(y)
 model.eval();vls=vc=vn=0
 with torch.no_grad():
  for x,p,y in vl:
   z=model(x,p); l=lossfn(z,y);vls+=l.item()*len(y);vc+=(z.argmax(1)==y).sum().item();vn+=len(y)
 rec={'epoch':ep,'train_loss':ls/n,'train_acc':correct/n,'val_loss':vls/vn,'val_acc':vc/vn};hist.append(rec);print(rec)
 if rec['val_loss']<best:
  best=rec['val_loss'];torch.save({'state_dict':model.state_dict(),'aa':AA,'maxlen':MAXLEN},MODEL/'masked_peptide_lm.pt')
(OUT/'masked_lm_metrics.json').write_text(json.dumps({'n_amp':len(amp),'n_collagen':len(coll),'n_total':len(seqs),'history':hist},indent=2))
