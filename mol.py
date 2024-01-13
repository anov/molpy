import re
import math

# mol object .atoms[] .bonds[]
# atom .label .x .y .bonds[]
# bond .atom1 .atom2 .btype (1,2,3, 5-up, 9-down)

class Atom:
 label="C"
 def __repr__(self):
  try: return "{0}: {1} {2} {3}".format(self.num,self.x,self.y,self.label)
  except: return "{0} {1} {2}".format(self.x,self.y,self.label)
 def needsHs(self):
  nb=0
  if self.label in ["O","S","Se"]: nb=2
  elif self.label in ["N","P","B"]: nb=3
  if not nb: return 0
  for b in self.bonds:
   if b.btype==2: nb=-1
   elif b.btype==3: nb=-2
  return nb-len(self.bonds) if nb>len(self.bonds) else 0

class Bond:
 btype=1
 def __repr__(self):
  return "a1: {0} a2: {1} {2}".format(self.atom1,self.atom2,self.btype)
 def toatom(a):
  if self.atom1==a: return self.atom2
  elif self.atom2==a: return self.atom1
  else: return None
 
class Molecule:
 def __init__(self):
  self.atoms=[]
  self.bonds=[]
  
 def __repr__(self):
  return "atoms: {0} bonds: {1}".format(self.atoms,self.bonds)

 def addatom(self,x,y,label):
  atom=Atom()
  atom.bonds=[]
  atom.x=x
  atom.y=y
  atom.label=label;
  self.atoms.append(atom)
  return len(self.atoms)-1;  

 def addbond(self,natom1,natom2,btype):
  bond=Bond()
  bond.atom1=self.atoms[natom1]
  bond.atom2=self.atoms[natom2]
  bond.btype=int(btype)
  self.bonds.append(bond)
  return len(self.bonds)-1

 def findbond(self,natom1,natom2):
  a1=self.atoms[natom1]
  a2=self.atoms[natom2]
  for i,b in enumerate(self.bonds):
   if (b.atom1==a1) and (b.atom2==a2): return i
   if (b.atom2==a1) and (b.atom1==a2): return i
  return -1

 def rebuild_bond_lists(self):
  for i,a in enumerate(self.atoms):
   a.bonds=[];
   a.num=i;
  for i,b in enumerate(self.bonds):
   b.num=i;
   if b.btype==2:
    b.atom1.bonds.append(b)
    b.atom2.bonds.append(b)
   else:
    b.atom1.bonds.insert(0,b);
    b.atom2.bonds.insert(0,b);

 def update_mol_data(self):
  self.rebuild_bond_lists();
  for b in self.bonds:
   dx=b.atom2.x-b.atom1.x;
   dy=b.atom2.y-b.atom1.y;
   b.blen=math.sqrt(dx*dx+dy*dy)
   if b.blen:
    b.nx=dx/b.blen; b.ny=dy/b.blen;
   else:
    b.nx=0
    b.ny=1

def stringtomol(s):
 a=s.split(";")
 if len(a)!=2: return None
 aline=a[0].split(" ")
 bline=a[1].split(" ")
 mol=nMolecule()
 if len(aline) % 3 or len(bline) % 3: return mol
 for i in range(0,len(aline),3): mol.addatom(aline[i],aline[i+1],aline[i+2])
 for i in range(0,len(bline),3): mol.addbond(bline[i]-1,bline[i+1]-1,bline[i+2])
 mol.update_mol_data()
 return mol
 
def moltostring(mol):
 mol.rebuild_bond_lists()
 s=''
 for i in range(len(mol.atoms)):
  a=mol.atoms[i]
  if i: s+=" "
  s+=str(a.x)+" "+str(a.y)+" "+a.label 

 s+=";"
 for i in range(len(mol.bonds)): 
  a=mol.bonds[i]
  if i: s+=" "
  s+=str(a.atom1.num+1)+" "+str(a.atom2.num+1)+" "+str(a.btype)
 return s


def loadmol(s):
 mol=Molecule()
 
 lines=s.splitlines()
  
 m=re.match("(^[\\d ]{3})([\\d ]{3})",lines[3]) 
 natoms=int(m.group(1))
 nbonds=int(m.group(2))
 
 rg=re.compile("^(.{10})(.{10}).{10}(.{3}).{2}.{3}(.{3})")
 
 for i in range(4,natoms+4):
  m=rg.match(lines[i])
  atom=mol.addatom(m.group(1)*50,m.group(2)*50,m.group(3).trim())
  a=int(m.group(4))
  if a==0: continue
  if a==3: mol.atoms[atom].charge=1
  elif a==4: mol.atoms[atom].label+="'"
  elif a==5: mol.atoms[atom].charge=-1
  elif a>5: mol.atoms[atom].charge=-(a-4) # 3=+1 4=radical, 5=-1, etc  
  else: mol.atoms[atom].charge=4-a

 rg=re.compile("^(.{3})(.{3})(.{3})(.{3})")
 for i in rage(natoms+4,nbonds+natoms+4):
  m=rg.match(lines[i])

  # m[3] - 1,2,3; 4-arom, more - queries
  # m[4] - 1 - up, 6 - down, 4-either
  
  a=int(m.group(4))
  b=1;
  if a==1: b=5
  elif a==6: b=9
  else:
   a=int(m.group(3))
   if a<4: b=a
  mol.addbond(int(m.group(1))-1,int(m.group2)-1,b)


# M XXX extensions, override charge in atom description 
# charge
#        1  atom#  charge
# M  CHG  1   4   1  
# radical  
#         1  atom# 1 - rad, 2-carbene
# M  RAD  1   4   0
 
 rg=re.compile("^M (.{3})(.{3})(.{3})?(.{3})")
 for i in range(natoms+nbonds+4,len(lines)):
  m=rg.match(lines[i])
  if not m: continue
  if m.group(1)=="CHG": pass
  elif m.group(1)=="RAD": pass
 return mol

def savemol():
 pass

def getlrud(a): # find left,right,up,down bonds
 if (len(a.bonds)<3) or (len(a.bonds)>4): return None
 lr=[];u=-1;d=-1
 for i in range(len(a.bonds)): 
  st=getstereo(a,i)
  if st==0: lr.push(i)
  elif st==1:
   if u==-1: u=i
   else: return None
  elif st==2:
   if d==-1: d=i
   else: return None

 if (u==-1)  and  (d==-1):  return None
 if (u!=-1)  and  (d!=-1) and (len(a.bonds)==3): return None
 if len(lr)==3:
  if u==-1: d=lr[0]; lr.splice(0,1)
  elif d==-1:  u=lr[0]; lr.splice(0,1)
 elif len(lr)!=2: return None

 b1=a.bonds[lr[0]]; b2=a.bonds[lr[1]]
 if b1.atom1==a: dx1=b1.nx;dy1=b1.ny
 else: dx1=-b1.nx;dy1=-b1.ny
 if b2.atom1==a: dx2=b2.nx;dy2=b2.ny
 else: dx2=-b2.nx;dy2=-b2.ny

 vx=dx1+dx2;vy=dy1+dy2;
 if abs(vx+vy)<0.001: return None # lr not at angle
 if vy*dx1-vx*dy1<0: l=lr[0]; r=lr[1];
 else: l=lr[1]; r=lr[0]

 if u!=-1: 
  b1=a.bonds[u]
  if b1.atom1==a: dx1=b1.nx;dy1=b1.ny
  else: dx1=-b1.nx;dy1=-b1.ny
  if vx*dx1+vy*dy1>0: return None
 if d!=-1:
  b1=a.bonds[d] 
  if b1.atom1==a: dx1=b1.nx;dy1=b1.ny
  else: dx1=-b1.nx;dy1=-b1.ny 
  if vx*dx1+vy*dy1>0: return None
 
 return [l,r,u,d]


def getstereo(a,n): # 0 - none, 1 - up, 2-down
 b=a.bonds[n]
 if b.atom1!=a: return 0
 if b.btype==5: return 1
 if b.btype==9: return 2
 return 0
 
def bdir(nx,ny,a,b):
 p=nx*b.nx+ny+b.ny
 return p if a==b.atom1 else -p

def countstereo(a):
 if len(a.bonds)<3: return 0
 n=0
 for i in range(a.bonds): 
  if getstereo(a,i): n+=1
 return n
 

def showmatch(m):
 print "as:",
 for i,a in enumerate(m.mol1.atoms):
  if a.mappedto: print i,a.mappedto.num,":",
  else: print i,"N :",
 print
 print "bs:",
 for i,a in enumerate(m.mol1.bonds): 
  if a.mappedto: print i,a.mappedto.num,":",
  else: print i,"N :",
 print
 
class molmatch:
 matchedatoms=[]
 matchedbonds=[]

 def matchlabels(self,a1,a2):
  if self.strictarom and ((not a1.arom)!=(not a2.arom)): return False
  if a1.label==a2.label: return True
  if a1.label=="*": return True; # any 
  if a1.label=="X": return a2.label in ["F","Cl","Br","I"]
  if "," in a1.label: return a2.label in a1.label.split(',')
  return False

 def btypematch(self,b1,b2):
  if b1.btype%4!=b2.btype%4: return 0
  if self.ez and b1.ez and not b2.ez: return 0
  return 1
 
 def unmatch(self,iatoms,ibonds): # clear all mappings staring from the specified indices and up
  for i in reversed(range(iatoms,len(self.matchedatoms))):
   a=self.matchedatoms[i] 
   a.mappedto.mappedto=None
   a.mappedto=None
  del self.matchedatoms[iatoms:]
 
  for i in reversed(range(ibonds,len(self.matchedbonds))):
   a=self.matchedbonds[i] 
   a.mappedto.mappedto=None;
   a.mappedto=None
  del self.matchedbonds[ibonds:]
  
 def atomsmatch(self,a1,a2,bybond=None): 
 # Save state markers
  inmatchedatoms=len(self.matchedatoms)
  inmatchedbonds=len(self.matchedbonds)
  statestack=[]
 
  if (a1.mappedto==a2) and  (a2.mappedto==a1): return True # atoms already matched to each other (ring detected), return match
  if a1.mappedto or a2.mappedto: return False # one of the atoms already matched to something else, return mismatch
  if not self.matchlabels(a1,a2): return False # label don't match, return mismatch
  cbond1=0; cbond2=0
  # match up atoms - set mappings, and add atoms to the list 
  a1.mappedto=a2; a2.mappedto=a1; 
  self.matchedatoms.append(a1)
  
  
  bondmatched=1
 
  while cbond1<len(a1.bonds):
   bondmatched=0
   b1=a1.bonds[cbond1]

   if b1.mappedto:
    try:   
     a2.bonds.index(b1.mappedto)
     bondmatched=1; cbond1+=1; cbond2=0; continue
    except: break
	 
   while cbond2<len(a2.bonds):
   
    b2=a2.bonds[cbond2]
    if b2.mappedto or not self.btypematch(b1,b2):  cbond2+=1; continue # b2 is taken, or bonds mismatch, next
    else:
	# remember the state variables for possible further use
      nmatchedatoms=len(self.matchedatoms)
      nmatchedbonds=len(self.matchedbonds)
	
	# match up bonds - set mappings, and add bonds to the list 
      b1.mappedto=b2; b2.mappedto=b1
      self.matchedbonds.append(b1)
		
      nexta1=b1.atom2 if (b1.atom1==a1) else b1.atom1
      nexta2=b2.atom2 if (b2.atom1==a2) else b2.atom1
	
      bondmatched=self.atomsmatch(nexta1,nexta2,b1)
      if bondmatched and bybond and (bybond.btype==2) and bybond.ez: #check EZ
	
       nb1=bybond.ezbs.index(b1)
       if nb1<2: nb1+=2
       else: nb1-=2
       nb2=bybond.mappedto.ezbs.indexOf(b2)
       if nb2<2: nb2+=2
       else: nb2-=2
       if not bybond.ezbs[nb1]: bondmatched=not not bybond.mappedto.ezbs[nb2]
       else: bondmatched=bybond.ezbs[nb1].mappedto==bybond.mappedto.ezbs[nb2]
 
      if bondmatched:	
   # the bond and the substituent it attaches match the source. Save the state variables in stack (for possible rollback) and keep checking
  	 
       statestack.append([cbond1,cbond2,nmatchedatoms,nmatchedbonds])
       break
      else: self.unmatch(nmatchedatoms,nmatchedbonds); cbond2+=1
     
  
   if not bondmatched: # No match found for one of the bonds and the substituent it attaches. Roll back to the last saved state and try to match in a different way
    if len(statestack)==0: break; # No saved state to roll back to. So, substructures do not match in any way, return mismatch
    a=statestack.pop()
    cbond1=a[0]
    cbond2=a[1]
    nmatchedatoms=a[2]
    nmatchedbonds=a[3]
    # Roll back all matches since the last saved state
    self.unmatch(nmatchedatoms,nmatchedbonds)
    cbond2+=1
   else: cbond1+=1; cbond2=0

  
  if bondmatched:
   if not self.chiral: return True
   else: 
    if issameRS(a1,a2): return True
 
 # Atoms don't match. Roll back to where it was at the start of the function
  self.unmatch(inmatchedatoms,inmatchedbonds)  
  return False

 
def  prepareEZ(b): # prepare EZ info; return: 1 - has ez, 2 - no ez; 0 - faulty 
 if b.btype!=2: return 2
 b.ez=0;
 if (len(b.atom1.bonds)<2) or (len(b.atom2.bonds)<2): return 2
 if len(b.atom1.bonds)>3 or len(b.atom2.bonds)>3: return 0

 ndx=-b.ny
 ndy=b.nx
 b11=None
 b12=None
 b21=None
 b22=None
 for b1 in b.atom1.bonds: 
  if b1!=b:
   if not b11: b11=b1
   else: b12=b1 
 for b1 in b.atom2:
  if b1!=b: 
   if not b21: b21=b1
   else: b22=b1

 pa=bdir(ndx,ndy,b.atom1,b11)
 pb=bdir(ndx,ndy,b.atom2,b21)
  
 if (abs(pa)<0.01) or (abs(pb)<0.01): return 0

 if (b12 and bdir(ndx,ndy,b.atom1,b12)>0==(pa>0)): return 0 
 if (b22 and bdir(ndx,ndy,b.atom2,b22)>0==(pb>0)): return 0
 b.ez=1
 b.ezbs=[]
 if ((pa>0)==(pb>0)):
  b.ezbs[0]=b11
  b.ezbs[1]=b12
 else:
  b.ezbs[0]=b12
  b.ezbs[1]=b11

 b.ezbs[2]=b21
 b.ezbs[3]=b22
 return 1
  
def  prepareRS(a): # prepare chirality info; return: 1- OK; 0 - faulty 
 a.cc=false
 if not countstereo(a): return 1

 c=getlrud(a)
 
 if not c: return 0
 # check if all 4 diffrent if not - return 1;
 
 a.lrud=[]
 for i in range(4):
  n1=c[i]
  if n1!=-1: a.lrud[i]=a.bonds[n1]
  else: a.lrud[i]=None;
 
 a.cc=true 
 return 1

def issameRS(a1,a2): # assumes substituents already matched. Return 1 - match, 0 - no match
 if not (a1.cc or a2.cc): return 1
 elif not (a1.cc and a2.cc): return 0 # both achiral - match, one is, other isn't - mismatch

 l=[]
 for i in range(4):
  b=a1.lrud[i]
  try: l[i]=a2.lrud.index(b.mappedto if b else None)
  except: l[i]=-1

 l1=l[1]-l[0]
 l2=l[3]-l[2]
 if l1% 2: return (l1<0)==(l2<0)
 else: return (l1<0)!=(l2<0)

 

def matchmolecules(amol1,amol2,how):
 amol1.update_mol_data()
 amol2.update_mol_data() 

 m=molmatch()
 m.exact='E' in how
 m.chiral='C' in how 
 m.ez='Z' in how
 m.strictarom='A' in how  

 if m.exact and ((len(amol1.atoms)!=len(amol2.atoms)) or (len(amol1.bonds)!=len(amol2.bonds))): return False

 if m.chiral: #  find chiral centers, parse them
  for a in amol1.atoms: prepareRS(a)
  for a in amol2.atoms: prepareRS(a)

 if m.ez: #  find double, parse them
  for b in amol1.bonds: prepareEZ(b)
  for b in amol2.bonds: prepareEZ(b)
 
 #initialize key variables
 m.matchedatoms=[]
 matchedbonds=[]
 
 m.mol1=amol1
 m.mol2=amol2 
 if len(amol1.atoms)==0: return True # empty molecule matches anything
 for a in amol1.atoms: a.mappedto=None
 for a in amol2.atoms: a.mappedto=None
 
 for b in amol1.bonds: b.mappedto=None
 for b in amol2.bonds: b.mappedto=None
 
 n=0
 if 'N' in how and amol1.lastmatchN!=undefined: n=amol1.lastmatchN+1

 for i in range(n,len(amol2.atoms)): 
  if m.atomsmatch(amol1.atoms[0],amol2.atoms[i]):
   amol1.lastmatchN=i
   return True

 return False
   
class inser:
 s=''
 insn=0
 c=0
 ins={}
 
 def addins(self,r,p,s,first=False):
  if p==-1: p=len(self.s)
  #var k="n"+("0000"+String(p)).slice(-4)+String((first?r.insn++:(r.insn++)+500));
 
  k=p*9000+self.insn+1000 if first else self.insn
  self.insn+=1
  self.ins[k]=s
  return k


 def delins(r,k):
  del self.ins[k]

def addatomsmile(r,a,bybond):
 subs=[]

 if a.mappedto!=-1: # ring
  r.c+=1
  #s=(r.c<10)?String(r.c):"%"+r.c
  s=r.c if r.c<10 else "%"+str(r.c)
  addins(r,bybond.atom1.mappedto,s,True)
  addins(r,bybond.atom2.mappedto,s,True)
  return 'r'
 
 s=a.label if a.label else "C"
 lrud=0; lastps=0; invert=0
 if countstereo(a):
  lrud=getlrud(a)
  if lrud:  
   kchir=addins(r,len(r.s)+len(s)+1,'') 
   s='['+s+']'
 elif s in ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"]: s='['+s+']';
 if a.arom: s=s.lower()
 
 r.s+=s; a.mappedto=len(r.s)
 
 for i,b in enumerate(a.bonds):
  b=a.bonds[i]
  if b.mappedto:
   if b.mappedto=='r':
    if lrud:
     try: subs.append(lrud.index(i))
     except: subs.append(-1)
    continue
  b.mappedto=1
  p=addins(r,-1,'(')
  if b.btype==3: r.s+='#'
  elif (b.btype==2) and not (b.atom1.arom and b.atom2.arom): r.s+='=';
  if addatomsmile(r,b.atom2 if (b.atom1==a) else b.atom1,b)=='r': 
   delins(r,p); b.mappedto='r'
   if lrud: subs.append(lrud.indexOf(i));
  else: lastps=[p,addins(r,-1,')')]


 if lastps: delins(r,lastps[0]); delins(r,lastps[1])
 
 if lrud:
  s=''
  if -1 in lrud: s="H"; subs.insert(0,lrud.indexOf(-1))
  for i,b in enumerate(a.bonds): 
   if b==bybond: subs.insert(0,lrud.indexOf(i))
   elif b.mappedto!='r': subs.append(lrud.indexOf(i)) 
  l1=subs[1]-subs[0]; l2=subs[3]-subs[2]	
  if l1 % 2: i=(l1<0)==(l2<0)
  else: i=(l1<0)!=(l2<0) 
  r.ins[kchir]=("@@" if i else "@")+s; 
 return ''


def getSmiles(mol):
 mol.rebuild_bond_lists()
 markaromatic(mol,True)
 for a in mol.atoms: a.mappedto=-1
 for b in mol.bonds: b.mappedto=0
 r=inser()
 
 done=False
 while not done:
  done=True
  for a in mol.atoms:
   if a.mappedto<0:
    if r.s: r.s+='.'
    addatomsmile(r,mol.atoms[i])
    done=False; break

 lp=0
 while 1:
  p=9999;i=-1;k=0
  for f in r.ins.keys():
   pos=f // 9000
   pr=f % 9000
   if (pos<p) or (pos==p and (pr<i)): k=f; p=pos; i=pr
 
  if p>9000: break
  s+=r.s[lp:p]+r.ins[k]; 
  lp=p
  del r.ins[k];
 
 s+=r.s[lp:];
 return s

def loadSmiles(s):
 m=Molecule()
 ra=re.compile("(Cl|Br|[BCNOSPIFbcnosp\*\.\>\d]|%\d\d|\[.*?\]|)")
 rb=re.compile("([\)\(]+|)")
 astack=[]
 bheads={}
 cbtype=1
 satom=-1 
 
 curpos=0
 
 while curpos<len(s):
  #print "pos:",curpos, "m:", m
  r=ra.match(s,curpos)
  l=r.group(0)
  #print "l:", l
  if not l: print "Bad character"; break # bad character
  c=l[0]
  if c=='.' or c=='>': satom=-1
  elif c in '0123456789%': 
   c=int(l) if c!='%' else int(l[1:])
   if c in bheads.keys(): m.addbond(satom,bheads[c],cbtype); del bheads[c]
   else: bheads[c]=satom
  else:
   if c=='[': c=l[1:-1]; a=m.addatom(0,0,c)
   elif l!='*' and l.lower()==l: a=m.addatom(0,0,l.upper()); m.atoms[a].arom=True
   else: a=m.addatom(0,0,l)

  
   
  if satom>=0: m.addbond(satom,a,cbtype)
  satom=a
  curpos=r.end()
  r=rb.match(s,curpos)

  for c in r.group(1):
   if c=='(': astack.append(satom)
   else: satom=astack.pop()
 
  curpos=r.end()
  if curpos>=len(s): break  
  
  c=s[curpos]; curpos+=1
  if c=='-' or c=='\\' or c=='/': cbtype=1
  elif c=='=' or c==':': cbtype=2
  elif c=='#': cbtype=3
  else:
   cbtype=1
   curpos-=1

 m.rebuild_bond_lists()
 return m

bzarom=loadSmiles('*1=**=**=*1')
furarom=loadSmiles('[O,S,N]1*=**=*1')
naphhalf=loadSmiles('c1c*=**=*1')

def markaromatic(m,clear=False):
 if clear:
  for a in m.atoms: a.arom=0
 f=matchmolecules(bzarom,m)
 while f:
  for a in bzarom.atoms: a.mappedto.arom=1
  f=matchmolecules(bzarom,m,"N")
 
 f=matchmolecules(furarom,m)
 while f:
  for a in furarom.atoms: a.mappedto.arom=1
  f=matchmolecules(furarom,m,"N")

 f=matchmolecules(naphhalf,m,"A")
 while f:
  for a in naphhalf.atoms: a.mappedto.arom=1
  f=matchmolecules(naphhalf,m,"AN")
  
def getmolbbox(m):
 if not len(m.atoms): return [0,0,0,0]
 minx=m.atoms[0].x
 maxx=minx 
 miny=m.atoms[0].y
 maxy=miny;
 
 for a in m.atoms:
  if a.x<minx: minx=a.x
  if a.x>maxx: maxx=a.x
  if a.y<miny: miny=a.y
  if a.y>maxy: maxy=a.y
 return [minx,miny,maxx,maxy]

def c2num(c):
 if c>94: return c-58
 elif c>63: return c-54
 else: return c-48

def num2c(c):
 if c<10: c+=48
 elif c<37: c+=54
 else: c+=58
 return chr(c)

def getpackedcoords(m):
 s=''
 bb=getmolbbox(m)
 for i in range(len(m.atoms)):
  n=(m.atoms[i].x-bb[0])
  s+=num2c(n/65)
  s+=num2c(n%65)
  n=(m.atoms[i].y-bb[1])
  s+=num2c(n/65)
  s+=num2c(n%65)
 return s

def loadpackedcoords(m,s):
 for i,a in enumerate(m.atoms):
  a.x=c2num(ord(s[i*4]))*65+c2num(ord(s[i*4+1]))
  a.y=c2num(ord(s[i*4+2]))*65+c2num(ord(s[i*4+3]))

def packedstringtomol(s):
 a=s.split(";")
 if len(a)!=3: return None
 aline=a[0].split(" ");
 mol=Molecule()
 for l in aline: mol.addatom(0,0,l)
 bline=a[1]
 for i in range(0,len(bline),4): 
  na1=c2num(ord(bline[i]));na2=c2num(ord(bline[i+2]))
  m=c2num(ord(bline[i+1]))
  na1=(m>>4)*65+na1; na2=(m & 7)*65+na2
  mol.addbond(na1,na2,c2num(ord(bline[i+3])))
 loadpackedcoords(mol,a[2])
 mol.update_mol_data()
 return mol

def moltopackedstring(mol):
 mol.rebuild_bond_lists()
 s=""
 for i,a in enumerate(mol.atoms): 
  if i: s+=" "
  s+=a.label
 
 s+=";";
 for i in range(len(mol.bonds)): 
  b=mol.bonds[i];
  s+=num2c(b.atom1.num%65)
  s+=num2c((b.atom1.num/65>>0)*8+(b.atom2.num/65>>0))
  s+=num2c(b.atom2.num%65)
  s+=num2c(b.btype)
 return s+';'+getpackedcoords(mol)


def loadBKC():
 pass

def mol2BKC():
 pass

def angled_line(length,angle): # dx ,dx of a line with length and angle
 return (length * math.cos(math.radians(angle)),length * math.sin(math.radians(angle)))

def shiftAtoms(al,dx,dy):
 for a in al:
  a.x+=dx
  a.y+=dy

def rotateAtoms(al,angle,cx=0,cy=0):
 sina=math.cos(math.radians(angle))
 cosaa=math.sin(math.radians(angle))
 for a in al:
  x=a.x-cx
  y=a.y-cy
  a.x=x * cosa - y * sina + cx
  a.y=y * sina + y * cosa + cy

def saveAtomPositions(al):
 xy=[]
 for a in al:
  xy.append([a.x,a.y])
 return xy

def setAtomPositions(al,xy,dxy=(0,0)):
 for i in range(len(al)):
  al[i].x=xy[i][0]+dxy[0]
  al[i].y=xy[i][1]+dxy[1]
 
class moleculeBuilder:
 def __init__(self,bondLength=40):
  self.m=Molecule()
  self.fragments=[]
  self.cfragment=-1
  self.bl=bondLength
 
 def startFragment(self):
  if self.cfragment>=0: self.endFragment()
  self.fragments.append([len(self.m.atoms),len(self.m.bonds),len(self.m.atoms),len(self.m.bonds)])
  self.cfragment=len(self.fragments)-1
  return self.cfragment

 def getMolecule(self):
  return self.m
  
 def endFragment(self):
  if self.cfragment>=0:
   cf=self.fragments[self.cfragment]
   cf[2]=len(self.m.atoms)
   cf[3]=len(self.m.bonds)
  self.cfragment=-1

 def getFragmentAtoms(self,n):
  try:
   cf=self.fragments[n]
   return [self.m.atoms[i] for i in range(cf[0],cf[2])]
  except:
   return []
   
 def getFragmentBonds(self,n):
  cf=self.fragments[n]
  return [self.m.bonds[i] for i in range(cf[1],cf[3])]

 def getXY(self,fromA, angle=None, dxy=(0,0), bondLength=None): # fromA = reference atom #, -1 if last added, or (x,y) coordinates; angle if bondLength away at this angle, dxy if shifted by its (x,y) 
  bl=bondLength if bondLength else self.bl
  if type(fromA) is int:
   if fromA==-1: fromA=len(self.m.atoms)-1
   x0,y0=(self.m.atoms[fromA].x,self.m.atoms[fromA].y)
  else:
   x0,y0=fromA
  if angle is None:
   dx,dy=0,0
  else:
   dx,dy=angled_line(bl,angle)
  return (x0+dx+dxy[0],y0+dy+dxy[1])

 def newAtom(self, label, fromA, angle=None, dxy=(0,0), bondLength=None):   
  x,y=self.getXY(fromA,angle,dxy,bondLength)
  return self.m.addatom(x,y,label)

 def getAtom(self, n):   
  return self.m.atoms[n]

 def getAtomXY(self, n):   
  return (self.m.atoms[n].x,self.m.atoms[n].y)

 def setAtomXY(self, n, xy):   
  a=self.m.atoms[n]
  a.x,a.y=xy

 def getFragmentXYs(self, n):   
  return saveAtomPositions(self.getFragmentAtoms(n))

 def setFragmentXYs(self, n,XYs,dxy=(0,0)):   
  setAtomPositions(self.getFragmentAtoms(n),XYs,dxy)

 def newBond(self,na1,na2,btype=1): # bond length away from fromA, with angle, and shifted by dxy  
  if na1==-1 and na2==-1:
   na1=len(self.m.atoms)-2
   na2=len(self.m.atoms)-1
  elif na1==-1:
   na1=len(self.m.atoms)-1
  elif na2==-1:
   na2=len(self.m.atoms)-1
  
  if btype=="up" or btype=="wedge":
   btype=5
  if btype=="down" or btype=="dash":
   btype=9
   
  return self.m.addbond(na1,na2,btype)

 def getBond(self,n):
  return self.m.bonds[n]
