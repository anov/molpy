import mol
import math
import cgi
import re

class bitStream:
 def __init__(self,s,p=0):
  self.s=s
  self.p=p
  self.cbits=0
  self.v=0
  
 def getBits(self,n):
  while self.cbits<n:
   self.v+=c2n(self.s[self.p])<<self.cbits
   self.cbits+=6
   self.p+=1
  r=self.v & (0xFFFF>>(16-n))
  self.v>>=n
  self.cbits-=n
  return r
 
 def addBits(self,n,bits):
  self.v+=bits<<self.cbits
  self.cbits+=n
  while self.cbits>=6:
   self.s+=urln2c(self.v & 0x3F)
   self.v>>=6
   self.cbits-=6

 def flush(self):
  if self.cbits: self.s+=n2c(self.v)

 def isEnd(self):
  return (not self.cbits) and (len(self.s)<=self.p)


def bsLength(bits):
 return (bits//6)+(bits%6!=0)

def bsBitsFor(v):
 i=1 
 while v>1<<i: i+=1 
 return i

def c2n(c):
 if c>='0' and c<='9':
  d=ord(c)-ord('0')
 elif c>='A' and c<='Z':
  d=ord(c)-ord('A')+10
 elif c>='a' and c<='z':
  d=ord(c)-ord('a')+26+10
 elif c=='-':
  d=26+26+10
 elif c=='_':
  d=26+26+10+1
 elif c=='~':
  d=26+26+10+2
 elif c=='!':
  d=26+26+10+3
 elif c=='*':
  d=26+26+10+4
 elif c=='.':
  d=26+26+10+5
 elif c=='(':
  d=26+26+10+6
 elif c==')':
  d=26+26+10+6

 return d


def loadmurl(s):
# add atom labels
 n=0
 m=mol.Molecule()
 while n<len(s):
  c=s[n]
  n+=1
  if c=='.': break
  if c>='0' and c<='9': # RLE
   if s[n]>='0' and s[n]<='9':
	c+=s[n]
	n+=1
   c=int(c)
   l=m.atoms[a].label
   for i in range(c-1):
    m.addatom(0,0,l)
   continue
  
  if c in 'BCDFHIKNOPSX':
   l=c
  else:
   try:
    l='LnMAscaTrfuZGebUhdgtmp'.index(c)
    l=['Li',"Na",'Mg','Al','Si','Cl',"Ca","Ti","Cr","Fe","Cu","Zn","Ge","Se",'Br',"Ru","Rh",'Pd',"Ag",'Sn',"Hg",'Pb'][l];
   except:
    if c in 'wxy':
     i=ord(c)-118
     l=s[n:n+i]
     n+=i
    elif c=='z':
     i=c2n(s[n])
     l=s[n+1,n+1+i]
     n+=i+1
    else:
	 l=c 
  a=m.addatom(0,0,l)

 bs=bitStream(s,n);
 xbits=bs.getBits(3)+6;
 ybits=bs.getBits(3)+6;
  
 if len(s)<n+bsLength((xbits+ybits)*len(m.atoms)): a=1 # not enough data
  
 for i in range(len(m.atoms)):
  a=m.atoms[i]  
  a.x=bs.getBits(xbits) 
  a.y=bs.getBits(ybits)

 n+=bsLength((xbits+ybits)*len(m.atoms)+6)

 bs=bitStream(s,n)
 bbits=bs.getBits(3)
 abits=bs.getBits(3)+2
 
 while not bs.isEnd(): # bonds
  c=bs.getBits(3)
  if c==0: break
  if c<4: bt=c
  elif c==4: bt=5
  elif c==5: bt=9
  m.addbond(bs.getBits(abits),bs.getBits(abits),bt)

 m.update_mol_data()
 return m


def loadsmurl(s):
 m=mol.Molecule()
 ra=re.compile("(Cl|Br|[BCNOSPIFbcnosp\*\.\>\d]| |)")
 rb=re.compile("([\)\(]+|)")
 astack=[]
 bheads={}
 satom=-1 
 
 curpos=0
 
 while curpos<len(s):
  #print "pos:",curpos, "m:", m
  r=ra.match(s,curpos)
  l=r.group(0)
  #print "l:", l
  if not l: print "Bad character"; break # bad character
  if l==' ': break
  c=l[0]
  if c=='.' or c=='>': satom=-1
  elif c in '0123456789%':
  # c=(c!='%')?Number(l):Number(l.slice(1)); 
   c=int(l)
   if c in bheads.keys(): m.addbond(satom,bheads[c],cbtype); del bheads[c]
   else: bheads[c]=satom
  else:
   a=m.addatom(0,0,l)
  
  if satom>=0: m.addbond(satom,a,1)
  satom=a
  curpos=r.end()
  r=rb.match(s,curpos)

  for c in r.group(1):
   if c=='(': astack.append(satom)
   else: satom=astack.pop()
 
  curpos=r.end()
  if curpos>=len(s): break  

 curpos+=1 # skip space
 
 for i in range(len(m.bonds)):
  c=c2n(s[curpos])
  negx=c & 16
  negy=c & 32
  nocoords=c & 8
  
  c=c // 8
  if c<3: b.btype=c+1
  elif c==3: b.btype=5
  elif c==4: b.btype=9

  if nocoords: continue
  b.atom2.x=b.atom1.x+c2n(s[curpos+1])
  b.atom2.y=b.atom1.y+c2n(s[curpos+2])
  curpos+=3  

 for i in range(len(m.atoms)):
  a=m.atoms[i]
  if a.label=='s': pass
  
 m.rebuild_bond_lists()
 return m

class svgCanvas:
 w=0
 h=0
 scale=1
 xoff=0
 yoff=0
 lab_rad=6
 def calcmolfit(self,m):
  self.scale=1
  bb=mol.getmolbbox(m)
  dx=20
  dy=20
  self.xoffs=-bb[0]+dx
  self.yoffs=bb[3]+dy
  #if (m.largefont) { m.xoffs+=20; m.yoffs+=20; }
  self.w=bb[2]-bb[0]+2*dx
  self.h=bb[3]-bb[1]+2*dy
  
 def writeheader(self):
  print '<svg xmlns="http://www.w3.org/2000/svg" height="{0}" width="{1}" font-size="14">'.format(self.h,self.w)

 def XYtoatom(self,x,y):
  return [x/self.scale-self.xoffs,-y/self.scale+self.yoffs]

 def atomtoXY(self,x,y):
  return [x*self.scale+self.xoffs,-y*self.scale+self.yoffs]

 def addline(self,p1,p2,ltype=0):
  print '<line stroke-width="1" style="stroke:black" x1="{0}" y1="{1}" x2="{2}" y2="{3}"></line>'.format(p1[0],p1[1],p2[0],p2[1])
  
  #'<line stroke-width="{0}" style="stroke:black" y2="110" x2="70" y1="98" x1="50.14"></line>'

 def addtriangle(self,p1,p2,p3):
  print '<polygon style="fill:black;stroke:black;stroke-width:1" points="{0},{1} {2},{3} {4},{5}"></polygon>'.format(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1])

 def addtext(self,p,s,ttype=0):
  print '<text dy="5.5" pointer-events="none" text-anchor="middle" fill="black" x="{0}" y="{1}">{2}</text>'.format(p[0],p[1],s)

 def close(self):
  print '</svg>'

#<text dy="5.5" pointer-events="none" text-anchor="middle" fill="black" y="50" x="70">O</text>
#<tspan baseline-shift=sub>S</tspan>

 #if (m.largefont) txt.setAttribute("font-size","30");
 #if (ttype=="P") txt.setAttribute("font-weight","bold"); else
 #if (ttype=="E") txt.setAttribute("fill","blue"); else
 
 #bb=txt.getBBox();
 #txt.setAttribute("dy",xy[1]-bb.y-(bb.height/2));

def Hshift(m,na):
 a=m.atoms[na]
 if len(a.bonds)==0: return 1 
 b=a.bonds[0]
 dx=b.atom2.x-b.atom1.x
 if a!=b.atom1: dx=-dx

 for i in range(1,len(a.bonds)):
  b=a.bonds[i]
  dx1=b.atom2.x-b.atom1.x
  if a!=b.atom1: dx1=-dx1
  dx+=dx1
  
 return 1 if dx<0 else -1

def padd(p,dx,dy):
 return (p[0]+dx,p[1]+dy)

def tosvg(m):
 svg=svgCanvas()
 svg.calcmolfit(m)
 svg.writeheader()
 if not len(m.atoms): return
 m.update_mol_data()
 
 for i,b in enumerate(m.bonds):
  a1=b.atom1
  a2=b.atom2 

  dx=b.nx
  dy=-b.ny
  p1=svg.atomtoXY(a1.x,a1.y)
  p2=svg.atomtoXY(a2.x,a2.y)
  
  showCs=0
  if a1.label and (showCs or a1.label!='C'):
   p1[0]+=dx*svg.lab_rad
   p1[1]+=dy*svg.lab_rad;
  if a2.label and (showCs or a2.label!='C'):
   p2[0]-=dx*svg.lab_rad
   p2[1]-=dy*svg.lab_rad
   
  bl=""
  if b.btype==1:
   svg.addline(p1,p2,bl)
  elif b.btype==2:
   svg.addline(padd(p1,3*dy,-3*dx),padd(p2,3*dy,-3*dx),bl)
   svg.addline(padd(p1,-3*dy,3*dx),padd(p2,-3*dy,3*dx),bl)
  elif b.btype==3:
   svg.addline(p1,p2,bl);
   svg.addline(padd(p1,4*dy,-4*dx),padd(p2,4*dy,-4*dx),bl)
   svg.addline(padd(p1,-4*dy,4*dx),padd(p2,-4*dy,4*dx),bl)
  elif b.btype==5 or b.btype==9:
   p3=padd(p2,-3*dy,3*dx)
   p2=padd(p2,3*dy,-3*dx)
     
   if (b.btype==5):
    svg.addtriangle(p1,p2,p3)
   else:
    a=int(b.blen*svg.scale // 6)   
    dx=(p2[0]-p1[0])/a
    dy=(p2[1]-p1[1])/a
    dx1=(p3[0]-p1[0])/a
    dy1=(p3[1]-p1[1])/a
    p2=[p1[0],p1[1]]
    p3=[p1[0],p1[1]]
    for j in range(a):
     svg.addline(p2,p3)
     p2[0]+=dx
     p2[1]+=dy
     p3[0]+=dx1
     p3[1]+=dy1
 
 for i,a in enumerate(m.atoms):
  if not a.label: a.label="C"
  p1=svg.atomtoXY(a.x,a.y);
  
  if showCs:
   l=a.label
  else:
   l=a.label if (a.label!="C") or (not len(a.bonds)) else ""  
   if l: svg.addtext(p1,l)
  
  l=a.needsHs()
  if l:
   if l==1:
    l="H"
   elif l==2:
	l="H2"
   elif l==3:
	l="H3"
   p1[0]+=a.vtext[0].getBBox().width*Hshift(mol,i);
   svg.addtext(p1,l)
   
 svg.close()

#m=loadsmurl('C(C)C 0AZ1ZA2NN')
#print m
#m=loadmurl(d="bxBzC10c.y1U09100p0710071P0U091U000j0p0j0P0M1Z1j091M1Z171P000y207UA06U807UC19U1U4U5U2U0UBU9")
#tosvg(m)

print "Content-type: image/svg\n"

#form = cgi.FieldStorage()
#d=form.getvalue("murl")

d="C10.10Qdz0u93_iscydRp0q4uy3J1a0J6AKboYbHkSH2COcW0"


if d:
 m=loadmurl(d)
 tosvg(m)
 exit(0)

d=form.getvalue("pstr")
#d="C C C C;001110211031;0O0g0O0E00000m00"
if d:
 m=mol.packedstringtomol(d)
 tosvg(m)
 exit(0)

print '<svg xmlns="http://www.w3.org/2000/svg" height="14" width="100" font-size="14">'
print '<text fill="red" x="0" y="14">No structure</text>'
print '</svg>'

