import cv2
import mol
import numpy as np

class cvCanvas:
 w=0
 h=0
 scale=2
 xoff=0
 yoff=0
 lab_rad=9
 img=None
 def calcmolfit(self,m):
  bb=mol.getmolbbox(m)
  dx=20
  dy=20
  self.xoffs=-bb[0]*self.scale+dx
  self.yoffs=bb[3]*self.scale+dy
  #if (m.largefont) { m.xoffs+=20; m.yoffs+=20; }
  self.w=(bb[2]-bb[0])*self.scale+2*dx
  self.h=(bb[3]-bb[1])*self.scale+2*dy
  
 def createImage(self):
  self.img=np.zeros((self.h,self.w,3), np.uint8)
  self.img[:]=(255,255,255)
  
 def XYtoatom(self,x,y):
  return [x/self.scale-self.xoffs,-y/self.scale+self.yoffs]

 def atomtoXY(self,x,y):
  return [x*self.scale+self.xoffs,-y*self.scale+self.yoffs]

 def getatomXY(self,a,toXY=None): # get of for atom, outside the atom label if direction is specified 
  return self.atomtoXY(a.x,a.y)

 def addline(self,p1,p2,color=(0,0,0)):
  cv2.line(self.img,(int(p1[0]),int(p1[1])),(int(p2[0]),int(p2[1])),color,1)  

 def addellipse(self,p,w,h,angle=0,color=(0,0,0)):
  cv2.ellipse(img,(p[0],p[1]),(w,h),angle,0,360,color,-1) # 

 def adddashedline(self,p1,p2,color=(0,0,0)):
  l=((p2[0]-p1[0])**2+(p2[1]-p1[1])**2)**0.5
  extra=(l-4)%8 # dash len=2 4=2*
  extra/=2
  dx=p2[0]-p1[0]
  dy=p2[1]-p1[1]
  i=0
  while i<l:    
   d=(extra+i)/l
   d1=(extra+i+2)/l
   cvc.addline(padd(p1,d*dx,d*dy),padd(p1,d1*dx,d1*dy),color)
   i+=8
   
 def addtriangle(self,p1,p2,p3):
  points = np.int32(np.array([(p1[0],p1[1]),(p2[0],p2[1]),(p3[0],p3[1])]))
  
  cv2.fillPoly(self.img, [points], (0,0,0))

 def addtext(self,p,s,ttype=0):
#  font = cv2.FONT_HERSHEY_COMPLEX_SMALL
  font = cv2.FONT_HERSHEY_PLAIN
  ldim,baseline=cv2.getTextSize(s,font,2,1)
  print ldim,baseline
  cv2.putText(self.img,s, (p[0]-ldim[0]/2,p[1]+ldim[1]/2), font, 2, (0,0,0), 1, cv2.LINE_AA)
  return max(ldim[0],ldim[1])/2


def padd(p,dx,dy):
 return (p[0]+dx,p[1]+dy)

def cv_drawmol(cvc,m):
 if not len(m.atoms): return
 m.update_mol_data()

 showCs=0

 for i,a in enumerate(m.atoms):
  if not a.label: a.label="C"
  p1=cvc.atomtoXY(a.x,a.y);
  
  if showCs:
   l=a.label
  else:
   l=a.label if (a.label!="C") or (not len(a.bonds)) else ""  
   if l: 
    a.lab_rad=cvc.addtext(p1,l)

 
 for i,b in enumerate(m.bonds):
  a1=b.atom1
  a2=b.atom2 

  dx=b.nx
  dy=-b.ny
  p1=cvc.atomtoXY(a1.x,a1.y)
  p2=cvc.atomtoXY(a2.x,a2.y)
  
  if a1.label and (showCs or a1.label!='C'):
   lab_rad=a1.lab_rad if a1.lab_rad else cvc.lab_rad
   p1[0]+=dx*lab_rad
   p1[1]+=dy*lab_rad
  if a2.label and (showCs or a2.label!='C'):
   lab_rad=a2.lab_rad if a2.lab_rad else cvc.lab_rad   
   p2[0]-=dx*lab_rad
   p2[1]-=dy*lab_rad
   
  bl=""
  if b.btype==1:
   cvc.addline(p1,p2,bl)
  elif b.btype==2:
   cvc.addline(padd(p1,3*dy,-3*dx),padd(p2,3*dy,-3*dx),bl)
   cvc.addline(padd(p1,-3*dy,3*dx),padd(p2,-3*dy,3*dx),bl)
  elif b.btype==3:
   cvc.addline(p1,p2,bl);
   cvc.addline(padd(p1,4*dy,-4*dx),padd(p2,4*dy,-4*dx),bl)
   cvc.addline(padd(p1,-4*dy,4*dx),padd(p2,-4*dy,4*dx),bl)
  elif b.btype==5 or b.btype==9:
   p3=padd(p2,-3*dy,3*dx)
   p2=padd(p2,3*dy,-3*dx)
     
   if (b.btype==5):
    cvc.addtriangle(p1,p2,p3)
   else:
    a=int(b.blen*cvc.scale // 6)   
    dx=(p2[0]-p1[0])/a
    dy=(p2[1]-p1[1])/a
    dx1=(p3[0]-p1[0])/a
    dy1=(p3[1]-p1[1])/a
    p2=[p1[0],p1[1]]
    p3=[p1[0],p1[1]]
    for j in range(a):
     cvc.addline(p2,p3)
     p2[0]+=dx
     p2[1]+=dy
     p3[0]+=dx1
     p3[1]+=dy1
        
 return cvc.img

# "180 -145 C 220 -145 Br 160 -110.35898384862247 H 145.35898384862247 -165 H 160 -179.64101615137753 H 79 -142 H 119 -142 O;1 2 1 1 3 1 1 4 9 1 5 5 6 7 1"
# 

m=mol.packedstringtomol("C Cl H H H H O;00110021003900455061;1Z0X2A0X1F14110D1F00000_0b0_")
cvc=cvCanvas()
cvc.calcmolfit(m)
cvc.createImage()
cv_drawmol(cvc,m)
cvc.adddashedline(cvc.getatomXY(m.atoms[0]),cvc.getatomXY(m.atoms[6]),(200,200,200))

m.bonds[0].btype=10


cv2.imshow('image',cvc.img)
cv2.waitKey(0)
cv2.destroyAllWindows()
