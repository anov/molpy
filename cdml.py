import xml.etree.ElementTree as ET
import mol

def coord(v):
 return float(v[:-2])

def btype(v):
 # w1  "h1"
 if v=="n2": return 2
 if v=="n3": return 3
 return 1


ns = {'c':'http://www.freesoftware.fsf.org/bkchem/cdml'}

#tree = ET.parse('sub.cdml')
tree = ET.parse('/home/chip/chemfun/samples/r1.cdml')

root = tree.getroot()

for s in root.findall('c:standard',ns):
 #print s.attrib
 bset=s.find('c:bond',ns).attrib  #length: 0.7 cm
 blen=coord(bset['length'])
 
for m in root.findall('c:molecule',ns):
 cm=mol.Molecule()
 ids={}
 for a in m.findall('c:atom',ns):
  xy=a.find('c:point',ns).attrib  #x,y  
  ids[a.attrib["id"]]=cm.addatom(coord(xy["x"]),coord(xy["y"]),a.attrib["name"])
  #print a.find('c:mark',ns)
 
 for a in m.findall('c:group',ns):
  xy=a.find('c:point',ns).attrib  #x,y  
  ids[a.attrib["id"]]=cm.addatom(coord(xy["x"]),coord(xy["y"]),a.attrib["name"])
  #print a.find('c:mark',ns)
  
 for b in m.findall('c:bond',ns):  
  cm.addbond(ids[b.attrib["start"]],ids[b.attrib["end"]],btype(b.attrib["type"]))
 
 for a in cm.atoms:
  a.x=int(a.x/0.7*40);
  a.y=int(-a.y/0.7*40);
 print mol.moltostring(cm)

arrows=[]
for a in root.findall('c:arrow',ns):
 #print a.attrib
 xys=a.findall('c:point',ns)
 xy1=xys[0].attrib
 xy2=xys[1].attrib
 arrows.append((coord(xy1["x"]),coord(xy1["y"]),coord(xy2["x"]),coord(xy2["y"])))

#print arrows

texts=[] 
for t in root.findall('c:text',ns):
 xy=t.find('c:point',ns).attrib  #x,y
 text=t.find('c:ftext',ns).text
 texts.append((coord(xy["x"]),coord(xy["y"]),text))

#print texts

marks=[]
for a in root.findall('c:plus',ns):
 #print a.attrib
 xy=a.find('c:point',ns).attrib  #x,y
 marks.append((coord(xy["x"]),coord(xy["y"])))
 
#print marks
