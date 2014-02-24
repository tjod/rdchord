import sys
import rdchord
from rdkit import Chem

def rdmol(imol, mol):
  if mol:
    print str(imol) + "\t" + "\\\\x" + rc.pickle(mol).encode("hex") + "\t" + rc.cansmiles(mol) + "\t" + rc.isosmiles(mol).replace("\\","\\\\") + "\t" + rc.fingerprint(mol,1024,2).ToBitString() + "\t" + mol.GetProp("_Name");
  else:
    print str(imol) + "\t\\N\t\\N\t\\N\t\\N"

def rdmolb(imol, molblock):
  m = Chem.MolFromMolBlock(molblock)
  #print m.GetProp("_Name")
  #for p in mol.GetPropNames():
  #  print p
  print str(imol) + "\t" + molblock.replace("\r","").replace("\n", "\\n");

def rdprop(imol, mol):
  if mol:
    for p in mol.GetPropNames():
      print str(imol) + "\t" + p + "\t" +  mol.GetProp(p)
  else:
    print str(imol) + "\t\\N\t\\N"

def newCopy(properties):
  if properties:
    print "Copy properties (id, name, value) From Stdin;"
  elif molblocks:
    print "Copy molblocks (id, molblock) From Stdin;"
  else:
    print "Copy structure (id, rdmol, cansmiles, isosmiles, fp, name) From Stdin;"

properties = False
molblocks = False
nchunk = 1000
if len(sys.argv) > 2:
  schema = sys.argv[1];
  sys.argv.pop(1);
  molfile = sys.argv[1]
  sys.argv.pop(1);
  while len(sys.argv) > 1:
    print sys.argv
    if sys.argv[1] == "--properties":
      properties = True
      nchunk = 10000
    elif sys.argv[1] == "--molblocks":
      molblocks = True
      nchunk = 10000
    sys.argv.pop(1)
else:
  print >> sys.stderr, "usage: python sdfloader.py schema sdfile [--properties][--molblocks]\n; writes to stdout intended for psql input"
  exit(0)

print """\\timing
--Drop Schema If Exists %s Cascade;
Create Schema %s;
Grant Usage On Schema %s To Public;
Set search_path=%s;""" % (schema, schema, schema, schema)

print """
Create Table structure (
 id Integer Primary Key,
 name Text,
 rdmol Bytea,
 cansmiles Text,
 isosmiles Text,
 fp Bit Varying,
 fpbits Int
 );
Create Table properties (
 id Integer,
 name Text,
 value Text);
Create Table molblocks (
 id Integer,
 molblock Text);
"""

newCopy(properties)

GD = dict()
rc = rdchord.rdchord(GD)
molblock = list()
imol = 0
stream = Chem.SDMolSupplier(molfile)
for idx,mol in enumerate(stream):
    imol += 1
    if properties:
      rdprop(imol, mol)
    elif molblocks:
      rdmolb(imol, stream.GetItemText(idx))
    else:
      rdmol(imol, mol)
    molblock = list()
    if 0 == (imol % nchunk):
      print >> sys.stderr, imol,
      print "\\."
      newCopy(properties)
