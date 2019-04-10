import rdchord
from plpy import plpy

GD = {}
rd = rdchord.rdchord(GD)

mol = rd.parse_smi("c1ccccc1C(=O)NC")
print rd.cansmiles(mol)
print rd.smilesBuffer(None)
plpy().notice("hello")