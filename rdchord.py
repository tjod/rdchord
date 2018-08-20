from rdkit.Chem import Mol,MolFromSmiles,MolToSmiles,MolFromMolBlock,MolToMolBlock,MolFromSmarts,Kekulize,SDMolSupplier,MolToSmarts,GetFormalCharge
import plpy

class rdchord:

  def __init__(self, Global):
    """use Global dictionary supplied by plpython to buffer parsed
       smiles and smarts for quicker lookup when reused.
    """
    # buffering smiles and smarts can eat up lots of memory,
    #  but speeds things up considerably
    #  certainly if an entire table can acutally fit
    # but also if this is the maximum size of a hit list, say from match
    if not Global.has_key("RDKIT"): Global["RDKIT"] = dict()
    self.GRD = Global["RDKIT"]
    self.maxsmi = 1000
    if not self.GRD.has_key("mol"): self.GRD["mol"] = dict()
    self.mol = self.GRD["mol"]
    # pick a reasonable number of smarts patterns you expect to use often,
    # say 166 public keys or even 1000 fragment keys
    self.maxsma = 1000
    if not self.GRD.has_key("pat"): self.GRD["pat"] = dict()
    self.pat = self.GRD["pat"]
    try:
        from rdkit.Chem import MolToInchi,InchiToInchiKey
        self.hasInchi = True
    except ImportError:
        self.hasInchi = False

  def smilesBuffer(self,smi):
    """one, or all keys(smiles) stored in global
    """
    if smi:
      return [self.mol[smi]]
    else:
      return self.mol.keys()
  
  def smartsBuffer(self,sma):
    """one, or all keys(smarts) stored in global
    """
    if sma:
      return [self.pat[sma]]
    else:
      return self.pat.keys()
  
  def parse_smi(self,smi):
    """parse smiles and return Mol after storing in global dict
       or return from global dict"""
    if self.mol.has_key(smi):
      # return copy is slower, but safer?
      #return Mol(self.mol[smi])
      #plpy.notice('found mol for %s' % smi)
      return self.mol[smi]

    newmol = MolFromSmiles(smi)
    if newmol:
      if len(self.mol) < self.maxsmi:
        #plpy.notice('new mol for %s' % smi)
        pass
      else:
        self.mol.popitem()
        #key,psmi = self.mol.popitem()
        #plpy.notice('mol reuse %s for %s' % (key,psmi))
      self.mol[smi] = newmol
      return newmol
    else:
      return None
  
  def parse_sma(self,sma):
    """parse smarts and return SmartsPattern after storing in global dict
       or return from global dict"""
    if self.pat.has_key(sma):
      return self.pat[sma]

    newpat = MolFromSmarts(sma)
    if newpat:
      if len(self.pat) < self.maxsma:
        #plpy.notice('new pat for "%s"' % sma)
        pass
      else:
        self.pat.popitem()
        #key,pat = self.pat.popitem()
        #plpy.notice('pattern reuse %s for %s' % (key,sma))
      self.pat[sma] = newpat
      return newpat
    else:
      #plpy.notice('pattern None')
      return None

  def parse_molblock(self,mb):
    """parse molblock and return mol"""
    #mol = MolFromMolBlock(mb)
    sd = SDMolSupplier()
    sd.SetData(mb)
    mol = sd.next()
    if mol:
      return mol
    else:
      return None

  def attempt_fix(self,mb):
    """try to fix molblock - typically an aromatic N without H"""
    return None

  def add_2Dcoords(self,m):
    """add computed 2D coords to input mol; return conformer number"""
    from rdkit.Chem import rdDepictor
    mcopy=Mol(m.ToBinary())
    rdDepictor.Compute2DCoords(mcopy)
    return mcopy

  def smiles(self,m):
    """make non-canoncial smiles from molecule: same as input smiles?"""
    return MolToSmiles(m, isomericSmiles=False, canonical=False)

  def cansmiles(self,m):
    """make canonical smiles from molecule"""
    return MolToSmiles(m, isomericSmiles=False)

  def isosmiles(self,m):
    """make isomeric smiles from molecule"""
    return MolToSmiles(m, isomericSmiles=True)

  def impsmiles(self,m,keepiso=False):
    """make implicit smiles from molecule"""
    import re
    #smi = MolToSmiles(m, canonical=True, isomericSmiles=True)
    smi = MolToSmiles(m, canonical=False, isomericSmiles=keepiso)
    if keepiso == False:
      smi=smi.replace('([H])','')
      smi=smi.replace('[H]','')
    mol = MolFromSmiles(smi)
    outsmi = ""
    iatom = 0
    atoms = mol.GetAtoms()
    #parts = re.split("([\d\(\)\+\[\]-=#:])",smi);
    parts = re.split("(\[.*?\]|Cl|Br|F|I|B|C|c|N|n|O|o|S|s|P|p)",smi);
    #return [a.GetSymbol() for a in atoms]
    for p in parts:
      if len(p) == 0:
        pass
      elif p.isalpha():
        hcount = atoms[iatom].GetImplicitValence()
        if hcount == 0: outsmi += p
        elif hcount == 1: outsmi += "[%sH]" % p
        else: outsmi += "[%sH%d]" % (p,hcount)
        iatom += 1
      elif p.startswith("["):
        hcount = atoms[iatom].GetNumImplicitHs()
        if hcount == 0: outsmi += p
        elif hcount == 1: outsmi += "[%sH%d]" % (atoms[iatom].GetSymbol(), atoms[iatom].GetFormalCharge())
        else: outsmi += "[%sH%d%+d]" % (atoms[iatom].GetSymbol(), hcount, atoms[iatom].GetFormalCharge())
        iatom += 1
      else:
        outsmi += p
    return outsmi

  def keksmiles(self,m):
    """make Kekule smiles from molecule"""
    mcopy=Mol(m.ToBinary())
    #mcopy=ROMol(m)
    Kekulize(mcopy)
    return MolToSmiles(mcopy, isomericSmiles=True, kekuleSmiles=True)

  def graph(self,m):
    from rdkit.Chem import EditableMol, Atom, rdchem
    hcount = m.GetNumAtoms(False) - m.GetNumAtoms(True)
    # create new molecule using single bonds only
    em = EditableMol(Mol())
    nbridx = [None] * m.GetNumAtoms()
    iatom = 0
    for atom in m.GetAtoms():
      atnum = atom.GetAtomicNum()
      if atnum == 1:
        #if atom.GetMass() > 1: pass
        hcount += 1
      else:
        newatom = Atom(atnum)
        #if atom.GetTotalDegree() == 0: newatom.SetNoImplicit(True) # otherwise [Na]. becomes [NaH].
        #newatom.SetFormalCharge(atom.GetFormalCharge())
        newatom.SetFormalCharge(0)
        em.AddAtom(newatom)
        aidx = atom.GetIdx()
        nbridx[aidx] = iatom
        iatom += 1
        for a2 in atom.GetNeighbors():
          a2idx = nbridx[a2.GetIdx()]
          if a2idx != None:
            em.AddBond(aidx, a2idx, rdchem.BondType.SINGLE)
    cansmi = self.cansmiles(em.GetMol())
    #cansmi = cansmi.replace('+','').replace('-','').replace('[N]','N').replace('[O]','O').replace('[C]','C').replace('[I]','I').replace('[S]','S').replace('[P]','P').replace('[B]','B').replace('[Br]','Br').replace('[Cl]','Cl')
    return "%s%s%d%+d" % (cansmi, '.H',  hcount, GetFormalCharge(m))

  def graph2(self,m):
    from rdkit.Chem import EditableMol, RemoveHs, Atom, rdchem, SanitizeMol, rdmolops
    natoms = m.GetNumAtoms()
    # create new molecule using single bonds only
    em = EditableMol(Mol())
    hcount = 0
    iatom = 0
    for atom in m.GetAtoms():
      atnum = atom.GetAtomicNum()
      hcount += atom.GetTotalNumHs(False)
      newatom = Atom(atnum)
      #newatom.SetFormalCharge(atom.GetFormalCharge())
      em.AddAtom(newatom)
    for bond in m.GetBonds():
      em.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), rdchem.BondType.SINGLE)
    try:
      mol = RemoveHs(em.GetMol())
    except:
      mol = em.GetMol()
    #mol = em.GetMol()
    #SanitizeMol(mol, SanitizeFlags.SANITIZE_ADJUSTHS)
    #Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS

    cansmi = self.cansmiles(mol)
    return "%s%s%d%+d" % (cansmi, '.H',  hcount, GetFormalCharge(m))

  def inchi(self,m):
    """make InChi from molecule"""
    if self.hasInchi:
        from rdkit.Chem import MolToInchi
        return MolToInchi(m)
    else:
        plpy.notice('InChi not available')
        return None

  def inchikey(self,m):
    """make InChi from molecule"""
    if self.hasInchi:
        from rdkit.Chem import MolToInchi,InchiToInchiKey
        return InchiToInchiKey(MolToInchi(m))
    else:
        plpy.notice('InChi not available')
        return None

  def molfile(self,m):
    """make molfile from molecule"""
    return MolToMolBlock(m)

  def to_binary(self,m):
    """return binary representation of mol"""
    return m.ToBinary()

  def from_binary(self,b):
    """reconstitute binary mol"""
    try:
      return Mol(b)
    except:
      return None

  def matches(self,m,p):
    """Tell whether molecule m matches pattern p"""
    if m.HasSubstructMatch(p):
      return True
    else:
      return False

  def get_matches(self,m,p):
    """Tell whether molecule m matches pattern p"""
    matches = m.GetSubstructMatches(p)
    if matches:
      return matches
    else:
      return None

  def fingerprint(self,m,fpsize=1024,bitsperhash=2,tgtDensity=0.3):
    """Compute bit fingerprint of molecule m"""
    from rdkit.Chem import RDKFingerprint
    #from rdkit import DataStructs
    #fp=RDKFingerprint(m, minPath=1, maxPath=7,fpSize=1024, bitsPerHash=2, useHs=False, tgtDensity=0.3)
    return RDKFingerprint(m, minPath=1, maxPath=7, fpSize=fpsize, nBitsPerHash=bitsperhash, tgtDensity=tgtDensity, minSize=fpsize)

  def maccskeys(self,m):
    """Compute maccs public166 key of molecule m"""
    from rdkit.Chem import MACCSkeys
    return MACCSkeys.GenMACCSKeys(m)

  def svg(self, m, width=250, height=250, matchmol=None, kekulize=True, adjust=False, orient=False, highlight=False):
    from rdkit.Chem.rdmolfiles import SDMolSupplier
    from rdkit.Chem import rdDepictor
    from rdkit.Chem import AllChem
    from rdkit.Chem import Kekulize
    from rdkit.Chem.Draw import rdMolDraw2D

    if adjust:
      # compute bounding box
      minWidth = width / 4
      minHeight = height / 4
      xmax = m.GetConformer().GetAtomPosition(0).x
      xmin = xmax
      ymax = m.GetConformer().GetAtomPosition(0).y
      ymin = ymax
      for i in range(0, m.GetNumAtoms()):
          pos = m.GetConformer().GetAtomPosition(i)
          xmax = max(xmax, pos.x)
          xmin = min(xmin, pos.x)
          ymax = max(ymax, pos.y)
          ymin = min(ymin, pos.y)
      # set pixels per Angstrom
      xscale = (xmax-xmin) / 10.0
      yscale = (ymax-ymin) / 10.0
      iwidth =  int(max(minWidth,   width * xscale))
      iheight = int(max(minHeight, height * yscale))
    else:
      iwidth = width
      iheight = height
    drawer = rdMolDraw2D.MolDraw2DSVG(iwidth, iheight)
    mcopy = rdMolDraw2D.PrepareMolForDrawing(m, kekulize=False, wedgeBonds=True, addChiralHs=True)
    matched_atoms=list(m.GetSubstructMatch(matchmol)) if matchmol else []
    if len(matched_atoms) > 0 and orient:
        AllChem.GenerateDepictionMatching2DStructure(mcopy, matchmol)
    if kekulize: Kekulize(mcopy)
    if len(matched_atoms) > 0 and highlight:
        drawer.DrawMolecule(mcopy, highlightAtoms=matched_atoms)
    else:
        drawer.DrawMolecule(mcopy)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
