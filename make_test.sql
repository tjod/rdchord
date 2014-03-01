\timing
/*
Create Table rdtest.structure (
 id int,
 rdmol bytea,
 fpbits integer[],
 fp bit varying,
 mol rdx.mol
);
*/

Drop Table rdtest.structure;
Create Table rdtest.structure As
 Select structure_id As id, mol_id as mid,
 rdkit.molfile_to_rdmol(molblock) As rdmol,
 rdx.mol_from_ctab(molblock::cstring,true) mol
 from chemistry.structure_mol
 Where import_id < 3
 Limit 100000
;

Create Index molidx On rdtest.structure using gist(mol);

Alter Table rdtest.structure Add Column fp Bit Varying;
Alter Table rdtest.structure Add Column fpbits Integer[];
Update rdtest.structure Set
 fp=rdkit.fp(rdmol,512), fpbits=rdkit.fpbits(rdmol,512);
Create Index fpidx On rdtest.structure using gin(fpbits);
