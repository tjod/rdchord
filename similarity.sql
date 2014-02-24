-- Similarity functions
Create or Replace FUNCTION tanimoto(bit, bit)
 Returns real AS '
Select nbits_set($1 & $2)::real /
      (nbits_set($1 & ~$2) + nbits_set($2 & ~$1) +
       nbits_set($1 & $2))::real;
' LANGUAGE sql IMMUTABLE;
Comment On FUNCTION tanimoto(bit, bit)
 Is 'standard tanimoto similarity of fingerprints';

/*
Create or Replace FUNCTION tani2(bit, bit)
 Returns real AS '
Select nbits_set($1 & $2)::real /
      (nbits_set($1) + nbits_set($2) - nbits_set($1 & $2))::real;
' LANGUAGE sql IMMUTABLE;
Comment On FUNCTION tani2(bit, bit)
 Is 'According to the Daylight Theory Manual.';

Create or Replace FUNCTION tani3(bit, bit)
 Returns real AS '
Select nbits_set($1 & $2)::real /
       nbits_set($1 | $2)::real;
' LANGUAGE sql IMMUTABLE;
Comment On FUNCTION tani3(bit, bit)
 Is 'According to Sayle';
*/

Create or Replace FUNCTION euclid(bit, bit)
 Returns real AS '
Select sqrt((nbits_set($1 & $2) + nbits_set(~$1 & ~$2))::real /
  length($1))::real;
' LANGUAGE sql IMMUTABLE;
Comment On FUNCTION euclid(bit, bit)
 Is 'standard euclid similarity of fingerprints';

Create or Replace FUNCTION hamming(bit, bit)
 Returns real AS '
Select ((nbits_set($1 & ~$2) + nbits_set(~$1 & $2))::real /
   length($1))::real;
' LANGUAGE sql IMMUTABLE;
Comment On FUNCTION hamming(bit, bit)
 Is 'standard hamming similarity of fingerprints';

Create or Replace FUNCTION tversky(bit, bit, real, real)
 Returns real AS '
Select nbits_set($1 & $2)::real /
  ($3*nbits_set($1 & ~$2) + $4*nbits_set($2 & ~$1) +
      nbits_set($1 & $2))::real;
' LANGUAGE sql IMMUTABLE;
Comment On FUNCTION tversky(bit, bit, real, real)
 Is 'standard tversky similarity of fingerprints';

Create or Replace FUNCTION similarity(bit, bit)
 Returns real AS '
Select tversky($1, $2, 0.5, 0.5);
' LANGUAGE sql IMMUTABLE;
Comment On FUNCTION similarity(bit, bit)
 Is 'standard dice similarity of fingerprints';
                                                                                
Create or Replace FUNCTION similarity(Bytea, Bytea)
 Returns real AS '
Select tversky(public166fp($1), public166fp($2), 0.5, 0.5);
' LANGUAGE sql IMMUTABLE;
Comment On FUNCTION similarity(Bytea, Bytea)
 Is 'standard dice similarity of rdmol using public166key';
