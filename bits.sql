-- aggregate functions

Create AGGREGATE andsum (
    BASETYPE = bit,
    SFUNC = bitand,
    STYPE = bit
);

Create AGGREGATE orsum (
    BASETYPE = bit,
    SFUNC = bitor,
    STYPE = bit
);

-- bit operator functions

Create or Replace FUNCTION contains(bit, bit)
 Returns boolean AS '
Select $2 = ($1 & $2);
' LANGUAGE sql IMMUTABLE;
Comment On FUNCTION contains(bit, bit)
 Is 'sql equivalent of bit_contains for fixed lengths';

Create or Replace FUNCTION bit_density(bit)
 Returns real AS '
Select (nbits_set($1)::real / length($1))::real;
' LANGUAGE sql IMMUTABLE;
Comment On FUNCTION bit_density(bit)
 Is 'useful for analysis of fingerprints';

Create or Replace FUNCTION hex(bit)
 Returns text AS $EOSQL$
Select 'X''' || substr(encode(bit_send($1),'hex'),9) ||
  '''::bit(' || length($1) || ')';
$EOSQL$ LANGUAGE sql IMMUTABLE;
Comment On FUNCTION hex(bit)
 Is 'hexadecimal representation of bit string';
