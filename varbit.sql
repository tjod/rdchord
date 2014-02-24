CREATE or REPLACE FUNCTION nbits_set(bit)
 RETURNS integer AS 'varbit', 'nbits_set'
 LANGUAGE c IMMUTABLE STRICT;
COMMENT on FUNCTION nbits_set(bit)
 IS 'number of bits set';

CREATE or REPLACE FUNCTION bit_set(bit, integer)
 RETURNS bit AS 'varbit', 'bit_set'
 LANGUAGE c IMMUTABLE STRICT;
COMMENT on FUNCTION bit_set(bit, integer)
 IS 'set Nth bit';

CREATE or REPLACE FUNCTION isbit_set(bit, integer)
 RETURNS boolean AS 'varbit', 'isbit_set'
 LANGUAGE c IMMUTABLE STRICT;
COMMENT on FUNCTION isbit_set(bit, integer)
 IS 'is Nth bit set?';

CREATE or REPLACE FUNCTION bit_contains(bit, bit)
 RETURNS boolean AS 'varbit', 'bit_contains'
 LANGUAGE c IMMUTABLE STRICT;
COMMENT on FUNCTION bit_contains(bit, bit)
 IS 'does first bitstring have all bits of second set?';

