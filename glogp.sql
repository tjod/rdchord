CREATE TABLE public.glogp (
    smarts character varying,
    train_freq integer,
    coefficient numeric(5,3),
    error numeric(5,3)
);

COPY public.glogp (smarts, train_freq, coefficient, error) FROM stdin;
[$([N+0v3H0D2](=C)(Cl))]	1	\N	\N
[$([S+0v2H0D2](C)(S))]	1	\N	\N
[$([C+0v4H0D2](#C)(c))]	2	\N	\N
[$([C+0v4H0D2](#C)(C))]	4	\N	\N
[$([C+0v4H0D2](#N)(S))]	3	\N	\N
[$([C+0v4H0D3](=C)(c)(c))]	3	\N	\N
[$([C+0v4H0D3](=C)(C)(c))]	3	\N	\N
[$([C+0v4H0D3](=C)(C)(Cl))]	2	\N	\N
[$([C+0v4H0D3](=C)(Cl)(Cl))]	3	\N	\N
[$([c+0v4H0D3](-c)(:c)(:n))]	7	\N	\N
[$([C+0v4H0D3](=C)(C)(N))]	5	\N	\N
[$([c+0v4H0D3](-c)(:c)(:o))]	2	\N	\N
[$([c+0v4H0D3](:c)(:c)(:o))]	6	\N	\N
[$([c+0v4H0D3](:c)(:c)(=O))]	5	\N	\N
[$([C+0v4H0D3](=C)(C)(O))]	3	\N	\N
[$([c+0v4H0D3](-c)(:c)(:s))]	2	\N	\N
[$([c+0v4H0D3](:c)(:c)(:s))]	3	\N	\N
[$([C+0v4H0D3](=C)(C)(S))]	2	\N	\N
[$([C+0v4H0D3](=C)(F)(F))]	1	\N	\N
[$([c+0v4H0D3](:c)(:n)(Br))]	5	\N	\N
[$([c+0v4H0D3](:c)(:n)(F))]	4	\N	\N
[$([c+0v4H0D3](:c)(:n)(I))]	2	\N	\N
[$([c+0v4H0D3](-c)(:n)(:n))]	2	\N	\N
[$([c+0v4H0D3](:c)(:n)(O))]	9	\N	\N
[$([c+0v4H0D3](-c)(:n)(:s))]	1	\N	\N
[$([c+0v4H0D3](:c)(:n)(S))]	2	\N	\N
[$([c+0v4H0D3](:c)(:o)(C))]	1	\N	\N
[$([c+0v4H0D3](:c)(:o)(N))]	2	\N	\N
[$([c+0v4H0D3](:c)(:o)(=O))]	2	\N	\N
[$([c+0v4H0D3](:c)(:s)(C))]	1	\N	\N
[$([C+0v4H0D3](=N)(Br)(N))]	1	\N	\N
[$([c+0v4H0D3](-n)(:c)(:c))]	7	\N	\N
[$([C+0v4H0D3](=N)(c)(c))]	7	\N	\N
[$([C+0v4H0D3](=N)(C)(c))]	2	\N	\N
[$([c+0v4H0D3](:n)(:n)(Br))]	1	\N	\N
[$([c+0v4H0D3](:n)(:n)(Cl))]	2	\N	\N
[$([c+0v4H0D3](:n)(:n)(F))]	2	\N	\N
[$([C+0v4H0D3](=N)(N)(N))]	1	\N	\N
[$([c+0v4H0D3](:n)(:n)(O))]	2	\N	\N
[$([c+0v4H0D3](:n)(:n)(=S))]	2	\N	\N
[$([c+0v4H0D3](:n)(:n)(S))]	8	\N	\N
[$([C+0v4H0D3](=N)(O)(c))]	1	\N	\N
[$([c+0v4H0D3](:n)(:o)(N))]	1	\N	\N
[$([c+0v4H0D3](:n)(:s)(N))]	3	\N	\N
[$([c+0v4H0D3](:n)(:s)(=O))]	1	\N	\N
[$([c+0v4H0D3](:n)(:s)(O))]	1	\N	\N
[$([c+0v4H0D3](:n)(:s)(S))]	1	\N	\N
[$([C+0v4H0D3](=O)(c)(c))]	6	\N	\N
[$([C+0v4H0D3](=O)(C)(S))]	1	\N	\N
[$([C+0v4H0D3](=S)(C)(N))]	4	\N	\N
[$([C+0v4H0D3](=S)(C)(O))]	1	\N	\N
[$([C+0v4H0D3](=S)(N)(c))]	1	\N	\N
[$([C+0v4H0D3](=S)(N)(N))]	6	\N	\N
[$([C+0v4H0D4](Br)(Br)(Br)(C))]	1	\N	\N
[$([C+0v4H0D4](Br)(C)(C)(C))]	1	\N	\N
[$([C+0v4H0D4](C)(C)(c)(c))]	2	\N	\N
[$([C+0v4H0D4](C)(C)(C)(N))]	6	\N	\N
[$([C+0v4H0D4](C)(C)(C)(S))]	5	\N	\N
[$([C+0v4H0D4](C)(C)(F)(F))]	2	\N	\N
[$([C+0v4H0D4](C)(Cl)(Cl)(Cl))]	7	\N	\N
[$([C+0v4H0D4](C)(Cl)(F)(F))]	1	\N	\N
[$([C+0v4H0D4](C)(C)(N)(c))]	1	\N	\N
[$([C+0v4H0D4](C)(C)(N)(N))]	1	\N	\N
[$([C+0v4H0D4](C)(C)(O)(c))]	2	\N	\N
[$([C+0v4H0D4](C)(C)(O)(O))]	1	\N	\N
[$([C+0v4H0D4](C)(F)(F)(O))]	2	\N	\N
[$([C+0v4H0D4](Cl)(Cl)(Cl)(c))]	2	\N	\N
[$([C+0v4H0D4](Cl)(Cl)(Cl)(Cl))]	1	\N	\N
[$([C+0v4H0D4](Cl)(Cl)(Cl)(S))]	1	\N	\N
[$([C+0v4H0D4](C)(N)(c)(c))]	1	\N	\N
[$([C+0v4H0D4](F)(F)(F)(O))]	2	\N	\N
[$([C+0v4H0D4](F)(F)(F)(S))]	8	\N	\N
[$([C+0v4H1D1](#C))]	4	\N	\N
[$([C+0v4H1D2](=C)(Br))]	1	\N	\N
[$([C+0v4H1D2](=C)(Cl))]	1	\N	\N
[$([c+0v4H1D2](:c)(:o))]	4	\N	\N
[$([C+0v4H1D2](=C)(O))]	2	\N	\N
[$([c+0v4H1D2](:c)(:s))]	7	\N	\N
[$([C+0v4H1D2](=N)(c))]	3	\N	\N
[$([C+0v4H1D2](=N)(N))]	3	\N	\N
[$([c+0v4H1D2](:n)(:s))]	2	\N	\N
[$([C+0v4H1D2](=O)(C))]	10	\N	\N
[$([C+0v4H1D2](=O)(n))]	1	\N	\N
[$([C+0v4H1D2](=O)(O))]	3	\N	\N
[$([C+0v4H1D3](Br)(C)(C))]	3	\N	\N
[$([C+0v4H1D3](Br)(C)(Cl))]	1	\N	\N
[$([C+0v4H1D3](Br)(N)(N))]	1	\N	\N
[$([C+0v4H1D3](c)(c)(c))]	1	\N	\N
[$([C+0v4H1D3](C)(c)(c))]	4	\N	\N
[$([C+0v4H1D3](C)(Cl)(Cl))]	3	\N	\N
[$([C+0v4H1D3](C)(C)(n))]	3	\N	\N
[$([C+0v4H1D3](C)(C)(S))]	1	\N	\N
[$([C+0v4H1D3](C)(F)(F))]	2	\N	\N
[$([C+0v4H1D3](Cl)(Cl)(Cl))]	1	\N	\N
[$([C+0v4H1D3](C)(N)(c))]	4	\N	\N
[$([C+0v4H1D3](C)(N)(N))]	2	\N	\N
[$([C+0v4H1D3](C)(N)(O))]	3	\N	\N
[$([C+0v4H1D3](C)(N)(S))]	5	\N	\N
[$([C+0v4H1D3](C)(O)(P))]	1	\N	\N
[$([C+0v4H1D3](F)(F)(c))]	1	\N	\N
[$([C+0v4H1D3](O)(c)(c))]	4	\N	\N
[$([C+0v4H2D2](Br)(c))]	1	\N	\N
[$([C+0v4H2D2](c)(c))]	9	\N	\N
[$([C+0v4H2D2](C)(F))]	4	\N	\N
[$([C+0v4H2D2](C)(I))]	3	\N	\N
[$([C+0v4H2D2](Cl)(c))]	1	\N	\N
[$([C+0v4H2D2](Cl)(Cl))]	1	\N	\N
[$([C+0v4H2D2](C)(n))]	9	\N	\N
[$([C+0v4H2D2](C)(P))]	1	\N	\N
[$([C+0v4H2D2](F)(F))]	1	\N	\N
[$([C+0v4H2D2](N)(N))]	2	\N	\N
[$([C+0v4H2D2](O)(O))]	5	\N	\N
[$([C+0v4H2D2](S)(c))]	1	\N	\N
[$([C+0v4H2D2](S)(n))]	1	\N	\N
[$([C+0v4H2D2](S)(S))]	1	\N	\N
[$([C+0v4H3D1](Br))]	1	\N	\N
[$([C+0v4H3D1](Cl))]	1	\N	\N
[$([C+0v4H3D1](F))]	1	\N	\N
[$([C+0v4H3D1](I))]	1	\N	\N
[$([C+0v4H3D1](P))]	3	\N	\N
[$([Cl+0v1H0D1](N))]	1	\N	\N
[$([F+0v1H0D1](S))]	2	\N	\N
[$([I+0v1H0D1](C))]	4	\N	\N
[$([n+0v3H0D2](:c)(:o))]	3	\N	\N
[$([N+0v3H0D2](=C)(O))]	4	\N	\N
[$([n+0v3H0D2](:c)(:s))]	1	\N	\N
[$([N+0v3H0D2](=N)(c))]	2	\N	\N
[$([N+0v3H0D2](=N)(C))]	2	\N	\N
[$([N+0v3H0D2](=N)(N))]	1	\N	\N
[$([N+0v3H0D2](=O)(c))]	4	\N	\N
[$([n+0v3H0D3](-c)(:c)(:c))]	1	\N	\N
[$([N+0v3H0D3](c)(c)(c))]	1	\N	\N
[$([n+0v3H0D3](-c)(:c)(:n))]	5	\N	\N
[$([n+0v3H0D3](:c)(:c)(N))]	10	\N	\N
[$([n+0v3H0D3](:c)(:c)(O))]	1	\N	\N
[$([N+0v3H0D3](C)(C)(O))]	1	\N	\N
[$([N+0v3H0D3](C)(C)(P))]	3	\N	\N
[$([N+0v3H0D3](C)(C)(S))]	4	\N	\N
[$([n+0v3H0D3](:c)(:n)(C))]	9	\N	\N
[$([N+0v3H0D3](C)(N)(c))]	2	\N	\N
[$([n+0v3H0D3](-c)(:n)(:n))]	1	\N	\N
[$([N+0v3H0D3](N)(c)(c))]	1	\N	\N
[$([N+0v3H1D1](=C))]	1	\N	\N
[$([N+0v3H1D2](c)(c))]	9	\N	\N
[$([N+0v3H1D2](C)(N))]	8	\N	\N
[$([N+0v3H1D2](C)(O))]	4	\N	\N
[$([N+0v3H1D2](C)(P))]	4	\N	\N
[$([N+0v3H1D2](O)(c))]	1	\N	\N
[$([N+0v3H2D1](n))]	10	\N	\N
[$([N+0v3H2D1](N))]	6	\N	\N
[$([N+0v3H2D1](P))]	1	\N	\N
[$([N+0v5H0D3](=O)(=O)(O))]	1	\N	\N
[$([o+0v2H0D2](:c)(:c))]	9	\N	\N
[$([O+0v2H0D2](c)(c))]	7	\N	\N
[$([o+0v2H0D2](:c)(:n))]	3	\N	\N
[$([O+0v2H0D2](C)(N))]	1	\N	\N
[$([O+0v2H0D2](C)(S))]	1	\N	\N
[$([O+0v2H1D1](n))]	1	\N	\N
[$([O+0v2H1D1](N))]	10	\N	\N
[$([P+0v5H0D4](=O)(C)(N)(N))]	1	\N	\N
[$([P+0v5H0D4](=O)(C)(O)(O))]	2	\N	\N
[$([P+0v5H0D4](=O)(C)(O)(S))]	2	\N	\N
[$([P+0v5H0D4](=O)(N)(N)(N))]	2	\N	\N
[$([P+0v5H0D4](=O)(N)(O)(O))]	2	\N	\N
[$([P+0v5H0D4](=O)(N)(O)(S))]	3	\N	\N
[$([P+0v5H0D4](=O)(O)(O)(O))]	6	\N	\N
[$([P+0v5H0D4](=O)(O)(O)(S))]	1	\N	\N
[$([P+0v5H0D4](=O)(O)(S)(S))]	1	\N	\N
[$([P+0v5H0D4](=S)(O)(O)(O))]	5	\N	\N
[$([P+0v5H0D4](=S)(O)(O)(S))]	4	\N	\N
[$([S+0v2H0D1](=c))]	2	\N	\N
[$([S+0v2H0D1](=P))]	9	\N	\N
[$([S+0v2H0D2](C)(N))]	1	\N	\N
[$([s+0v2H0D2](:n)(:n))]	1	\N	\N
[$([S+0v2H1D1](c))]	2	\N	\N
[$([S+0v2H1D1](C))]	3	\N	\N
[$([S+0v4H0D3](=O)(c)(c))]	2	\N	\N
[$([S+0v4H0D3](=O)(C)(c))]	2	\N	\N
[$([S+0v4H0D3](=O)(C)(C))]	1	\N	\N
[$([S+0v6H0D4](=O)(=O)(c)(c))]	2	\N	\N
[$([S+0v6H0D4](=O)(=O)(C)(c))]	9	\N	\N
[$([S+0v6H0D4](=O)(=O)(C)(N))]	3	\N	\N
[$([S+0v6H0D4](=O)(=O)(C)(O))]	1	\N	\N
[$([S+0v6H0D4](=O)(=O)(F)(c))]	2	\N	\N
[$([S+0v6H0D4](=O)(=O)(N)(N))]	1	\N	\N
[$([Br+0v1H0D1](c))]	65	0.733	0.140
[$([Br+0v1H0D1](C))]	22	0.999	0.108
[$([C+0v4H0D2](#N)(c))]	31	0.694	0.309
[$([C+0v4H0D2](#N)(C))]	46	0.418	0.323
[$([C+0v4H0D2](=N)(=S))]	23	1.735	0.255
[$([c+0v4H0D3](:c)(:c)(Br))]	60	0.412	0.147
[$([c+0v4H0D3](-c)(:c)(:c))]	39	0.272	0.042
[$([c+0v4H0D3](:c)(:c)(:c))]	139	0.287	0.023
[$([c+0v4H0D3](:c)(:c)(C))]	734	0.126	0.045
[$([C+0v4H0D3](=C)(C)(C))]	24	-0.032	0.105
[$([c+0v4H0D3](:c)(:c)(Cl))]	176	0.357	0.376
[$([c+0v4H0D3](:c)(:c)(F))]	43	0.551	0.229
[$([c+0v4H0D3](:c)(:c)(I))]	34	1.031	0.369
[$([c+0v4H0D3](:c)(:c)(:n))]	168	0.129	0.069
[$([c+0v4H0D3](:c)(:c)(N))]	456	-0.012	0.064
[$([c+0v4H0D3](:c)(:c)(O))]	509	-0.107	0.081
[$([c+0v4H0D3](:c)(:c)(S))]	110	-0.327	0.125
[$([c+0v4H0D3](:c)(:n)(C))]	57	-0.378	0.091
[$([c+0v4H0D3](:c)(:n)(Cl))]	11	-0.082	0.374
[$([c+0v4H0D3](:c)(:n)(:n))]	38	0.305	0.211
[$([c+0v4H0D3](:c)(:n)(N))]	53	0.384	0.117
[$([c+0v4H0D3](:c)(:n)(=O))]	58	0.053	0.210
[$([C+0v4H0D3](=N)(C)(C))]	39	0.665	0.371
[$([C+0v4H0D3](=N)(C)(N))]	11	-0.360	0.242
[$([c+0v4H0D3](:n)(:n)(C))]	29	-0.445	0.169
[$([c+0v4H0D3](:n)(:n)(N))]	28	0.491	0.161
[$([c+0v4H0D3](:n)(:n)(=O))]	39	-0.001	0.249
[$([C+0v4H0D3](=O)(C)(c))]	52	-0.071	0.118
[$([C+0v4H0D3](=O)(C)(C))]	46	-0.708	0.130
[$([C+0v4H0D3](=O)(C)(N))]	145	-0.290	0.124
[$([C+0v4H0D3](=O)(C)(O))]	238	-0.092	0.128
[$([C+0v4H0D3](=O)(N)(c))]	60	0.064	0.132
[$([C+0v4H0D3](=O)(N)(N))]	58	0.460	0.159
[$([C+0v4H0D3](=O)(N)(O))]	111	0.155	0.150
[$([C+0v4H0D3](=O)(O)(c))]	111	0.844	0.131
[$([C+0v4H0D4](C)(C)(C)(c))]	18	0.040	0.126
[$([C+0v4H0D4](C)(C)(C)(C))]	38	0.102	0.093
[$([C+0v4H0D4](C)(C)(C)(O))]	14	-0.531	0.144
[$([C+0v4H0D4](C)(F)(F)(F))]	16	-0.206	0.169
[$([C+0v4H0D4](F)(F)(F)(c))]	51	-0.005	0.153
[$([c+0v4H1D2](:c)(:c))]	1385	0.328	0.010
[$([C+0v4H1D2](=C)(c))]	41	0.552	0.090
[$([C+0v4H1D2](=C)(C))]	85	0.347	0.039
[$([c+0v4H1D2](:c)(:n))]	254	-0.087	0.065
[$([C+0v4H1D2](=C)(N))]	18	-0.020	0.180
[$([c+0v4H1D2](:n)(:n))]	47	-0.100	0.149
[$([C+0v4H1D2](=O)(c))]	11	0.176	0.179
[$([C+0v4H1D2](=O)(N))]	11	0.049	0.189
[$([C+0v4H1D3](C)(C)(c))]	41	0.327	0.089
[$([C+0v4H1D3](C)(C)(C))]	64	-0.003	0.046
[$([C+0v4H1D3](C)(C)(Cl))]	12	0.033	0.047
[$([C+0v4H1D3](C)(C)(N))]	61	-0.524	0.075
[$([C+0v4H1D3](C)(C)(O))]	117	-0.083	0.083
[$([C+0v4H1D3](C)(O)(c))]	12	0.049	0.169
[$([C+0v4H1D3](C)(O)(n))]	40	0.076	0.204
[$([C+0v4H1D3](C)(O)(O))]	16	-0.216	0.245
[$([C+0v4H2D1](=C))]	41	0.437	0.068
[$([C+0v4H2D2](Br)(C))]	12	-0.053	0.175
[$([C+0v4H2D2](C)(c))]	175	0.368	0.054
[$([C+0v4H2D2](C)(C))]	378	0.446	0.012
[$([C+0v4H2D2](C)(Cl))]	18	0.160	0.090
[$([C+0v4H2D2](C)(N))]	185	-0.155	0.046
[$([C+0v4H2D2](C)(O))]	308	0.004	0.075
[$([C+0v4H2D2](C)(S))]	24	-0.405	0.105
[$([C+0v4H2D2](N)(c))]	13	0.024	0.145
[$([C+0v4H2D2](O)(c))]	28	-0.146	0.123
[$([C+0v4H3D1](c))]	181	0.587	0.054
[$([C+0v4H3D1](C))]	556	0.563	0.023
[$([C+0v4H3D1](n))]	16	0.082	0.175
[$([C+0v4H3D1](N))]	175	-0.003	0.055
[$([C+0v4H3D1](O))]	156	0.161	0.082
[$([C+0v4H3D1](S))]	31	-0.431	0.117
[$([Cl+0v1H0D1](c))]	184	0.589	0.372
[$([Cl+0v1H0D1](C))]	55	0.603	0.037
[$([F+0v1H0D1](c))]	49	-0.085	0.225
[$([F+0v1H0D1](C))]	89	0.382	0.044
[$([I+0v1H0D1](c))]	36	0.421	0.360
[$([N+0v3H0D1](#C))]	79	-0.829	0.299
[$([n+0v3H0D2](:c)(:c))]	308	-0.228	0.120
[$([N+0v3H0D2](=C)(c))]	32	-0.166	0.212
[$([N+0v3H0D2](=C)(C))]	16	-0.249	0.158
[$([n+0v3H0D2](:c)(:n))]	37	-0.234	0.105
[$([N+0v3H0D2](=C)(N))]	37	1.594	0.469
[$([n+0v3H0D2](:n)(:n))]	13	-0.478	0.131
[$([N+0v3H0D2](=O)(N))]	33	-0.297	0.378
[$([n+0v3H0D3](:c)(:c)(C))]	61	0.108	0.199
[$([N+0v3H0D3](C)(c)(c))]	13	0.663	0.225
[$([N+0v3H0D3](C)(C)(c))]	35	0.265	0.130
[$([N+0v3H0D3](C)(C)(C))]	81	-0.116	0.135
[$([N+0v3H0D3](C)(C)(N))]	31	-0.050	0.357
[$([n+0v3H1D2](:c)(:c))]	95	0.027	0.129
[$([N+0v3H1D2](C)(c))]	108	-0.307	0.086
[$([N+0v3H1D2](C)(C))]	181	-0.310	0.094
[$([n+0v3H1D2](:c)(:n))]	11	0.702	0.171
[$([N+0v3H1D2](C)(S))]	19	-0.301	0.255
[$([N+0v3H1D2](N)(c))]	38	-0.231	0.243
[$([N+0v3H1D2](S)(c))]	19	-0.590	0.197
[$([N+0v3H2D1](c))]	176	-0.594	0.075
[$([N+0v3H2D1](C))]	154	-0.812	0.062
[$([N+0v3H2D1](S))]	33	-0.772	0.205
[$([N+0v5H0D3](=O)(=O)(c))]	114	-0.085	0.342
[$([N+0v5H0D3](=O)(=O)(C))]	27	-0.589	0.372
[$([O+0v2H0D1](=c))]	78	-0.843	0.177
[$([O+0v2H0D1](=C))]	787	-0.324	0.110
[$([O+0v2H0D1](=N))]	175	0.159	0.173
[$([O+0v2H0D1](=P))]	20	-1.204	0.168
[$([O+0v2H0D1](=S))]	86	0.027	0.094
[$([O+0v2H0D2](C)(c))]	341	0.219	0.108
[$([O+0v2H0D2](C)(C))]	199	-0.132	0.141
[$([O+0v2H0D2](C)(P))]	26	0.107	0.117
[$([O+0v2H0D2](P)(c))]	11	0.547	0.228
[$([O+0v2H1D1](c))]	200	-0.018	0.087
[$([O+0v2H1D1](C))]	411	-0.457	0.076
[$([S+0v2H0D1](=C))]	35	-0.115	0.157
[$([s+0v2H0D2](:c)(:c))]	12	0.935	0.150
[$([S+0v2H0D2](c)(c))]	12	1.528	0.324
[$([S+0v2H0D2](C)(c))]	28	1.464	0.176
[$([S+0v2H0D2](C)(C))]	17	1.504	0.181
[$([S+0v2H0D2](C)(P))]	11	1.199	0.199
[$([S+0v6H0D4](=O)(=O)(N)(c))]	63	-0.210	0.193
\.


Create OR Replace FUNCTION public.glogp(Bytea)
 Returns numeric AS '
 
Select sum(coefficient*rdkit.count_matches($1,smarts))
   from public.glogp;
' LANGUAGE sql IMMUTABLE;
Comment on FUNCTION public.glogp(Bytea)
 Is 'glogP using gNova fragments';
