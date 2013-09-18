##SEM for the lichen co-occurrence network data
#geno = genotype, composite variable, 1
#T6, RL6, 1017, WC5, H10, 1005, 996, 1023, 999,1008,11,1012,1007,10 = indicators of geno, 11 - 114
####net.lat = latitude of the trees, 2
#net.rough = bark roughness, 3
#pc1 = principal component axis of network structure, 4

#composite variable construction (geno)
T6 -> geno, l.11.1, NA
RL6 -> geno, l.12.1, NA
1017 -> geno, l.13.1, NA
WC5 -> geno, l.14.1, NA
H10 -> geno, l.15.1, NA
1005 -> geno, l.16.1, NA
996 -> geno, l.17.1, NA
1023 -> geno, l.18.1, NA
999 -> geno, l.19.1, NA
1008 -> geno, l.110.1, NA
11 -> geno, l.111.1, NA
1012 -> geno, l.112.1, NA
1007 -> geno, l.113.1, NA
10 -> geno, l.114.1, 1

#composite paths
geno -> net.rough, g.1.3, NA
geno -> pc1, g.1.4, NA

#observed exogenous paths
#net.lat -> pc1, g.2.4, NA

#observed endogenous paths
net.rough -> pc1, g.3.4, NA

#errors
geno <-> geno, p.1.1, NA
T6 <-> T6, p.11.11, 1
RL6 <-> RL6, p.12.12, 1
1017 <-> 1017, p.13.13, 1
WC5 <-> WC5, p.14.14, 1
H10 <-> H10, p.15.15, 1
1005 <-> 1005, p.16.16, 1
996 <-> 996, p.17.17, 1
1023 <-> 1023, p.18.18, 1
999 <-> 999, p.19.19, 1
1008 <-> 1008, p.110.110, 1
11 <-> 11, p.111.111, 1
1012 <-> 1012, p.112.112, 1
1007 <-> 1007, p.113.113, 1
10 <-> 10, p.114.114, 1

#net.lat <-> net.lat, e.2.2, 1
net.rough <-> net.rough, e.3.3, NA
pc1 <-> pc1, e.4.4, NA
