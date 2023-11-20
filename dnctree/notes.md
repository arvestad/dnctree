Computed 5820774 distances for 56771 taxa. A distance matrix contains 1611444835 pairs. Savings: 99.63878540092873 %
 dnctree  --base_case_size 100 --verbose  -f fasta sample_data/PF00001_full.txt


dnctree --verbose -f fasta sample_data/PF00001_rp15.txt
Computed 1412289 distances for 16338 taxa. A distance matrix contains 133456953 pairs. Savings: 98.94176439049976 %
Takes 45 mins

dnctree  --verbose  -w sample_data/PF00001_rp15.txt > /dev/null
Computed 1396279 distances for 16338 taxa. A full distance matrix would contain 133456953 pairs. Savings: 99.0 %
real	41m4.461s
user	40m37.568s
sys	0m6.762s
