def  fourthnight '
{
shopen
plotselect keith
michkon
harm1
gscan_on u35d

umv energy E0
ascan energy E0-0.003 E0+0.003 60 1

for (etrans=-5;etrans<=40;etrans=etrans+0.5) {
umv energy E0+etrans/1000
umv schi2 -A[mono]
dscan schi2 -0.015 0.015 60 .1
umv schi2 CEN
dscan schi2 -0.005 0.005 20 7
dscan schi2 -0.5 0.5 1 7
}

for (etrans=42;etrans<=120;etrans=etrans+2) {
umv energy E0+etrans/1000
umv schi2 -A[mono]
dscan schi2 -0.015 0.015 60 .1
umv schi2 CEN
dscan schi2 -0.005 0.005 20 7
dscan schi2 -0.5 0.5 1 7
}


umv energy E0
ascan energy E0-0.003 E0+0.003 60 1

for (etrans=-5;etrans<=40;etrans=etrans+0.5) {
umv energy E0+etrans/1000
umv schi2 -A[mono]
dscan schi2 -0.015 0.015 60 .1
umv schi2 CEN
dscan schi2 -0.005 0.005 20 7
dscan schi2 -0.5 0.5 1 7
}

}
'







