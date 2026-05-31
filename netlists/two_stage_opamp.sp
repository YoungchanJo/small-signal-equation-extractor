* Two stage OPAMP
.include 45nm_bulk.txt
.include @PARAM_PATH@

mp1 n4 n4 VDD VDD pmos l={L0} w={W0} m={M0}
mp2 n5 n4 VDD VDD pmos l={L0} w={W0} m={M0}

mn1 n4 n2 n3 n3 nmos l={L1} w={W1} m={M1}
mn2 n5 n1 n3 n3 nmos l={L1} w={W1} m={M1}

mn3 n3 n6 VSS VSS nmos l={L2} w={W2} m={M2}

mn4 n6 n6 VSS VSS nmos l={L3} w={W3} m={M3}

mp3 vout n5 VDD VDD pmos l={L4} w={W4} m={M4}

mn5 vout n6 VSS VSS nmos l={L5} w={W5} m={M5}

cc n5 vout {CL}
ibias VDD n6 {IREF}

vin in 0 dc=0 ac=1.0
ein1 n1 cm in 0 0.5
ein2 n2 cm in 0 -0.5
vcm cm 0 dc=0.6
vdd VDD 0 dc=1.2
vss 0 VSS dc=0

.control
option numdgt=4
set wr_singlescale
op
wrdata @DC_PATH@ i(vdd) v(n3) v(n4) v(n5) v(n6)
reset
set units=degrees
ac dec 10 1k 10G
save v(vout)
run
wrdata @AC_PATH@ vdb(vout) vp(vout)

.endc

.end