ccc =  6668;  		<!-- ok     noshow	-->
ccc;			<!-- ok     show 6668	-->
aaa =  8;		<!-- ok     noshow	-->	
aaa + 7	;		<!-- ok     show 15	-->
bbb;			<!-- ok     uninitialised: should show bbb -->
aaa+ccc	;		<!-- ok     show 6676	-->
bbb = -66666;		<!-- ok     noshow	-->
aaa+ccc+bbb;		<!-- ok     show -59990 -->
