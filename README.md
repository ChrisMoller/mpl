# mpl
Yet another array/math oriented language w/o APL's weird symbols and J's general weirdness and wonky installation.


**NOT EVEN REMOTELY READY FOR PRIME TIME -- PLEASE DON'T GRAB IT YET!**


But just as a sample of what works so far, here's a test log:  

#tests/test2.mpl

*ccc			<!-- ok     show 6668	-->ccc			<!-- ok     show 6668	-->

*	\>\> 6668

*aaa + 7			<!-- ok     show 15	-->aaa + 7			<!-- ok     show 15	-->

*	\>\> 15

*bbb			<!-- ok     uninitialised: should show bbb --> >> bbb

*aaa+ccc			<!-- ok     show 6676	-->aaa+ccc			<!-- ok     show 6676	-->

*	\>\> 6676

*aaa+ccc+bbb		<!-- ok     show -59990 -->aaa+ccc+bbb		<!-- ok     show -59990 -->

*	\>\> -59990

#tests/test3.mpl

*ccc			<!-- ok     show 6668	-->ccc			<!-- ok     show 6668	-->

*	\>\> 6668

#tests/testarray.mpl

*1 2 3 4 5 * 6>> 6 12 18 24 30 
*7 * 8 9 10 11>> 56 63 70 77 
*1 2 3 4 5 * 6 7 8 9 10>> 6 14 24 36 50 
*(::20)[2 3#4 7 3 9 1 11]
*>> 4 7 3 
*>> 9 1 11 
*(6*::6)[4 1]>> 24 6 
*(::6)[ 3](::6)[ 3]
*	>> 3
*a = ?5#10; b => a;  a[b]>> 2.34495 2.36975 2.4625 5.28537 7.27973 
*>a>> 2 1 0 4 3 
*tests/test.mpl
*4-74-7
*	>> -3
*4+84+8
*	>> 12
*5^35^3
*	>> 125
*10 log 310 log 3
*	>> 0.477121
*2 root 22 root 2
*	>> 1.41421
*ln 2.7ln 2.7
*	>> 0.993252
*log 2log 2
*	>> 0.30103
*exp 2exp 2
*	>> 7.38906
*root 25root 25
*	>> 5
*(root 25)^2(root 25)^2
*	>> 25
*((sin 0.5)^2) + (cos 0.5)^2 ((sin 0.5)^2) + (cos 0.5)^2 
*	>> 1
*atan (sin .5) / cos .5atan (sin .5) / cos .5
*	>> 0.5
*.5 atan .5.5 atan .5
*	>> 0.785398
*tan .5 atan .5tan .5 atan .5
*	>> 1
*|-8|-8
*	>> 8
*!5!5
*	>> 120
*ceil 5.5ceil 5.5
*	>> 6
*floor 5.5floor 5.5
*	>> 5
*round 5.4round 5.4
*	>> 5
*round 5.6round 5.6
*	>> 6
*sind 60sind 60
*	>> 0.866025
*cosd 60cosd 60
*	>> 0.5
*tand 45tand 45
*	>> 1
*asind sind 60asind sind 60
*	>> 60
*acosd cosd 60acosd cosd 60
*	>> 60
*atand tand 45atand tand 45
*	>> 45
*1 atand 11 atand 1
*	>> 45
*