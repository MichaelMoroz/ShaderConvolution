//The first 6 columns of the U matrix: 
float U[546] = float[](-1.2516E-07, -1.2715E-06, -3.7479E-06, -7.6960E-06, -1.3266E-05, -2.0624E-05, -2.9956E-05, -4.1471E-05, -5.5407E-05, -7.2033E-05, -9.1658E-05, -1.1464E-04, -1.4137E-04, -1.7234E-04, -2.0809E-04, -2.4925E-04, -2.9658E-04, -3.5094E-04, -4.1340E-04, -4.8519E-04, -5.6781E-04, -6.6306E-04, -7.7315E-04, -9.0076E-04, -1.0492E-03, -1.2228E-03, -1.4266E-03, -1.6675E-03, -1.9541E-03, -2.2979E-03, -2.7140E-03, -3.2233E-03, -3.8541E-03, -4.6468E-03, -5.6603E-03, -6.9822E-03, -8.7489E-03, -1.1181E-02, -1.4654E-02, -1.9853E-02, -2.8130E-02, -4.2461E-02, -7.0466E-02, -1.3634E-01, -3.4306E-01, -8.4285E-01, -3.4306E-01, -1.3634E-01, -7.0466E-02, -4.2461E-02, -2.8130E-02, -1.9853E-02, -1.4654E-02, -1.1181E-02, -8.7489E-03, -6.9822E-03, -5.6603E-03, -4.6468E-03, -3.8541E-03, -3.2233E-03, -2.7140E-03, -2.2979E-03, -1.9541E-03, -1.6675E-03, -1.4266E-03, -1.2228E-03, -1.0492E-03, -9.0076E-04, -7.7315E-04, -6.6306E-04, -5.6781E-04, -4.8519E-04, -4.1340E-04, -3.5094E-04, -2.9658E-04, -2.4925E-04, -2.0809E-04, -1.7234E-04, -1.4137E-04, -1.1464E-04, -9.1658E-05, -7.2033E-05, -5.5407E-05, -4.1471E-05, -2.9956E-05, -2.0624E-05, -1.3266E-05, -7.6960E-06, -3.7479E-06, -1.2715E-06, -1.2516E-07, 4.5871E-07, 6.4334E-06, 2.0853E-05, 4.4966E-05, 7.9899E-05, 1.2682E-04, 1.8699E-04, 2.6182E-04, 3.5285E-04, 4.6186E-04, 5.9085E-04, 7.4208E-04, 9.1816E-04, 1.1221E-03, 1.3573E-03, 1.6277E-03, 1.9381E-03, 2.2939E-03, 2.7013E-03, 3.1681E-03, 3.7031E-03, 4.3171E-03, 5.0230E-03, 5.8365E-03, 6.7766E-03, 7.8670E-03, 9.1369E-03, 1.0623E-02, 1.2372E-02, 1.4444E-02, 1.6916E-02, 1.9891E-02, 2.3507E-02, 2.7952E-02, 3.3485E-02, 4.0481E-02, 4.9482E-02, 6.1305E-02, 7.7220E-02, 9.9271E-02, 1.3086E-01, 1.7779E-01, 2.4931E-01, 3.5075E-01, 3.5119E-01, -4.7991E-01, 3.5119E-01, 3.5075E-01, 2.4931E-01, 1.7779E-01, 1.3086E-01, 9.9271E-02, 7.7220E-02, 6.1305E-02, 4.9482E-02, 4.0481E-02, 3.3485E-02, 2.7952E-02, 2.3507E-02, 1.9891E-02, 1.6916E-02, 1.4444E-02, 1.2372E-02, 1.0623E-02, 9.1369E-03, 7.8670E-03, 6.7766E-03, 5.8365E-03, 5.0230E-03, 4.3171E-03, 3.7031E-03, 3.1681E-03, 2.7013E-03, 2.2939E-03, 1.9381E-03, 1.6277E-03, 1.3573E-03, 1.1221E-03, 9.1816E-04, 7.4208E-04, 5.9085E-04, 4.6186E-04, 3.5285E-04, 2.6182E-04, 1.8699E-04, 1.2682E-04, 7.9899E-05, 4.4966E-05, 2.0853E-05, 6.4334E-06, 4.5871E-07, 5.4940E-08, 1.1276E-05, 4.8633E-05, 1.1934E-04, 2.2875E-04, 3.8176E-04, 5.8329E-04, 8.3852E-04, 1.1531E-03, 1.5331E-03, 1.9854E-03, 2.5178E-03, 3.1388E-03, 3.8584E-03, 4.6875E-03, 5.6389E-03, 6.7272E-03, 7.9691E-03, 9.3838E-03, 1.0994E-02, 1.2825E-02, 1.4908E-02, 1.7277E-02, 1.9977E-02, 2.3054E-02, 2.6570E-02, 3.0595E-02, 3.5212E-02, 4.0526E-02, 4.6658E-02, 5.3760E-02, 6.2013E-02, 7.1640E-02, 8.2908E-02, 9.6135E-02, 1.1168E-01, 1.2994E-01, 1.5121E-01, 1.7551E-01, 2.0187E-01, 2.2661E-01, 2.3823E-01, 2.0252E-01, 2.2027E-02, -4.0653E-01, 2.2198E-01, -4.0653E-01, 2.2027E-02, 2.0252E-01, 2.3823E-01, 2.2661E-01, 2.0187E-01, 1.7551E-01, 1.5121E-01, 1.2994E-01, 1.1168E-01, 9.6135E-02, 8.2908E-02, 7.1640E-02, 6.2013E-02, 5.3760E-02, 4.6658E-02, 4.0526E-02, 3.5212E-02, 3.0595E-02, 2.6570E-02, 2.3054E-02, 1.9977E-02, 1.7277E-02, 1.4908E-02, 1.2825E-02, 1.0994E-02, 9.3838E-03, 7.9691E-03, 6.7272E-03, 5.6389E-03, 4.6875E-03, 3.8584E-03, 3.1388E-03, 2.5178E-03, 1.9854E-03, 1.5331E-03, 1.1531E-03, 8.3852E-04, 5.8329E-04, 3.8176E-04, 2.2875E-04, 1.1934E-04, 4.8633E-05, 1.1276E-05, 5.4940E-08, 2.2927E-06, 1.2485E-05, -1.2711E-05, -1.0779E-04, -2.9609E-04, -5.9598E-04, -1.0236E-03, -1.5941E-03, -2.3229E-03, -3.2254E-03, -4.3184E-03, -5.6194E-03, -7.1475E-03, -8.9238E-03, -1.0971E-02, -1.3313E-02, -1.5979E-02, -1.8998E-02, -2.2403E-02, -2.6230E-02, -3.0519E-02, -3.5310E-02, -4.0650E-02, -4.6587E-02, -5.3170E-02, -6.0448E-02, -6.8469E-02, -7.7271E-02, -8.6880E-02, -9.7294E-02, -1.0847E-01, -1.2028E-01, -1.3249E-01, -1.4467E-01, -1.5606E-01, -1.6541E-01, -1.7060E-01, -1.6819E-01, -1.5254E-01, -1.1467E-01, -4.1217E-02, 8.3146E-02, 2.5268E-01, 3.2006E-01, -2.6119E-01, 9.1322E-02, -2.6119E-01, 3.2006E-01, 2.5268E-01, 8.3146E-02, -4.1217E-02, -1.1467E-01, -1.5254E-01, -1.6819E-01, -1.7060E-01, -1.6541E-01, -1.5606E-01, -1.4467E-01, -1.3249E-01, -1.2028E-01, -1.0847E-01, -9.7294E-02, -8.6880E-02, -7.7271E-02, -6.8469E-02, -6.0448E-02, -5.3170E-02, -4.6587E-02, -4.0650E-02, -3.5310E-02, -3.0519E-02, -2.6230E-02, -2.2403E-02, -1.8998E-02, -1.5979E-02, -1.3313E-02, -1.0971E-02, -8.9238E-03, -7.1475E-03, -5.6194E-03, -4.3184E-03, -3.2254E-03, -2.3229E-03, -1.5941E-03, -1.0236E-03, -5.9598E-04, -2.9609E-04, -1.0779E-04, -1.2711E-05, 1.2485E-05, 2.2927E-06, -1.0138E-06, 8.1916E-05, 2.4808E-04, 3.8300E-04, 3.9185E-04, 2.0259E-04, -2.4147E-04, -9.8701E-04, -2.0740E-03, -3.5381E-03, -5.4125E-03, -7.7283E-03, -1.0516E-02, -1.3803E-02, -1.7619E-02, -2.1990E-02, -2.6938E-02, -3.2486E-02, -3.8648E-02, -4.5434E-02, -5.2842E-02, -6.0857E-02, -6.9443E-02, -7.8541E-02, -8.8053E-02, -9.7834E-02, -1.0767E-01, -1.1727E-01, -1.2621E-01, -1.3393E-01, -1.3965E-01, -1.4232E-01, -1.4058E-01, -1.3262E-01, -1.1619E-01, -8.8552E-02, -4.6747E-02, 1.1662E-02, 8.6578E-02, 1.7053E-01, 2.3800E-01, 2.2568E-01, 2.0588E-02, -3.6559E-01, 1.4046E-01, -3.7472E-02, 1.4046E-01, -3.6559E-01, 2.0588E-02, 2.2568E-01, 2.3800E-01, 1.7053E-01, 8.6578E-02, 1.1662E-02, -4.6747E-02, -8.8552E-02, -1.1619E-01, -1.3262E-01, -1.4058E-01, -1.4232E-01, -1.3965E-01, -1.3393E-01, -1.2621E-01, -1.1727E-01, -1.0767E-01, -9.7834E-02, -8.8053E-02, -7.8541E-02, -6.9443E-02, -6.0857E-02, -5.2842E-02, -4.5434E-02, -3.8648E-02, -3.2486E-02, -2.6938E-02, -2.1990E-02, -1.7619E-02, -1.3803E-02, -1.0516E-02, -7.7283E-03, -5.4125E-03, -3.5381E-03, -2.0740E-03, -9.8701E-04, -2.4147E-04, 2.0259E-04, 3.9185E-04, 3.8300E-04, 2.4808E-04, 8.1916E-05, -1.0138E-06, -1.4852E-05, 2.6268E-05, 6.7113E-04, 1.8106E-03, 3.1240E-03, 4.3138E-03, 5.1407E-03, 5.4200E-03, 5.0125E-03, 3.8129E-03, 1.7434E-03, -1.2535E-03, -5.2163E-03, -1.0167E-02, -1.6111E-02, -2.3037E-02, -3.0911E-02, -3.9676E-02, -4.9245E-02, -5.9493E-02, -7.0251E-02, -8.1291E-02, -9.2318E-02, -1.0295E-01, -1.1270E-01, -1.2098E-01, -1.2702E-01, -1.2993E-01, -1.2864E-01, -1.2191E-01, -1.0837E-01, -8.6647E-02, -5.5503E-02, -1.4259E-02, 3.6561E-02, 9.4109E-02, 1.5145E-01, 1.9505E-01, 2.0217E-01, 1.4108E-01, -1.5200E-02, -2.3213E-01, -2.7172E-01, 2.8423E-01, -6.8474E-02, 1.5249E-02, -6.8474E-02, 2.8423E-01, -2.7172E-01, -2.3213E-01, -1.5200E-02, 1.4108E-01, 2.0217E-01, 1.9505E-01, 1.5145E-01, 9.4109E-02, 3.6561E-02, -1.4259E-02, -5.5503E-02, -8.6647E-02, -1.0837E-01, -1.2191E-01, -1.2864E-01, -1.2993E-01, -1.2702E-01, -1.2098E-01, -1.1270E-01, -1.0295E-01, -9.2318E-02, -8.1291E-02, -7.0251E-02, -5.9493E-02, -4.9245E-02, -3.9676E-02, -3.0911E-02, -2.3037E-02, -1.6111E-02, -1.0167E-02, -5.2163E-03, -1.2535E-03, 1.7434E-03, 3.8129E-03, 5.0125E-03, 5.4200E-03, 5.1407E-03, 4.3138E-03, 3.1240E-03, 1.8106E-03, 6.7113E-04, 2.6268E-05, -1.4852E-05); 
//The first 6 rows of the V matrix: 
float V[546] = float[](-1.2516E-07, -1.2715E-06, -3.7479E-06, -7.6960E-06, -1.3266E-05, -2.0624E-05, -2.9956E-05, -4.1471E-05, -5.5407E-05, -7.2033E-05, -9.1658E-05, -1.1464E-04, -1.4137E-04, -1.7234E-04, -2.0809E-04, -2.4925E-04, -2.9658E-04, -3.5094E-04, -4.1340E-04, -4.8519E-04, -5.6781E-04, -6.6306E-04, -7.7315E-04, -9.0076E-04, -1.0492E-03, -1.2228E-03, -1.4266E-03, -1.6675E-03, -1.9541E-03, -2.2979E-03, -2.7140E-03, -3.2233E-03, -3.8541E-03, -4.6468E-03, -5.6603E-03, -6.9822E-03, -8.7489E-03, -1.1181E-02, -1.4654E-02, -1.9853E-02, -2.8130E-02, -4.2461E-02, -7.0466E-02, -1.3634E-01, -3.4306E-01, -8.4285E-01, -3.4306E-01, -1.3634E-01, -7.0466E-02, -4.2461E-02, -2.8130E-02, -1.9853E-02, -1.4654E-02, -1.1181E-02, -8.7489E-03, -6.9822E-03, -5.6603E-03, -4.6468E-03, -3.8541E-03, -3.2233E-03, -2.7140E-03, -2.2979E-03, -1.9541E-03, -1.6675E-03, -1.4266E-03, -1.2228E-03, -1.0492E-03, -9.0076E-04, -7.7315E-04, -6.6306E-04, -5.6781E-04, -4.8519E-04, -4.1340E-04, -3.5094E-04, -2.9658E-04, -2.4925E-04, -2.0809E-04, -1.7234E-04, -1.4137E-04, -1.1464E-04, -9.1658E-05, -7.2033E-05, -5.5407E-05, -4.1471E-05, -2.9956E-05, -2.0624E-05, -1.3266E-05, -7.6960E-06, -3.7479E-06, -1.2715E-06, -1.2516E-07, 4.5871E-07, 6.4334E-06, 2.0853E-05, 4.4966E-05, 7.9899E-05, 1.2682E-04, 1.8699E-04, 2.6182E-04, 3.5285E-04, 4.6186E-04, 5.9085E-04, 7.4208E-04, 9.1816E-04, 1.1221E-03, 1.3573E-03, 1.6277E-03, 1.9381E-03, 2.2939E-03, 2.7013E-03, 3.1681E-03, 3.7031E-03, 4.3171E-03, 5.0230E-03, 5.8365E-03, 6.7766E-03, 7.8670E-03, 9.1369E-03, 1.0623E-02, 1.2372E-02, 1.4444E-02, 1.6916E-02, 1.9891E-02, 2.3507E-02, 2.7952E-02, 3.3485E-02, 4.0481E-02, 4.9482E-02, 6.1305E-02, 7.7220E-02, 9.9271E-02, 1.3086E-01, 1.7779E-01, 2.4931E-01, 3.5075E-01, 3.5119E-01, -4.7991E-01, 3.5119E-01, 3.5075E-01, 2.4931E-01, 1.7779E-01, 1.3086E-01, 9.9271E-02, 7.7220E-02, 6.1305E-02, 4.9482E-02, 4.0481E-02, 3.3485E-02, 2.7952E-02, 2.3507E-02, 1.9891E-02, 1.6916E-02, 1.4444E-02, 1.2372E-02, 1.0623E-02, 9.1369E-03, 7.8670E-03, 6.7766E-03, 5.8365E-03, 5.0230E-03, 4.3171E-03, 3.7031E-03, 3.1681E-03, 2.7013E-03, 2.2939E-03, 1.9381E-03, 1.6277E-03, 1.3573E-03, 1.1221E-03, 9.1816E-04, 7.4208E-04, 5.9085E-04, 4.6186E-04, 3.5285E-04, 2.6182E-04, 1.8699E-04, 1.2682E-04, 7.9899E-05, 4.4966E-05, 2.0853E-05, 6.4334E-06, 4.5871E-07, 5.4940E-08, 1.1276E-05, 4.8633E-05, 1.1934E-04, 2.2875E-04, 3.8176E-04, 5.8329E-04, 8.3852E-04, 1.1531E-03, 1.5331E-03, 1.9854E-03, 2.5178E-03, 3.1388E-03, 3.8584E-03, 4.6875E-03, 5.6389E-03, 6.7272E-03, 7.9691E-03, 9.3838E-03, 1.0994E-02, 1.2825E-02, 1.4908E-02, 1.7277E-02, 1.9977E-02, 2.3054E-02, 2.6570E-02, 3.0595E-02, 3.5212E-02, 4.0526E-02, 4.6658E-02, 5.3760E-02, 6.2013E-02, 7.1640E-02, 8.2908E-02, 9.6135E-02, 1.1168E-01, 1.2994E-01, 1.5121E-01, 1.7551E-01, 2.0187E-01, 2.2661E-01, 2.3823E-01, 2.0252E-01, 2.2027E-02, -4.0653E-01, 2.2198E-01, -4.0653E-01, 2.2027E-02, 2.0252E-01, 2.3823E-01, 2.2661E-01, 2.0187E-01, 1.7551E-01, 1.5121E-01, 1.2994E-01, 1.1168E-01, 9.6135E-02, 8.2908E-02, 7.1640E-02, 6.2013E-02, 5.3760E-02, 4.6658E-02, 4.0526E-02, 3.5212E-02, 3.0595E-02, 2.6570E-02, 2.3054E-02, 1.9977E-02, 1.7277E-02, 1.4908E-02, 1.2825E-02, 1.0994E-02, 9.3838E-03, 7.9691E-03, 6.7272E-03, 5.6389E-03, 4.6875E-03, 3.8584E-03, 3.1388E-03, 2.5178E-03, 1.9854E-03, 1.5331E-03, 1.1531E-03, 8.3852E-04, 5.8329E-04, 3.8176E-04, 2.2875E-04, 1.1934E-04, 4.8633E-05, 1.1276E-05, 5.4940E-08, 2.2927E-06, 1.2485E-05, -1.2711E-05, -1.0779E-04, -2.9609E-04, -5.9598E-04, -1.0236E-03, -1.5941E-03, -2.3229E-03, -3.2254E-03, -4.3184E-03, -5.6194E-03, -7.1475E-03, -8.9238E-03, -1.0971E-02, -1.3313E-02, -1.5979E-02, -1.8998E-02, -2.2403E-02, -2.6230E-02, -3.0519E-02, -3.5310E-02, -4.0650E-02, -4.6587E-02, -5.3170E-02, -6.0448E-02, -6.8469E-02, -7.7271E-02, -8.6880E-02, -9.7294E-02, -1.0847E-01, -1.2028E-01, -1.3249E-01, -1.4467E-01, -1.5606E-01, -1.6541E-01, -1.7060E-01, -1.6819E-01, -1.5254E-01, -1.1467E-01, -4.1217E-02, 8.3146E-02, 2.5268E-01, 3.2006E-01, -2.6119E-01, 9.1322E-02, -2.6119E-01, 3.2006E-01, 2.5268E-01, 8.3146E-02, -4.1217E-02, -1.1467E-01, -1.5254E-01, -1.6819E-01, -1.7060E-01, -1.6541E-01, -1.5606E-01, -1.4467E-01, -1.3249E-01, -1.2028E-01, -1.0847E-01, -9.7294E-02, -8.6880E-02, -7.7271E-02, -6.8469E-02, -6.0448E-02, -5.3170E-02, -4.6587E-02, -4.0650E-02, -3.5310E-02, -3.0519E-02, -2.6230E-02, -2.2403E-02, -1.8998E-02, -1.5979E-02, -1.3313E-02, -1.0971E-02, -8.9238E-03, -7.1475E-03, -5.6194E-03, -4.3184E-03, -3.2254E-03, -2.3229E-03, -1.5941E-03, -1.0236E-03, -5.9598E-04, -2.9609E-04, -1.0779E-04, -1.2711E-05, 1.2485E-05, 2.2927E-06, -1.0138E-06, 8.1916E-05, 2.4808E-04, 3.8300E-04, 3.9185E-04, 2.0259E-04, -2.4147E-04, -9.8701E-04, -2.0740E-03, -3.5381E-03, -5.4125E-03, -7.7283E-03, -1.0516E-02, -1.3803E-02, -1.7619E-02, -2.1990E-02, -2.6938E-02, -3.2486E-02, -3.8648E-02, -4.5434E-02, -5.2842E-02, -6.0857E-02, -6.9443E-02, -7.8541E-02, -8.8053E-02, -9.7834E-02, -1.0767E-01, -1.1727E-01, -1.2621E-01, -1.3393E-01, -1.3965E-01, -1.4232E-01, -1.4058E-01, -1.3262E-01, -1.1619E-01, -8.8552E-02, -4.6747E-02, 1.1662E-02, 8.6578E-02, 1.7053E-01, 2.3800E-01, 2.2568E-01, 2.0588E-02, -3.6559E-01, 1.4046E-01, -3.7472E-02, 1.4046E-01, -3.6559E-01, 2.0588E-02, 2.2568E-01, 2.3800E-01, 1.7053E-01, 8.6578E-02, 1.1662E-02, -4.6747E-02, -8.8552E-02, -1.1619E-01, -1.3262E-01, -1.4058E-01, -1.4232E-01, -1.3965E-01, -1.3393E-01, -1.2621E-01, -1.1727E-01, -1.0767E-01, -9.7834E-02, -8.8053E-02, -7.8541E-02, -6.9443E-02, -6.0857E-02, -5.2842E-02, -4.5434E-02, -3.8648E-02, -3.2486E-02, -2.6938E-02, -2.1990E-02, -1.7619E-02, -1.3803E-02, -1.0516E-02, -7.7283E-03, -5.4125E-03, -3.5381E-03, -2.0740E-03, -9.8701E-04, -2.4147E-04, 2.0259E-04, 3.9185E-04, 3.8300E-04, 2.4808E-04, 8.1916E-05, -1.0138E-06, -1.4852E-05, 2.6268E-05, 6.7113E-04, 1.8106E-03, 3.1240E-03, 4.3138E-03, 5.1407E-03, 5.4200E-03, 5.0125E-03, 3.8129E-03, 1.7434E-03, -1.2535E-03, -5.2163E-03, -1.0167E-02, -1.6111E-02, -2.3037E-02, -3.0911E-02, -3.9676E-02, -4.9245E-02, -5.9493E-02, -7.0251E-02, -8.1291E-02, -9.2318E-02, -1.0295E-01, -1.1270E-01, -1.2098E-01, -1.2702E-01, -1.2993E-01, -1.2864E-01, -1.2191E-01, -1.0837E-01, -8.6647E-02, -5.5503E-02, -1.4259E-02, 3.6561E-02, 9.4109E-02, 1.5145E-01, 1.9505E-01, 2.0217E-01, 1.4108E-01, -1.5200E-02, -2.3213E-01, -2.7172E-01, 2.8423E-01, -6.8474E-02, 1.5249E-02, -6.8474E-02, 2.8423E-01, -2.7172E-01, -2.3213E-01, -1.5200E-02, 1.4108E-01, 2.0217E-01, 1.9505E-01, 1.5145E-01, 9.4109E-02, 3.6561E-02, -1.4259E-02, -5.5503E-02, -8.6647E-02, -1.0837E-01, -1.2191E-01, -1.2864E-01, -1.2993E-01, -1.2702E-01, -1.2098E-01, -1.1270E-01, -1.0295E-01, -9.2318E-02, -8.1291E-02, -7.0251E-02, -5.9493E-02, -4.9245E-02, -3.9676E-02, -3.0911E-02, -2.3037E-02, -1.6111E-02, -1.0167E-02, -5.2163E-03, -1.2535E-03, 1.7434E-03, 3.8129E-03, 5.0125E-03, 5.4200E-03, 5.1407E-03, 4.3138E-03, 3.1240E-03, 1.8106E-03, 6.7113E-04, 2.6268E-05, -1.4852E-05); 
//center of the convolution 
int Nc = 45; 
//number of ranks 
int Nr = 6; 
//convolution size 
int N = 91; 

float pack2(vec2 a)
{
    return uintBitsToFloat(packHalf2x16(a));
}

vec2 unpack2(float packed)
{
    return unpackHalf2x16(floatBitsToUint(packed));
}

vec3 pack2vec3(vec3 a, vec3 b) 
{
    uvec3 packed = uvec3(packHalf2x16(vec2(a.x,b.x)), packHalf2x16(vec2(a.y,b.y)), packHalf2x16(vec2(a.z,b.z)));
    return vec3(uintBitsToFloat(packed.x), uintBitsToFloat(packed.y), uintBitsToFloat(packed.z));
}

void unpack2vec3(vec3 packed, out vec3 a, out vec3 b) 
{
    uvec3 unpacked = uvec3(floatBitsToUint(packed.x), floatBitsToUint(packed.y), floatBitsToUint(packed.z));
    vec2 unpackedA = unpackHalf2x16(unpacked.x);
    vec2 unpackedB = unpackHalf2x16(unpacked.y);
    vec2 unpackedC = unpackHalf2x16(unpacked.z);
    a = vec3(unpackedA.x, unpackedB.x, unpackedC.x);
    b = vec3(unpackedA.y, unpackedB.y, unpackedC.y);
}

//first pass of separable convolution
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{    
    //do a convolution in the x direction on iChannel0
    //since we have 6 ranks and we can only store 2 in a vec3 we need to do the partial sums in 1/2 resolution chunks in a single buffer(which means max 8 ranks)
    ivec2 coord = ivec2(fragCoord.xy) % ivec2(iResolution.xy/2.0);
    ivec2 block = ivec2(fragCoord.xy) / ivec2(iResolution.xy/2.0);
    int rank_offset = 2 * (block.x + block.y * 2);

    if(rank_offset >= Nr)
        discard;

    vec3 sum0 = vec3(0.0);
    vec3 sum1 = vec3(0.0);
    for (int i = -Nc; i < Nc; i++) 
    {
        //use texelFetch to get the pixel at the current index
        //since we are doing it at 1/2 resolution we need to use 1 mip level higher
        vec4 pixel = texelFetch(iChannel0, coord + ivec2(i, 0), 1);
        vec3 pcolor = pixel.xyz/(pixel.w+1e-6);
        
        //get the rank index
        int offset = i + Nc;
        sum0 += pcolor * U[offset + rank_offset * N];
        sum1 += pcolor * U[offset + (rank_offset + 1) * N];
    }

    //write the sum packed to the output
    fragColor = vec4(pack2vec3(sum0, sum1), 1.0);
}

//second pass: do the convolution in the y direction
//and write the result to the output    
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    //do a convolution in the y direction on iChannel0
    ivec2 resol = ivec2(iResolution.xy/2.0);
    ivec2 coord = ivec2(fragCoord.xy) % resol;
    ivec2 block = ivec2(fragCoord.xy) / resol;

    //we only need 1 block, since we can sum all the ranks here
    if(block.x + block.y * 2 >= 1)
        discard;

    vec3 sum = vec3(0.0);
    for (int i = -Nc; i < Nc; i++) 
    {
        //use texelFetch to get the pixel at the current index
        ivec2 pos = ivec2(coord.x, coord.y + i);

        //skip if the coordinate is outside the image
        if(pos.x < 0 || pos.x >= iResolution.x/2 || pos.y < 0 || pos.y >= iResolution.y/2)
            continue;
        
        //loop over pairs of ranks
        for(int j = 0; j < Nr; j += 2)
        {
            int block_id = j / 2;
            ivec2 block_offset = ivec2(block_id % 2, block_id / 2);
            vec4 pixel = texelFetch(iChannel0, pos + block_offset * resol, 0);

            //unpack the values
            vec3 val1, val2;
            unpack2vec3(pixel.xyz, val1, val2);

            //add the values to the sum
            sum += val1 * V[i + Nc + j * N];
            sum += val2 * V[i + Nc + (j + 1) * N];
        }
    }

    //return the sum
    fragColor = vec4(sum, 1.0);
}


//single component convolution pass x
//first pass of separable convolution
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{    
    //do a convolution in the x direction on iChannel0
    ivec2 resol = ivec2(iResolution.xy);
    ivec2 coord = ivec2(fragCoord.xy);

    float sums[Nr];
    for(int i = 0; i < Nr; i++)
        sums[i] = 0.0;

    for (int i = -Nc; i < Nc; i++) 
    {
        ivec2 pos = coord + ivec2(i, 0);

        if(pos.x < 0 || pos.x >= resol.x || pos.y < 0 || pos.y >= resol.y)
            continue;

        float density = texelFetch(iChannel0, pos, 0).z;
        
        //loop over ranks
        for(int j = 0; j < Nr; j++)
        {
            sums[j] += density * U[i + Nc + j * N];
        }
    }

    //write the sum packed to the output
    fragColor = vec4(0.0);

    //loop over pairs of ranks
    for(int j = 0; j < Nr; j += 2)
    {
        fragColor[j / 2] = pack2(sums[j], sums[j + 1]);
    }
}

//second pass: do the convolution in the y direction
//and write the result to the output
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    //do a convolution in the y direction on iChannel0
    ivec2 resol = ivec2(iResolution.xy);
    ivec2 coord = ivec2(fragCoord.xy);

    float sum = 0.0;
    for (int i = -Nc; i < Nc; i++) 
    {
        //use texelFetch to get the pixel at the current index
        ivec2 pos = coord + ivec2(0, i);

        //skip if the coordinate is outside the image
        if(pos.x < 0 || pos.x >= resol.x || pos.y < 0 || pos.y >= resol.y)
            continue;
        
        vec4 data = texelFetch(iChannel0, pos, 0);

        //loop over pairs of ranks
        for(int j = 0; j < Nr; j += 2)
        {
            vec2 val = unpack2(data[j / 2]);
            sum += val.x * V[i + Nc + j * N];
            sum += val.y * V[i + Nc + (j + 1) * N];
        }
    }

    //return the sum
    fragColor = vec4(sum);
}

//The complex filter: 
const vec2 K[516] = vec2[](vec2(7.8033e-02, -3.5162e-03), vec2(-5.9367e-03, 1.3257e-01), vec2(1.5799e-01, 5.7222e-02), vec2(1.3545e-01, 1.3120e-01), vec2(1.6784e-01, 1.0130e-01), vec2(1.9344e-01, 1.8925e-02), vec2(1.3579e-01, 1.3113e-01), vec2(1.7159e-01, 5.5419e-02), vec2(1.3454e-01, 1.0573e-01), vec2(1.0828e-01, 1.2367e-01), vec2(1.2028e-01, 1.0419e-01), vec2(9.4762e-02, 1.2295e-01), vec2(9.8577e-02, 1.0607e-01), vec2(6.4038e-02, 9.1733e-02), vec2(8.9744e-02, 2.2344e-02), vec2(9.5298e-02, -1.9267e-02), vec2(9.8076e-02, -4.8947e-02), vec2(8.9666e-02, -8.2757e-02), vec2(1.1143e-01, -7.3276e-02), vec2(1.2142e-01, -7.6871e-02), vec2(1.5251e-01, -1.8297e-02), vec2(1.4652e-01, -7.2444e-02), vec2(1.4124e-01, -1.0131e-01), vec2(1.5848e-01, -9.3226e-02), vec2(1.6593e-01, -9.9065e-02), vec2(1.6098e-01, -1.1981e-01), vec2(1.3104e-01, -1.5756e-01), vec2(1.7208e-01, -1.1987e-01), vec2(1.2853e-01, -1.6426e-01), vec2(1.8529e-01, -8.5550e-02), vec2(2.0123e-01, -1.0622e-02), vec2(1.4510e-01, -1.2972e-01), vec2(1.7940e-01, -5.2185e-02), vec2(1.6196e-01, -8.2630e-02), vec2(1.6838e-01, -5.1045e-02), vec2(1.6795e-01, 4.2197e-02), vec2(1.6518e-01, 4.3514e-02), vec2(1.7038e-01, 2.4463e-02), vec2(1.7220e-01, 3.5690e-02), vec2(1.2862e-01, 1.3097e-01), vec2(1.5943e-01, 1.0965e-01), vec2(1.7570e-01, 1.1449e-01), vec2(1.6027e-01, 1.5934e-01), vec2(2.0694e-01, 1.2571e-01), vec2(1.9833e-01, 1.7570e-01), vec2(9.7105e-02, 2.6780e-01), vec2(1.2592e-01, 2.7657e-01), vec2(1.4264e-01, 2.9430e-01), vec2(3.0130e-01, 1.7379e-01), vec2(1.8407e-01, 3.1838e-01), vec2(2.0875e-01, 3.2697e-01), vec2(2.5311e-01, 3.2340e-01), vec2(2.3617e-01, 3.6194e-01), vec2(3.2403e-01, 3.1648e-01), vec2(2.1782e-01, 4.1992e-01), vec2(2.6578e-01, 4.1463e-01), vec2(2.5658e-01, 4.4186e-01), vec2(3.4195e-01, 4.0065e-01), vec2(1.2538e-01, 5.2456e-01), vec2(2.7152e-01, 4.7683e-01), vec2(3.0216e-01, 4.6548e-01), vec2(2.9388e-01, 4.7391e-01), vec2(4.7297e-01, 2.9449e-01), vec2(4.7037e-01, 3.0652e-01), vec2(6.8051e-01, 1.2736e-01), vec2(4.0317e-01, 3.9071e-01), vec2(3.8534e-01, 4.0241e-01), vec2(4.5883e-01, 3.1691e-01), vec2(4.4414e-01, 3.3273e-01), vec2(2.8986e-02, 5.4795e-01), vec2(3.3097e-01, 4.2584e-01), vec2(3.2473e-01, 4.1474e-01), vec2(3.0188e-01, 4.1222e-01), vec2(2.6935e-01, 4.1232e-01), vec2(1.9527e-01, 4.3087e-01), vec2(2.8211e-01, 3.5435e-01), vec2(2.2912e-01, 3.6644e-01), vec2(1.8369e-01, 3.6729e-01), vec2(1.2437e-01, 3.6745e-01), vec2(5.7786e-02, 3.6319e-01), vec2(1.6186e-01, 3.0788e-01), vec2(1.4968e-01, 2.9078e-01), vec2(1.8790e-01, 2.3881e-01), vec2(1.6292e-01, 2.3367e-01), vec2(1.4529e-01, 2.2159e-01), vec2(1.4621e-01, 1.9301e-01), vec2(1.2384e-01, 1.8905e-01), vec2(1.4237e-01, 1.5397e-01), vec2(1.6145e-01, 1.0665e-01), vec2(1.4723e-01, 1.0963e-01), vec2(1.7521e-01, 1.5090e-02), vec2(1.6225e-01, 5.7465e-02), vec2(1.6454e-01, 4.5861e-02), vec2(1.6948e-01, -3.5517e-02), vec2(1.7555e-01, 1.1838e-02), vec2(1.6676e-01, -7.2450e-02), vec2(1.7524e-01, -6.4791e-02), vec2(1.6524e-01, -1.0284e-01), vec2(1.2241e-01, -1.6009e-01), vec2(1.8727e-01, -8.1149e-02), vec2(1.2124e-01, -1.6969e-01), vec2(1.1312e-01, -1.7658e-01), vec2(1.5828e-01, -1.3016e-01), vec2(1.0815e-01, -1.6903e-01), vec2(1.6469e-01, -1.0110e-01), vec2(1.5007e-01, -1.0620e-01), vec2(1.4407e-01, -9.7251e-02), vec2(1.3834e-01, -8.7056e-02), vec2(1.3279e-01, -7.7225e-02), vec2(1.3007e-01, -6.1122e-02), vec2(1.3135e-01, -2.3068e-02), vec2(1.0602e-01, -6.0398e-02), vec2(9.5810e-02, -5.3245e-02), vec2(8.6674e-02, -4.4055e-02), vec2(9.2471e-02, -1.6855e-03), vec2(1.0416e-01, 4.0697e-02), vec2(8.9984e-02, 1.1348e-01), vec2(9.4781e-02, 1.2293e-01), vec2(9.7718e-02, 1.2581e-01), vec2(1.2576e-01, 1.0587e-01), vec2(1.5932e-01, 6.2413e-02), vec2(1.1941e-01, 1.3512e-01), vec2(1.5058e-01, 1.1385e-01), vec2(1.5160e-01, 1.2165e-01), vec2(1.7379e-01, 9.0742e-02), vec2(1.8448e-01, 3.9117e-02), vec2(2.3122e-02, 1.6642e-01), vec2(1.2933e-01, 2.9864e-02), vec2(2.4712e-02, 2.8660e-02), vec2(9.3864e-03, 2.6798e-02), vec2(-3.2414e-02, -1.2782e-02), vec2(-2.2546e-02, 1.8302e-02), vec2(-1.5880e-03, 3.4607e-03), vec2(4.0220e-02, -7.2688e-03), vec2(5.7986e-02, -5.1965e-02), vec2(1.0983e-01, -1.6674e-02), vec2(1.1149e-01, -8.3110e-02), vec2(1.5254e-01, -4.6029e-02), vec2(1.7426e-01, -2.2705e-02), vec2(1.6914e-01, -4.5279e-02), vec2(1.2200e-01, 7.5358e-03), vec2(1.0446e-01, 2.1118e-02), vec2(6.1637e-02, 9.8972e-02), vec2(6.9651e-02, 1.3046e-01), vec2(8.6983e-02, 1.4604e-01), vec2(1.0844e-01, 1.4842e-01), vec2(1.4570e-01, 1.2651e-01), vec2(1.2768e-01, 1.5276e-01), vec2(1.3231e-01, 1.5607e-01), vec2(5.7881e-02, 2.0328e-01), vec2(1.2876e-01, 1.7811e-01), vec2(1.6121e-01, 1.6398e-01), vec2(1.5364e-01, 1.8464e-01), vec2(1.6442e-01, 1.8929e-01), vec2(1.9291e-01, 1.7537e-01), vec2(2.3812e-01, 1.2619e-01), vec2(2.0394e-01, 1.9027e-01), vec2(2.5876e-01, 1.2202e-01), vec2(1.8234e-01, 2.2825e-01), vec2(9.2028e-02, 2.8454e-01), vec2(2.5946e-01, 1.5737e-01), vec2(1.8098e-01, 2.4711e-01), vec2(2.3769e-01, 1.9820e-01), vec2(2.1764e-01, 2.2084e-01), vec2(9.6459e-02, 2.9552e-01), vec2(1.2417e-01, 2.8236e-01), vec2(1.8032e-01, 2.4772e-01), vec2(1.9555e-01, 2.2763e-01), vec2(6.6748e-02, 2.8744e-01), vec2(1.3959e-01, 2.5552e-01), vec2(1.7198e-01, 2.2051e-01), vec2(1.3796e-01, 2.3492e-01), vec2(1.9425e-01, 1.8360e-01), vec2(1.6710e-01, 1.9392e-01), vec2(6.1975e-02, 2.4165e-01), vec2(8.5163e-02, 2.2980e-01), vec2(9.6310e-02, 2.1735e-01), vec2(2.0144e-01, 1.1939e-01), vec2(1.1774e-01, 2.0099e-01), vec2(1.3016e-01, 1.9347e-01), vec2(1.5293e-01, 1.7562e-01), vec2(1.4397e-01, 1.8709e-01), vec2(1.9021e-01, 1.5080e-01), vec2(1.4483e-01, 2.0781e-01), vec2(1.7875e-01, 2.0024e-01), vec2(1.8971e-01, 2.1848e-01), vec2(2.5385e-01, 1.9246e-01), vec2(1.6262e-01, 3.1851e-01), vec2(2.8386e-01, 2.9214e-01), vec2(3.4706e-01, 3.1343e-01), vec2(3.8895e-01, 3.7177e-01), vec2(5.8227e-01, 2.0172e-01), vec2(6.2961e-01, 2.7422e-01), vec2(1.6895e-01, 4.0820e-01), vec2(5.6562e-01, 3.8943e-01), vec2(5.1419e-01, 3.3960e-01), vec2(5.0515e-01, 1.8527e-01), vec2(4.3386e-01, 1.7452e-01), vec2(1.2280e-01, 3.8839e-01), vec2(2.7972e-01, 2.2282e-01), vec2(2.4550e-01, 2.0300e-01), vec2(2.1176e-01, 1.9719e-01), vec2(1.8048e-01, 1.9869e-01), vec2(1.3361e-01, 2.1519e-01), vec2(1.6996e-01, 1.7331e-01), vec2(1.4032e-01, 1.8984e-01), vec2(1.1493e-01, 2.0256e-01), vec2(8.0042e-02, 2.1900e-01), vec2(3.7909e-02, 2.2983e-01), vec2(1.0651e-01, 2.0854e-01), vec2(1.0152e-01, 2.1497e-01), vec2(1.3723e-01, 2.0307e-01), vec2(1.2224e-01, 2.1746e-01), vec2(1.1037e-01, 2.3096e-01), vec2(1.1313e-01, 2.4216e-01), vec2(8.6390e-02, 2.5837e-01), vec2(1.1283e-01, 2.5589e-01), vec2(1.4434e-01, 2.5287e-01), vec2(1.1014e-01, 2.7376e-01), vec2(2.2107e-01, 2.0295e-01), vec2(1.2814e-01, 2.7831e-01), vec2(1.2014e-01, 2.8410e-01), vec2(2.1603e-01, 2.2354e-01), vec2(1.2539e-01, 2.8357e-01), vec2(2.2497e-01, 2.1252e-01), vec2(1.9807e-01, 2.3362e-01), vec2(2.2854e-01, 1.9964e-01), vec2(2.7629e-01, 1.1446e-01), vec2(1.7688e-01, 2.3251e-01), vec2(2.6383e-01, 1.1062e-01), vec2(2.6122e-01, 9.7767e-02), vec2(2.1019e-01, 1.6865e-01), vec2(2.4249e-01, 9.5734e-02), vec2(1.6674e-01, 1.8725e-01), vec2(1.6861e-01, 1.7107e-01), vec2(1.5647e-01, 1.6852e-01), vec2(1.4633e-01, 1.6400e-01), vec2(1.3367e-01, 1.6372e-01), vec2(1.1177e-01, 1.7139e-01), vec2(5.6610e-02, 1.9087e-01), vec2(1.1340e-01, 1.5613e-01), vec2(1.1488e-01, 1.4352e-01), vec2(1.2288e-01, 1.1745e-01), vec2(1.0109e-01, 1.0795e-01), vec2(1.0620e-01, 4.8183e-02), vec2(1.0252e-01, 2.9193e-02), vec2(1.2235e-01, 7.4439e-03), vec2(1.6980e-01, -9.9199e-03), vec2(1.6876e-01, -4.8809e-02), vec2(1.3276e-01, -8.8047e-02), vec2(1.3827e-01, -1.4697e-02), vec2(1.0702e-01, -2.9767e-02), vec2(7.6974e-02, -1.1847e-02), vec2(3.9688e-02, -9.7502e-03), vec2(4.8521e-04, 3.7722e-03), vec2(-2.6705e-02, -1.1381e-02), vec2(-1.8497e-02, 2.9551e-02), vec2(-6.2024e-03, 4.9840e-02), vec2(-2.7878e-02, -1.6046e-02), vec2(1.8511e-02, -3.0806e-02), vec2(-2.7836e-02, -1.9845e-02), vec2(-3.2493e-02, -2.0923e-02), vec2(-5.2068e-02, 3.6690e-03), vec2(-3.7953e-02, 3.0688e-02), vec2(-2.6219e-02, 4.3786e-03), vec2(6.1092e-03, 7.6991e-03), vec2(3.2029e-02, 1.2180e-02), vec2(5.6061e-02, 2.3768e-02), vec2(9.2551e-02, 1.1001e-02), vec2(1.6358e-01, 1.6288e-02), vec2(1.9310e-01, 8.9260e-03), vec2(2.0121e-01, 1.0889e-01), vec2(2.5219e-01, 8.6307e-02), vec2(2.7329e-01, 5.9290e-02), vec2(2.7451e-01, 2.7632e-02), vec2(2.6359e-01, -3.1746e-02), vec2(2.4931e-01, 1.6139e-02), vec2(2.3089e-01, 2.3064e-02), vec2(1.7980e-01, 1.1813e-01), vec2(1.9144e-01, 5.9333e-02), vec2(1.8657e-01, 4.0101e-02), vec2(1.6515e-01, 6.7862e-02), vec2(1.4725e-01, 7.8381e-02), vec2(1.3879e-01, 7.6263e-02), vec2(1.4489e-01, 5.6267e-02), vec2(1.1548e-01, 1.0476e-01), vec2(1.4301e-01, 7.7030e-02), vec2(9.6349e-02, 1.4381e-01), vec2(3.7010e-02, 1.7958e-01), vec2(1.5803e-01, 1.1801e-01), vec2(1.1381e-01, 1.7920e-01), vec2(1.6672e-01, 1.5215e-01), vec2(1.6814e-01, 1.7248e-01), vec2(8.5091e-02, 2.3985e-01), vec2(1.2587e-01, 2.3916e-01), vec2(1.9086e-01, 2.1090e-01), vec2(2.2916e-01, 1.9630e-01), vec2(1.2816e-01, 2.8944e-01), vec2(2.1682e-01, 2.4798e-01), vec2(2.7741e-01, 2.1019e-01), vec2(2.6480e-01, 2.4670e-01), vec2(3.3724e-01, 1.6024e-01), vec2(3.3806e-01, 1.9225e-01), vec2(2.3027e-01, 3.2706e-01), vec2(2.7321e-01, 3.0423e-01), vec2(3.0286e-01, 2.8991e-01), vec2(4.2219e-01, 5.8910e-02), vec2(3.4458e-01, 2.5899e-01), vec2(3.6127e-01, 2.4084e-01), vec2(3.8869e-01, 1.9969e-01), vec2(3.7148e-01, 2.2963e-01), vec2(4.0927e-01, 1.4280e-01), vec2(3.3649e-01, 2.6276e-01), vec2(3.4674e-01, 2.3110e-01), vec2(3.1799e-01, 2.4598e-01), vec2(3.3036e-01, 1.9071e-01), vec2(1.7086e-01, 3.1042e-01), vec2(2.0081e-01, 2.5015e-01), vec2(1.5159e-01, 2.3784e-01), vec2(6.6222e-02, 2.3250e-01), vec2(5.8810e-02, 2.0414e-01), vec2(-8.8347e-02, 2.0675e-01), vec2(-7.1377e-02, -2.0455e-01), vec2(-1.2621e-01, 1.8606e-01), vec2(6.4406e-03, 2.1235e-01), vec2(1.5365e-01, 1.8664e-01), vec2(2.2431e-01, 1.7098e-01), vec2(6.7375e-02, 3.1363e-01), vec2(2.8388e-01, 2.1206e-01), vec2(3.2203e-01, 2.0448e-01), vec2(3.4222e-01, 2.1096e-01), vec2(3.4872e-01, 2.2810e-01), vec2(3.2210e-01, 2.8021e-01), vec2(3.8830e-01, 1.9264e-01), vec2(3.6696e-01, 2.3678e-01), vec2(3.4116e-01, 2.7305e-01), vec2(2.9310e-01, 3.2035e-01), vec2(2.2891e-01, 3.6526e-01), vec2(3.2542e-01, 2.7534e-01), vec2(3.0975e-01, 2.8253e-01), vec2(3.3754e-01, 2.3078e-01), vec2(3.0688e-01, 2.5656e-01), vec2(2.7580e-01, 2.7420e-01), vec2(2.5471e-01, 2.7302e-01), vec2(2.0806e-01, 2.9614e-01), vec2(2.1756e-01, 2.7165e-01), vec2(2.2141e-01, 2.4389e-01), vec2(1.7113e-01, 2.6630e-01), vec2(2.5073e-01, 1.6786e-01), vec2(1.4571e-01, 2.4428e-01), vec2(1.2245e-01, 2.4093e-01), vec2(1.8145e-01, 1.7846e-01), vec2(9.6222e-02, 2.2082e-01), vec2(1.5700e-01, 1.6217e-01), vec2(1.2625e-01, 1.7066e-01), vec2(1.3539e-01, 1.4342e-01), vec2(1.6069e-01, 8.8286e-02), vec2(9.2909e-02, 1.4605e-01), vec2(1.4625e-01, 7.0714e-02), vec2(1.4678e-01, 5.2612e-02), vec2(1.3176e-01, 8.2468e-02), vec2(1.5679e-01, 2.2208e-02), vec2(1.4820e-01, 7.6558e-02), vec2(1.7028e-01, 5.3733e-02), vec2(1.8536e-01, 4.5399e-02), vec2(1.9648e-01, 3.9445e-02), vec2(2.1189e-01, 3.7233e-02), vec2(2.2623e-01, 5.1718e-02), vec2(2.2247e-01, 1.1370e-01), vec2(2.6397e-01, 2.8539e-02), vec2(2.7549e-01, 1.5499e-02), vec2(2.7919e-01, -1.5977e-02), vec2(2.6594e-01, 1.7990e-02), vec2(2.2787e-01, -2.1081e-02), vec2(1.9177e-01, 2.3994e-02), vec2(1.6327e-01, 1.6339e-02), vec2(9.3779e-02, 2.7057e-02), vec2(5.9091e-02, 1.4963e-02), vec2(3.4185e-02, 2.4138e-03), vec2(1.3352e-03, 9.7323e-03), vec2(-2.5501e-02, 7.4968e-03), vec2(-4.8507e-02, 4.9073e-03), vec2(-5.1735e-02, 6.8650e-03), vec2(-3.8654e-02, -4.4996e-04), vec2(4.5480e-03, -3.3879e-02), vec2(-2.6949e-02, -2.3830e-02), vec2(4.1352e-02, -1.1384e-02), vec2(-1.5993e-03, 2.8030e-03), vec2(-1.8635e-02, -1.0771e-03), vec2(-1.2721e-02, 3.6453e-02), vec2(-4.6442e-02, 4.1117e-02), vec2(-5.5826e-02, 6.8541e-02), vec2(-1.6294e-02, 1.2263e-01), vec2(-1.0851e-01, 1.2571e-01), vec2(-3.9954e-02, 1.9389e-01), vec2(-1.0381e-01, 1.8934e-01), vec2(-1.3049e-01, 1.8295e-01), vec2(-8.5747e-02, 2.0859e-01), vec2(-1.0265e-01, 2.0179e-01), vec2(-7.3032e-02, 2.1568e-01), vec2(-1.7631e-01, 1.7876e-01), vec2(-1.5288e-01, 2.1235e-01), vec2(-1.1970e-01, 2.4036e-01), vec2(-7.7720e-02, 2.7128e-01), vec2(4.5249e-03, 2.9849e-01), vec2(-2.7657e-02, 3.1424e-01), vec2(-1.0724e-02, 3.3068e-01), vec2(-1.3991e-01, 3.1457e-01), vec2(-1.4552e-02, 3.5621e-01), vec2(4.6618e-02, 3.6618e-01), vec2(2.0576e-02, 3.7988e-01), vec2(3.2950e-02, 3.8856e-01), vec2(8.7975e-02, 3.8622e-01), vec2(1.9714e-01, 3.4527e-01), vec2(1.1487e-01, 3.8328e-01), vec2(2.4659e-01, 3.0840e-01), vec2(1.1595e-01, 3.6853e-01), vec2(2.9666e-03, 3.8121e-01), vec2(2.6820e-01, 2.5793e-01), vec2(1.7242e-01, 3.1930e-01), vec2(2.5692e-01, 2.4770e-01), vec2(2.4525e-01, 2.4937e-01), vec2(1.2549e-01, 3.2155e-01), vec2(1.7255e-01, 2.9190e-01), vec2(2.4068e-01, 2.3334e-01), vec2(2.6475e-01, 1.9488e-01), vec2(1.5555e-01, 2.8544e-01), vec2(2.3385e-01, 2.2330e-01), vec2(2.6805e-01, 1.6808e-01), vec2(2.4783e-01, 1.9233e-01), vec2(2.9360e-01, 1.0742e-01), vec2(2.7975e-01, 1.2809e-01), vec2(1.9687e-01, 2.3375e-01), vec2(2.2201e-01, 2.0794e-01), vec2(2.3319e-01, 1.8939e-01), vec2(2.9721e-01, 1.7685e-02), vec2(2.4895e-01, 1.5831e-01), vec2(2.5433e-01, 1.4259e-01), vec2(2.6335e-01, 1.1190e-01), vec2(2.4828e-01, 1.3027e-01), vec2(2.6419e-01, 7.1948e-02), vec2(2.2086e-01, 1.4847e-01), vec2(2.2396e-01, 1.2616e-01), vec2(2.0658e-01, 1.3412e-01), vec2(2.1354e-01, 9.6431e-02), vec2(1.3090e-01, 1.7807e-01), vec2(1.5685e-01, 1.3474e-01), vec2(1.4442e-01, 1.2845e-01), vec2(1.2195e-01, 1.3892e-01), vec2(1.6446e-01, 9.6410e-02), vec2(1.7567e-01, 1.5357e-01), vec2(4.2450e-01, 7.2618e-01), vec2(1.4311e-01, 1.8431e-01), vec2(1.3548e-01, 1.3413e-01), vec2(1.6732e-01, 7.8581e-02), vec2(1.7986e-01, 7.0770e-02), vec2(7.9812e-02, 1.9075e-01), vec2(1.9279e-01, 1.0806e-01), vec2(2.0928e-01, 1.0536e-01), vec2(2.1961e-01, 1.1151e-01), vec2(2.2504e-01, 1.2422e-01), vec2(2.1269e-01, 1.5996e-01), vec2(2.5318e-01, 1.0428e-01), vec2(2.4571e-01, 1.3505e-01), vec2(2.3581e-01, 1.6208e-01), vec2(2.1279e-01, 1.9934e-01), vec2(1.7571e-01, 2.3699e-01), vec2(2.4180e-01, 1.7373e-01), vec2(2.3768e-01, 1.8373e-01), vec2(2.6499e-01, 1.4938e-01), vec2(2.5052e-01, 1.7502e-01), vec2(2.3636e-01, 1.9698e-01), vec2(2.3352e-01, 2.0784e-01), vec2(2.0269e-01, 2.3942e-01), vec2(2.1880e-01, 2.2854e-01), vec2(2.3798e-01, 2.1889e-01), vec2(1.9758e-01, 2.5815e-01), vec2(2.8592e-01, 1.6225e-01), vec2(1.9015e-01, 2.7607e-01), vec2(1.6838e-01, 2.9433e-01), vec2(2.5351e-01, 2.3424e-01), vec2(1.4113e-01, 3.2004e-01), vec2(2.4111e-01, 2.6311e-01), vec2(1.9467e-01, 3.0623e-01), vec2(2.1986e-01, 3.0021e-01), vec2(2.9211e-01, 2.4497e-01), vec2(1.0718e-01, 3.7119e-01), vec2(2.5982e-01, 2.9731e-01), vec2(2.5277e-01, 3.1016e-01), vec2(1.2883e-01, 3.7612e-01), vec2(2.1898e-01, 3.3009e-01), vec2(3.7738e-02, 3.8812e-01), vec2(5.2490e-02, 3.7677e-01), vec2(3.6166e-02, 3.6736e-01), vec2(2.1966e-02, 3.5583e-01), vec2(-3.8568e-03, 3.4427e-01), vec2(-5.1958e-02, 3.2675e-01), vec2(-1.5001e-01, 2.7750e-01), vec2(-6.2927e-02, 2.9180e-01), vec2(-6.5690e-02, 2.7443e-01), vec2(-5.1047e-02, 2.6362e-01), vec2(-9.2611e-02, 2.4472e-01), vec2(-4.7519e-02, 2.4657e-01), vec2(-8.9646e-02, 2.0930e-01), vec2(-1.0267e-01, 2.0178e-01), vec2(-1.2481e-01, 1.8879e-01), vec2(-1.0130e-01, 2.0061e-01), vec2(-4.4770e-02, 2.1125e-01), vec2(-1.3318e-01, 1.4648e-01), vec2(-9.2603e-02, 1.3785e-01), vec2(-8.0729e-02, 9.3810e-02), vec2(-5.1488e-02, 7.1875e-02), vec2(-1.7479e-02, 5.9520e-02), vec2(-3.8167e-02, 5.7530e-03), vec2(-4.4117e-03, 1.8182e-02), vec2(-5.3527e-02, -4.2471e-03));//center of the convolution 
//center of the convolution
const int Nc = 64; 
//number of ranks 
const int Nr = 4; 
//convolution size 
const int N = 129; 

vec2 cmul(vec2 a, vec2 b) 
{ 
    return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); 
}

//single component convolution pass x
//first pass of separable complex convolution
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{    
    //do a convolution in the x direction on iChannel0
    ivec2 resol = ivec2(iResolution.xy);
    ivec2 coord = ivec2(fragCoord.xy);

    vec2 s0 = vec2(0.0);
    vec2 s1 = vec2(0.0);
    vec2 s2 = vec2(0.0);
    vec2 s3 = vec2(0.0);

    for (int i = -Nc; i < Nc; i++) 
    {
        ivec2 pos = coord + ivec2(i, 0);

        if(pos.x < 0 || pos.x >= resol.x || pos.y < 0 || pos.y >= resol.y)
            continue;

        float density = texelFetch(iChannel0, pos, 0).z;
        
        int id = i + Nc;
        s0 += cmul(K[id + 0*N], vec2(density, 0.0));
        s1 += cmul(K[id + 1*N], vec2(density, 0.0));
        s2 += cmul(K[id + 2*N], vec2(density, 0.0));
        s3 += cmul(K[id + 3*N], vec2(density, 0.0));
    }

    //write the sum packed to the output
    fragColor = vec4(pack2(s0), pack2(s1), pack2(s2), pack2(s3));
}

//single component convolution pass y
//second pass of separable complex convolution
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{    
    //do a convolution in the y direction on iChannel0
    ivec2 resol = ivec2(iResolution.xy);
    ivec2 coord = ivec2(fragCoord.xy);

    vec2 s = vec2(0.0);

    for (int i = -Nc; i < Nc; i++) 
    {
        ivec2 pos = coord + ivec2(0, i);

        if(pos.x < 0 || pos.x >= resol.x || pos.y < 0 || pos.y >= resol.y)
            continue;

        vec4 data = texelFetch(iChannel0, pos, 0);
        
        int id = i + Nc;
        s += cmul(K[id + 0*N], unpack2(data.x));
        s += cmul(K[id + 1*N], unpack2(data.y));
        s += cmul(K[id + 2*N], unpack2(data.z));
        s += cmul(K[id + 3*N], unpack2(data.w));
    }

    //write the sum packed to the output
    float result = dot(s,s);
    fragColor = vec4(result);
}
