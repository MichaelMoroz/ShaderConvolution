//The first 8 columns of the U matrix: 
const float U[512] = float[](-4.8700E-04, -1.2393E-03, -2.1780E-03, -3.2455E-03, -4.4319E-03, -5.7229E-03, -7.0822E-03, -8.5081E-03, -9.9956E-03, -1.1524E-02, -1.3115E-02, -1.4747E-02, -1.6437E-02, -1.8165E-02, -1.9935E-02, -2.1748E-02, -2.3615E-02, -2.5521E-02, -2.7467E-02, -2.9454E-02, -3.1483E-02, -3.3552E-02, -3.5656E-02, -3.7809E-02, -4.0008E-02, -4.2253E-02, -4.4544E-02, -4.6877E-02, -4.9262E-02, -5.1697E-02, -5.4179E-02, -5.6717E-02, -5.9309E-02, -6.1955E-02, -6.4660E-02, -6.7423E-02, -7.0250E-02, -7.3139E-02, -7.6096E-02, -7.9121E-02, -8.2219E-02, -8.5393E-02, -8.8645E-02, -9.1979E-02, -9.5400E-02, -9.8911E-02, -1.0252E-01, -1.0622E-01, -1.1003E-01, -1.1395E-01, -1.1799E-01, -1.2215E-01, -1.2644E-01, -1.3087E-01, -1.3545E-01, -1.4019E-01, -1.4509E-01, -1.5019E-01, -1.5549E-01, -1.6101E-01, -1.6677E-01, -1.7282E-01, -1.7919E-01, -1.8594E-01, 3.8441E-03, 7.0405E-03, 9.1937E-03, 1.0459E-02, 1.0810E-02, 1.0314E-02, 9.3210E-03, 7.8464E-03, 5.9459E-03, 3.8550E-03, 1.3471E-03, -1.3269E-03, -4.3436E-03, -7.4768E-03, -1.0742E-02, -1.4145E-02, -1.7780E-02, -2.1471E-02, -2.5218E-02, -2.9011E-02, -3.2834E-02, -3.6669E-02, -4.0472E-02, -4.4304E-02, -4.8134E-02, -5.1941E-02, -5.5708E-02, -5.9415E-02, -6.3064E-02, -6.6634E-02, -7.0117E-02, -7.3489E-02, -7.6735E-02, -7.9851E-02, -8.2796E-02, -8.5573E-02, -8.8129E-02, -9.0468E-02, -9.2529E-02, -9.4310E-02, -9.5762E-02, -9.6811E-02, -9.7448E-02, -9.7602E-02, -9.7180E-02, -9.6145E-02, -9.4400E-02, -9.1818E-02, -8.8317E-02, -8.3752E-02, -7.7957E-02, -7.0722E-02, -6.1851E-02, -5.1066E-02, -3.8035E-02, -2.2350E-02, -3.4984E-03, 1.9136E-02, 4.6374E-02, 7.9306E-02, 1.1943E-01, 1.6891E-01, 2.3113E-01, 3.1208E-01, 3.8780E-03, 1.3323E-02, 2.3249E-02, 3.1961E-02, 3.8579E-02, 4.2719E-02, 4.5218E-02, 4.6011E-02, 4.5227E-02, 4.3743E-02, 4.0658E-02, 3.6968E-02, 3.1993E-02, 2.6639E-02, 2.0864E-02, 1.4672E-02, 7.8253E-03, 9.2183E-04, -5.9968E-03, -1.2866E-02, -1.9612E-02, -2.6177E-02, -3.2618E-02, -3.8747E-02, -4.4527E-02, -4.9918E-02, -5.4911E-02, -5.9772E-02, -6.3961E-02, -6.7619E-02, -7.1102E-02, -7.3763E-02, -7.5807E-02, -7.7626E-02, -7.8602E-02, -7.9264E-02, -7.8991E-02, -7.8320E-02, -7.6672E-02, -7.4539E-02, -7.1693E-02, -6.7631E-02, -6.2852E-02, -5.7139E-02, -5.0033E-02, -4.1961E-02, -3.2719E-02, -2.1960E-02, -9.9851E-03, 3.3321E-03, 1.8084E-02, 3.4415E-02, 5.2011E-02, 7.0739E-02, 9.0328E-02, 1.1027E-01, 1.2973E-01, 1.4719E-01, 1.6028E-01, 1.6515E-01, 1.5511E-01, 1.1801E-01, 2.9671E-02, -1.6494E-01, -3.5540E-04, 2.0967E-03, 8.6209E-03, 1.8563E-02, 3.1646E-02, 4.6963E-02, 6.3099E-02, 7.9456E-02, 9.5368E-02, 1.1087E-01, 1.2472E-01, 1.3765E-01, 1.4785E-01, 1.5687E-01, 1.6430E-01, 1.6980E-01, 1.7103E-01, 1.7090E-01, 1.6914E-01, 1.6570E-01, 1.6068E-01, 1.5446E-01, 1.4901E-01, 1.4077E-01, 1.3071E-01, 1.1933E-01, 1.0723E-01, 9.6199E-02, 8.3432E-02, 7.0140E-02, 5.7980E-02, 4.4693E-02, 3.1409E-02, 1.9397E-02, 7.2246E-03, -3.7034E-03, -1.4528E-02, -2.4071E-02, -3.3149E-02, -4.0878E-02, -4.7477E-02, -5.3472E-02, -5.8147E-02, -6.1644E-02, -6.4292E-02, -6.5635E-02, -6.5764E-02, -6.4860E-02, -6.2719E-02, -5.9421E-02, -5.5069E-02, -4.9810E-02, -4.3674E-02, -3.6819E-02, -2.9452E-02, -2.1831E-02, -1.4254E-02, -7.0495E-03, -6.1328E-04, 4.5980E-03, 8.0935E-03, 9.4065E-03, 8.2309E-03, 4.7655E-03, -1.0751E-02, -5.0997E-03, 1.1055E-02, 3.1040E-02, 5.1374E-02, 6.9059E-02, 8.3718E-02, 9.4425E-02, 1.0077E-01, 1.0489E-01, 1.0369E-01, 1.0027E-01, 9.2030E-02, 8.2261E-02, 7.0767E-02, 5.7584E-02, 4.1982E-02, 2.6253E-02, 1.0631E-02, -4.5384E-03, -1.8931E-02, -3.2333E-02, -4.5708E-02, -5.7000E-02, -6.6373E-02, -7.3788E-02, -7.9487E-02, -8.5550E-02, -8.8285E-02, -8.8960E-02, -9.0129E-02, -8.7818E-02, -8.3704E-02, -8.0207E-02, -7.4144E-02, -6.8539E-02, -6.0532E-02, -5.2927E-02, -4.3433E-02, -3.4362E-02, -2.5041E-02, -1.3821E-02, -2.8284E-03, 8.4234E-03, 2.0931E-02, 3.3216E-02, 4.5423E-02, 5.7869E-02, 6.9653E-02, 8.0558E-02, 9.0237E-02, 9.8149E-02, 1.0354E-01, 1.0546E-01, 1.0271E-01, 9.3924E-02, 7.7119E-02, 5.0622E-02, 1.2212E-02, -3.9707E-02, -1.0456E-01, -1.7460E-01, -2.1946E-01, -1.3176E-01, -1.1608E-02, -4.1497E-02, -5.8424E-02, -5.8730E-02, -4.2288E-02, -1.3921E-02, 1.8135E-02, 5.0399E-02, 7.9162E-02, 1.0472E-01, 1.2078E-01, 1.3175E-01, 1.3044E-01, 1.2448E-01, 1.1298E-01, 9.5729E-02, 6.9944E-02, 4.3317E-02, 1.6612E-02, -9.0493E-03, -3.2492E-02, -5.3057E-02, -7.4192E-02, -8.8730E-02, -9.7519E-02, -1.0090E-01, -1.0004E-01, -1.0137E-01, -9.4199E-02, -8.3009E-02, -7.4215E-02, -5.8870E-02, -4.1686E-02, -2.6940E-02, -1.0183E-02, 4.7403E-03, 2.0031E-02, 3.3876E-02, 4.6223E-02, 5.7492E-02, 6.7493E-02, 7.4221E-02, 8.0422E-02, 8.5075E-02, 8.4909E-02, 8.4108E-02, 8.1289E-02, 7.2857E-02, 6.3279E-02, 5.1226E-02, 3.5946E-02, 1.4927E-02, -7.8408E-03, -3.2714E-02, -5.9234E-02, -8.6278E-02, -1.1213E-01, -1.3105E-01, -1.3678E-01, -1.2019E-01, -6.8706E-02, 3.0595E-02, 1.7353E-01, 2.5477E-01, -5.2558E-03, 5.7421E-04, 2.8679E-02, 7.1520E-02, 1.2013E-01, 1.6229E-01, 1.9293E-01, 2.0637E-01, 2.0044E-01, 1.8468E-01, 1.4803E-01, 1.0465E-01, 5.1659E-02, -8.6808E-04, -5.0693E-02, -9.4777E-02, -1.2117E-01, -1.3936E-01, -1.4732E-01, -1.4474E-01, -1.3277E-01, -1.1422E-01, -9.7506E-02, -6.9057E-02, -3.5777E-02, -2.2474E-03, 2.8151E-02, 5.9320E-02, 8.1875E-02, 9.6887E-02, 1.1571E-01, 1.1865E-01, 1.1298E-01, 1.1389E-01, 1.0079E-01, 9.4466E-02, 7.5018E-02, 6.2853E-02, 4.1418E-02, 2.7969E-02, 1.7640E-02, -1.4046E-03, -1.3750E-02, -2.3031E-02, -3.7013E-02, -4.5500E-02, -5.0654E-02, -5.6221E-02, -5.7448E-02, -5.5244E-02, -5.0687E-02, -4.4902E-02, -3.7034E-02, -2.7487E-02, -1.7137E-02, -7.0352E-03, 2.1131E-03, 9.2634E-03, 1.3798E-02, 1.5012E-02, 1.2549E-02, 6.9539E-03, 3.3534E-04, -3.2169E-03, -3.9837E-02, -3.7920E-02, 9.3629E-03, 6.9318E-02, 1.1881E-01, 1.4046E-01, 1.3889E-01, 1.1331E-01, 6.9891E-02, 2.1323E-02, -3.0613E-02, -7.8570E-02, -1.0936E-01, -1.3006E-01, -1.3699E-01, -1.2800E-01, -9.2327E-02, -5.3098E-02, -1.2135E-02, 2.7024E-02, 6.0740E-02, 8.6744E-02, 1.1452E-01, 1.2451E-01, 1.2063E-01, 1.0564E-01, 8.4358E-02, 6.7550E-02, 3.9803E-02, 1.0865E-02, -1.5696E-02, -4.0295E-02, -5.8812E-02, -7.8366E-02, -8.5841E-02, -9.6379E-02, -9.2088E-02, -9.2549E-02, -7.9259E-02, -7.2315E-02, -6.6029E-02, -4.5276E-02, -3.0431E-02, -1.6668E-02, 6.3617E-03, 2.4718E-02, 4.1191E-02, 5.8795E-02, 7.3061E-02, 8.3723E-02, 9.0979E-02, 9.1335E-02, 8.7853E-02, 7.7802E-02, 5.9924E-02, 3.3282E-02, -5.8637E-03, -4.9998E-02, -9.5989E-02, -1.3414E-01, -1.4457E-01, -9.4312E-02, 5.3627E-02, -3.9837E-02); 
//The first 8 rows of the V matrix: 
const float V[512] = float[](-4.8700E-04, -1.2393E-03, -2.1780E-03, -3.2455E-03, -4.4319E-03, -5.7229E-03, -7.0822E-03, -8.5081E-03, -9.9956E-03, -1.1524E-02, -1.3115E-02, -1.4747E-02, -1.6437E-02, -1.8165E-02, -1.9935E-02, -2.1748E-02, -2.3615E-02, -2.5521E-02, -2.7467E-02, -2.9454E-02, -3.1483E-02, -3.3552E-02, -3.5656E-02, -3.7809E-02, -4.0008E-02, -4.2253E-02, -4.4544E-02, -4.6877E-02, -4.9262E-02, -5.1697E-02, -5.4179E-02, -5.6717E-02, -5.9309E-02, -6.1955E-02, -6.4660E-02, -6.7423E-02, -7.0250E-02, -7.3139E-02, -7.6096E-02, -7.9121E-02, -8.2219E-02, -8.5393E-02, -8.8645E-02, -9.1979E-02, -9.5400E-02, -9.8911E-02, -1.0252E-01, -1.0622E-01, -1.1003E-01, -1.1395E-01, -1.1799E-01, -1.2215E-01, -1.2644E-01, -1.3087E-01, -1.3545E-01, -1.4019E-01, -1.4509E-01, -1.5019E-01, -1.5549E-01, -1.6101E-01, -1.6677E-01, -1.7282E-01, -1.7919E-01, -1.8594E-01, 3.8441E-03, 7.0405E-03, 9.1937E-03, 1.0459E-02, 1.0810E-02, 1.0314E-02, 9.3210E-03, 7.8464E-03, 5.9459E-03, 3.8550E-03, 1.3471E-03, -1.3269E-03, -4.3436E-03, -7.4768E-03, -1.0742E-02, -1.4145E-02, -1.7780E-02, -2.1471E-02, -2.5218E-02, -2.9011E-02, -3.2834E-02, -3.6669E-02, -4.0472E-02, -4.4304E-02, -4.8134E-02, -5.1941E-02, -5.5708E-02, -5.9415E-02, -6.3064E-02, -6.6634E-02, -7.0117E-02, -7.3489E-02, -7.6735E-02, -7.9851E-02, -8.2796E-02, -8.5573E-02, -8.8129E-02, -9.0468E-02, -9.2529E-02, -9.4310E-02, -9.5762E-02, -9.6811E-02, -9.7448E-02, -9.7602E-02, -9.7180E-02, -9.6145E-02, -9.4400E-02, -9.1818E-02, -8.8317E-02, -8.3752E-02, -7.7957E-02, -7.0722E-02, -6.1851E-02, -5.1066E-02, -3.8035E-02, -2.2350E-02, -3.4984E-03, 1.9136E-02, 4.6374E-02, 7.9306E-02, 1.1943E-01, 1.6891E-01, 2.3113E-01, 3.1208E-01, 3.8780E-03, 1.3323E-02, 2.3249E-02, 3.1961E-02, 3.8579E-02, 4.2719E-02, 4.5218E-02, 4.6011E-02, 4.5227E-02, 4.3743E-02, 4.0658E-02, 3.6968E-02, 3.1993E-02, 2.6639E-02, 2.0864E-02, 1.4672E-02, 7.8253E-03, 9.2183E-04, -5.9968E-03, -1.2866E-02, -1.9612E-02, -2.6177E-02, -3.2618E-02, -3.8747E-02, -4.4527E-02, -4.9918E-02, -5.4911E-02, -5.9772E-02, -6.3961E-02, -6.7619E-02, -7.1102E-02, -7.3763E-02, -7.5807E-02, -7.7626E-02, -7.8602E-02, -7.9264E-02, -7.8991E-02, -7.8320E-02, -7.6672E-02, -7.4539E-02, -7.1693E-02, -6.7631E-02, -6.2852E-02, -5.7139E-02, -5.0033E-02, -4.1961E-02, -3.2719E-02, -2.1960E-02, -9.9851E-03, 3.3321E-03, 1.8084E-02, 3.4415E-02, 5.2011E-02, 7.0739E-02, 9.0328E-02, 1.1027E-01, 1.2973E-01, 1.4719E-01, 1.6028E-01, 1.6515E-01, 1.5511E-01, 1.1801E-01, 2.9671E-02, -1.6494E-01, 3.5540E-04, -2.0967E-03, -8.6209E-03, -1.8563E-02, -3.1646E-02, -4.6963E-02, -6.3099E-02, -7.9456E-02, -9.5368E-02, -1.1087E-01, -1.2472E-01, -1.3765E-01, -1.4785E-01, -1.5687E-01, -1.6430E-01, -1.6980E-01, -1.7103E-01, -1.7090E-01, -1.6914E-01, -1.6570E-01, -1.6068E-01, -1.5446E-01, -1.4901E-01, -1.4077E-01, -1.3071E-01, -1.1933E-01, -1.0723E-01, -9.6199E-02, -8.3432E-02, -7.0140E-02, -5.7980E-02, -4.4693E-02, -3.1409E-02, -1.9397E-02, -7.2246E-03, 3.7034E-03, 1.4528E-02, 2.4071E-02, 3.3149E-02, 4.0878E-02, 4.7477E-02, 5.3472E-02, 5.8147E-02, 6.1644E-02, 6.4292E-02, 6.5635E-02, 6.5764E-02, 6.4860E-02, 6.2719E-02, 5.9421E-02, 5.5069E-02, 4.9810E-02, 4.3674E-02, 3.6819E-02, 2.9452E-02, 2.1831E-02, 1.4254E-02, 7.0495E-03, 6.1328E-04, -4.5980E-03, -8.0935E-03, -9.4065E-03, -8.2309E-03, -4.7655E-03, -1.0751E-02, -5.0997E-03, 1.1055E-02, 3.1040E-02, 5.1374E-02, 6.9059E-02, 8.3718E-02, 9.4425E-02, 1.0077E-01, 1.0489E-01, 1.0369E-01, 1.0027E-01, 9.2030E-02, 8.2261E-02, 7.0767E-02, 5.7584E-02, 4.1982E-02, 2.6253E-02, 1.0631E-02, -4.5384E-03, -1.8931E-02, -3.2333E-02, -4.5708E-02, -5.7000E-02, -6.6373E-02, -7.3788E-02, -7.9487E-02, -8.5550E-02, -8.8285E-02, -8.8960E-02, -9.0129E-02, -8.7818E-02, -8.3704E-02, -8.0207E-02, -7.4144E-02, -6.8539E-02, -6.0532E-02, -5.2927E-02, -4.3433E-02, -3.4362E-02, -2.5041E-02, -1.3821E-02, -2.8284E-03, 8.4234E-03, 2.0931E-02, 3.3216E-02, 4.5423E-02, 5.7869E-02, 6.9653E-02, 8.0558E-02, 9.0237E-02, 9.8149E-02, 1.0354E-01, 1.0546E-01, 1.0271E-01, 9.3924E-02, 7.7119E-02, 5.0622E-02, 1.2212E-02, -3.9707E-02, -1.0456E-01, -1.7460E-01, -2.1946E-01, -1.3176E-01, -1.1608E-02, -4.1497E-02, -5.8424E-02, -5.8730E-02, -4.2288E-02, -1.3921E-02, 1.8135E-02, 5.0399E-02, 7.9162E-02, 1.0472E-01, 1.2078E-01, 1.3175E-01, 1.3044E-01, 1.2448E-01, 1.1298E-01, 9.5729E-02, 6.9944E-02, 4.3317E-02, 1.6612E-02, -9.0493E-03, -3.2492E-02, -5.3057E-02, -7.4192E-02, -8.8730E-02, -9.7519E-02, -1.0090E-01, -1.0004E-01, -1.0137E-01, -9.4199E-02, -8.3009E-02, -7.4215E-02, -5.8870E-02, -4.1686E-02, -2.6940E-02, -1.0183E-02, 4.7403E-03, 2.0031E-02, 3.3876E-02, 4.6223E-02, 5.7492E-02, 6.7493E-02, 7.4221E-02, 8.0422E-02, 8.5075E-02, 8.4909E-02, 8.4108E-02, 8.1289E-02, 7.2857E-02, 6.3279E-02, 5.1226E-02, 3.5946E-02, 1.4927E-02, -7.8408E-03, -3.2714E-02, -5.9234E-02, -8.6278E-02, -1.1213E-01, -1.3105E-01, -1.3678E-01, -1.2019E-01, -6.8706E-02, 3.0595E-02, 1.7353E-01, 2.5477E-01, 5.2558E-03, -5.7421E-04, -2.8679E-02, -7.1520E-02, -1.2013E-01, -1.6229E-01, -1.9293E-01, -2.0637E-01, -2.0044E-01, -1.8468E-01, -1.4803E-01, -1.0465E-01, -5.1659E-02, 8.6808E-04, 5.0693E-02, 9.4777E-02, 1.2117E-01, 1.3936E-01, 1.4732E-01, 1.4474E-01, 1.3277E-01, 1.1422E-01, 9.7506E-02, 6.9057E-02, 3.5777E-02, 2.2474E-03, -2.8151E-02, -5.9320E-02, -8.1875E-02, -9.6887E-02, -1.1571E-01, -1.1865E-01, -1.1298E-01, -1.1389E-01, -1.0079E-01, -9.4466E-02, -7.5018E-02, -6.2853E-02, -4.1418E-02, -2.7969E-02, -1.7640E-02, 1.4046E-03, 1.3750E-02, 2.3031E-02, 3.7013E-02, 4.5500E-02, 5.0654E-02, 5.6221E-02, 5.7448E-02, 5.5244E-02, 5.0687E-02, 4.4902E-02, 3.7034E-02, 2.7487E-02, 1.7137E-02, 7.0352E-03, -2.1131E-03, -9.2634E-03, -1.3798E-02, -1.5012E-02, -1.2549E-02, -6.9539E-03, -3.3534E-04, 3.2169E-03, -3.9837E-02, -3.7920E-02, 9.3629E-03, 6.9318E-02, 1.1881E-01, 1.4046E-01, 1.3889E-01, 1.1331E-01, 6.9891E-02, 2.1323E-02, -3.0613E-02, -7.8570E-02, -1.0936E-01, -1.3006E-01, -1.3699E-01, -1.2800E-01, -9.2327E-02, -5.3098E-02, -1.2135E-02, 2.7024E-02, 6.0740E-02, 8.6744E-02, 1.1452E-01, 1.2451E-01, 1.2063E-01, 1.0564E-01, 8.4358E-02, 6.7550E-02, 3.9803E-02, 1.0865E-02, -1.5696E-02, -4.0295E-02, -5.8812E-02, -7.8366E-02, -8.5841E-02, -9.6379E-02, -9.2088E-02, -9.2549E-02, -7.9259E-02, -7.2315E-02, -6.6029E-02, -4.5276E-02, -3.0431E-02, -1.6668E-02, 6.3617E-03, 2.4718E-02, 4.1191E-02, 5.8795E-02, 7.3061E-02, 8.3723E-02, 9.0979E-02, 9.1335E-02, 8.7853E-02, 7.7802E-02, 5.9924E-02, 3.3282E-02, -5.8637E-03, -4.9998E-02, -9.5989E-02, -1.3414E-01, -1.4457E-01, -9.4312E-02, 5.3627E-02, -3.9837E-02); 
//center of the convolution 
const int Nc = 64; 
//number of ranks 
const int Nr = 8; 
//convolution size 
const int N = 129; 
//symmetric convolution 
#define SYMMETRIC

float pack2(vec2 a)
{
    return uintBitsToFloat(packHalf2x16(a));
}

vec2 unpack2(float packed)
{
    return unpackHalf2x16(floatBitsToUint(packed));
}

float getDensity(ivec2 p)
{
    return texelFetch(iChannel0, p, 0).z;
}

float getEstimatedDensity(ivec2 p)
{
    //the laplacian of the density field is the source term for the poisson equation
    return texelFetch(iChannel1, p + ivec2(1,0), 0).x + texelFetch(iChannel1, p - ivec2(1,0), 0).x + texelFetch(iChannel1, p + ivec2(0,1), 0).x + texelFetch(iChannel1, p - ivec2(0,1), 0).x - 4.0*texelFetch(iChannel1, p, 0).x;
}

float getDensityError(ivec2 p)
{
    //find how far off the estimated density is from the actual density and use that as the source term
    return getDensity(p) - getEstimatedDensity(p);
}

//single component convolution pass x
//first pass of separable poisson filter convolution
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{    
    //do a convolution in the x direction on iChannel0
    ivec2 resol = ivec2(iResolution.xy);
    ivec2 coord = ivec2(fragCoord.xy);

    vec4 sum0 = vec4(0.0);
    vec4 sum1 = vec4(0.0);

    for (int i = -Nc; i <= Nc; i++) 
    {
        ivec2 pos = coord + ivec2(i, 0);

        //wrap around
        pos.x = (pos.x + resol.x) % resol.x;
        pos.y = (pos.y + resol.y) % resol.y;
       
        float density = getDensityError(pos);

        int id = (i > Nc) ? (2 * Nc - i) : (i + Nc);
        //loop over ranks
        for(int j = 0; j < 4; j++)
        {
            sum0[j] += density * U[id + 2 * j * Nc];
            sum1[j] += density * U[id + (2 * j + 1) * Nc];
        }
    }

    //loop over pairs of ranks
    for(int j = 0; j < 4; j++)
    {
        fragColor[j] = pack2(vec2(sum0[j], sum1[j]));
    }
}

//single component convolution pass y
//second pass of separable poisson filter convolution
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{    
    //do a convolution in the y direction on iChannel0
    ivec2 resol = ivec2(iResolution.xy);
    ivec2 coord = ivec2(fragCoord.xy);

    float sum = 0.0;
    for (int i = -Nc; i <= Nc; i++) 
    {
        //use texelFetch to get the pixel at the current index
        ivec2 pos = coord + ivec2(0, i);
        pos.x = pos.x % resol.x;
        pos.y = pos.y % resol.y;
        
        vec4 data = texelFetch(iChannel0, pos, 0);

        int id = (i > Nc) ? (2 * Nc - i) : (i + Nc);
        //loop over pairs of ranks
        for(int j = 0; j < 4; j++)
        {
            vec2 val = unpack2(data[j / 2]);
            sum += val.x * V[id + 2 * j * Nc];
            sum += val.y * V[id + (2 * j + 1) * Nc];
        }
    }

    //get previous poisson solution
    float prev = texelFetch(iChannel1, coord, 0).x;

    //add the previous poisson solution to the convolution
    sum += prev;

    //return the final poisson solution
    fragColor = vec4(sum);
}