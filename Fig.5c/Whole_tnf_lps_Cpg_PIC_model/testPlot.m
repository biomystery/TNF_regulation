figure('position',[10 50 800 800])
subplot 441
plot(t1,r1(:,48))
title('c1')

subplot 442
plot(t1,r1(:,49))
title('c1off')

subplot 443
plot(t1,r1(:,44))
title('tnf')

subplot 444
plot(t1,r1(:,45))
title('tnfrm')

subplot 445
plot(t1,r1(:,46))
title('TNFR')

subplot 446
plot(t1,r1(:,47))
title('TNFRtnf')

subplot 447
plot(t1,r1(:,50))
title('C1tnf')

subplot 448
plot(t1,r1(:,51))
title('C1tnfoff')

subplot 449
plot(t1,r1(:,52))
title('TTR')

subplot(4,4,10)
plot(t1,r1(:,53))
title('A20')

subplot(4,4,11)
plot(t1,r1(:,54))
title('A20t')
