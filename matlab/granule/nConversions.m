function dy = nConversions(t, y)

KNH4 = 0.0024;
KNO2 = 0.000238;

C1 = 0.035;
C2 = 10;

dy = zeros(3,1);    % a column vector

dy(1) = -C1 * y(1) / (y(1) + KNH4);
dy(3) = C2 * y(2) / (y(2) + KNO2)* KNH4 / (y(1) + KNH4);
dy(2) = - dy(1) - dy(3) ;
