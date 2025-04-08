function dh=ellipcor(latitude)
%the function calculate ellipsodial correction bewteen two ellipsoids by
%giving the latitude of point

% continue later

% t/p
a1=6378136.3;
e1=0.081819221456;

% GRS80
a2=6378137.0;
e2=0.081819190842621;


h1=(a1*(1-e1^2))/sqrt(1-(e1^2)*((sin(latitude))^2));
h2=(a2*(1-e2^2))/sqrt(1-(e2^2)*((sin(latitude))^2));

dh=h1-h2;
end
