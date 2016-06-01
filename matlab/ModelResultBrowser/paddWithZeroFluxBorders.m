function up = paddWithZeroFluxBorders(u)

up = zeros(size(u) + 2);
%center
up(2:end-1, 2:end-1) = u;

%sides
up(2:end-1, 1) = u(1:end, 1); 
up(2:end-1, end) = u(1:end, end); 

%top+bottom
up(1, 2:end-1) = u(1, 1:end); 
up(end, 2:end-1) = u(end, 1:end); 

%corners
up(1,1)=up(2,1);
up(1,end)=up(2,end);
up(end,1)=up(end-1,1);
up(end,end)=up(end-1,end);