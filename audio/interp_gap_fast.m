function y_intrp = interp_gap_fast(x,y,gap)

x=reshape(x,1,[]);
y=reshape(y,1,[]);

D = diff(x);
E = find(D > gap);
S = [1,E(1:end-1)+1];

ne = length(S);

Z = cell(1,ne);

for ie=1:ne
    y_s = find(y>x(S(ie)),1);
    y_e = find(y>x(E(ie)),1)-1;
    Z{ie} = y_s : y_e;
end
y_intrp = ismember(1:length(y),cell2mat(Z));



    