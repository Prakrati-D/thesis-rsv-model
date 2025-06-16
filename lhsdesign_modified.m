
lb = [0 0 0 0 6 0 0 4 0 0 0 0];
ub = [100 100 100 100 12 100 100 8 100 100 100 1];
[x,y] = lhsdesign_modifie(10,lb,ub);

function [X_normalized, X_scaled] = lhsdesign_modifie(n,min_ranges_p,max_ranges_p)
p=length(min_ranges_p);
[M,N]=size(min_ranges_p);
if M<N
    min_ranges_p=min_ranges_p';
end
    
[M,N]=size(max_ranges_p);
if M<N
    max_ranges_p=max_ranges_p';
end

slope=max_ranges_p-min_ranges_p;
offset=min_ranges_p;

SLOPE=ones(p,n);
OFFSET=ones(p,n);

for i=1:p
    SLOPE(i,:)=ones(1,n).*slope(i);
    OFFSET(i,:)=ones(1,n).*offset(i);
end
X_normalized = lhsdesign(p,n);

X_scaled=SLOPE.*X_normalized+OFFSET;
end