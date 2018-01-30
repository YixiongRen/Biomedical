function out=Random_Choose(max,n)
out=randi([1,max],[1,n]);
out=unique(out);
while length(out)<n
    out=randi([1,max],[1,n]);
    out=unique(out);
end
end