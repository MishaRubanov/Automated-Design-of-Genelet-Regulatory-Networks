function P = topology_plotter_SM2(rtop)
% close all;
% global P;
B = rtop;
A = rtop;
A(rtop<0) = 0;
B(rtop>0) = 0;
T = A+B;
G = digraph(T');
% figure;
P = plot(G, 'Layout', 'Layered', 'Direction', 'right','Sources', 1, 'Sinks', 3);

[Ax, Ay] = find(A);
[Bx, By] = find(B);
highlight(P, Ay', Ax', 'EdgeColor', 'r')
highlight(P, By', Bx','EdgeColor', 'b')

rCtext = cell(1,length(Ax));
rRtext = cell(1,length(Bx));

for i = 1:length(Ax)
    rCtext{i} = ['rC', num2str(Ax(i))];
end
for i = 1:length(Bx)
    rRtext{i} = ['rR', num2str(Bx(i))];
end
if any(Ax)
    labeledge(P, Ay', Ax', rCtext);
end
if any(Bx)
    labeledge(P, By', Bx', rRtext);
end

P.ArrowSize = 30;
P.ArrowPosition = 0.98;
P.MarkerSize = 5;
P.LineWidth = 2;
P.NodeFontSize = 30;
P.EdgeFontSize = 15;
end