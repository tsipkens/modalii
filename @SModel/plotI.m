
function [] = plotI(smodel,T,C,data,n)

Jmod = @(l) C(n).*(smodel.blackbody(T(n),l).*smodel.prop.Em(l,smodel.prop.dp0)./...
    (l.*1e-9)./smodel.data_sc)';
J = squeeze(data(n,1,:));

if length(smodel.l)<4
    plotSpec = 'o';
else
    plotSpec = '.';
end

eval_vec = floor(min(smodel.l)/50)*50:1:ceil(max(smodel.l)/50)*50;
plot(eval_vec,Jmod(eval_vec));
hold on;
plot(smodel.l,J,plotSpec,'markers',4);
hold off;

end

