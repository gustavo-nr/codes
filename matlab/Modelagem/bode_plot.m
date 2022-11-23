function bode_plot(A,B,C,D)
%%
%Definição do sistema
res = ss(A,B,C,D);
%%
%Diagramas de Bode
input={'on','off','off'};
output={'on','off','off','off','off','off'};
for in=1:3
    if in>1
        input{in-1}='off';
        input{in}='on';
    end
    for ou=1:6
        if ou>1
            output{ou-1}='off';
            output{ou}='on';
        elseif ou==1
            output{end}='off';
            output{1}='on';
        end
        figure(10*in+ou)
        optbode=bodeoptions('cstprefs');
        optbode.InputVisible=input;
        optbode.OutputVisible=output;
        optbode.YLabel.String={'Magnitude','Fase'};
        optbode.XLabel.String='Frequência';
        optbode.Title.String='Diagrama de Bode';
        bodao=bodeplot(res,{10^-30,10^10},optbode);
        grid on
    end
end