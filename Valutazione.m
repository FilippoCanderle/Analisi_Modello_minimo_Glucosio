function[t,y] =Valutazione(parameters,time,G0,x0,glucose,insuline)

[t,y] = ode45(@(t,y) odefcn(t,y,insuline,time,parameters), [time(1), time(end)],[G0,x0]);

%plot plasmatic insuline
figure;
subplot(2,1,1);
time_for_interp=linspace(time(1),time(end),length(time));
I_inter=interp1(2*time,insuline,time_for_interp);
plot(time,I_inter,'-')
hold on
title('Insulina')
xlabel('Time[min]')
ylabel('Insulina [{\mu}U/ml]')

subplot(2,1,2);
plot(t,y(:,1),'-')
hold on
title('Glucosio')
xlabel('Time[min]')
ylabel('Glucosio [mg/dl]')

%Controllo la percentuale di valori al di fuori dei valori clinicamente
%normali
s=0;
for p=1:length(y)
    if y(p)<60 | y(p)>130
        s=s+1;
    end
end
p=s.*100/length(y);
disp(['Valori fuori dal fisiologico ',num2str(p),' %'])

%Valutazione dell'errore rispetto alle metriche cliniche in letteratura.
sol = ode45(@(t,y) odefcn(t,y,insuline,time,parameters), [time(1), time(end)],[G0,x0]);
evaluated_sol=deval(sol,time);
error=100*abs((glucose-evaluated_sol(1,:)')./glucose);
error_value=mean(error(2:end));
disp(['Errore medio: ',num2str(error_value),' %'])
end