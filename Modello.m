% Questo capitolo, dedicato all'analisi del modello selezionato, si pone il
% semplice ma valido obiettivo di risolvere il sistema di equazioni differenziali
% con l'ausilio del Software MATLAB, per poi trarre alcune semplicissime
% conclusioni.
% Avendo creato il modello come un compartimentale in cui l'input è il
% flusso di insulina i(t) e tenendo a mente che:
% -L'insulina ha il compito di assorbire il glucosio plasmatico ed evitarne
% il suo accumulo a livello plasmatico.
% -I pazienti diabetici devono sempre assumere una dose intravenosa di
% insulina prima dei pasti.
% Valuteremo come reagisce la concentrazione di glucosio plasmatico a
% diversi tipi di input.
% 
% Per avere un riferimento clinicamente veritiero si utilizzeranno dei dati
% presenti in letteratura dallo studio di Pacini et Al.[](glucosio e
% insulina). 
% Dapprima si analizzeranno questi dati discreti .
% Successivamente, si costruirà la funzione descrittiva di i(t) attraverso
% un metodo di interpolazione lineare e si risolverà g(t) utilizzando la
% funzione built-in di matlab ode45(), che implementa il metodo di Runge-Kutta.
% 
% Oltre ai dati di letteratura, si caricheranno altri vettori che
% simuleranno degli input artificiali di insulina. Seguirà una risoluzione
% analoga a quella presentata per il data-set, da cui si potrà valutare la
% reazione di g(t) ai diversi input.
% 
% Per fornire alcune banalissime metriche, infine, si valuterà l'aderenza
% tra le varie funzioni di g(t) i dati relativi clinici di glucosio
% andando a valutarne l'errore percentuale e valutando la percentuale di
% valori, per ogni curva, entro il range fisiologicamente accettabile.
% 
% Per iniziare, si procede caricando i dati e le costanti nel file
% modello.m:
% Load data
clear all
close all
load("experimental_data.mat")

%Valori da letteratura da Pacini et Al.
time=tgi(:,1);
glucose=tgi(:,2);
insuline=tgi(:,3);

G0=279;
x0=0;
Gb=93; 
Ib=11;
Sg=2.6E-2;
k=0.025;
Si=5.0E-4;

parameters=[Sg,Gb,k,Ib,Si];

% Si forniscono ora i parametri e dati alla funzione ode45(), che in
% sinergia con la funzione ode(In appendice a questo lavoro) risolverà g(t).
% Segue il plot dell'insulina a livello plasmatico e del livello di glucosio
[t,y] = ode45(@(t,y) odefcn(t,y,insuline,time,parameters), [time(1), time(end)],[G0,x0]);

% plot dell'insulina
% figure;
% subplot(2,1,1);
% time_for_interp=linspace(time(1),time(end),length(time));
% I_inter=interp1(2*time,insuline,time_for_interp);
% plot(time,I_inter,'-*')
% hold on
% plot(time,insuline,'or')
% title('Plasmatic levels of insuline')
% legend({'Interpolation','Samples'})
% xlabel('Time[min]')
% ylabel('Insuline [{\mu}U/ml]')
% 
% %plot livello di glucosio
% subplot(2,1,2);
% plot(t,y(:,1),'-')
% hold on
% plot(time,glucose,'o')
% legend({'Solved', 'Measured'})
% title('Glucose')
% xlabel('Time[min]')
% ylabel('Glucose [mg/dl]')
% legend({'Model','Samples'})

%Controllo la percentuale di valori al di fuori dei valori clinicamente
%normali
% s=0;
% for p=(length(y)+1)./2:length(y)
%     if y(p)<60 | y(p)>130
%         s=s+1;
%     end
% end
% p=s.*100/length(y)/2;
% disp(['Valori fuori dal fisiologico ',num2str(p),' %'])
% 
% %Valutazione rispetto ai dati clinici
% sol = ode45(@(t,y) odefcn(t,y,insuline,time,parameters), [time(1), time(end)],[G0,x0]);
% evaluated_sol=deval(sol,time);
% error=100*abs((glucose-evaluated_sol(1,:)')./glucose);
% error_value=mean(error(2:end));
%disp(['Errore medio: ',num2str(error_value),' %'])

% %%%%%%%%%%%%%%%%%%%%%%%%
% %%VALUTAZIONI COMPLETAMENTE RANDOM (STANDARD)
%%valutazione - rect iniziale
insuline=100.*[1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);

%%valutazione - Funzione exp decrescente
%insuline=100.^(-100.*time);
%[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);

%Step
insuline=100.*[1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2];
[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);

%valutazione - Flusso costante di insulina basale
insuline=100.*ones(24,1);
[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);

%valutazione - impulso secco di insulina e stop
insuline=200.*[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);

%valutazione - nessun flusso di insulina
insuline=zeros(24,1);
[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);

%%VALUTAZIONI CERCANDO DI AVERE L'ERRORE MINORE POSSIBILE
%valutazione random 1
insuline=[90 86 83 78 73 65 61 50 45 39 33 28 21 18 11 7 8 7 5 4 3 1 1 1 ];
[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);

%valutazione random 2
insuline=[15 32 47 82 90 87 81 77 69 62 55 49 41 36 26 21 15 11 9 6 5 4 2 1 ];
[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);

%valutazione random 3
insuline=[15 32 47 82 90 87 81 77 69 60 51 49 41 36 26 21 15 11 9 6 5 4 2 1 ];
[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);

%valutazione random 4
insuline=[15 32 47 82 90 87 81 77 69 60 51 49 41 34 28 21 15 11 9 6 5 4 2 1 ];
[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);
% 
%%VALUTAZIONI IN CUI IL GLUCOSIO E' UNA FUNZIONE CRESCENTE 
%%(MI ASPETTO GLUCOSIO MOLTO BASSO)
%Valutazione-Funzione lineare crescente
insuline=7.*time;
[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);

%Valutazione-Funzione quadratica
insuline=2.*(time).^2;
[t,y]=Valutazione(parameters,time,G0,x0,glucose,insuline);