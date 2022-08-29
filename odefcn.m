
function dydt=odefcn(t,y,insuline,time,parameters)
% exctract parameters from the vector

Sg=parameters(1);
Gb=parameters(2);
k=parameters(3);
Ib=parameters(4);
Si=parameters(5);

I_inter=interp1(time,insuline,t);

%I(t)=I_interp
dydt(1)=Sg*(Gb-y(1))-y(2)*y(1);
dydt(2)=k*(Si*(I_inter-Ib)-y(2));
dydt=dydt'; % transpose for ODE45 sintax
end

