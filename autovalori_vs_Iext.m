clear
close all
clc

a = 0.7; b = 0.8; c = 12.5;

% Definisco il range in cui varia Iext
Iext = linspace(0, 3, 1000);

Lx = [0 3]; Ly = [-2 1];
ticksX = Lx(1):0.2:Lx(2);
ticksY = Ly(1):0.2:Ly(2); 

% Inizializzo i vettori che conterranno le parti reali di λ1 e λ2
Rel1 = zeros(1,1000);
Rel2 = zeros(1,1000);

% Riempio i vettori Rel1 e Rel2 calcolando gli autovalori della matrice
% Jacobiana per ogni valore di v0 ricavato da un'equazione di terzo grado in funzione di Iext
for j = 1:1000
    v = [-1/3 0 (1 - 1/b) (Iext(j) - a/b)];
    v0 = roots(v);
    for i = 1:3 %le radici del polinomio saranno una reale e due complesse coniugate
        if imag(v0(i)) == 0  % considero solo la radice reale, che corrisponde al caso fisico studiato
            J = [1 - v0(i)^2 -1; 1/c -b/c];
            autovalori = eig(J);
        end
    end
    Rel1(j) = real(autovalori(1));
    Rel2(j) = real(autovalori(2));
end

FS = 9;  % fontsize 

% Grafico
figure;
plot(Iext, Rel1, '-r', 'LineWidth', 2);  % Grafico della prima funzione in rosso
hold on;
plot(Iext, Rel2, '-b', 'LineWidth', 2); % Grafico della seconda funzione in blu
hold off;

set(gca,'fontsize',FS)
xlim(Lx)
ylim(Ly)
set(gca,'xtick',ticksX);
set(gca,'ytick',ticksY);

title('Parte reale degli autovalori al variare della corrente esterna');
xlabel('Corrente esterna Iext');
ylabel('Parte reale di λ1 e λ2');
legend('Re{λ1}', 'Re{λ2}');
grid on;