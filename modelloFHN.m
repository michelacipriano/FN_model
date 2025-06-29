clear 
close all
clc

% INPUT ========================================================
% Corrente esterna  [Iext = 1.0]
Iext = 1.0;
% Vettore delle variabili y = [v, w]
% Condizioni iniziali: y0 = [ v0 = -2.8, w0 = -1.8] 
y0 = [-2.8; -1.8];
% Tempo di simulazione [tMin = 0 tMax = 200]
tSpan = [0 200];
% Dimensione del plot per il piano delle fasi: asse X [-3 3] / asse Y [-2 3]
Lx = [-3 3]; Ly = [-2 3];
% Spaziatura assi X e Y
ticksX = Lx(1):1:Lx(2);
ticksY = Ly(1):1:Ly(2); 
% Numero di vettori per il campo delle velocita' nX [20] 
nX = 20;
% Parametri di controllo a, b, c = tau
a = 0.7; b = 0.8; c = 12.5;
% Vettore dei parametri di controllo K
K = [Iext; a; b; c];

% SOLUZIONE tramite ode45 ===================================================
[t,y] = ode45(@(t,y) FNode(t,y,K), tSpan,y0);

% Campo vettoriale
x1 = linspace(Lx(1),Lx(2),nX);
x2 = linspace(Ly(1),Ly(2),nX);
[xx, yy] = meshgrid(x1,x2);
f = (xx - xx.^3/3 - yy + Iext);   
g = (1/c).*(xx + a - b.*yy); 
fs = f./sqrt(f.^2 + g.^2);    % normalizzazione per avere vettori unitari
gs = g./sqrt(f.^2  +g.^2);
  
% Calcolo del punto critico
syms p
Sp = vpasolve(p-p^3/3-(p+a)/b + Iext == 0,p,[-3 3]);  % risolvo per una variabile simbolica p
Sq = (Sp+a)/b;
Sp = double(Sp); Sq = double(Sq);
disp('Critical point');
fprintf('   v_C =  %2.2f\n', Sp);
disp('   ')
fprintf('   w_C =  %2.2f\n', Sq);


  
% GRAFICI =======================================================  

FS = 14;  % fontsize

% Plot piano delle fasi: v vs w   ------------------------------------------
figure(1)
   pos = [0.35 0.05 0.29 0.39];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   hold on
   box on
  
% Campo delle velocita'
   hq = quiver(xx,yy,fs,gs);
   set(hq,'color',[0.2 0.2 0.2],'AutoScaleFactor',0.6);
   set(gca,'fontsize',FS)
   xlim(Lx)
   ylim(Ly)
   set(gca,'xtick',ticksX);
   set(gca,'ytick',ticksY);
   grid on   
   xlabel('potenziale di membrana v'); ylabel('variabile di recupero  w');

% v nullcline
     v = linspace(Lx(1),Lx(2),200);
     w = (v - v.^3/3 + Iext);
        xP = v; yP = w;
          plot(xP,yP,'r','linewidth',1.5)
% w nullcline 
    w = (v + a)/b;
        xP = v; yP = w;
          plot(xP,yP,'m','linewidth',1.5)

% Phase portrait
    xP = y(:,1); yP = y(:,2);
      plot(xP,yP,'b','linewidth',2)
    xP = y(1,1); yP = y(1,2);    % condizioni iniziali: inizio della traiettoria
      Hplot = plot(xP,yP,'o');
      set(Hplot,'markersize',8,'markerfacecolor',[0 1 0],'markeredgecolor',[0 1 0])
    xP = y(end,1); yP = y(end,2);   % fine della traiettoria
      Hplot = plot(xP,yP,'o');
      set(Hplot,'markersize',8,'markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0])
   
   tm1 = 'I_{ext} = ';
   tm2 = num2str(Iext,'%3.3f');
   tm3 = '  v_C = ';
   tm4 = num2str(Sp,'%3.2f');
   tm5 = '  w_C = ';
   tm6 = num2str(Sq,'%3.2f');
   tm = [tm1 tm2 tm3 tm4 tm5 tm6];   
   hT = title(tm,'FontName','Courier');  

   
% Evoluzione temporale di v e w   
figure(2)
  pos = [0.05 0.05 0.29 0.29];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  xP = t; yP = y(:,1);  % v
    plot(xP,yP,'b','linewidth',2)
  hold on
    yP = y(:,2);        % w
  plot(xP,yP,'r','linewidth',2)

  legend('v','w','location','south','orientation','horizontal')
  xlabel('t')
  ylabel('v & w')
  title(tm,'fontName','Courier')
  grid on
  set(gca,'fontsize',FS)
  box on 

  disp('  ')

% FUNZIONE ===========================================================

function dydt = FNode(t,y,K)
   a = K(2); b = K(3); c = K(4); Iext = K(1);
   dydt = [(y(1) - y(1)^3/3 - y(2) + Iext); (1/c)*(y(1) + a - b*y(2))];

end

