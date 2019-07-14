% indici degli elementi di a entro b
function [v]=posizione(a,b)

% a e b sono vettori
% il vettore v ha la stessa lunghezza di a
% v(j) è l'indice dell'elemento a(j) dentro il vettore b
% ciascun elemento di a deve essere contenuto esattamente una volta in b, 
% altrimenti si ottiene un errore

% dimensiona v, fittiziamente pari ad a
v=a;
% cerca le posizioni di ciascun elemento di a in b
for i=1:length(a)
    v(i)=find(b==a(i));
end

return
