function s=pu_states(Ne,modelo)
%
% El estado del Tx se modela como una cadena de Markov con dos estados: 
% 0:inactivo 1:activo.
%
% Por defecto, la funcion matlab "hmmgenerate" genera una secuencia de
% estados empezando siempre pòr el primer estado.
% Para conseguir que el estado inicial sea aleatorio y de acuerdo a las
% probabilidades de estado, se considera una cadena de Markov aumentada con
% 3 estados. El segundo y tercer estados son los correspondientes a 
% canal inactivo y canal activo. El primer estado es transitorio y la 
% primera transicion se produce a los estados 2 o 3 con probabilidades de 
% transicion dadas por las probabilidades de estado, que se calculan a 
% partir de las probabilidades de estado de la cadena de Markov original
% (con dos estados)
% 
% Variables de entrada:
% - Ne: no. de realizaciones del estado de los Tx's.
% - modelo: modelo y parametros del proceso estocastico s(n): 
%       'ML'
%       'iid': p0 (probabilidad de s(n)=0)
%       'markov': t00,t11 (probabilidades de transicion)
%
% Variables de salida:
% - s: realizaciones del estado del Tx en binario (Ne x 1)

% Modelo del proceso estocastico s(n)
switch modelo{1}
   case 'iid'
      p0 = modelo{2};
      s=1*(rand(Ne,1) > p0);
   case 'markov'
      t00 = modelo{2};
      t11 = modelo{3};  
      % Probabilidades estados
      p0=(1-t11)/(2-t00-t11);
      p1=(1-t00)/(2-t00-t11);
      % Matriz de transiciones aumentada
      Ta=zeros(3,3);
      Ta(1,2)=p0;
      Ta(1,3)=p1;
      Ta(2,2)=t00;
      Ta(2,3)=1-t00;
      Ta(3,2)=1-t11;
      Ta(3,3)=t11;
      % Genera una secuencia (realizacion)
      [~,estados] = hmmgenerate(Ne,Ta,eye(3));
      s=estados'-2;
   case 'ML'
      s=randi([0 1],Ne,1);
   otherwise
      error('modelo erroneo');
end



