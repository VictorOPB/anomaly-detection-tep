function [u_hat, gamma_out] = dos_gate_bernoulli(u, P_DEL, PER_CHANNEL, HOLD_MODE)
%#codegen
% Parâmetros vindos da máscara (escalares):
%   P_DEL       : double 0..1
%   PER_CHANNEL : logical (0/1)
%   HOLD_MODE   : uint8/double {1=Hold last, 2=Zeros, 3=NaNs}

% Declaração dos parâmetros (o Simulink liga pela máscara)
coder.extrinsic('coder'); %#ok<EXTRIN>  % só pra suprimir warnings do editor
% (Os símbolos P_DEL, PER_CHANNEL, HOLD_MODE são "Parameter" no Edit Data)

% --- Inicializações para ajudar o analisador de tamanhos ---
u_hat     = u;
gamma_out = zeros(size(u));

% Persistente para memória do hold
persistent init u_prev
if isempty(init)
    u_prev = u;    % mesmo tamanho de u
    init   = true;
end

% --- Gera máscara lógica gamma: 1=ENTREGA, 0=DROP ---
if PER_CHANNEL
    gamma = rand(size(u)) > P_DEL;   % logical, mesmo tamanho de u
else
    g     = rand > P_DEL;
    gamma = repmat(g, size(u));      % logical
end

% --- Aplica o modo de retenção ---
mode = uint8(HOLD_MODE);
if mode==1
    % Hold last (sample-and-hold)
    u_hat = u;
    m = ~gamma;               % onde caiu o pacote
    u_hat(m) = u_prev(m);     % usa último valor válido
    % atualiza memória apenas quando entregou
    u_prev(gamma) = u(gamma);

elseif mode==2
    % Zeros
    u_hat = u .* double(gamma);

else
    % NaNs
    u_hat = u;
    m = ~gamma;
    u_hat(m) = NaN;
end

gamma_out = double(gamma);   % útil para log