%% === Configuração ===
mdl    = 'MultiLoop_mode1';
blkDos = [mdl '/TE Plant/dos_gate_bernoulli'];  % caminho do seu bloco
Ts     = 0.1;             % seu passo (s)
Tsim_h = 12;              % duração (horas) da simulação sob DoS
rng(42,'twister');        % semente global reprodutível

% parâmetros do Bernoulli:
p_del       = 0.30;       % prob. de DROP (0->entrega, 1->drop)
per_channel = 'on';       % 'on' = cada canal tem sorteio independente
hold_mode   = 1.0;        % 1=hold-last, 2=zeros, 3=NaNs (mesmo que na máscara)

% aplica na MÁSCARA do bloco
set_param(blkDos, 'P_DEL', num2str(p_del), ...
                  'PER_CHANNEL', per_channel, ...
                  'HOLD_MODE', hold_mode);

%% === Carregamento Gx ===

% carregando Gx salvo
load('/Users/insider/TCC/tep-attack-simulator/temexd_mod/Gx_lti.mat','Gx');
if Gx.Ts==0, Gx = d2d(Gx,Ts); end

% (opcional) nomes coerentes (se necessário)
T = readtable('TEP_IO_map.csv');
Gx.InputName  = T.Name(strcmp(T.Type,'xmv'));
Gx.OutputName = T.Name(strcmp(T.Type,'xmeas'));
%% === Run ===
open_system(mdl);                         %#ok<*UNRCH>
set_param(mdl,'StopTime', num2str(Tsim_h));   % StopTime em horas

% (Opcional) deixe consistente o solver do seu pipeline
set_param(mdl,'Solver','ode23tb','MaxStep', num2str(Ts));

simOut_dos = sim(mdl, 'ReturnWorkspaceOutputs','on');

%% === Coleta ===
t_s   = simOut_dos.get('tout');      % tempo em segundos na prática / horas no POV do modelo
xmv   = simOut_dos.get('xmv');       % [N x 12] comandos dos controladores não contaminados
xmeas = simOut_dos.get('simout');    % [N x 41] medições da planta
uhat  = simOut_dos.get('u_hat');     % [N x 12] atuador após DoS
gamma = simOut_dos.get('gamma_out'); % máscara 0/1 por canal - atualmente dims [N x 1 x 12] pois fiz um vetor coluna, e não vetor linha, tal qual as outras saídas.

if ndims(gamma) == 3, gamma = squeeze(permute(gamma, [3 1 2])); end % -> [N x 12]
disp(size(gamma));

%% === Obtenção y_lin ===

% Alinhamento DC (mesma lógica do benigno)
N0 = max(1, round(200/Ts));
u0 = mean(xmv(1:N0,:),1);
y0 = mean(xmeas(1:N0,:),1);
du = xmv - u0;
tt = (0:size(xmv,1)-1)'*Ts;

% Previsão LTI
ydev = lsim(Gx, du,tt);
yhat_lin = y0 + ydev;
res = xmeas - yhat_lin;

%% === Monta tabela e salva ===
VN = ['t_s', ...
      strcat('xmv_',string(1:size(xmv,2))), ... % entradas não contaminadas
      strcat('xmeas_',string(1:size(xmeas,2))), ... % saída da planta contaminada
      strcat('yhat_', string(1:size(yhat_lin,2))), ... % o que o modelo linear preveu
      strcat('res_', string(1:size(res,2))), ... % resíduos planta - modelo linear, agora contaminados
      strcat('uhat_',string(1:size(uhat,2))), ... % entradas contaminadas
      strcat('gamma_',string(1:size(gamma,2)))]; % se houve ou não perda de pacote

T = array2table([t_s, xmv, xmeas, yhat_lin, res, uhat, gamma], 'VariableNames', VN);

% anexa metadados do ensaio (útil p/ reprodutibilidade)
T.meta_p_del(1)       = p_del;
T.meta_per_channel(1) = strcmpi(per_channel,'on');
T.meta_hold_mode(1)   = str2double(hold_mode);
T.meta_Ts_s(1)        = Ts;
T.meta_Tsim_h(1)      = Tsim_h;

outdir = '/Users/insider/TCC/tep-attack-simulator/temexd_mod';
fname  = sprintf('dos_bern_p%.2f_per%s_hold%s_%dh.csv', ...
                  p_del, per_channel, hold_mode, Tsim_h);
writetable(T, fullfile(outdir, fname));
disp("Salvo: " + fullfile(outdir, fname));
%% === Rodar vários cenários (grade de p_bernoulli) ===
p_list     = [0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90];
sd=1;
per_opt    = {'on','off'};  % por-canal ou pacote único
hold_list  = {1.0};     % 1=hold-last, 2=zeros mas sim não roda com 2!

for per = 1:numel(per_opt)
  for h = 1:numel(hold_list)
    for p = 1:numel(p_list)
        rng(sd,'twister');
        fname = sprintf('dosB_seed%02d_p%.2f_per%s_hold%s_%dh.csv', ...
                         sd, p_list(p), per_opt{per}, string(hold_list{h}), Tsim_h);
        disp(fname);
        set_param(blkDos, 'P_DEL', num2str(p_list(p)), ...
                          'PER_CHANNEL', per_opt{per}, ...
                          'HOLD_MODE', hold_list{h});
    
        simOut_dos = sim(mdl,'ReturnWorkspaceOutputs','on');
    
        t_s   = simOut_dos.get('tout');
        xmv   = simOut_dos.get('xmv');
        xmeas = simOut_dos.get('simout');
        uhat  = simOut_dos.get('u_hat');
        gamma = simOut_dos.get('gamma_out');

        if ndims(gamma) == 3, gamma = squeeze(permute(gamma, [3 1 2])); end % -> [N x 12]

        u0 = mean(xmv(1:N0,:),1);
        y0 = mean(xmeas(1:N0,:),1);
        du = xmv - u0;
        tt = (0:size(xmv,1)-1)'*Ts;
        
        % Previsão LTI
        ydev = lsim(Gx, du,tt);
        yhat_lin = y0 + ydev;
        res = xmeas - yhat_lin;
    
        VN = ['t_s', ...
              strcat('xmv_',string(1:size(xmv,2))), ... % entradas não contaminadas
              strcat('xmeas_',string(1:size(xmeas,2))), ... % saída da planta contaminada
              strcat('yhat_', string(1:size(yhat_lin,2))), ... % o que o modelo linear preveu
              strcat('res_', string(1:size(res,2))), ... % resíduos planta - modelo linear, agora contaminados
              strcat('uhat_',string(1:size(uhat,2))), ... % entradas contaminadas
              strcat('gamma_',string(1:size(gamma,2)))]; % se houve ou não perda de pacote
        T = array2table([t_s, xmv, xmeas, yhat_lin, res, uhat, gamma], 'VariableNames', VN);
        T.meta_seed(1)        = sd;
        T.meta_p_del(1)       = p_list(p);
        T.meta_per_channel(1) = strcmpi(per_opt{per},'on');
        T.meta_hold_mode(1)   = hold_list{h};
        writetable(T, fullfile(outdir,fname));
    end
  end
end
