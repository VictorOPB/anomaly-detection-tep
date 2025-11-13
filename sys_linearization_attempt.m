%% 0) Abrir/Inicializar
clc; close all

Mode_1_Init
mdl = 'MultiLoop_mode1';
bdclose('all');              
load_system(mdl);            
open_system(mdl);

%% 1) Encontrar blocos únicos de xmv/xmeas
opts = {'LookUnderMasks','all','FollowLinks','on'};

candXMEAS = find_system(mdl, opts{:}, 'RegExp','on','Name','^xmeas$');
candXMV   = find_system(mdl, opts{:}, 'RegExp','on','Name','^xmv$');
candIDV   = find_system(mdl, opts{:}, 'RegExp', 'on','Name', '^idv$');

if isempty(candXMEAS) || isempty(candXMV) || isempty(candIDV)
    error('Não achei blocos chamados xmeas/xmv/idv. Confira o modelo aberto.');
end
blkXMEAS = candXMEAS{2};
blkXMV   = candXMV{2};
blkIDV   = candIDV{1};
blkPlant = [mdl '/TE Plant'];
blkDist = [mdl '/TE Plant/Disturbances'];

disp('Blocos detectados:');
disp(blkXMEAS);
disp(blkXMV);
disp(blkIDV);
disp(blkDist);
%% 2) Ponto de operação robusto (snapshot pela simulação)

opspec = operspec(mdl);
op = findop(mdl, opspec);
%% 3) Linearizar
linOpts = linearizeOptions('StoreAdvisor', 'on');
sysAll  = linearize(mdl, blkPlant, op, linOpts);  % deve ser 41 saídas x 12 entradas

% Verificar o tamanho
[ny, nu] = size(sysAll);
fprintf('Sistema linearizado: %d saídas x %d entradas\n', ny, nu);

if ny==0 || nu==0
    error('Linearização resultou vazia. Verifique "io" e o operating point.');
end

% As primeiras 12 entradas correspondem a xmv,
% as primeiras 41 saídas correspondem a xmeas.
Gx = sysAll(1:41, 1:12);
T = readtable('TEP_IO_map.csv');
Gx.InputName  = T.Name(strcmp(T.Type,'xmv'));
Gx.OutputName = T.Name(strcmp(T.Type,'xmeas'));
Gx.TimeUnit = 'hours';
Gx = d2d(Gx, 0.1);
% Extrair matrizes
[A,B,C,D] = ssdata(Gx);

%% 4) Validação rápida com degrau pequeno num subconjunto MIMO
% Entradas e saídas selecionadas: vide TEP_IO_map.csv

Gsub = Gx([7 12 36],[4 6]);
figure; step(Gsub); grid on;
title('LTI Validation (small step on selected inputs)');

%% 5) PRBS de baixa amplitude - excitar a dinâmica em várias frequências (NÃO USADO, PULAR)
% Gerar um sinal PRBS de baixa amplitude
amp = 0.005;
Tsim = 9600; % prbs por toda extensão dos dados de sim() abaixo.
Ts = 0.1000;
t_prbs = 0:Ts:Tsim;

U = zeros(numel(t_prbs),12);
for k=1:12
    % alterna a cada 0.1 / 0.02 = 5s
    U(:,k) = idinput(numel(t_prbs), 'prbs', [0 0.02], [-amp amp]);
end

% suavização do sinal
alpha = 0.9; % bem suave; 0 -> sem filtro
for k=1:12
    for n=2:numel(t_prbs)
        U(n,k) = alpha*U(n-1,k) + (1-alpha)*U(n,k);
    end
end
u_prbs = timeseries(U,t_prbs);
assignin("base","u_prbs",u_prbs); % Bloco "From Workspace"

%% 6) Comparação dados reais dos sinais extraídos com o sistema linearizado:
% a) Solver estável multirate
set_param('MultiLoop_mode1', 'Solver', 'ode23tb', 'MaxStep', '0.1000')
simOut = sim(mdl,'ReturnWorkspaceOutputs','on');

% b) Extrair sinais
u  = simOut.get('xmv');    % entradas manipuladas [N×12]
y  = simOut.get('simout');  % saídas medidas [N×41]
Ts = Gx.Ts; 
N0 = round(200/Ts); 
u0 = mean(u(1:N0, :), 1);    % [1x12]
y0 = mean(y(1:N0, :), 1);    % [1x41]

du = u - u0;                 % desvios de entrada
t  = (0:size(u,1)-1)'*Ts;

y_lin_dev = lsim(Gx, du, t); % saída em desvios
y_lin     = y0 + y_lin_dev;  % reconstrói valor absoluto para comparar

% vamos redimensionar t para a dimensão de amostragem interna da planta:
% 0.0005s
t = t * (0.0005 / Ts);

% (c) Plotar comparação
isY = strcmp(T.Type, 'xmeas');
labelsY = T.Name(isY);
unitsY = T.Unit(isY);
% idxY = [7 8 9];    % Exemplo: pressão, nível e temperatura do reator
idxY = [7 12 36];
figure; 
tiledlayout(numel(idxY),1)
for k = 1:numel(idxY)
    nexttile
    plot(t, y(:,idxY(k)),'-', t, y_lin(:,idxY(k)),'--','LineWidth',1)
    grid on
    ylabel(sprintf('%s (%s)', labelsY{idxY(k)}, unitsY{idxY(k)}), 'Interpreter','none');
end
xlabel('Time (hours)');
legend('TEP Plant','Linear Model');
sgtitle('Comparison: Simulated TEP vs Linearized Model');

%% 7) Qualidade do fit:

% Métricas globais
err    = y - y_lin;                     % [N x 41]
rmsAbs = sqrt(mean(err.^2, 1));         % RMSE absoluto por saída

% NRMSE por faixa (robusto quando std(y) é pequeno)
yrange = max(y,[],1) - min(y,[],1) + 1e-12;
NRMSE_range = 100*(rmsAbs ./ yrange);

namesY = Gx.OutputName(:);
Tmetrics = table((1:numel(namesY))', string(namesY), rmsAbs', NRMSE_range', ...
    'VariableNames', {'y_idx','y_name','RMSE_abs','RMSE_rel_percent'});
disp(Tmetrics(1:10,:))  % espiada nas 10 primeiras

writetable(Tmetrics, 'TEP_LTI_metrics_all_outputs_.csv');

% RMSE móvel (janela 300s)
Wsec = 300;                % janela de 300s = 1.5h na amostragem interna
W    = max(1, round(Wsec / Ts));     % em amostras

idxPlot = [7 12 36];        % pressão, nível do separador de produto, taxa de produto
err2     = err.^2;         % erro ao quadrado
rmsMov   = sqrt(movmean(err2, W, 1)); % [N x 41], média móvel centrada
y_mu    = movmean(y,W,1);
yhat_mu  = movmean(y_lin,W,1);

SS_res = movsum((y-y_lin).^2,W,1);
SS_tot = movsum((y-y_mu).^2,W,1)+1e-12;
R2_mov = 1- SS_res./SS_tot;

figure; tiledlayout(numel(idxPlot),1)
for k = 1:numel(idxPlot)
    j = idxPlot(k);
    nexttile
    plot(t, rmsMov(:,j),'LineWidth',1);
    grid on
    ylabel(sprintf('%s', Gx.OutputName{j}),'Interpreter','none')
end
xlabel('Time (hours)')
sgtitle(sprintf('Rolling RMSE (window = %sh)', '1,5'))

figure; tiledlayout(numel(idxPlot),1)
for k=1:numel(idxPlot)
  nexttile
  plot(t, R2_mov(:,idxPlot(k))); grid on
  ylim([-1 1]); ylabel(sprintf('%s', Gx.OutputName{idxPlot(k)}),'Interpreter','none')
end
xlabel('Time (hours)'); sgtitle(sprintf('Rolling R^2 (window = %sh)', '1,5'));

%% 8) extraindo dados dos controladores de interesse (DEPRECADO; PULAR)
% controlador de velocidade, fórmula:
% du(k) = Kc[e(k) - e(k-1)] + K_{i,inc}e(k)
% u(k) = u(k-1) + du(k)
% Kc = ganho do /Gain (prop. incremental)
% K_{i,inc} = ganho do /Gain1 (integral incremental)
% K_{i,cont} = K_{i,inc} / Ts
% Ti = (KcTs)/(T_{i,inc})

function Tpi = read_te_pi_params(mdl)
    % Lê parâmetros dos controladores PI incrementais do TEP.
    % Retorna tabela com Kc, Ki_inc, Ts, Ki_cont e Ti.
    
    opts = {'LookUnderMasks','all','FollowLinks','on'};
    % Candidatos a "raiz" de controlador
    piRoots = [ ...
        find_system([mdl '/TE Plant'], opts{:}, 'Name','Vel PI'); ...
        find_system([mdl '/TE Plant'], opts{:}, 'Name','yA control'); ...
        find_system([mdl '/TE Plant'], opts{:}, 'Name','yAC control') ...
    ];
    
    % Compila o diagrama para obter SampleTime "Compiled"
    set_param(mdl,'SimulationCommand','update');
    
    mw = get_param(mdl,'ModelWorkspace');
    
    rows = {};
    for i = 1:numel(piRoots)
        root = piRoots{i};
    
        % pegue os blocos do tipo Gain diretamente (evita colisão por nome)
        gP = find_system(root, 'SearchDepth',1, 'BlockType','Gain');    % pode retornar 1 ou mais; o "P" costuma ser o sem sufixo
        gI = find_system(root, 'SearchDepth',1, 'BlockType','Gain');    % idem — vamos filtrar por nome
    
        % escolha por nome
        gP = gP(contains(get_param(gP,'Name'),'Gain') & ~contains(get_param(gP,'Name'),'Gain1'));
        gI = gI(contains(get_param(gI,'Name'),'Gain1'));
    
        % paths esperados do atraso e do ZOH (podem não existir)
        udel = find_system(root, 'SearchDepth',1, 'BlockType','UnitDelay');
        zoh  = find_system(root, 'SearchDepth',1, 'BlockType','ZeroOrderHold');
    
        % ---- valores dos ganhos (podem ser expressões)
        Kc     = eval_param_num(gP);
        Ki_inc = eval_param_num(gI);
        disp('Ki_inc:'); disp(Ki_inc);
        disp('Kc:'); disp(Kc);
        % ---- Ts compilado
        Ts = NaN;
        if ~isempty(udel)
            try
                st = get_param(udel{1}, 'CompiledSampleTime');  % [offset Ts]
                Ts = st(2);
            catch
            end
        end
        if (isnan(Ts) || Ts<=0) && ~isempty(zoh)
            % tenta ZOH (parâmetro explícito)
            try
                Ts = str2double(get_param(zoh{1}, 'SampleTime'));
            catch
            end
            % se herdado, tenta compilado
            if (isnan(Ts) || Ts<=0)
                try
                    st = get_param(zoh{1}, 'CompiledSampleTime');
                    Ts = st(2);
                catch
                end
            end
        end
    
        % conversões contínuas
        if ~isnan(Ts) && Ts>0 && ~isnan(Ki_inc) && Ki_inc~=0
            Ki_cont = Ki_inc / Ts;
            Ti      = (Kc * Ts) / Ki_inc;
        else
            Ki_cont = NaN;
            Ti      = Inf;
        end
    
        % em qual malha estamos? (nome do subsistema pai)
        parent = get_param(get_param(root,'Parent'),'Name');
    
        rows(end+1,:) = {root, parent, Kc, Ki_inc, Ts, Ki_cont, Ti}; %#ok<AGROW>
    end

    Tpi = cell2table(rows, 'VariableNames', ...
        {'ControllerPath','LoopName','Kc','Ki_inc','Ts','Ki_cont','Ti_seconds'});
    
    % ordena para leitura
    Tpi = sortrows(Tpi, 'ControllerPath');
    
    % salva
    writetable(Tpi, 'TEP_PI_parameters.csv');

    % ---- helper interno ----
    function val = eval_param_num(blks)
        % devolve número resolvendo expressão do parâmetro 'Gain'
        val = NaN;
        if isempty(blks), return; end
        expr = get_param(blks{1}, 'Gain');     % string (p.ex. 'Kc_SL' ou '0.2/60')
        num  = str2double(expr);
        if ~isnan(num)
            val = num; return;
        end
        % tenta no Model Workspace
        try
            val = mw.evalin(expr);
            if isnumeric(val) && isscalar(val), val = double(val); return; end
            catch 
        end
        % tenta no Base Workspace
        try
            val = evalin('base', expr);
            if isnumeric(val) && isscalar(val), val = double(val); return; end
            catch
        end
        % se ainda não deu, deixa NaN (o usuário pode investigar a expressão)
    end
end

Tpi = read_te_pi_params(mdl);
% controladores ligados às saídas medidas em Gsub:
T_strip  = Tpi(contains(Tpi.LoopName,'Stripper level'), :);   % afeta xmv(10)
T_recyc  = Tpi(contains(Tpi.LoopName,'Recycle') | ...
               contains(Tpi.LoopName,'Compressor') | ...
               contains(Tpi.LoopName,'Production rate'), :); % próximo de xmv(5) no Gsub
T_purge  = Tpi(contains(Tpi.LoopName,'Purge'), :);
disp("T_strip:"); disp(T_strip);
disp("T_purge:"); disp(T_purge);
disp("T_recyc:"); disp(T_recyc);

%% 9) Salvar processo benigno

Ttbl = array2table([t(1:round(.75*end)) , u(1:round(.75*end),:) , y(1:round(.75*end),:) , y_lin(1:round(.75*end),:), err(1:round(.75*end),:)], ...
                    'VariableNames', ['t_s', strcat('xmv_',string(1:12)), strcat('xmeas_',string(1:41)), ...
                    strcat('yhat_', string(1:41)), strcat('res_', string(1:41))]);
writetable(Ttbl, '/Users/insider/TCC/tep-attack-simulator/temexd_mod/benign_36h.csv');

%% 10) Salvar Gx para carregar no .m de ataque
save('/Users/insider/TCC/tep-attack-simulator/temexd_mod/Gx_lti.mat', 'Gx')