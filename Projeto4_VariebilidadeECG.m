%Projeto Temático 4 - Variabilidade da Frequência Cardíaca
%Beatriz Silva, 56135
%Ines Gonçalves, 55067

%Carregamento dos dados
data = load('Dados67.mat');

%Transforma a tabela com os dados do ECG num vetor com esses dados
ECG_data = table2array(data.ECG_EDA11(:,3));
% ECG_data_2 = table2array(data.ECG_EDA21(:,3));
% ECG_data_3 = table2array(data.ECG_EDA31(:,3));
% ECG_data_5 = table2array(data.ECG_EDA51(:,3));
% ECG_data_6 = table2array(data.ECG_EDA61(:,3));

%Size dos vetores para cada indivíduo
size = length(ECG_data);
% size_2 = length(ECG_data_2);
% size_3 = length(ECG_data_3);
% size_5 = length(ECG_data_5);
% size_6 = length(ECG_data_6);

%% Indivíduo 1
%taxa de aquisição de 500 Hz => fs = 500/2;
fa = 500;
fs_max = fa / 2;
time_max = size/fa;
t = 1/fa: 1/fa: time_max;

% Visualização dos dados
figure(1);
subplot(1,2,1); plot(t, ECG_data); title("Indivíduo 1"); 
subplot(1,2,2); plot(t, ECG_data); title("Indivíduo 1 30 segundos"); xlim([0 30]);

% Encontrar picos QRS
% [PKS,LOCS] = findpeaks(ECG_data_1, fa, "MinPeakHeight", 42500);
% figure(??????);
% plot(t, ECG_data_1, 'blue'); hold on; 
% plot(LOCS, PKS, 'r*'); 

% Situacao 1 - - - - - - - - - - - - - - - - - - - - 
% 25min é o tempo
n_amostras_s1 = fa*(60*25);
t_s1 = 1/fa: 1/fa: (60*25);
highpass = High_Pass();
ECG_data_filtered_s1 = filter(highpass, ECG_data(1:n_amostras_s1, 1));

%tiramos a informação para os t = 25min
%t1_data_1 = ECG_data_1_filtered(1:n_amostras_1_1, 1);

[PKS_s1,LOCS_s1] = findpeaks(ECG_data_filtered_s1, fa, "MinPeakHeight", 9500, "MinPeakDistance", 0.4);

figure(2);
plot(t_s1, ECG_data_filtered_s1, 'blue'); hold on; 
plot(LOCS_s1, PKS_s1, 'r*');

%figure(4);
%spectrogram(t1_data_1, 4096, 512, fa, 'yaxis');

figure(3); 
[Pxx_s1, F_s1] = pwelch(ECG_data_filtered_s1, 4096, 2048, 4096, fa);
plot(F_s1, Pxx_s1);
xlim([0 4]);
xlabel('Frequência'); ylabel('Potência');

% Calculo da frequência cardíaca para os t = 25min
a1 = 1;
b1 = 1;
for i_s1 = 1:1:2248
    calc_fc1 = 60/(LOCS_s1(i_s1 + 1) - LOCS_s1(i_s1));
    if calc_fc1 > 70 && calc_fc1 < 110
        fc_s1(a1) = calc_fc1;
        t_fc_s1(b1) = LOCS_s1(i_s1);
        a1 = a1 + 1;
        b1 = b1 + 1;
    end
end

soma = sum(fc_s1(:, 10:(a1 -1))); %Aplicação de um filtro (contagem a partir dos 10 picos)
tamanho = a1-10;
BPM_mean_s1 = soma/ tamanho;
figure(4);
plot(t_fc_s1, fc_s1);

%% Situacao 2 - - - - - - - - - - - - - - - - - - - - 
% 25min-28min (3min) é o tempo
n_amostras_s2 = fa*(60*3);
t_s2 = (60*25) + 1/fa:1/fa:(60*28); 
data_s2 = ECG_data(n_amostras_s1 + 1: n_amostras_s1 + n_amostras_s2, 1);

figure(5);
plot(t_s2, data_s2);


figure(6); 
[Pxx_s2, F_s2] = pwelch(data_s2, 4096, 2048, 4096, fa);
plot(F_s2, Pxx_s2);
xlim([0 4]);
xlabel('Frequência'); ylabel('Potência');

% highpass = High_Pass();
% ECG_data_1_filtered = filter(highpass, ECG_data_1);

%tiramos a informação para os t = 3min (25-28min)
% t2_data_1 = ECG_data_1_filtered_s1(data_1_2, 1);
% figure(5);
% plot(t2(:,2:90001), data_1_2);

% Find Peaks
[PKS_s2,LOCS_s2] = findpeaks(data_s2, fa, "MinPeakHeight", 40000, "MinPeakDistance", 0.4);

for index_LOCS = 1:1:269
    LOCS_s2_uptodate(index_LOCS) = LOCS_s2(index_LOCS) + 1500; 
end
figure(7);
plot(t_s2, data_s2, 'blue'); hold on; 
plot(LOCS_s2_uptodate, PKS_s2, 'r*');

% Calculo da frequência cardíaca
a = 1;
b = 1;
for i_s2 = 1:1:268
    fc = 60/(LOCS_s2_uptodate(i_s2 + 1) - LOCS_s2_uptodate(i_s2));
    if fc > 70 
        fc_s2(a) = fc;
        t_fc_s2(b) = LOCS_s2_uptodate(i_s2);
        a = a + 1;
        b = b + 1;
    end
end

soma_s2 = sum(fc_s2); 
tamanho_s2 = a - 1;
BPM_mean_s2 = soma_s2/ tamanho_s2;
figure(8);
plot(t_fc_s2, fc_s2, 'red'); hold on;
plot(t_fc_s1, fc_s1, 'blue');

%% Situacao 3 - - - - - - - - - - - - - - - - - - - - 
% 28min-38min (10min) é o tempo
t_min_s3 = 28*60;
n_amostras_s3 = fa*(60*10);
t_s3 = (60*28) + 1/fa:1/fa:(60*38); 
data_s3 = ECG_data(n_amostras_s1 + n_amostras_s2 + 1: n_amostras_s3 + n_amostras_s2 + n_amostras_s1, 1);

figure(9);
plot(t_s3, data_s3);


figure(10); 
[Pxx_s3, F_s3] = pwelch(data_s3, 4096, 2048, 4096, fa);
plot(F_s3, Pxx_s3);
xlim([0 4]);
xlabel('Frequência'); ylabel('Potência');

% Find Peaks
[PKS_s3,LOCS_s3] = findpeaks(data_s3, fa, "MinPeakHeight", 40000, "MinPeakDistance", 0.4);
for index_LOCS_s3 = 1:1:800
    LOCS_s3(index_LOCS_s3) = LOCS_s3(index_LOCS_s3) + t_min_s3; 
end
figure(11);
plot(t_s3, data_s3, 'blue'); hold on; 
plot(LOCS_s3, PKS_s3, 'r*');

% Calculo da frequência cardíaca
a3 = 1;
b3 = 1;
for i_s3 = 1:1:799
    calc_fc_s3(i_s3) = 60/(LOCS_s3(i_s3 + 1) - LOCS_s3(i_s3));
%     if fc > 70 
%         fc_s3(a3) = calc_fc_s3;
%         t_fc_s3(b3) = LOCS_s3(i_s3);
%         a3 = a3 + 1;
%         b3 = b3 + 1;
%     end
end

soma_s3 = sum(fc_s3); 
tamanho_s3 = a3 - 1;
BPM_mean_s3 = soma_s3/ tamanho_s3;
figure(12);
plot(LOCS_s3(1:779,1), calc_fc_s3, 'green'); hold on;
plot(t_fc_s2, fc_s2, 'red'); hold on;
plot(t_fc_s1, fc_s1, 'blue');
