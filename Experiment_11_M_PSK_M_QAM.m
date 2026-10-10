function Experiment_11_M_PSK_M_QAM
% Experiment 11 - M-PSK and M-QAM (QPSK, 8-PSK, 16-QAM, 64-QAM)
% Gray vs binary mapping, energy normalization, min-distance detection,
% BER/SER and spectral efficiency. Works in MATLAB R2015a (base MATLAB only).
% Noise is specified as Eb/N0 (dB).
close all; clc; rng(1);

Mlist   = [4 8 16 64];
names   = {'QPSK','8-PSK','16-QAM','64-QAM'};
isPSK   = [true true false false];
Nsym    = 2e5;
EbN0dB  = 0:1:20;
testSNR = [0 5 10 20];
Q       = @(x) 0.5*erfc(x/sqrt(2));
col     = lines(4);
nM = numel(Mlist); nS = numel(EbN0dB);

%% Expected effects (written BEFORE simulation)
fprintf('================ EXPECTED EFFECTS (theory) ================\n');
fprintf('1. Higher M -> more bits/symbol (spectral efficiency log2(M)) but smaller\n');
fprintf('   minimum distance at fixed Es=1 -> higher BER/SER at the same Eb/N0.\n');
fprintf('2. Increasing Eb/N0 shrinks noise clouds -> BER/SER fall rapidly.\n');
fprintf('3. At 0 dB clouds overlap heavily (esp. 64-QAM); at 20 dB nearly separated.\n');
fprintf('4. Gray mapping: neighbours differ by 1 bit -> BER ~ SER/log2(M).\n');
fprintf('   Binary mapping: some neighbours differ by >1 bit -> BER higher than Gray.\n');
fprintf('5. Normalization fixes Es=1: it changes only the scale, not SER/BER.\n\n');

%% Mandatory validation: average energy before/after normalization
fprintf('=========== AVERAGE CONSTELLATION ENERGY ===========\n');
fprintf('%-8s %-18s %-18s\n','Mod','Es before norm.','Es after norm.');
Cg = cell(1,nM); Cb = cell(1,nM);
for m = 1:nM
    [Craw,Cn] = make_const(Mlist(m), isPSK(m), 'gray');
    Cg{m} = Cn;
    [~,Cn2] = make_const(Mlist(m), isPSK(m), 'binary');
    Cb{m} = Cn2;
    fprintf('%-8s %-18.6f %-18.6f\n', names{m}, mean(abs(Craw).^2), mean(abs(Cn).^2));
end
fprintf('(Raw: PSK radius 1; QAM levels +-1,+-3,... -> Es = 2(M-1)/3.)\n\n');

%% Diagnostic tests
fprintf('================ DIAGNOSTIC TESTS ================\n');
fprintf('%-8s %-14s %-14s %-14s\n','Mod','MaxHam Gray','MaxHam Binary','Noiseless SER');
for m = 1:nM
    hG = max_neighbor_hamming(Cg{m}, Mlist(m));
    hB = max_neighbor_hamming(Cb{m}, Mlist(m));
    lab = (0:Mlist(m)-1).';
    s0  = mean(mindist_detect(Cg{m}(lab+1), Cg{m}) ~= lab);
    fprintf('%-8s %-14d %-14d %-14g\n', names{m}, hG, hB, s0);
end
fprintf('(Gray => 1 ; Noiseless SER must be 0 => mapper/demapper consistent.)\n\n');
fprintf('Starting simulation...\n');

%% BER/SER sweep
berG = zeros(nM,nS); serG = berG; berB = berG; serB = berG;
berT = berG; serT = berG;
for m = 1:nM
    M = Mlist(m); k = log2(M);
    for s = 1:nS
        [berG(m,s),serG(m,s)] = simulate(Cg{m}, M, EbN0dB(s), Nsym);
        [berB(m,s),serB(m,s)] = simulate(Cb{m}, M, EbN0dB(s), Nsym);
        g = 10^(EbN0dB(s)/10);
        if isPSK(m)
            serT(m,s) = min(1, 2*Q(sqrt(2*k*g)*sin(pi/M)));
        else
            L = sqrt(M); x = sqrt(3*k*g/(M-1));
            pe = 2*(1-1/L)*Q(x);  serT(m,s) = 1-(1-pe)^2;
        end
        berT(m,s) = serT(m,s)/k;
    end
end

%% Results table
fprintf('=========== BER / SER at test points (Eb/N0) ===========\n');
for t = testSNR
    s = find(EbN0dB==t);
    fprintf('\nEb/N0 = %2d dB\n', t);
    fprintf('%-8s %-4s %-11s %-11s %-11s %-11s %-11s\n','Mod','SE','BER Gray','BER theory','BER Binary','SER Gray','SER theory');
    for m = 1:nM
        fprintf('%-8s %-4d %-11.3e %-11.3e %-11.3e %-11.3e %-11.3e\n', names{m}, log2(Mlist(m)), ...
            berG(m,s), berT(m,s), berB(m,s), serG(m,s), serT(m,s));
    end
end
fprintf('\nNote: BER = 0 only means below the resolution 1/(Nsym*log2M) (finite-sample effect).\n');
fprintf('Theory uses the nearest-neighbour approximation, so a small mismatch at 0 dB is expected.\n\n');

%% FIGURE 1: Ideal / received constellations
figure(1); set(gcf,'Name','Ideal and received constellations','Color','w','Units','normalized','Position',[0.02 0.05 0.95 0.85]);
Nplot = 1500;
for m = 1:nM
    for s = 1:numel(testSNR)
        [rx,~] = channel(Cg{m}, Mlist(m), testSNR(s), Nplot);
        subplot(nM,numel(testSNR),(m-1)*numel(testSNR)+s);
        plot(real(rx),imag(rx),'.','MarkerSize',3,'Color',[0.25 0.55 0.9]); hold on;
        plot(real(Cg{m}),imag(Cg{m}),'r+','MarkerSize',5,'LineWidth',1);
        axis equal; grid on; xlim([-2 2]); ylim([-2 2]);
        title(sprintf('%s, Eb/N0 = %d dB',names{m},testSNR(s)),'FontSize',8);
        if s==1, ylabel('Q'); end
        if m==nM, xlabel('I'); end
    end
end
annotation('textbox',[0 0.95 1 0.05],'String','Ideal (red +) and received (blue) constellations - Gray mapping, Es = 1', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');

%% FIGURE 2: Decision regions
figure(2); set(gcf,'Name','Decision regions','Color','w','Units','normalized','Position',[0.1 0.1 0.8 0.8]);
lim = 1.5; gg = linspace(-lim,lim,300); [X,Y] = meshgrid(gg,gg);
for m = 1:nM
    M = Mlist(m); k = log2(M);
    lab = mindist_detect(X(:)+1j*Y(:), Cg{m});
    Lmap = reshape(lab,size(X));
    subplot(2,2,m);
    image(gg,gg,ind2rgb(Lmap+1,hsv(M))); set(gca,'YDir','normal');
    hold on; plot(real(Cg{m}),imag(Cg{m}),'k.','MarkerSize',12);
    if M <= 16
        for q = 0:M-1
            text(real(Cg{m}(q+1)),imag(Cg{m}(q+1))+0.09,dec2bin(q,k), ...
                'HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
        end
    end
    axis equal tight; xlabel('I'); ylabel('Q');
    title([names{m} ' decision regions (Gray labels)']);
end

%% FIGURE 3: BER / SER
figure(3); set(gcf,'Name','BER and SER','Color','w','Units','normalized','Position',[0.02 0.15 0.95 0.6]);
subplot(1,3,1);
for m = 1:nM
    b = berG(m,:); b(b==0) = NaN;
    semilogy(EbN0dB,b,'o-','Color',col(m,:),'LineWidth',1.3); hold on;
end
for m = 1:nM, semilogy(EbN0dB,berT(m,:),'--','Color',col(m,:)); end
grid on; ylim([1e-6 1]); xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('BER (Gray): sim (o) vs theory (--)');
legend(names{1},names{2},names{3},names{4},'Location','southwest');

subplot(1,3,2);
for m = 1:nM
    v = serG(m,:); v(v==0) = NaN;
    semilogy(EbN0dB,v,'o-','Color',col(m,:),'LineWidth',1.3); hold on;
end
for m = 1:nM, semilogy(EbN0dB,serT(m,:),'--','Color',col(m,:)); end
grid on; ylim([1e-6 1]); xlabel('E_b/N_0 (dB)'); ylabel('SER');
title('SER (Gray): sim (o) vs theory (--)');
legend(names{1},names{2},names{3},names{4},'Location','southwest');

subplot(1,3,3);
for m = 1:nM
    a = berG(m,:); a(a==0) = NaN; b = berB(m,:); b(b==0) = NaN;
    semilogy(EbN0dB,a,'o-','Color',col(m,:),'LineWidth',1.3); hold on;
    semilogy(EbN0dB,b,'s--','Color',col(m,:),'LineWidth',1.0);
end
grid on; ylim([1e-6 1]); xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('Gray (solid, o) vs Binary (dashed, s)');
h1 = plot(NaN,NaN,'k-o'); h2 = plot(NaN,NaN,'k--s');
legend([h1 h2],{'Gray','Binary'},'Location','southwest');

%% FIGURE 4: BER vs spectral efficiency
figure(4); set(gcf,'Name','BER vs spectral efficiency','Color','w','Units','normalized','Position',[0.2 0.15 0.6 0.65]);
eta = log2(Mlist);
mk = {'o-','s-','d-','^-'};
for t = 1:numel(testSNR)
    s = find(EbN0dB==testSNR(t));
    v = berG(:,s).'; v(v==0) = NaN;
    semilogy(eta,v,mk{t},'LineWidth',1.4,'MarkerSize',7); hold on;
end
grid on;
set(gca,'XTick',eta,'XTickLabel',{'2 (QPSK)','3 (8-PSK)','4 (16-QAM)','6 (64-QAM)'});
xlabel('Spectral efficiency log_2(M)  (bit/s/Hz)'); ylabel('BER (Gray)'); ylim([1e-6 1]);
title('BER versus spectral efficiency at fixed E_b/N_0');
legend('E_b/N_0 = 0 dB','E_b/N_0 = 5 dB','E_b/N_0 = 10 dB','E_b/N_0 = 20 dB','Location','southeast');

fprintf('Done. Four figures generated (Figure 1-4).\n');


%% ======================= SUBFUNCTIONS =======================
function lab = pos2label(p, mapping)
    if strcmp(mapping,'gray')
        lab = bitxor(p, bitshift(p,-1));
    else
        lab = p;
    end

function [Craw,Cn] = make_const(M, ispsk, mapping)
    Craw = zeros(M,1);
    if ispsk
        p   = (0:M-1).';
        off = (M==4)*pi/4;
        pts = exp(1j*(2*pi*p/M + off));
        Craw(pos2label(p,mapping)+1) = pts;
    else
        L  = sqrt(M);
        l1 = pos2label((0:L-1).', mapping);
        for pI = 0:L-1
            for pQ = 0:L-1
                label = l1(pI+1)*L + l1(pQ+1);
                Craw(label+1) = (2*pI-(L-1)) + 1j*(2*pQ-(L-1));
            end
        end
    end
    Cn = Craw / sqrt(mean(abs(Craw).^2));

function lab = mindist_detect(r, C)
    dmin = inf(size(r)); lab = zeros(size(r));
    for m = 1:numel(C)
        d = abs(r - C(m)).^2;
        idx = d < dmin;
        dmin(idx) = d(idx); lab(idx) = m-1;
    end

function [rx,labels] = channel(Cn, M, EbN0dB, N)
    k = log2(M); N0 = 1/(k*10^(EbN0dB/10));
    labels = randi([0 M-1], N, 1);
    rx = Cn(labels+1) + sqrt(N0/2)*(randn(N,1)+1j*randn(N,1));

function [ber,ser] = simulate(Cn, M, EbN0dB, N)
    k = log2(M);
    [rx,labels] = channel(Cn, M, EbN0dB, N);
    det = mindist_detect(rx, Cn);
    ser = mean(det ~= labels);
    x = bitxor(labels, det); nb = 0;
    for i = 1:k, nb = nb + sum(bitget(x,i)); end
    ber = nb/(N*k);

function h = max_neighbor_hamming(C, M)
    k = log2(M);
    D = abs(bsxfun(@minus, C, C.')); D(D==0) = inf;
    dmin = min(D(:)); h = 0;
    for a = 1:M
        nb = find(D(a,:) < dmin*1.001);
        for b = nb
            x = bitxor(a-1,b-1); c = 0;
            for i = 1:k, c = c + bitget(x,i); end
            h = max(h,c);
        end
    end