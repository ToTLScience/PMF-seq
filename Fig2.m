load oxphos_annot.mat

%% Import raw counts
[numData textData rawData] = xlsread('raw counts.xlsx');

%% Names of sgRNAs
for i = 2:15272;
sgnames{i-1,1} = rawData{i,2};
end
clear i

%% Transform num data
numData2 = numData(:,6:end);
nsamples = min(size(numData2));

for i = 1:nsamples;
    totalct(i) = sum(numData2(:,i));
    numData3(:,i) = log2(numData2(:,i)/totalct(i)*1e6+1);
    totalct2(i) =  sum(numData3(:,i));
end

%% From sgRNA to gene
genenames = unique(sgnames);

for  i = 1:length(genenames); 
    ind = find(ismember(sgnames,genenames{i}));
    if length(ind)>1;
        geneData(i,:) = mean(numData3(ind,:));
        catog(i,:) = mean(numData(ind,1:5));
    else if length(ind) == 1;
            geneData(i,:) = numData3(ind,:);
            catog(i,:) = mean(numData(ind,1:5));
        end
    end
end

%% Calculate L2FC (high tail - low tail)

GM_1 = geneData(:,2)-geneData(:,1);
GM_2 = geneData(:,4)-geneData(:,3);

PSucc_1 = geneData(:,6)-geneData(:,5);
PSucc_2 = geneData(:,8)-geneData(:,7);

AAscTMPD_1 = geneData(:,10)-geneData(:,9);
AAscTMPD_2 = geneData(:,12)-geneData(:,11);

GM = (GM_1+GM_2)/2;
PSucc = (PSucc_1+PSucc_2)/2;
AAscTMPD = (AAscTMPD_1+AAscTMPD_2)/2;

m_LFC = [GM_1 GM_2 PSucc_1 PSucc_2 AAscTMPD_1 AAscTMPD_2...
    GM PSucc AAscTMPD];
    

%% Numbering the sets 

%1 Glu/Mal Rep 1
%2 Glu/Mal Rep 2
%3 Succinate(+Pier) Rep 1
%4 Succinate(+Pier) Rep 2
%5 Ascorbate+TMPD(+Anti) Rep 1
%6 Ascorbate+TMPD(+Anti) Rep 2
%7 Mean Glu/Mal
%8 Mean Succinate(+Pier)
%9 Mean Ascorbate+TMPD(+Anti)
boxlabel2 = {'Glu/Mal Rep 1','Glu/Mal Rep 2'...
    'Succinate(+Pier) Rep 1','Succinate(+Pier) Rep 2',...
    'Ascorbate+TMPD(+Anti) Rep 1','Ascorbate+TMPD(+Anti) Rep 2'...
    'Glu/Mal','Succinate(+Pier)','Ascorbate+TMPD(+Anti)'};   

%% null distribution (non-expressed genes) statistics
ind = find(catog(:,2)==1);
null_d = m_LFC(ind,:);
mean_null = mean(null_d);
std_null = std(null_d);
var_null = var(null_d);

%% Dummy distribution (non-cutting guides)
ind = find(catog(:,4)==max(catog(:,4)));
dummy_LFC = m_LFC(ind,:);
dummy_LFC = m_LFC(ind,:);

%% Core essential genes
ind = find(catog(:,5)==1);
core_ess_LFC = m_LFC(ind,:);
core_ess_LFC = m_LFC(ind,:);

%% Calculate z score

for i = 1:9;
    for j = 1:length(m_LFC);
        Z(j,i) = (m_LFC(j,i)-mean_null(i))/std_null(i);
    end
end


%% Z-scores for OXPHOS genes
%CI CII CIII CIV CV
Clear ind
ind = find(catog2(:,2)==1);
Z_CI = Z(ind,:);

clear ind
ind = find(catog2(:,3)==1);
Z_CII = Z(ind,:);

clear ind
ind = find(catog2(:,4)==1);
Z_CIII = Z(ind,:);

clear ind
ind = find(catog2(:,5)==1);
Z_CIV = Z(ind,:);

clear ind
ind = find(catog2(:,6)==1);
Z_CV = Z(ind,:);

GOI = 'CYCS';
clear ind Z_GOI;
ind = find(ismember(genenames,GOI));
Z_GOI = Z(ind,:);


%% Plotting Oxphos cocktail - Fig 2a
hFig = figure(1);
set(hFig, 'Position', [0 0 250 1500]);
clf

% For example: plotting G/M Rep 1 and G/M Rep 2
Z_bg = Z_nonexp; % non-expressed genes as the negative controls
b=1;
f=2;

sxlim = [-15 10];
sylim = [-15 10];


    subplot(5,1,1);
    plot(Z_bg(:,b),Z_bg(:,f),'.','markers',15,'color',[0.9,0.9,0.9]);hold on;
    plot(Z_CI(:,b),Z_CI(:,f),'.','markers',15,'color',[1,0,0]);hold on;
    plot([-20 20],[0 0],'k:');hold on;
    plot([0 0],[-20 20],'k:');hold off;
    box off
    set(gca,'xscale','linear');
    set(gca,'yscale','linear');
    set(0,'DefaultAxesTitleFontWeight','normal');
    xlim(sxlim);
    ylim(sylim);
    set(gca,'FontSize',24);
  
    title('complex I');
    h=ylabel(boxlabel2{f});
    set(h,'FontSize',24,'FontWeight','normal');

    h=xlabel(boxlabel2{b});
    set(h,'FontSize',24,'FontWeight','normal');


    subplot(5,1,2);
    plot(Z_bg(:,b),Z_bg(:,f),'.','markers',15,'color',[0.9,0.9,0.9]);hold on;
    plot(Z_CII(:,b),Z_CII(:,f),'.','markers',15,'color',[0.5,0.2,0.5]);hold on;
    plot([-20 20],[0 0],'k:');hold on;
    plot([0 0],[-20 20],'k:');hold off;
    box off
    set(gca,'xscale','linear');
    set(gca,'yscale','linear');
    set(0,'DefaultAxesTitleFontWeight','normal');
    xlim(sxlim);
    ylim(sylim);
    set(gca,'FontSize',24);
  
    title('complex II');
    h=ylabel(boxlabel2{f});
    set(h,'FontSize',24,'FontWeight','normal');

    h=xlabel(boxlabel2{b});
    set(h,'FontSize',24,'FontWeight','normal');


    subplot(5,1,3);
    plot(Z_bg(:,b),Z_bg(:,f),'.','markers',15,'color',[0.9,0.9,0.9]);hold on;
    plot(Z_CIII(:,b),Z_CIII(:,f),'.','markers',15,'color',[0.5,0.7,0.2]);hold on;
    plot([-20 20],[0 0],'k:');hold on;
    plot([0 0],[-20 20],'k:');hold off;
    box off
    set(gca,'xscale','linear');
    set(gca,'yscale','linear');
    set(0,'DefaultAxesTitleFontWeight','normal');
    xlim(sxlim);
    ylim(sylim);
    set(gca,'FontSize',24);
  
    title('complex III');
    h=ylabel(boxlabel2{f});
    set(h,'FontSize',24,'FontWeight','normal');

    h=xlabel(boxlabel2{b});
    set(h,'FontSize',24,'FontWeight','normal');


    subplot(5,1,4);
    plot(Z_bg(:,b),Z_bg(:,f),'.','markers',15,'color',[0.9,0.9,0.9]);hold on;
    plot(Z_GOI(:,b),Z_GOI(:,f),'.','markers',15,'color',[0.5,0.7,1]);hold on;
    plot([-20 20],[0 0],'k:');hold on;
    plot([0 0],[-20 20],'k:');hold off;
    box off
    set(gca,'xscale','linear');
    set(gca,'yscale','linear');
    set(0,'DefaultAxesTitleFontWeight','normal');
    xlim(sxlim);
    ylim(sylim);
    set(gca,'FontSize',24);
  
    title(GOI);
    h=ylabel(boxlabel2{f});
    set(h,'FontSize',24,'FontWeight','normal');

    h=xlabel(boxlabel2{b});
    set(h,'FontSize',24,'FontWeight','normal');

      
    subplot(5,1,5);
    plot(Z_bg(:,b),Z_bg(:,f),'.','markers',15,'color',[0.9,0.9,0.9]);hold on;
    plot(Z_CIV(:,b),Z_CIV(:,f),'.','markers',15,'color',[0,0,1]);hold on;
    plot([-20 20],[0 0],'k:');hold on;
    plot([0 0],[-20 20],'k:');hold off;
    box off
    set(gca,'xscale','linear');
    set(gca,'yscale','linear');
    set(0,'DefaultAxesTitleFontWeight','normal');
    xlim(sxlim);
    ylim(sylim);
    set(gca,'FontSize',24);
  
    title('complex IV');
    h=ylabel(boxlabel2{f});
    set(h,'FontSize',24,'FontWeight','normal');

    h=xlabel(boxlabel2{b});
    set(h,'FontSize',24,'FontWeight','normal');
