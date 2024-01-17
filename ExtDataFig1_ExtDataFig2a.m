%% Import data for MitoPathways (including sublocation and custom sets)
[MnumData MtextData MrawData] = xlsread('MitoPathways.xlsx','MitoPathways');
for i = 1:157;
    pwname{i,1} = MtextData{i,1};
    delimiter = ', ';
    pwset{i} = strsplit(MtextData{i,2}, delimiter);
end


%% Import the Z-scores (no dummies)
[numData textData rawData] = xlsread('All_Z_scores.xlsx');


%% Names of genes
for i = 2:length(rawData);
    genename{i-1,1} = rawData{i,1};
end


%% Samples
%1'G/M'	
%2'Succ+Pier'	
%3'Asc/TMPD+Anti'

for i = 1:3;
    setname{i} = rawData{1,i+1};
end

ngene = max(size(genename));


%% For the G/M set
x = 1;
CurrentSet = numData(:,x);

for i = 1:157;
    k = 0;
    for j = 1:length(pwset{i});
        if ~isempty(find(ismember(genename,pwset{i}{1,j}), 1));
               ind = find(ismember(genename,pwset{i}{1,j}));
               k = k+1;
               Z_pw{i}(k) = CurrentSet(ind(1));
        end
    end
end


%% Plotting selected multiple CDFs vs. all proteins - Ext Data Fig 1
hFig = figure(2);
set(hFig, 'Position', [0 0 300 1200]);

sorted_all = sort(numData(:,x));
n1 = numel(sorted_all);
cdf1 = (1:n1)/n1;

for i = 1:156;
    setsize(i) = length(Z_pw{i});
    sorted_pw{i} = sort(Z_pw{i});
    npw(i) = numel(sorted_pw{i});
    cdfpw{i} = (1:npw(i))/npw(i);
end

% Gene sets of interest (OXPHOS subunits and assembly factor)
pind = [1 39 40 42 43 45 46 54 48 49];

for i = 1:10;
    subplot(10,1,i);
    plot(sorted_all,cdf1,'color',[0,0,0]);hold on;
    plot(sorted_pw{pind(i)}, cdfpw{pind(i)},'color',[1,0,0]);hold on;
    set(gca,'FontSize',14);
  
    xlim([-8 5]);

    h=ylabel('proportion');
    set(h,'FontSize',14,'FontWeight','normal');

    h=xlabel(strcat('\itZ\rm (',setname{x},')'));
    set(h,'FontSize',14,'FontWeight','normal');

if i==1;
    h=legend('MitoPlus library (n=1864)', strcat(pwname{pind(i)},' (n=',string(setsize(pind(i))),')'),'Location','northwest');
else
    h=legend('', strcat(pwname{pind(i)},' (n=',string(setsize(pind(i))),')'),'Location','northwest');
end
    legend boxoff
    set(h,'FontSize',14,'FontWeight','normal');
end



%% Calculate the Culmulative hypergeometric P values - Ext Data Fig 2a
for x = 1:3;
clear Z_pw;
CurrentSet = numData(:,x);
cutoff = -2;
N = sum(CurrentSet < cutoff);
Ln = length(CurrentSet);

for i = 1:156;
    k = 0;
    for j = 1:length(pwset{i});
        if ~isempty(find(ismember(genename,pwset{i}{1,j}), 1));
               ind = find(ismember(genename,pwset{i}{1,j}));
               k = k+1;
               Z_pw{i}(k) = CurrentSet(ind(1));
        end
    end
end

for i = 1:156;
        if ~isempty(Z_pw{i});
            M = sum(Z_pw{i} < cutoff);
            Lm = length(Z_pw{i});
            chosPV(i,x) = 1-hygecdf(M,Ln,N,Lm);
        else
            chosPV(i,x) = 1;
        end
end

end


%% Numbering of the pathways
%1Mitochondrial central dogma
%2mtDNA maintenance
%3mtDNA replication
%4mtDNA nucleoid
%5mtDNA repair
%6mtDNA modifications
%7mtDNA stability and decay
%8mtRNA metabolism
%9Transcription
%10mtRNA granules
%11Polycistronic mtRNA processing
%12mt-tRNA modifications
%13mt-rRNA modifications
%14mt-mRNA modifications
%15mtRNA stability and decay
%16Translation
%17Mitochondrial ribosome
%18Mitochondrial ribosome assembly
%19Translation factors
%20mt-tRNA synthetases
%21fMet processing
%22Protein import, sorting and homeostasis
%23Protein import and sorting
%24TOM
%25SAM
%26MIA40
%27TIM22 carrier pathway
%28TIM23 presequence pathway
%29Import motor
%30OXA
%31Preprotein cleavage
%32Protein homeostasis
%33Proteases
%34Chaperones
%35OXPHOS
%36OXPHOS subunits
%37OXPHOS assembly factors
%38Complex I
%39CI subunits
%40CI assembly factors
%41Complex II
%42CII subunits
%43CII assembly factors
%44Complex III
%45CIII subunits
%46CIII assembly factors
%47Complex IV
%48CIV subunits
%49CIV assembly factors
%50Complex V
%51CV subunits
%52CV assembly factors
%53Respirasome assembly
%54Cytochrome C
%55Metabolism
%56Carbohydrate metabolism
%57Gluconeogenesis
%58Malate-aspartate shuttle
%59Glycerol phosphate shuttle
%60Pyruvate metabolism
%61TCA cycle
%62TCA-associated
%63Itaconate metabolism
%64Ketone metabolism
%65Propanoate metabolism
%66Lipid metabolism
%67Fatty acid oxidation
%68Type II fatty acid synthesis
%69Lipoate insertion
%70Cardiolipin synthesis
%71Cholesterol, bile acid, steroid synthesis
%72Cholesterol-associated
%73Phospholipid metabolism
%74Eicosanoid metabolism
%75Amino acid metabolism
%76Branched-chain amino acid metabolism
%77Branched-chain amino acid dehydrogenase complex
%78Lysine metabolism
%79Serine metabolism
%80Glycine metabolism
%81Glycine cleavage system
%82Glutamate metabolism
%83Proline metabolism
%84Glyoxylate metabolism
%85GABA metabolism
%86Catechol metabolism
%87Kynurenine metabolism
%88Urea cycle
%89Nucleotide metabolism
%90Nucleotide import
%91Nucleotide synthesis and processing
%92Creatine metabolism
%93Metals and cofactors
%94Carnitine synthesis and transport
%95Carnitine shuttle
%96Coenzyme A metabolism
%97Coenzyme Q metabolism
%98Copper metabolism
%99Iron homeostasis
%100Heme synthesis and processing
%101Fe-S cluster biosynthesis
%102Molybdenum cofactor synthesis and proteins
%103NAD biosynthesis and metabolism
%104Tetrahydrobiopterin synthesis
%105Fe-S-containing proteins
%106Heme-containing proteins
%107Vitamin metabolism
%108Choline and betaine metabolism
%109Biotin utilizing proteins
%110Vitamin A metabolism
%111Vitamin B1 metabolism
%112Vitamin B2 metabolism
%113Vitamin B6 metabolism
%114Vitamin B12 metabolism
%115Vitamin C metabolism
%116Vitamin D metabolism
%117Folate and 1-C metabolism
%118Detoxification
%119Xenobiotic metabolism
%120ROS and glutathione metabolism
%121Amidoxime reducing complex
%122Selenoproteins
%123Electron carriers
%124Cytochromes
%125Q-linked reactions, other
%126Sulfur metabolism
%127Small molecule transport
%128SLC25A family
%129ABC transporters
%130Sideroflexins
%131Calcium uniporter
%132Signaling
%133Calcium homeostasis
%134Calcium cycle
%135Mitochondrial permeability transition pore
%136EF hand proteins
%137Immune response
%138cAMP-PKA signaling
%139Mitochondrial dynamics and surveillance
%140Fusion
%141Fission
%142Organelle contact sites
%143Intramitochondrial membrane interactions
%144Trafficking
%145Mitophagy
%146Autophagy
%147Apoptosis
%148Cristae formation
%149MICOS complex
%150Complex V
%151Matrix
%152MIM
%153IMS
%154MOM
%155MitoCarta3
%156OXPHOS