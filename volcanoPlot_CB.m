% volcanoPlot_CB.m
%
% This will produce a volcano plot of the protein spectra identified by the
% proteomics core at Princeton University from 10 cerebellar tissue samples. The data
% were exported from Scaffold 5 using the export -> protein information
% option to an excel file. Then the relevant columns (p-value, normalized
% spectra count, and protein identification probability) were converted into
% matlab variables and saved as mass_spec_data.mat. This will load that
% matlab file and produce a volcano plot which has on the x-axis Log2(fold change (Crus2 / lobule1))
% and on the y-axis the p-value of the lobular sample difference. This was
% written and tested on Matlab R2021a.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change directory to location of this code and mass_spec_data.mat file
% before running the code below
dataDir = pwd;
load(fullfile(dataDir,'mass_spec_data.mat'))

% Let's add a pseudocount of 1 to avoid infinity when calculating
% log2(fold-change)
mass_spec_norm_total_spectra = mass_spec_norm_total_spectra + 1;

% Calculate fold change, first averaging samples, and finally log2
% First 5 columns are Crus 2
avgSpectra = [mean(mass_spec_norm_total_spectra(:,1:5),2) mean(mass_spec_norm_total_spectra(:,6:10),2) ];
foldChange = avgSpectra(:,1) ./ avgSpectra(:,2);
logSpectra = log2(foldChange);

% Volcano plot
fsize = [9 8];
f = figure('units','inches','position',[0 0 fsize],'color','w','PaperSize',fsize);
scatter(logSpectra, -log10(mass_spec_pvals), 35, 'markerfacecolor', [112,224,187]./255, 'markeredgecolor', [112,224,187]./255, 'markerfacealpha',0.2, 'markeredgealpha',0.2);
hold on; 
scatter(logSpectra(mass_spec_pvals<0.05), -log10(mass_spec_pvals(mass_spec_pvals<0.05)), 35, 'markerfacecolor', [55,142,142]./255, 'markeredgecolor', [55,142,142]./255, 'markerfacealpha',1, 'markeredgealpha',0.3);
protsOfInterest = {'PVALB','MBP','ATP1A3','NFASC','NCAM2','BSN','NLGN3','CALB','CALML3','LMNB1','HNRNPA1','VIM'};
for p = 1:length(protsOfInterest)
    ind = find(matches(names,protsOfInterest{p}));
    text(logSpectra(ind),-log10(mass_spec_pvals(ind)),['\leftarrow ' names{ind}],'fontsize',12,'fontweight','bold')
end
grid on;
set(gca,'xlim',[-3 3],'tickdir','out','fontsize',18,'linewidth',2,'GridColor',[0.65 0.65 0.65])
xlabel({'log2(Fold Change)',['Lob I-IV ' '\leftarrow' '      ' '\rightarrow' 'Crus 2']})
xlabel(['Lob I-IV ' '\leftarrow' '      ' 'log2(Fold Change)' '      ' '\rightarrow' 'Crus 2'])
ylabel('-log10(p-value)')

% Uncomment the line below if you want to save figure out
%saveas(gcf,'volcanoPlot','pdf')


