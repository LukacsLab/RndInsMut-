%% Plotting cluster results

addpath('C:\Users\phyjga\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis\CHO6-11LIB\Sequencing2analysis\04-cluster2020')
%addpath('C:\Users\julia\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis\CHO6-11LIB\Sequencing2analysis\04-cluster2020')

clc
clear all
close all

cd('C:\Users\phyjga\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis\CHO6-11LIB\Sequencing2analysis\03-merged2020')
%cd('C:\Users\julia\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis\CHO6-11LIB\Sequencing2analysis\03-merged2020')
FileClusters=fopen('CHO6-11LIB-Tol2_CPM.tsv');
% Currently at the begining of the file 
fseek(FileClusters,0,'eof') % move to the end
fileSize=ftell(FileClusters);
fseek(FileClusters,0,'bof') % to move to the beginning of the file
ClusterDataTemp={};
count=0;
rowcount=0;
while count<fileSize
    lineString = fgets(FileClusters);
    clusterNumber=strfind(lineString,'CLUSTER');
    if clusterNumber > 0
    else
        ClusterDataTemp = [ClusterDataTemp; textscan(lineString,'%s %f %f %s %s %f %f')];
    end
    count=ftell(FileClusters);
    if mod(rowcount,1000)==0
        disp(['On line number ' num2str(rowcount)]);
    end
    rowcount=rowcount+1;
end
disp('End')
% Where the columns correspond to:
% #1: Chromosome reference
% #2: Position 1
% #3: Position 2
% #4: Sample name
% #5: Raw counts
% #6: Counts per million
% #7: Polarity
ClusterData={};
% Change from array of cells to array of strings
ClusterData(1)={string(ClusterDataTemp(:,1))};
ClusterData(2)={ClusterDataTemp(:,2)};
ClusterData(3)={ClusterDataTemp(:,3)};
ClusterData(4)={string(ClusterDataTemp(:,4))};
ClusterData(5)={ClusterDataTemp(:,6)};
ClusterData(6)={ClusterDataTemp(:,7)};
ClusterData(7)={ClusterDataTemp(:,5)};

clear ClusterDataTemp

%% Plot counts per insertion location hystogram

figure
Increment=100000;
Start=0;
End=1.1e6;
binEdges=Start:Increment:End;
binCenters=Start+Increment/2:Increment:End-Increment/2;
n=histogram(cell2mat(ClusterData{6}),'BinEdges',binEdges);
bar(binCenters,n.Values)
set(gca, 'YScale', 'log')
xlabel('CPM per location')
ylabel('Number of locations')

figure
Increment=1;%500;
Start=0;
End=50;%7e3;
binEdges=Start:Increment:End;
binCenters=Start+Increment/2:Increment:End-Increment/2;
n=histogram(cell2mat(ClusterData{5}),'BinEdges',binEdges);
bar(binCenters,n.Values)
xlim([0.5 End])
%set(gca, 'YScale', 'log')
xlabel('Raw reads per location')
ylabel('Number of Locations')

%% Filter by number of raw reads and CPM

CPMThreshold=1;
rawReadsThreshold=2;

ClusterData{1}(cell2mat(ClusterData{6})<CPMThreshold)=[];
ClusterData{2}(cell2mat(ClusterData{6})<CPMThreshold)=[];
ClusterData{3}(cell2mat(ClusterData{6})<CPMThreshold)=[];
ClusterData{4}(cell2mat(ClusterData{6})<CPMThreshold)=[];
ClusterData{5}(cell2mat(ClusterData{6})<CPMThreshold)=[];
ClusterData{7}(cell2mat(ClusterData{6})<CPMThreshold)=[];
ClusterData{6}(cell2mat(ClusterData{6})<CPMThreshold)=[];

ClusterData{1}(cell2mat(ClusterData{5})<rawReadsThreshold)=[];
ClusterData{2}(cell2mat(ClusterData{5})<rawReadsThreshold)=[];
ClusterData{3}(cell2mat(ClusterData{5})<rawReadsThreshold)=[];
ClusterData{4}(cell2mat(ClusterData{5})<rawReadsThreshold)=[];
ClusterData{6}(cell2mat(ClusterData{5})<rawReadsThreshold)=[];
ClusterData{7}(cell2mat(ClusterData{5})<rawReadsThreshold)=[];
ClusterData{5}(cell2mat(ClusterData{5})<rawReadsThreshold)=[];

%% Create structure array

chromNames = unique(ClusterData{1}); 
ChromClusterData=struct;
ChromClusterData.Chromosomes=ClusterData{1};
ChromClusterData.StartPosition=ClusterData{2};
ChromClusterData.EndPosition=ClusterData{3};
ChromClusterData.SampleName=ClusterData{4};
ChromClusterData.CountsPerMillion=ClusterData{6};
ChromClusterData.RawCounts=ClusterData{5};
ChromClusterData.Polarity=ClusterData{7};

%% Clone labelling

% Create a cell with the labels for each of the sample names

% column 1 is sample name
% column 2 is phenotype
% column 3 is numeric certainty: 1 is certain, 0 is uncertain

% For the first sequencing
% SampleLabels={...
%     'A1_S1',    'RT?',      0,      'P2G4',     'A1';...
%     'A2_S7',    'AITC+',    1,      'P4F7',     'A2';...
%     'A3_S13',   'C+M+',     1,      'P5D6',     'A3';...
%     'A4_S19',   'C+M+(h)',  1,      'P1C11',    'A4';...
%     'A5_S25',   'WS12+',    1,      'P1D9',     'A5';...
%     'A6_S31',   'WS12+',    1,      'P1F2',     'A6';...
%     'B1_S2',    'C+  M-',   0,      'P5G6',     'B1';... % Read depth <10%
%     'B2_S8',    'AITC+',    1,      'P2C6',     'B2';...
%     'B3_S14',   'AITC+',    1,      'P2F7',     'B3';...
%     'B4_S20',   'AITC+',    1,      'P2G9',     'B4';...
%     'B5_S26',   'AITC+',    1,      'P1C3',     'B5';...
%     'B6_S32',   'N',        0,      'P5G3',     'B6';... % Read depth <10%
%     'C1_S3',    'C+  M-',   0,      'P1B8',     'C1';...
%     'C2_S9',    'N',        0,      'P2G11',    'C2';...
%     'C3_S15',   'N',        0,      'P3E3',     'C3';... % Read depth <10%
%     'C4_S21',   'N(*)',     0,      'P4G9',     'C4';...
%     'C5_S27',   'N',        0,      'P5C7',     'C5';...
%     'C6_S33',   'N',        0,      'P3C3',     'C6';... % Read depth <10%
%     'D1_S4',    'N',        0,      'P1C9',     'D1';... % Read depth <10%
%     'D2_S10',   'C+++M-',   1,      'P2B2',     'D2';...
%     'D3_S16',   'C+++M-',   1,      'P2C2',     'D3';...
%     'D4_S22',   'C+  M-',   0,      'P2D9',     'D4';... % Read depth <10%
%     'D5_S28',   'C+++M-(*)',1,      'P2D11',    'D5';...
%     'D6_S34',   'C+  M-',   0,      'P3D7',     'D6';... % Read depth <10%
%     'E1_S5',    'N',        0,      'P3E6',     'E1';...
%     'E2_S11',   'N',        1,      'P3E9',     'E2';...
%     'E3_S17',   'C+++M-(*)',0,      'P3F5',     'E3';...
%     'E4_S23',   'N',        1,      'P4B9',     'E4';...
%     'E5_S29',   'C++ M-',   1,      'P6F8',     'E5';...
%     'E6_S35',   'C+++M-',   1,      'P6E11',    'E6';...
%     'F1_S6',    'C++ M-',   0,      'P6D9',     'F1';... % Read depth <10%
%     'F2_S12',   'C++ M-',   1,      'P5E6',     'F2';...
%     'F3_S18',   'C+  M-',   0,      'P5E5',     'F3';...
%     'F4_S24',   'AITC+',    1,      'P5B2',     'F4';...
%     'F5_S30',   'C++ M-',   0,      'P4F8',     'F5';...
%     'F6_S36',   'N',        0,      'P4B11',    'F6'};

% For the second sequencing
SampleLabels={...   %Well   %Pheno  %Certainty  %Name   %WellShort  %GroupInterest  %GroupPheno
'G11_S46',  'C+M-',     1,	'P1B8',             'G11',      30,     3;...
'G3_S5',	'WS12+',    1,  'P1C11',            'G3',       0,      1;...
'G7_S19',   'N',        1,  'P1C11 neg pheno',  'G7',       0,      0.9;...
'G10_S39',  'AITC+',    1,  'P1C3',             'G10',      0,      2;...
'D7_S16',   'N',        0,  'P1C9',             'D7',       0,      0;... % low read depth
'A7_S13',   'C++M-',    1,  'P1D7',             'A7',       30,     0;...
'G4_S7',    'WS12+',    1,  'P1D9',             'G4',       0,      1;...
'A9_S27',   'N',        1,  'P1E2',             'A9',       40,     1;...
'A10_S33',  'unknown',  1,  'P1E4',             'A10',      0,      4;...
'G5_S9',    'WS12+',    1,  'P1F2',             'G5',       0,      1;...
'D8_S23',   'C+++M-',   1,  'P2B2',             'D8',       30,     3;...
'D9_S30',   'C+++M-',   1,  'P2C2',             'D9',       30,     3;...
'G6_S11',   'AITC+',    1,  'P2C6',             'G6',       0,      2;...
'D11_S43',  'C+++M-',   1,  'P2D11',            'D11',      20,     3;...
'D10_S36',  'C+  M-',   0,  'P2D9',             'D10',      20,     3;... % low read depth
'C7_S15',   'N',        1,  'P2F11',            'C7',       10,     0;...
'G8_S26',   'AITC+',    1,  'P2F7',             'G8',       0,      2;...
'G1_S1',    'AITC+',    1,  'P2G4',             'G1',       0,      2;...
'G9_S32',   'AITC+',    1,  'P2G9',             'G9',       0,      2;...
'C12_S49',  'N',        0,  'P3C3',             'C12',      40,     0;... % low read depth
'C10_S35',  'unknown',  1,  'P3C4',             'C10',      0,      2;...
'D12_S50',  'C+  M-',   0,  'P3D7',             'D12',      0,      3;... % low read depth
'C9_S29',   'N',        0,  'P3E3',             'C9',       0,      0;... % low read depth
'H2_S4',    'N',        1,  'P3E6',             'H2',       40,     0;...
'H3_S6',    'N',        1,  'P3E9',             'H3',       40,     0;...
'C11_S42',  'unknown',  1,  'P3F10',            'C11',      0,      2;...
'B9_S28',   'C++M-',    1,  'P3F4',             'B9',       20,     3;...
'E9_S31',   'C+++M-',   1,  'P3F5',             'E9',       40,     4;...
'A11_S40',  'AITC+',    1,  'P3F9',             'A11',      0,      2;...
'F10_S38',  'unknown',  0,  'P3G3',             'F10',      0,      4;... % low read depth
'H7_S20',   'C+M-?',    1,  'P4B11',            'H7',       10,     1;...
'H4_S8',    'N',        1,  'P4B9',             'H4',       40,     0;...
'F12_S52',  'N',        1,  'P4C5',             'F12',      40,     0;...
'B10_S34',  'C+M-?',    1,  'P4F3',             'B10',      10,     4;...
'A8_S21',   'AITC+',    1,  'P4F7',             'A8',       0,      2;...
'F11_S45',  'C++ M-',   1,  'P4F8',             'F11',      10,     3;...
'G12_S53',  'N',        1,  'P4G9',             'G12',      0,      0;...
'B11_S41',  'C+M-?',    1,  'P5B10',            'B11',      10,     1;...
'H6_S12',   'AITC+',    1,  'P5B2',             'H6',       0,      2;...
'H1_S2',    'N',        1,  'P5C7',             'H1',       0,      0;...
'E7_S17',   'C+M-',     1,  'P5D2',             'E7',       10,     4;...
'E8_S24',   'C++M-',    1,  'P5D4',             'E8',       30,     3;...
'G2_S3',    'WS12+',    1,  'P5D6',             'G2',       0,      1;...
'F8_S25',   'C++ M-',   1,  'P5E6',             'F8',       20,     4;...
'A12_S47',  'AITC+',    1,  'P5F9',             'A12',      0,      2;...
'B12_S48',  'unknown',  0,  'P5G3',             'B12',      0,      4;... % low read depth
'B7_S14',   'C+  M-',   0,  'P5G6',             'B7',       0,      3;... % low read depth
'E10_S37',  'C++M-',    1,  'P6C3',             'E10',      30,     0;...
'F7_S18',   'C++ M-',   0,  'P6D9',             'F7',       0,      3;... % low read depth
'H5_S10',   'C+++M-',   1,  'P6E11',            'H5',       30,     3;...
'B8_S22',   'N',        1,  'P6E2',             'B8',       30,     4;...
'E12_S51',  'AITC+',    1,  'P6F11',            'E12',      0,      2;...
'E11_S44',  'C++ M-',   1,  'P6F8',             'E11',      20,     3};

% Groups of Interest
% Group 10: Related to P4F8
% Group 20: Cold group of P2D11
% Group 30: Cold group of P2B2
% Group 40: Related to P3F5

Temp=sortrows(SampleLabels,[6,7,2]);
SampleLabels=Temp;

% To set all certainties to 1
%SampleLabels(:,3)={1};

%% Adding the phenotype label to ChromClusterData

ChromClusterData.Phenotype={};
for i=1:length(ChromClusterData.SampleName)
    for j=1:length(SampleLabels(:,1))
        if strcmp(ChromClusterData.SampleName(i),SampleLabels(j,1))
            ChromClusterData.Phenotype(i)=SampleLabels(j,2);
            ChromClusterData.Certainty(i)=SampleLabels(j,3);
            ChromClusterData.CloneName(i)=SampleLabels(j,4);
        end
    end
end

%% Once we have the correlation matrix, we can sort the samples to minimize
% distance between them, where correlated samples have low distance and
% non-correlated samples have high distance.
% I make up a metric called corrdist
% corrdist(a,b) = 1 - corr(a,b); % It is zero for perfectly equal samples.
% Calculate the distance matrix

% if exist('CorrMatAll')
%     DistMat = ones(size(CorrMatAll))-CorrMatAll; % Distances of clone i with clone j in position i,j or j,i
%     for i=1:size(DistMat,1) % Set diagonals as infinites
%         DistMat(i,i)=Inf;
%     end    
%     % Indices
%     idx = 1; % starting index
%     out_idx = zeros(size(DistMat,1),1); % output order
%     out_idx(1) = idx; % first index
%     DistMatTemp=DistMat;
%     for k = 2:size(DistMatTemp,1)
%         start_ind = idx;
%         [~,idx] = min(DistMatTemp(start_ind,:));
%         DistMatTemp(:,start_ind) = Inf;
%         out_idx(k) = idx;
%     end
%     count=1;
%     for i=1:length(out_idx)
%         DistMat(count,:)=DistMat(out_idx(i),out_idx);
%         count=count+1;
%     end
%     SampleLabels=SampleLabels(out_idx,:);
% end   


%% Create joint label for plots
SampleLabelsPlot={};
for i=1:length(SampleLabels)
    SampleLabelsPlot{i}=[SampleLabels{i,5} '-' SampleLabels{i,4} ' (' SampleLabels{i,2} ')']; 
end

ResponsesNumeric=cell2mat(SampleLabels(:,3));

%% Create global coordinates to be able to plot all insertions in the same plot

% Import contigs lengths for each of them

cd('C:\Users\phyjga\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis\CHO6-11LIB')
%cd('C:\Users\julia\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis\CHO6-11LIB')

FileContigs=fopen('CriGri_PICRH-1_0_contigs.csv');
% Currently at the begining of the file 
fseek(FileContigs,0,'eof'); % move to the end
fileSize=ftell(FileContigs);
fseek(FileContigs,0,'bof'); % to move to the beginning of the file
DataTemp=[];
Headers=fgets(FileContigs);
count=0;
rowCount=1;
ContigsLength=struct;
ContigsLength.Chromosomes={};
ContigsLength.TotalLength=[];
ContigsLength.Name={};
while count<fileSize
    lineString = fgets(FileContigs);
    DataTemp = strsplit(lineString,',');
    try
    ContigsLength.Chromosomes{rowCount}=DataTemp{1};
    ContigsLength.TotalLength(rowCount)=str2num(DataTemp{2});
    ContigsLength.Name(rowCount)=DataTemp(3);
    catch
        disp(['It didnt work for row' num2str(rowCount)])
    end
    count=ftell(FileContigs);
    rowCount=rowCount+1;
end

% Create global coordinates for start and end for each of the insertions

chromEndsGlobal=[];

for i=1:length(chromNames)
    LogicalIndex=ChromClusterData.Chromosomes==chromNames(i);
    StartGlobalTemp=[]; EndGlobalTemp=[]; StartTemp=[]; EndTemp=[];
    if i==1
        ChromClusterData.StartPositionGlobal(LogicalIndex)=...
            ChromClusterData.StartPosition(LogicalIndex);
        ChromClusterData.EndPositionGlobal(LogicalIndex)=...
            ChromClusterData.EndPosition(LogicalIndex);
        CumLength=double(ContigsLength.TotalLength(ContigsLength.Chromosomes==chromNames(i)));
    else
        StartTemp=[ChromClusterData.StartPosition{LogicalIndex}];
        EndTemp=[ChromClusterData.EndPosition{LogicalIndex}];
        for j=1:length(StartTemp)
            StartGlobalTemp(j)=StartTemp(j)+CumLength;
            EndGlobalTemp(j)=EndTemp(j)+CumLength;
        end
        count=1;
        for j=1:length(LogicalIndex)
            if LogicalIndex(j)==1
                ChromClusterData.StartPositionGlobal{j}=StartGlobalTemp(count);
                ChromClusterData.EndPositionGlobal{j}=EndGlobalTemp(count);
                count=count+1;
            end
        end
        CumLength=CumLength+ContigsLength.TotalLength(ContigsLength.Chromosomes==chromNames(i));
        clear *Temp
    end
    chromEndsGlobal(i)=CumLength;
end

%% Calculate the number of insertions per chromosome

ChromCount=[];
for i=1:length(chromNames)
    CountTemp=0;
    for j=1:length(ChromClusterData.Chromosomes)
        count=strcmp(chromNames(i),ChromClusterData.Chromosomes(j));
        CountTemp=CountTemp+count;
    end
    for j=1:length(ContigsLength.Chromosomes)
        if strcmp(chromNames(i),ContigsLength.Chromosomes(j))
            Length=ContigsLength.TotalLength(j);
        end
    end
    ChromCount(i,1)=CountTemp;
    disp(ChromCount(i))
    ChromCount(i,2)=ChromCount(i)/Length;
    disp(ChromCount(i,2))
    Length=[];
end

figure
bar(ChromCount(:,1))
set(gca,'XTick',1:length(chromNames))
set(gca,'XTickLabel',chromNames)
set(gca,'XTickLabelRotation',90)
xlabel('Chromosome')
ylabel('Number of Insertions (#)')

figure
bar(ChromCount(:,2))
set(gca,'XTick',1:length(chromNames))
set(gca,'XTickLabel',chromNames)
set(gca,'XTickLabelRotation',90)
xlabel('Chromosome')
ylabel('Insertion density (#/bp)')
set(gca,'YScale','log')



%% Quantify number of insertions per clone

NumInsertions.SampleName=SampleLabels(:,1);
NumInsertions.SampleName2=SampleLabels(:,4);
NumInsertions.Certainty=SampleLabels(:,3);
for i=1:length(NumInsertions.SampleName)
    NumInsertions.Count(i)=sum(ChromClusterData.SampleName==NumInsertions.SampleName(i));
    NumInsertions.Chromosomes{i}=ChromClusterData.Chromosomes(ChromClusterData.SampleName==NumInsertions.SampleName(i));
    NumInsertions.StartPosition{i}=ChromClusterData.StartPosition(ChromClusterData.SampleName==NumInsertions.SampleName(i));
    NumInsertions.StartPositionGlobal{i}=ChromClusterData.StartPositionGlobal(ChromClusterData.SampleName==NumInsertions.SampleName(i));
    NumInsertions.EndPosition{i}=ChromClusterData.EndPosition(ChromClusterData.SampleName==NumInsertions.SampleName(i));
    NumInsertions.EndPositionGlobal{i}=ChromClusterData.EndPositionGlobal(ChromClusterData.SampleName==NumInsertions.SampleName(i));
    NumInsertions.InsertionLength{i}=cell2mat(NumInsertions.StartPosition{i})-cell2mat(NumInsertions.EndPosition{i});
    NumInsertions.CountsPerMillion{i}=ChromClusterData.CountsPerMillion(ChromClusterData.SampleName==NumInsertions.SampleName(i));
    NumInsertions.RawCounts{i}=ChromClusterData.RawCounts(ChromClusterData.SampleName==NumInsertions.SampleName(i));
    NumInsertions.Polarity{i}=ChromClusterData.Polarity(ChromClusterData.SampleName==NumInsertions.SampleName(i));
    NumInsertions.Phenotype{i}=ChromClusterData.Phenotype(ChromClusterData.SampleName==NumInsertions.SampleName(i));
    NumInsertions.CloneName{i}=ChromClusterData.CloneName(ChromClusterData.SampleName==NumInsertions.SampleName(i));
end
NumInsertions.chromEndsGlobal=chromEndsGlobal;

figure
set(groot, 'defaultAxesTickLabelInterpreter','none')
bar(1:length(NumInsertions.Count),NumInsertions.Count)
set(gca,'XTick',1:length(NumInsertions.Count))
set(gca,'XTickLabel',SampleLabelsPlot)
set(gca,'XTickLabelRotation',90)
ylabel('Number of insertions')

figure
set(groot, 'defaultAxesTickLabelInterpreter','none')
histogram(NumInsertions.Count,10)
ylabel('Number of clones')
xlabel('Number of insertions')
ylim([0, 15])

%close all

%% Find local position and contig for a global position
% Only for start locations with an insertion
LocationQuerry=896001400;
IndexCandidate=find(abs(cell2mat(ChromClusterData.StartPositionGlobal)-LocationQuerry)<100000,1);
        
% Show data in that location
IndexStart=IndexCandidate-15;
IndexEnd=IndexCandidate+26;

for i=IndexStart:1:IndexEnd
disp([ChromClusterData.StartPositionGlobal(i)...
    ChromClusterData.SampleName{i}...
    ChromClusterData.CloneName{i}...
    ChromClusterData.Chromosomes(i)...
    ChromClusterData.StartPosition(i)...
    ChromClusterData.EndPosition(i)...
    ChromClusterData.Polarity(i)...
    ChromClusterData.CountsPerMillion(i)...
    ChromClusterData.RawCounts(i)])
end

%% Add data for certain genes of interest
% For the 2020 genome:

Genes=[];

Genes.TRPM8.Name='TRPM8';
Genes.TRPM8.Chromosome='NC_048595.1';
Genes.TRPM8.StartPosition=439889997;
Genes.TRPM8.EndPosition=439959602;
% Calculate global position
ChromIndex=find(strcmp(Genes.TRPM8.Chromosome,chromNames));
if ChromIndex==1
    Genes.TRPM8.StartPositionGlobal=Genes.TRPM8.StartPosition;
    Genes.TRPM8.EndPositionGlobal=Genes.TRPM8.EndPosition;
else
    Genes.TRPM8.StartPositionGlobal=Genes.TRPM8.StartPosition+chromEndsGlobal(ChromIndex-1);
    Genes.TRPM8.EndPositionGlobal=Genes.TRPM8.EndPosition+chromEndsGlobal(ChromIndex-1);
end    

Genes.TRPA1.Name='TRPA1';
Genes.TRPA1.Chromosome='NC_048595.1';
Genes.TRPA1.StartPosition=145849616;
Genes.TRPA1.EndPosition=145893304;
% Calculate global position
ChromIndex=find(strcmp(Genes.TRPA1.Chromosome,chromNames));
if ChromIndex==1
    Genes.TRPA1.StartPositionGlobal=Genes.TRPA1.StartPosition;
    Genes.TRPA1.EndPositionGlobal=Genes.TRPA1.EndPosition;
else
    Genes.TRPA1.StartPositionGlobal=Genes.TRPA1.StartPosition+chromEndsGlobal(ChromIndex-1);
    Genes.TRPA1.EndPositionGlobal=Genes.TRPA1.EndPosition+chromEndsGlobal(ChromIndex-1);
end

NumInsertions.Genes=Genes;

%% Plot the global coordinates for each of the clones
addpath('C:\Users\phyjga\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis')
%addpath('C:\Users\julia\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis')
cd('C:\Users\phyjga\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis\CHO6-11LIB\Sequencing2analysis\04-cluster2020')

close all
%GlobalInsertionsPlot(NumInsertions,'Raw reads',SampleLabelsPlot,1)
%GlobalInsertionsPlot_rectangle(NumInsertions,'Raw reads',SampleLabelsPlot,1)
GlobalInsertionsPlot_marker(NumInsertions,'Raw reads',SampleLabelsPlot,1)

%% Plot the global insertion landscape for only a few clones
% For WS12+ clones
%ChosenIndexes=round([SampleLabels{:,7}])==1;
% For AITC+ clones
%ChosenIndexes=round([SampleLabels{:,7}])==2;
% For C+M- clones
%ChosenIndexes=round([SampleLabels{:,7}])==3;
% For P4F8 group
% ChosenIndexes=[SampleLabels{:,6}]==10;
% For P2D11 group
% ChosenIndexes=[SampleLabels{:,6}]==20;
% For P2B2 group
% ChosenIndexes=[SampleLabels{:,6}]==30;
% For P3F5 group
% ChosenIndexes=[SampleLabels{:,6}]==40;
% For P1C11 and P1C11 neg pheno
%ChosenIndexes=zeros(length(SampleLabels(:,1)),1);
%ChosenIndexes(7)=1;
%ChosenIndexes(8)=1;

NumInsertionsTemp=NumInsertions;

% Edit the NumInsertions structure array
NumInsertionsTemp.SampleName(not(ChosenIndexes))=[];
NumInsertionsTemp.SampleName2(not(ChosenIndexes))=[];
NumInsertionsTemp.Certainty(not(ChosenIndexes))=[];
NumInsertionsTemp.Count(not(ChosenIndexes))=[];
NumInsertionsTemp.Chromosomes(not(ChosenIndexes))=[];
NumInsertionsTemp.StartPosition(not(ChosenIndexes))=[];
NumInsertionsTemp.StartPositionGlobal(not(ChosenIndexes))=[];
NumInsertionsTemp.EndPosition(not(ChosenIndexes))=[];
NumInsertionsTemp.EndPositionGlobal(not(ChosenIndexes))=[];
NumInsertionsTemp.InsertionLength(not(ChosenIndexes))=[];
NumInsertionsTemp.CountsPerMillion(not(ChosenIndexes))=[];
NumInsertionsTemp.RawCounts(not(ChosenIndexes))=[];
NumInsertionsTemp.Polarity(not(ChosenIndexes))=[];
NumInsertionsTemp.Phenotype(not(ChosenIndexes))=[];
NumInsertionsTemp.CloneName(not(ChosenIndexes))=[];

SampleLabelsPlotTemp=SampleLabelsPlot;
SampleLabelsPlotTemp(not(ChosenIndexes))=[];

GlobalInsertionsPlot_marker(NumInsertionsTemp,'Raw reads',SampleLabelsPlotTemp,1)

%% Create a locations matrix
% As this is unclustered data, we need to create a Matrix that contains all
% the possible positions in the genome with insertions in this dataset
% Use the global locations for this

% First list of locations
LocationsList=unique(cell2mat(ChromClusterData.StartPositionGlobal));

% Now lets reduce the list with a window
windowSize=2;
LocationsList2=unique(windowSize.*round(LocationsList/windowSize));
DataMatrix=zeros(length(NumInsertions.SampleName),length(LocationsList2));
for i=1:length(NumInsertions.SampleName)
    for j=1:length(NumInsertions.StartPositionGlobal{i})
        IndexTemp=find(abs(cell2mat(NumInsertions.StartPositionGlobal{i}(j))-LocationsList2)<windowSize);
        if cell2mat(NumInsertions.Polarity{i}{j})=='+'
            DataMatrix(i,IndexTemp)=1;
        else
            DataMatrix(i,IndexTemp)=-1;
        end
    end
end


%% Calculate clone correlation matrix

% For the samples
CorrMatAll=corr(DataMatrix.');
CorrMatAll2=corr(DataMatrix);

set(groot,'defaultAxesTickLabelInterpreter','none');  
figure
imagesc(CorrMatAll)
colormap('hot')
axis tight
set(gca,'XTick',1:size(CorrMatAll,1))
set(gca,'XTickLabel',SampleLabelsPlot)
set(gca,'XTickLabelRotation',90)
set(gca,'YTick',1:size(CorrMatAll,1))
set(gca,'YTickLabel',SampleLabelsPlot)
view(2)
h=colorbar;
ylabel(h,'Correlation')

figure
imagesc(CorrMatAll2)
colormap('hot')
axis tight
set(gca,'XTick',1:length(ClusterData))
set(gca,'XTickLabel',1:length(ClusterData))
set(gca,'XTickLabelRotation',90)
set(gca,'YTick',1:length(ClusterData))
set(gca,'YTickLabel',1:length(ClusterData))
xlabel('Cluster number')
ylabel('Cluster number')
view(2)
h=colorbar;
ylabel(h,'Correlation')

%% Build a phylogenetic tree with the correlation information
% In order to do this, I want to first calculate the distance between each
% of the clones

DistancesPhylo=pdist2(DataMatrix,DataMatrix,"squaredeuclidean");
TreeElement=seqlinkage(DistancesPhylo,'single',SampleLabelsPlot);
view(TreeElement)


%% Finding the gene!
% We need to identify the minimum number of common clusters between the
% clones of interest that can explain the phenotype of interest.

% We need to define a window size to cluster the insertions:

windowSize = 200000;

% Now lets find clusters in all the data, defined as locations with more
% than 1 insertion

clusterCount=0;
chromNames=unique(ChromClusterData.Chromosomes);
Cluster={};

clc
for i=1:length(chromNames)
    disp(['Contig number ' num2str(i) ' out of ' num2str(length(chromNames))])
    ChromIndex=find(ChromClusterData.Chromosomes==chromNames{i});
    if length(ChromIndex)>2
        count=1;
        Step=1;
        while count+Step<=length(ChromIndex)
            Dist=abs(ChromClusterData.StartPosition{ChromIndex(count)}-...
                ChromClusterData.StartPosition{ChromIndex(count+Step)});
%            disp(['Distance: ' num2str(Dist) ' Count: ' num2str(count) ' Step: ' num2str(Step)])
            disp(ChromClusterData.StartPosition{ChromIndex(count)})
            if Dist<windowSize
                Step=Step+1;
                disp(['Step value ' num2str(Step)])
            elseif Step==1
                count=count+Step;
                disp(['Resetting count to ' num2str(count)])
            else %if Step>1
                clusterCount=clusterCount+1;
                disp(['Cluster number ' num2str(clusterCount)])                
                ClusterTemp=struct;
                ClusterTemp.Chromosome=chromNames{i};
                ClusterTemp.StartPosition=cell2mat(ChromClusterData.StartPosition(ChromIndex(count:count+Step-1)));
                ClusterTemp.EndPosition=cell2mat(ChromClusterData.EndPosition(ChromIndex(count:count+Step-1)));
                ClusterTemp.SampleName=ChromClusterData.SampleName(ChromIndex(count:count+Step-1));
                ClusterTemp.CloneName=ChromClusterData.CloneName(ChromIndex(count:count+Step-1));
                ClusterTemp.RawCounts=cell2mat(ChromClusterData.RawCounts(ChromIndex(count:count+Step-1)));
                ClusterTemp.CountsPerMillion=cell2mat(ChromClusterData.CountsPerMillion(ChromIndex(count:count+Step-1)));
                ClusterTemp.Polarity=ChromClusterData.Polarity(ChromIndex(count:count+Step-1));
                ClusterTemp.Phenotype=ChromClusterData.Phenotype(ChromIndex(count:count+Step-1));
                ClusterTemp.Certainty=cell2mat(ChromClusterData.Certainty(ChromIndex(count:count+Step-1)));
                Cluster{clusterCount}=ClusterTemp;
                count=count+Step;
                Step=1;
            end
        end
        if Step>1 % in case the contig ends with a cluster
                clusterCount=clusterCount+1;
                disp(['Cluster number ' num2str(clusterCount)])                
                ClusterTemp=struct;
                ClusterTemp.Chromosome=chromNames{i};
                ClusterTemp.StartPosition=cell2mat(ChromClusterData.StartPosition(ChromIndex(count:count+Step-1)));
                ClusterTemp.EndPosition=cell2mat(ChromClusterData.EndPosition(ChromIndex(count:count+Step-1)));
                ClusterTemp.SampleName=ChromClusterData.SampleName(ChromIndex(count:count+Step-1));
                ClusterTemp.CloneName=ChromClusterData.CloneName(ChromIndex(count:count+Step-1));
                ClusterTemp.RawCounts=cell2mat(ChromClusterData.RawCounts(ChromIndex(count:count+Step-1)));
                ClusterTemp.CountsPerMillion=cell2mat(ChromClusterData.CountsPerMillion(ChromIndex(count:count+Step-1)));
                ClusterTemp.Polarity=ChromClusterData.Polarity(ChromIndex(count:count+Step-1));
                ClusterTemp.Phenotype=ChromClusterData.Phenotype(ChromIndex(count:count+Step-1));
                ClusterTemp.Certainty=cell2mat(ChromClusterData.Certainty(ChromIndex(count:count+Step-1)));
                Cluster{clusterCount}=ClusterTemp;
                count=count+Step;
                Step=1;
        end
    end
end
clear clusterCount
ClusterFields={'Chromosome','StartPosition','EndPosition','SampleName','CloneName','RawCounts','CountsPerMillion','Polarity','Phenotype','Certainty'};

clc
% Displaying the resulting clusters
for i=1:length(Cluster)
    disp(['Cluster number ' num2str(i) ' in contig ' Cluster{i}.Chromosome])
    for j=1:length(Cluster{i}.SampleName)
        disp(strjoin([Cluster{i}.SampleName{j},' ',Cluster{i}.Phenotype{j},num2str(Cluster{i}.StartPosition(j)),...
            num2str(Cluster{i}.EndPosition(j)),Cluster{i}.Polarity{j},num2str(Cluster{i}.RawCounts(j)),...
            num2str(Cluster{i}.CountsPerMillion(j))],'\t'))
    end
    disp([' '])
end

% Removing cluster that contain certain negatives

% for i=1:length(Cluster)
% 	hasNegatives=0;
% 	for j=1:length(Cluster{i}.Phenotype)
%         if and(Cluster{i}.Certainty(j)==1,strcmp('N',Cluster{i}.Phenotype(j)))
%             hasNegatives=hasNegatives+1;
%         end
%     end
%     if hasNegatives>0
%         Cluster{i}=[];
%         disp(['Removing cluster ' num2str(i)])
%     end
% end
% 
% Cluster=Cluster(~cellfun('isempty',Cluster));

% Filtering clusters
% 1) Considering polarity
% If a cluster has mixed polarity, it should be
% a) If size is 2, it should be removed
% b) If size if size > 2 and any of the polarities is single, it should be removed from the
% cluster
% c) If size >2 and the size of the polarity subclusters is >2, it should
% be split into two clusters.

% clc
% n=length(Cluster);
% for i=1:n
%     disp(['Cluster number ' num2str(i)])
%     NumForward=0;
%     NumReverse=0;
%     ForwardIndex=[];
%     ReverseIndex=[];
%     for j=1:length(Cluster{i}.Polarity)
%         if strcmp(Cluster{i}.Polarity{j},'+')
%             NumForward=NumForward+1;
%             ForwardIndex=[ForwardIndex j];
%         elseif strcmp(Cluster{i}.Polarity{j},'-')
%             NumReverse=NumReverse+1;
%             ReverseIndex=[ReverseIndex j];
%         end
%     end
%     disp(['Forward: ' num2str(NumForward) ' Reverse: ' num2str(NumReverse)])
%     if or(NumForward==0,NumReverse==0)
%     elseif length(Cluster{i}.Polarity)==2
%         if and(NumForward==1,NumReverse==1)
%             disp(['Removing cluster number ' num2str(i)])
%             Cluster{i}=[]; % Remove cluster
%         end
%     elseif NumForward==1
%         disp(['Removing forward insertion from cluster number ' num2str(i)])
%         for z=2:length(ClusterFields)
%             Cluster{i}.(ClusterFields{z})(ForwardIndex)=[]; % Remove the odd polarity
%         end
%     elseif NumReverse==1
%         disp(['Removing reverse insertion from cluster number ' num2str(i)])
%         for z=2:length(ClusterFields)
%             Cluster{i}.(ClusterFields{z})(ReverseIndex)=[]; % Remove the odd polarity
%         end
%     else % Remove the Reverse Polarity (for example) and put it into a new Cluster
%         disp(['Splitting cluster number ' num2str(i)])
%         m=length(Cluster);
%         ClusterTemp.Chromosome=Cluster{i}.Chromosome;
%         for z=2:length(ClusterFields)
%             ClusterTemp.(ClusterFields{z})=Cluster{i}.(ClusterFields{z})(ReverseIndex);
%         end
%         Cluster{m+1}=ClusterTemp;
%     end
% end

% Remove empty clusters
for i=1:length(Cluster)
    if isempty(Cluster{i})
        disp(['Removing Cluster ' num2str(i)])
        Cluster{i}=[];        
    elseif length(Cluster{i}.SampleName)<2
        disp(['Removing Cluster ' num2str(i)])
        Cluster{i}=[];
    end
end

Cluster=Cluster(~cellfun('isempty',Cluster));

%% Display number of insertions, number of clusters, CPM limits and raw reads limits
disp(['Raw reads threshold: ' num2str(rawReadsThreshold)])
disp(['CPM threshold: ' num2str(CPMThreshold)])
disp(['Number of insertions: ' num2str(sum(NumInsertions.Count))])
disp(['Number of clusters: ' num2str(length(Cluster))])

%% Give the clusters a correlation score

for i=1:length(Cluster)
    ClusterScore=0;
    Count=0;
    for j=1:length(Cluster{i}.SampleName)
        index1=find(strcmp(Cluster{i}.SampleName{j},SampleLabels(:,1))==1);
        for z=j+1:length(Cluster{i}.SampleName)
            index2=find(strcmp(Cluster{i}.SampleName{z},SampleLabels(:,1))==1);
            ClusterScore=ClusterScore+CorrMatAll(index1,index2);
            Count=Count+1;
        end
    end
    Cluster{i}.ClusterScore=ClusterScore/Count;
end

%% Ordering clusters by number of members
A=[];
for i=1:length(Cluster)
    A=[A length(Cluster{i}.SampleName)];
end
[~,B]=sort(A);
Cluster=Cluster(B(end:-1:1));


% %% Separating clusters with AITC+ members
% 
% AITCIndex=[];
% WS12Index=[];
% for i=1:length(Cluster)
%    if sum(strcmp('AITC+',Cluster{i}.Phenotype))>0
%        AITCIndex=[AITCIndex i];
%    end       
%    if sum(strcmp('WS12+',Cluster{i}.Phenotype))>0
%        WS12Index=[WS12Index i];
%    end
% end
% 
% %% Displaying AITC clusters
% disp('AITC+ Clusters')
% clc
% for i=AITCIndex
%     disp(['Cluster number ' num2str(i) ' in contig ' Cluster{i}.Chromosome '. CORRELATION SCORE = ' num2str(Cluster{i}.ClusterScore)])
%     for j=1:length(Cluster{i}.SampleName)
%         disp(strjoin([Cluster{i}.SampleName{j},' ',Cluster{i}.CloneName{j},' ',Cluster{i}.Phenotype{j},num2str(Cluster{i}.StartPosition(j)),...
%             num2str(Cluster{i}.EndPosition(j)),Cluster{i}.Polarity{j},num2str(Cluster{i}.RawCounts(j)),...
%             num2str(Cluster{i}.CountsPerMillion(j))],'\t'))
%     end
%     disp([' '])
% end
% 
% %% Displaying WS12 clusters
% disp('WS12+ Clusters')
% clc
% for i=WS12Index
%     disp(['Cluster number ' num2str(i) ' in contig ' Cluster{i}.Chromosome '. CORRELATION SCORE = ' num2str(Cluster{i}.ClusterScore)])
%     for j=1:length(Cluster{i}.SampleName)
%         disp(strjoin([Cluster{i}.SampleName{j},' ',Cluster{i}.CloneName{j},' ',Cluster{i}.Phenotype{j},num2str(Cluster{i}.StartPosition(j)),...
%             num2str(Cluster{i}.EndPosition(j)),Cluster{i}.Polarity{j},num2str(Cluster{i}.RawCounts(j)),...
%             num2str(Cluster{i}.CountsPerMillion(j))],'\t'))
%     end
%     disp([' '])
% end

%% Displaying the cluster results

clc
for i=1:length(Cluster)
    disp(['Cluster number ' num2str(i) ' in contig ' Cluster{i}.Chromosome '. CORRELATION SCORE = ' num2str(Cluster{i}.ClusterScore) '. N = ' num2str(length(Cluster{i}.SampleName))])
    for j=1:length(Cluster{i}.SampleName)
        disp(strjoin([Cluster{i}.SampleName{j},' ',Cluster{i}.CloneName{j},' ',Cluster{i}.Phenotype{j},num2str(Cluster{i}.StartPosition(j)),...
            num2str(Cluster{i}.EndPosition(j)),Cluster{i}.Polarity{j},num2str(Cluster{i}.RawCounts(j)),...
            num2str(Cluster{i}.CountsPerMillion(j))],'\t'))
    end
    disp([' '])
end

%% Plotting the insertions for a chosen cluster
ClusterPlot(Cluster,25,'Raw reads',SampleLabelsPlot,SampleLabels(:,4))

%% Finding all clusters for a certain group:

% WS12+ cells
%GroupIndex=find([SampleLabels{:,7}]==1);
% AITC+ cells
%GroupIndex=find([SampleLabels{:,7}]==2);
% P4F8 group members
%GroupIndex=find([SampleLabels{:,6}]==10);
% P2D11 group members
%GroupIndex=find([SampleLabels{:,6}]==20);
% P2B2 group members
%GroupIndex=find([SampleLabels{:,6}]==30);
% P3F5 group members
GroupIndex=find([SampleLabels{:,6}]==40);
% For P1C11 and P1C11 neg pheno
%GroupIndex=[7,8];

% Obtain list of clusters that contain at least "threshold" of the members
ClusterGroup=[];
Threshold=3;
for j=1:length(Cluster)
    count=0;
    for i=1:length(GroupIndex)
        count=count+sum(strcmp(SampleLabels(GroupIndex(i),4),Cluster{j}.CloneName));
    end
    if count>=Threshold
        ClusterGroup=[ClusterGroup j];
    end
end

% Plotting a few clusters together

ChosenClusters=ClusterGroup;
ClusterPlot(Cluster,ChosenClusters,'Raw reads',SampleLabelsPlot,SampleLabels(:,4))

% Displaying the insertion information

clc
for i=ChosenClusters
    disp(['Cluster number ' num2str(i) ' in contig ' Cluster{i}.Chromosome '. CORRELATION SCORE = ' num2str(Cluster{i}.ClusterScore) '. N = ' num2str(length(Cluster{i}.SampleName))])
    for j=1:length(Cluster{i}.SampleName)
        disp(strjoin([Cluster{i}.SampleName{j},' ',Cluster{i}.CloneName{j},' ',Cluster{i}.Phenotype{j},num2str(Cluster{i}.StartPosition(j)),...
            num2str(Cluster{i}.EndPosition(j)),Cluster{i}.Polarity{j},num2str(Cluster{i}.RawCounts(j)),...
            num2str(Cluster{i}.CountsPerMillion(j))],'\t'))
    end
    disp([' '])
end

%% Finding all insertions for a certain clone (and seeing who shares them!)

% P4F8
CloneName="P4F8";
CloneIndex=find(strcmp({SampleLabels{:,4}},CloneName));

% Obtain list of insertions that contain this clone
Chrom=NumInsertions.Chromosomes{NumInsertions.SampleName2==CloneName};
StartPos=cell2mat(NumInsertions.StartPosition{NumInsertions.SampleName2==CloneName});
StartPosGlobal=cell2mat(NumInsertions.StartPositionGlobal{NumInsertions.SampleName2==CloneName});
EndPos=cell2mat(NumInsertions.EndPosition{NumInsertions.SampleName2==CloneName});

% Obtaining the list of clones that also contain these insertions
ClusterClone={};
ClusterTemp=[];
CloneIndex=find(NumInsertions.SampleName2==CloneName);

for i=1:length(StartPosGlobal)
    ClusterTemp.Chromosome=Chrom(i);
    ClusterTemp.StartPosition=StartPos(i);
    ClusterTemp.EndPosition=EndPos(i);
    ClusterTemp.SampleName=NumInsertions.SampleName{CloneIndex};
    ClusterTemp.CloneName={CloneName};
    ClusterTemp.RawCounts=NumInsertions.RawCounts{CloneIndex}{i};
    ClusterTemp.CountsPerMillion=NumInsertions.CountsPerMillion{CloneIndex}{i};
    ClusterTemp.Polarity=NumInsertions.Polarity{CloneIndex}{i};
    IndexTemp=[];
    for j=1:length(NumInsertions.StartPositionGlobal)
        IndexTemp=find(abs(cell2mat(NumInsertions.StartPositionGlobal{j})-StartPosGlobal(i))<windowSize);
        if not(isempty(IndexTemp))
            for z=1:length(IndexTemp)
                ClusterTemp.StartPosition=[ClusterTemp.StartPosition NumInsertions.StartPosition{j}{IndexTemp(z)}];
                ClusterTemp.EndPosition=[ClusterTemp.EndPosition NumInsertions.EndPosition{j}{IndexTemp(z)}];
                ClusterTemp.SampleName=[ClusterTemp.SampleName [NumInsertions.SampleName{j}]];
                ClusterTemp.CloneName=[ClusterTemp.CloneName NumInsertions.SampleName2{j}];
                ClusterTemp.RawCounts=[ClusterTemp.RawCounts NumInsertions.RawCounts{j}{IndexTemp(z)}];
                ClusterTemp.CountsPerMillion=[ClusterTemp.CountsPerMillion NumInsertions.CountsPerMillion{j}{IndexTemp(z)}];
                ClusterTemp.Polarity=[ClusterTemp.Polarity NumInsertions.Polarity{j}{IndexTemp(z)}];
            end
        end
    end
    ClusterClone{i}=ClusterTemp;
    clear ClusterTemp
end

% Plotting together

ChosenClusters=1:length(ClusterClone);
ClusterPlot(ClusterClone,ChosenClusters,'Raw reads',SampleLabelsPlot,SampleLabels(:,4))

%% Obtaining the list of genes for the insertion of interest
addpath('C:\Users\phyjga\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis\CHO6-11LIB\Sequencing2analysis\03-merged2020')
Chromosome='NC_048600.1';
StartPos=   108344025;
EndPos=     108344149;
Polarity=   '+';

cd('C:\Users\phyjga\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis\CHO6-11LIB\Sequencing2analysis\03-merged2020\GCF_003668045.3_annotations\ncbi_dataset\data\GCF_003668045.3')
AnnotationsFile='genomic.gtf';
savingPath='C:\Users\phyjga\OneDrive - University of Leeds\Lukacs Lab\LibraryAnalysis\CHO6-11LIB\Sequencing2analysis\04-cluster2020\P3F5 group\';
GeneList=geneList(Chromosome,StartPos,EndPos,Polarity,AnnotationsFile,300000,[savingPath,Chromosome,'_',num2str(StartPos),'.csv']);

cd(savingPath)
