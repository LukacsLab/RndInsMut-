function GlobalInsertionsPlot_rectangle(NumInsertions,Counts,SampleLabelsPlot,Genes)
    % GlobalInsertionsPlot takes insertions Start and End coordinates for certain
    % Samples and returns a plot of those coordinates. It color codes the
    % raw counts for each insertion and shows the polarity of the
    % insertion.
    %
    % Eg: GlobalInsertionsPlot(NumInsertions,'RawReads',SampleLabelsPlot,1)
    %
    % Input arguments:
    % Structure array organized by sample name with following fields:
    % - SampleName: list of the name for the sample organization of the
    % data.
    % - RawCounts: array of cells, where each cell contains the raw counts 
    % of each of the insertions for each sample.
    % - CPM: array of cells, where each cell contains the CPM of each of 
    % the insertions for each sample.
    % - StartPositionGlobal: array of cells, where each cell contains the 
    % list with the start position for each of the insertions for each 
    % sample. If using multiple chromosomes, a global coordinate should
    % be computed beforehand to be able to add all insertions in the same
    % plot.
    % - EndPositionGlobal: array of cells, where each cell contains the 
    % list with the start position for each of the insertions for each 
    % sample.
    % - Polarity: array of cells where each cell contains a list with 
    % either '+' or '-' for each of the insertions of each sample.
    % - chromEndsGlobal list of the ends of the global positions calculated
    % for each of the chromosomes when plotting multiple chromosomes. If
    % plotting only one chromosome, leave empty.
    % - Genes: only needed if choosing Genes to 1. Containing all genes to 
    % be plotted with structure fields:
    % * Name
    % * Chromosome
    % * StartPosition
    % * EndPosition
    % * StartPositionGlobal
    % * EndPositionGlobal
    % 
    % Argument 2: Counts
    % String with value 'CPM' or 'Raw reads'
    %
    % Argument 3: SampleLabelsPlot
    % Cell with a list of strings that will label the samples on the Y-axis.
    %
    % Argument 4: Genes
    % 1 for yes, 0 for no, it adds TRPM8 and TRPA1 genes plot lines
    %
    if strcmp(Counts,'CPM')
        disp('Plotting counts per million')
        CPM=1; RawReads=0;
    elseif strcmp(Counts,'Raw reads')
        disp('Plotting raw reads')
        CPM=0; RawReads=1;
    end

    count=1;
    figure('Name','Summary plot','Position',[100,150,800,400],'Units','pixels')
    if CPM==1 % Currently assuming the maximum number of CPM or Raw reads
        colorIndexRef=1e7;
    elseif RawReads==1
        colorIndexRef=1e5;
    end
    c=cool(ceil(log10(colorIndexRef)));
%     c=bone(ceil(log10(colorIndexRef)));
%     c = flipud(c);
    count=1;
    polarity=1;

    if Genes==1
        GeneList=fieldnames(NumInsertions.Genes);
        for i=1:length(GeneList)
            GeneTemp=getfield(NumInsertions.Genes,GeneList{i});
            rectangle('Position',[GeneTemp.StartPositionGlobal,...
                0,GeneTemp.EndPositionGlobal-...
                GeneTemp.StartPositionGlobal,length(SampleLabelsPlot)],...
               'EdgeColor','k','FaceColor','k','LineWidth',1);
            text(GeneTemp.StartPositionGlobal,length(SampleLabelsPlot)+1,GeneTemp.Name,'Rotation',90)
        end
    end

    for z=1:length(NumInsertions.SampleName)
        StartingIndexes=cell2mat(NumInsertions.StartPositionGlobal{z});
        EndingIndexes=cell2mat(NumInsertions.EndPositionGlobal{z});
        if CPM==1
            colorIndex=floor(log10(double(cell2mat(NumInsertions.CountsPerMillion{z}))))+1; % For CPM color coding 
        elseif RawReads==1
            colorIndex=round(log10(double(cell2mat(NumInsertions.RawCounts{z}))),0)+1; % For Raw reads color coding / +1 to avoid the 0 in index
        end
        for j=1:length(StartingIndexes)
            h=rectangle('Position',[StartingIndexes(j),count-1,EndingIndexes(j)-StartingIndexes(j),1],...
               'EdgeColor',c(colorIndex(j),:),'FaceColor',c(colorIndex(j),:),'LineWidth',2);
            hold on
%             if polarity==1
%                 if cell2mat(NumInsertions.Polarity{z}{j})=='+'
%                     scatter(EndingIndexes(j),count-0.5,'Marker','>','Color',c(colorIndex(j),:),'MarkerFaceColor',c(colorIndex(j),:),'MarkerEdgeColor',c(colorIndex(j),:),'LineWidth',1)
%                 else
%                     scatter(EndingIndexes(j),count-0.5,'Marker','<','Color',c(colorIndex(j),:),'MarkerFaceColor',c(colorIndex(j),:),'MarkerEdgeColor',c(colorIndex(j),:),'LineWidth',1)
%                 end
%             end
        end
        count=count+1;
    end
    try
        for i=1:length(NumInsertions.chromEndsGlobal)
            line([NumInsertions.chromEndsGlobal(i) NumInsertions.chromEndsGlobal(i)],[0 z],'Color',[0.8 0.8 0.8],'LineStyle','-')
        end
    catch
        disp('Single Chromosome plot')
    end
    count=count-1;
    xlim([0 NumInsertions.chromEndsGlobal(end)])
    ylim([0,count])
    xl=xlim();
    for z=1:length(NumInsertions.SampleName)
        line([xl(1) xl(2)],[z z],'Color',[0.95 0.95 0.95],'LineStyle','-')
    end
    box on
    h=gca; 
    h.LineWidth=1;
    h.FontSize=9;
    h.YAxis.TickLength = [0 0];
    h.YAxis.TickValues=[0.5:1:count-0.5];
    h.YAxis.TickLabels=SampleLabelsPlot;
    h.YAxis.TickLabelInterpreter='none';
    colormap(c)
    cb=colorbar;
    cb.Ticks=[0:1/(ceil(log10(colorIndexRef))):1]; % If we are plotting from 1 to 10^N, we will have N+1 ticks in the plot.
    cb.TickLabels=num2cell(10.^[0:1:(ceil(log10(colorIndexRef)))]);
    cb.LineWidth=1.5;
    if CPM==1
        ylabel(cb,'Counts per million'); % For CPM color code
    elseif RawReads==1
        ylabel(cb,'Raw counts'); % For Raw counts color code
    end
    xlabel('Global coordinate/ bp')
    set(gca,'TickDir','out');
    return