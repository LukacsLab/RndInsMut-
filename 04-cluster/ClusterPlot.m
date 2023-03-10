function ClusterPlot(Cluster,ChosenCluster,Counts,SampleLabelsPlot,SampleLabels)
    % ClusterPlot takes insertions Start and End coordinates for certain
    % Samples in a cluster and returns a plot of those coordinates. 
    % It color codes the raw counts for each insertion and shows the polarity of the
    % insertion.
    %
    % Eg: ClusterPlot(Cluster,ChosenCluster,Counts,SampleLabelsPlot)
    %
    % Input arguments:
    % Structure array organized by cluster number with following fields:
    % - Chromosome
    % - StartPosition
    % - EndPosition
    % - SampleName: list of the samples in the cluster (PCR)
    % - CloneName: list of the samples in the cluster (plate position)
    % - RawCounts: array of raw counts of each of the insertions for each sample.
    % - CountsPerMillion: array of CPMs of each of the insertions for each sample.
    % - Polarity: list with either '+' or '-' for the insertions in each sample.
    % - Phenotype: of each sample
    % - Certainty
    % - ClusterScore
    %
    % Argument 2: ChosenCluster
    % Array with the selected clusters to plot
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
    % Argument 5: SampleLabels
    % Cell array with all the clone names.
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
    hold on
    if CPM==1 % Currently assuming the maximum number of CPM or Raw reads
        colorIndexRef=1e6;
    end
    if RawReads==1
        colorIndexRef=1e5;
    end
    c=cool(ceil(log10(colorIndexRef)));
    count=1;
    polarity=1;
    
    if length(ChosenCluster)==1
        if CPM==1
            colorIndex=floor(log10(double(Cluster{ChosenCluster}.CountsPerMillion)));
        elseif RawReads==1
            colorIndex=round(log10(double(Cluster{ChosenCluster}.RawCounts)),0)+1;
        end
        for j=1:length(Cluster{ChosenCluster}.StartPosition)
            Height=find(strcmp(Cluster{ChosenCluster}.CloneName{j},SampleLabels),1); % Assuming one insertion per sample in the area
            h=rectangle('Position',[Cluster{ChosenCluster}.StartPosition(j),Height-1,Cluster{ChosenCluster}.EndPosition(j)-Cluster{ChosenCluster}.StartPosition(j),1],...
                'EdgeColor',c(colorIndex(j),:),'FaceColor',c(colorIndex(j),:));
            if polarity==1
                if strcmp(Cluster{ChosenCluster}.Polarity{j},'+')
                    scatter(Cluster{ChosenCluster}.EndPosition(j),Height-0.5,'Marker','>','Color',c(colorIndex(j),:),'MarkerFaceColor',c(colorIndex(j),:),'MarkerEdgeColor',c(colorIndex(j),:),'LineWidth',1)
                else
                    scatter(Cluster{ChosenCluster}.StartPosition(j),Height-0.5,'Marker','<','Color',c(colorIndex(j),:),'MarkerFaceColor',c(colorIndex(j),:),'MarkerEdgeColor',c(colorIndex(j),:),'LineWidth',1)
                end
            end
            count=count+1;
            title([Cluster{ChosenCluster}.Chromosome],'Interpreter','none')
            Start=min(Cluster{ChosenCluster}.StartPosition)-0.5e2;
            End=max(Cluster{ChosenCluster}.EndPosition)+0.5e2;
            xlim([Start End])
            ylim([0,length(SampleLabelsPlot)])
            xl=xlim();
            xlabel('Coordinate/ bp')
            set(gca,'TickDir','out');            
            for z=1:length(SampleLabelsPlot)
                line([xl(1) xl(2)],[z z],'Color',[0.95 0.95 0.95],'LineStyle','-')
            end 
            box on
            h=gca; 
            h.LineWidth=1;
            h.FontSize=9;
            h.YAxis.TickLength = [0 0];
            h.YAxis.TickValues=[0.5:1:length(SampleLabelsPlot)-0.5];
            h.YAxis.TickLabels=SampleLabelsPlot;
            h.YAxis.TickLabelInterpreter='none';
            colormap(c)
        end
    elseif length(ChosenCluster)>1
        for ClustInd=1:length(ChosenCluster)
            if CPM==1
                colorIndex=floor(log10(double(Cluster{ChosenCluster(ClustInd)}.CountsPerMillion)));
            elseif RawReads==1
                colorIndex=round(log10(double(Cluster{ChosenCluster(ClustInd)}.RawCounts)),0)+1;
            end
            if ClustInd==ChosenCluster(end)
            subplot('Position',...
                [0.15+0.8/length(ChosenCluster)*(ClustInd-1),... % Left
                0.1,... % bottom
                1/length(ChosenCluster),... % width
                0.8]) % height                
            else
            subplot('Position',...
                [0.15+0.8/length(ChosenCluster)*(ClustInd-1),... % Left
                0.1,... % bottom
                0.8/length(ChosenCluster),... % width
                0.8]) % height
            end
            hold on
            for j=1:length(Cluster{ChosenCluster(ClustInd)}.StartPosition)
                Height=find(strcmp(Cluster{ChosenCluster(ClustInd)}.CloneName{j},SampleLabels),1); % Assuming one insertion per sample in the area
                h=rectangle('Position',[Cluster{ChosenCluster(ClustInd)}.StartPosition(j),Height-1,...
                    Cluster{ChosenCluster(ClustInd)}.EndPosition(j)-Cluster{ChosenCluster(ClustInd)}.StartPosition(j),1],...
                    'EdgeColor',c(colorIndex(j),:),'FaceColor',c(colorIndex(j),:));
                if polarity==1
                    if strcmp(Cluster{ChosenCluster(ClustInd)}.Polarity{j},'+')
                        scatter(Cluster{ChosenCluster(ClustInd)}.EndPosition(j),...
                            Height-0.5,'Marker','>','Color',c(colorIndex(j),:),...
                            'MarkerFaceColor',c(colorIndex(j),:),'MarkerEdgeColor',c(colorIndex(j),:),'LineWidth',1)
                    else
                        scatter(Cluster{ChosenCluster(ClustInd)}.StartPosition(j),...
                            Height-0.5,'Marker','<','Color',c(colorIndex(j),:),...
                            'MarkerFaceColor',c(colorIndex(j),:),'MarkerEdgeColor',c(colorIndex(j),:),'LineWidth',1)
                    end
                end
                count=count+1;
            end
            title([Cluster{ChosenCluster(ClustInd)}.Chromosome],'Interpreter','none')
            Start=min(Cluster{ChosenCluster(ClustInd)}.StartPosition);
            End=max(Cluster{ChosenCluster(ClustInd)}.EndPosition);
            xlim([Start-0.1*(End-Start) End+0.1*(End-Start)])
            ylim([0,length(SampleLabelsPlot)])
            xl=xlim();
            set(gca,'TickDir','out');
            for z=1:length(SampleLabelsPlot)
                line([xl(1) xl(2)],[z z],'Color',[0.95 0.95 0.95],'LineStyle','-')
            end 
            box on
            h=gca; 
            h.LineWidth=1;
            h.FontSize=9;
            colormap(c)
            h.YAxis.TickLabels=[];
            if ClustInd==1
                h.YAxis.TickLength = [0 0];      
                h.YAxis.TickValues=[0.5:1:length(SampleLabels(:,1))-0.5];
                h.YAxis.TickLabels=SampleLabelsPlot;
                h.YAxis.TickLabelInterpreter='none';
            end
            if ClustInd==floor(length(ChosenCluster)/2)
                xlabel('Coordinate/ bp')
            end
        end
    end
    cb=colorbar;
    cb.Ticks=[0:1/(ceil(log10(colorIndexRef))):1]; % If we are plotting from 1 to 10^N, we will have N+1 ticks in the plot.
    cb.TickLabels=num2cell(10.^[0:1:(ceil(log10(colorIndexRef)))]);
    cb.LineWidth=1.5;
    if CPM==1
        ylabel(cb,'Counts per million'); % For CPM color code
    elseif RawReads==1
        ylabel(cb,'Raw counts'); % For Raw counts color code
    end
    return