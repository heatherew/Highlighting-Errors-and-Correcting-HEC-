function [allEMG, allIMU, allSinusoids] = plotSegmentedData(allEMG, ...
    allIMU, allSinusoids, movements, positions, participantID)
% This function plots the segmented EMG data (and sinusoids when
% applicable) in a user input figure window, and calls the fixdata function
% if the user identifies any problems.
% Inputs:
%   - allEMG: nested cell array containing all EMG data. Each cell is a
%             different "trial" of data collection. Within each trial, each
%             cell is a different muscle contraction type
%   - allIMU: nested cell array (like allEMG) of accelerometer data
%   - allSinusoids: nested cell array (like allEMG) of the sinuoids when
%                   applicable. When not applicable, cells associated with 
%                   each muscle contraction type will be empty.
%   - movements: nested cell array, where cells associated with each muscle
%                contraction indicate the type of muscle contraction.
%   - positions: nested cell array (like movements) with limb positions.
%   - participantID: string with participant ID, for figure title.
% Outputs: same as inputs, but updated if problems were corrected


% create the figure window
fig = uifigure('HandleVisibility', 'on');
set(fig,'Units','normalized')
set(fig,'Position',[0 0 1 1])

colorArray = [0 0.5 0; 0.85 0.325 0.098; 0.929 0.694 0.125; ...
    0.494 0.184 0.556; 0.466 0.674 0.188; 0.301 0.745 0.933; ...
    0.635 0.078 0.184; 0 0.447 0.741];

% loop through each trial
for currentTrial = 1:numel(allEMG)
    undo = true;

    trialEMG = allEMG{currentTrial};
    trialIMU = allIMU{currentTrial};
    trialSinusoids = allSinusoids{currentTrial};
    trialMovements = movements{currentTrial};
    trialPositions = positions{currentTrial};
    
    while undo
        % loop through muscle contractions and plot
        fig.Name = ['HEC: ' participantID ' Trial ' num2str(currentTrial)];
        tiledlayout(fig,numel(trialPositions),numel(trialMovements))
        for currentPosition = find(any(~cellfun(@isempty,trialEMG),2))'
            for currentMovement = 1:numel(trialMovements)
                axes_all(currentPosition,currentMovement) = nexttile;
                hold on
                plot(trialEMG{currentPosition,currentMovement}); 
                colororder(colorArray)
                ylim([-1.01 1.01])
                xlim([1 length(trialEMG{currentPosition,currentMovement})])

                % plot sinusoids if necessary
                if contains(trialMovements{currentMovement},{'FE','OC','PS','UD'})
                    plot(trialSinusoids{currentPosition, ...
                        currentMovement},'k','LineWidth',2)
                end

                % plot position and movement labels
                if currentMovement == 1
                    ylabel(trialPositions{currentPosition-sum(any(cellfun(@isempty,trialEMG),2))}, ...
                        FontWeight='bold')
                end
                if currentPosition-sum(any(cellfun(@isempty,trialEMG),2)) == numel(trialPositions)
                    xlabel(trialMovements{currentMovement}, ...
                        FontWeight='bold')
                end
            end % currMovement
        end % currPosition

        buttonGroup = uibuttongroup(fig);
        set(buttonGroup,'Units','normalized')
        set(buttonGroup,'Position',[0.25 0.94 0.5 0.04])
        set(buttonGroup,'Units','pixels');

        toggleButton_1 = uitogglebutton(buttonGroup,Text='all good', ...
            Position=[1 1 buttonGroup.Position(3)/3 ...
            buttonGroup.Position(4)-1],Value=0);
        toggleButton_0 = uitogglebutton(buttonGroup, ...
            Text='(possibly) some problems', ...
            Position=[1+buttonGroup.Position(3)/3 1 ...
            buttonGroup.Position(3)/3-1 ...
            buttonGroup.Position(4)-1],Value=0);
        toggleButton_continue = uibutton(buttonGroup,'state', ...
            Text='continue',Position=[buttonGroup.Position(3)/3*2 1 ...
            buttonGroup.Position(3)/3-1 buttonGroup.Position(4)-1]);
        
        waitfor(toggleButton_continue,'Value')

        if toggleButton_1.Value
            clf
            undo = false;
            currEMG_new = trialEMG; 
            currSinusoids_new = trialSinusoids; 
            currIMU_new = trialIMU;
        elseif toggleButton_0.Value
            delete(buttonGroup)
            [fig, currEMG_new, currSinusoids_new, currIMU_new,undo] = ...
                fixdata(fig,axes_all,trialEMG,trialSinusoids,trialIMU);
        end
    end
    allEMG{currentTrial} = currEMG_new; 
    allSinusoids{currentTrial} = currSinusoids_new; 
    allIMU{currentTrial} = currIMU_new;

end % currentTrial

close all