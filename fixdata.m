function [fig, EMG_all, sinusoid_all, IMU_all, undo] = ...
    fixdata(fig,axes_all,EMG_all,sinusoid_all,IMU_all)
% Hightlight and Correct Errors
% This code walks the user through the process of correcting sinusoids and
% cutting out erroneous data.
% Inputs:
%   - fig:          handle of user interface figure (matlab uifigure)
%   - axes_all:     array of handles of all existing axes in fig,
%   - EMG_all:      cell array of EMG data corresponding with each axis
%   - sinusoid_all: cell array of sinusoids corresponding with each axis.
%                   If static training data, just input as [].
%   - IMU_all:      cell array of accelerometer data corresponding with
%                   each axis.

% Initialize some settings
cutDataXAxisLimit = 3000;
fontSize = 16;
colorArray = [0 0.5 0; 0.85 0.325 0.098; 0.929 0.694 0.125; ...
    0.494 0.184 0.556; 0.466 0.674 0.188; 0.301 0.745 0.933; ...
    0.635 0.078 0.184; 0 0.447 0.741];
buttonWidth = 100;
buttonHeight = 22;
undo = false;

% If sinusoid_all == [], then remake it as an array of empty cells
if isempty(sinusoid_all)
    sinusoid_all = repmat({[]},size(axes_all,1),size(axes_all,2));
end

% plot checkboxes for each axis
for currentRow = 1:size(axes_all,1)
    for currentColumn = 1:size(axes_all,2)
        set(axes_all(currentRow,currentColumn),'Units','pixels')

        % Position: [left bottom width height]
        checkBoxes(currentRow,currentColumn) = uicheckbox(fig, ...
            'Text','Problem?','Position', ...
            [axes_all(currentRow,currentColumn).Position(1) + ...
            axes_all(currentRow,currentColumn).Position(3)/2.5 ...
            axes_all(currentRow,currentColumn).Position(2) + ...
            axes_all(currentRow,currentColumn).Position(4)*1.01  ...
            84 22]);
    end
end

% determine size of computer screen
set(fig,'Units','Pixels');
screenSize = fig.Position(3:4);

% plot continue button
continueButton = createStateButton(fig, 'Continue', ...
    screenSize(1)/2-buttonWidth/2, axes_all(end).Position(2)/3, ...
    buttonWidth, buttonHeight);

% wait for continue button to be clicked
waitfor(continueButton,'Value')

% pull out checkbox values, & delete the continue button and checkboxes
checkBoxValues = arrayfun(@(x) x.Value, checkBoxes);
delete([continueButton reshape(checkBoxes,1,[])])

% plot an undo button
undoButton = createStateButton(fig, 'Undo – Go Back', ...
    10, 10, buttonWidth, buttonHeight);

if any(checkBoxValues,'all') % if none, will go back w/ no data change
    delete(axes_all)

    % set up indices to represent each axis (for mapping checkBoxValues)
    axisIndices = reshape(1:numel(checkBoxValues), ...
        [size(checkBoxValues,1), size(checkBoxValues,2)]);

    for currPlot = find(checkBoxValues)'
        [x,y] = find(axisIndices==currPlot);
        fig.Name = ['Position ' num2str(x) ' Movement ' num2str(y)];
        % fix sinusoid peaks
        if isempty(sinusoid_all{currPlot})
            cutDataValue = true;

%% Fix Sinusoids!
        else % aka if there's a sinusoid here
            [continueButton,cutDataButton,updateButton, ...
                flipSinusoidButton] = createStateButton(fig, ...
                {'Continue', 'Cut Some Data', 'Update Peaks', ...
                'Flip Sinusoid'}, screenSize(1)*0.7778, ...
                [150 200 250 300], buttonWidth, buttonHeight);
            
            tiledlayout(fig,3,3)

            % plot EMG and sinusoid
            axis1 = nexttile([1 2]);
            hold on
            plot(EMG_all{currPlot})
            colororder(colorArray)
            sinusoid_axis1=plotSinusoid(axis1,sinusoid_all{currPlot},0);
            ylim([-1.01 1.01])
            xlim([1 length(EMG_all{currPlot})])
            title('EMG')

            % set up some space for adjusting things
            textArea = nexttile([3 1]);

            % smooth EMG, then plot with sinusoid (rescaled to [0 1])
            EMG_smoothed = rescale(smoothdata(abs(EMG_all{currPlot}),...
                1,'gaussian',300));
            axis2 = nexttile([1 2]); 
            hold on
            plot(EMG_smoothed)
            colororder(colorArray)
            sinusoid_axis2=plotSinusoid(axis2,sinusoid_all{currPlot},1);
            ylim([-0.005 1.005])
            xlim([1 length(EMG_all{currPlot})])
            title('Smoothed EMG')

            % plot accelerometer and sinusoid
            axis3 = nexttile([1 2]); 
            hold on
            plot(IMU_all{currPlot}) % default colours for accelerometer
            sinusoid_axis3=plotSinusoid(axis3,sinusoid_all{currPlot},0);
            ylim([-1 1])
            xlim([1 length(EMG_all{currPlot})])
            title('Accelerometer')

            % go back to text area and set it up
            % the current peaks, boxes to input new peaks, a button to
            % refresh the sinusoids and the current peaks, and a button
            % to confirm that it's all good
            axes(textArea)
            hold on
            removeTicks(gca)
            plotBigText(textArea,0.1,0.9,'Peaks:',1,[]);
            [~,peaks_maxima] = findpeaks(sinusoid_all{currPlot});
            [~,peaks_minima] = findpeaks(-sinusoid_all{currPlot});
            peaks = sort([peaks_maxima'; peaks_minima']);
            for currPeak = 1:10
                axes(textArea); hold on
                plotBigText(textArea,0.15,0.9-0.05*currPeak, ...
                    [num2str(currPeak) ':'],0,'right');
                peakTextObjects(currPeak) = plotBigText(textArea,0.2, ...
                    0.9-0.05*currPeak,num2str(peaks(currPeak)),0,[]);
                textInput(currPeak) = createTextInput(fig, ...
                    screenSize(1)*0.7778, ...
                    (0.8264*(0.9-0.05*currPeak)+0.0885)*screenSize(2), ...
                    screenSize(1)*0.0521, screenSize(2)*0.0295);

                peakLabel(:,currPeak) = labelPeaks( ...
                    [axis1, axis2, axis3], peaks(currPeak), ...
                    sinusoid_all{currPlot}(peaks(currPeak)), ...
                    currPeak, [0 1 0]);
            end

            while ~continueButton.Value %until the confirm peaks is pressed
                while ~updateButton.Value && ~continueButton.Value
                    pause(1)
                end

                % delete sinusoids and current peaks
                delete([sinusoid_axis1 sinusoid_axis2 sinusoid_axis3 ...
                    peakTextObjects reshape(peakLabel,1,[])])

                if exist('errorText','var') % if i've posted an error msg
                    delete(errorText); clear errorText
                end

                % calculate new sinusoid
                textInputValues = cellfun(@str2double, ...
                    arrayfun(@(x) x.Value, textInput));
                oldpeaks = peaks;
                peaks(~isnan(textInputValues)) = ...
                    textInputValues(~isnan(textInputValues));
                if all(peaks == sort(peaks)) %if the peak order makes sense
                    startFrame = peaks(1) - round(mean(diff(peaks)))/2;
                    endFrame = peaks(end) + round(mean(diff(peaks)))/2;
                    newPeaks = [startFrame peaks' endFrame];
                    newSinusoid = interp1(newPeaks', ...
                        [0 -1 1 -1 1 -1 1 -1 1 -1 1 0], ...
                        1:length(EMG_all{currPlot}),'linear','extrap');
                    newSinusoid = sin(pi()/2*newSinusoid)';

                    if flipSinusoidButton.Value
                        newSinusoid = -newSinusoid; 
                    end

                    % update sinusoid_all with new sinusoid
                    sinusoid_all{currPlot} = newSinusoid;

                else % else -- don't change data + display an error message
                    axes(textArea); hold on; 
                    errorText = plotBigText(textArea,0.5,0.95, ...
                        'Those numbers do not work, try again!',0, ...
                        'center');
                    peaks = oldpeaks;
                end
                % update the peaks onscreen
                for currPeak = 1:10
                    axes(textArea); hold on
                    peakTextObjects(currPeak) = plotBigText(textArea, ...
                        0.2,0.9-0.05*currPeak,num2str(peaks(currPeak)), ...
                        0,[]);

                    peakLabel(:,currPeak) = labelPeaks( ...
                        [axis1, axis2, axis3], peaks(currPeak), ...
                        sinusoid_all{currPlot}(peaks(currPeak)), ...
                        currPeak, [0 1 0]);
                end

                % clear data tips
                delete([findobj(axis1,'Type','datatip') ...
                    findobj(axis2,'Type','datatip') ...
                    findobj(axis3,'Type','datatip')])

                % plot new sinusoids
                axes(axis1); hold on; 
                sinusoid_axis1=plotSinusoid(axis1,sinusoid_all{currPlot},0);
                
                axes(axis2); hold on; 
                sinusoid_axis2=plotSinusoid(axis2,sinusoid_all{currPlot},1);
                
                axes(axis3); hold on; 
                sinusoid_axis3=plotSinusoid(axis3,sinusoid_all{currPlot},0);

                % reset update button
                updateButton.Value = false;

                % get cut button value
                cutDataValue = cutDataButton.Value;

                % clear the textInput boxes
                for currPeak = 1:10
                    textInput(currPeak).Value = {''};
                end

                if undoButton.Value % if the undo button was pressed
                    undo = true;
                    delete([continueButton updateButton textInput ...
                        flipSinusoidButton cutDataButton undoButton])
                    return % leave this whole function and start over
                end
            end

            delete([continueButton updateButton textInput ...
                flipSinusoidButton cutDataButton])
            clear textInput
        end

%% Cut out Bad Data!
        if isempty(sinusoid_all{currPlot}) || cutDataValue 
            tiledlayout(fig,2,1)

            % setup tile for EMG (and sinusoid if it exists)
            axis1 = nexttile; 

            % setup tile for buttons and input
            axis2 = nexttile; hold on
            removeTicks(gca)

            % 3 possible cuts: cut off start, cut off end, cut out a chunk
            % buttons first:
            [continueButton,updateButton] = createStateButton(fig, ...
                {'Continue', 'Update Data'}, ...
                [screenSize(1)/2+10 screenSize(1)/2-10-100], ...
                screenSize(2)/7, buttonWidth, buttonHeight);

            % text boxes now:
            textInput = createTextInput(fig, ...
                [screenSize(1)/3.85 ...
                screenSize(1)/2-screenSize(1)*0.0521/2 ...
                screenSize(1)-screenSize(1)/3.2 ...
                screenSize(1)-screenSize(1)/3.2], ...
                [screenSize(2)/7+150 screenSize(2)/7+150 ...
                screenSize(2)/7+150 screenSize(2)/7+100], ...
                screenSize(1)*0.0521, screenSize(2)*0.0295);
            plotBigText(axis2,[0.25 0.5 0.75],0.75,{'New Start:', ...
                'New End:','Remove Section:'},1,'center');
            plotBigText(axis2,0.7,[0.65 0.475],{'From:','To:'},0,'right');

            if isempty(sinusoid_all{currPlot})
                currentColorArray = colorArray;
            else
                currentColorArray = [colorArray; [0 0 0]];
            end

            while ~continueButton.Value
                axes(axis1); hold on
                EMG_axis1 = plot(EMG_all{currPlot});
                colororder(currentColorArray)
                try
                    sinusoid_axis1 = plotSinusoid(axis1, ...
                        sinusoid_all{currPlot},0);
                catch
                    sinusoid_axis1 = [];
                end
                ylim([-1.01 1.01])
                xlim([1 ...
                    max([length(EMG_all{currPlot}) cutDataXAxisLimit])])

                while ~updateButton.Value
                    if ~continueButton.Value
                        pause(0)
                    else 
                        break % break out of this while loop
                    end
                end

                % update EMG, IMU, and sinusoid
                toCutValues = cellfun(@str2double,arrayfun(@(x) x.Value,...
                    textInput));
                try
                    [EMG_all{currPlot},IMU_all{currPlot}, ...
                        sinusoid_all{currPlot}] = cutData(toCutValues, ...
                        EMG_all{currPlot},IMU_all{currPlot}, ...
                        sinusoid_all{currPlot});

                    % get rid of error text if it's on the screen
                    if any(~isnan(toCutValues)) && exist('errorText','var')
                        delete(errorText)
                        clear errorText
                    end
                catch % if anything went wrong above^
                    axes(axis1); 
                    errorText = plotBigText(axis1, ...
                        max([length(EMG_all{currPlot}) ...
                        cutDataXAxisLimit])/2,0.5, ...
                        'Those numbers do not work, try again!',0, ...
                        'center');
                end

                % error message if from/to catch is in wrong order
                if ~exist('errorText','var') && ...
                    (toCutValues(3) > toCutValues(4))
                    axes(axis1); 
                    errorText = plotBigText(axis1, ...
                        max([length(EMG_all{currPlot}) ...
                        cutDataXAxisLimit])/2,0.5, ...
                        'Those numbers do not work, try again!',0, ...
                        'center');
                end

                delete(findobj(axis1,'Type','datatip'))

                % now delete plots (will update when looping around)
                % delete(e1); delete(s1)

                updateButton.Value = false;
                for currentTextBox = 1:4
                    textInput(currentTextBox).Value = {''};
                end
                delete([EMG_axis1' sinusoid_axis1])

                if undoButton.Value % if the undo button was pressed
                    undo = true;
                    % get rid of buttons and text boxes
                    delete([textInput continueButton ...
                        updateButton undoButton])
                    return
                end
            end

            delete([textInput continueButton updateButton])
        end
    end

    % plot all data and say "Updated" for axes where data was corrected
    tiledlayout(fig,size(axes_all,1),size(axes_all,2))
    for currentPosition = 1:size(axes_all,1)
        for currentMovement = 1:size(axes_all,2)
            axes_all(currentPosition,currentMovement) = nexttile; hold on

            plot(EMG_all{currentPosition,currentMovement})
            colororder(colorArray)
            ylim([-1.01 1.01])
            xlim([1 length(EMG_all{currentPosition,currentMovement})])
            plotSinusoid(gca,sinusoid_all{currentPosition, ...
                currentMovement},0);
            if checkBoxValues(currentPosition,currentMovement)
                title('Updated','FontSize',fontSize)
            end
        end % currentMovement
    end % currentPosition

    set(axes_all(end),'Units','pixels')
    continueButton = createStateButton(fig, 'Continue', ...
                screenSize(1)/2-buttonWidth/2, ...
                axes_all(end).Position(2)/2, buttonWidth, buttonHeight);

    waitfor(continueButton,'Value')
    delete(continueButton)
    if undoButton.Value % if the undo button was pressed
        undo = true;
        delete(undoButton)
        return % leave this whole function and start over
    end
end
delete(undoButton)

end % main function end

%% Helper functions nested within fixdata:

function handle = plotSinusoid(targetAxis,sinusoid,rescaleNeeded)
% plots sinusoid on in the target axis
% also rescales sinusoid to [0 1] if needed
axes(targetAxis)
if rescaleNeeded
    sinusoid = rescale(sinusoid);
end
handle = plot(sinusoid, 'k','LineWidth',2);
end


function varargout = createStateButton(figureHandle,buttonText, ...
    left,bottom,width,height)
% creates state button(s)
% text, left, and bottom can have multiple elements or just one
%   if multiple elements, then will loop through and make
%   multiple buttons!
% currently assuming same figureHandle and same width and height

% if only one button, make buttonText a cell
if ~iscell(buttonText)
    buttonText = {buttonText};
end

% if inputted just one left but multiple bottoms, repeat left
if length(left) < length(bottom) && length(left) == 1
    left = repmat(left,size(bottom));
end

% if inputted just one bottom but multiple lefts, repeat bottom
if length(bottom) < length(left) && length(bottom) == 1
    bottom = repmat(bottom,size(left));
end

for currentButton = 1:length(left)
    varargout{currentButton} = uibutton(figureHandle,'state', ...
        Text=buttonText{currentButton}, ...
        Position=[left(currentButton) bottom(currentButton) ...
        width height]);
end
end

function handle = createTextInput(figureHandle, ...
    left,bottom,width,height)
% creates text input area
for currentBox = 1:length(left)
    handle(currentBox) = uitextarea(figureHandle, ...
        Position=[left(currentBox) bottom(currentBox) ...
        width height]);
end
end

function removeTicks(targetAxis)
set(targetAxis,'xtick',[],'xticklabel',[],'ytick',[], ...
    'yticklabel',[])
end

function handle = plotBigText(targetAxis,x,y,textToPlot,bold,alignment)
% this will plot text on a target axis with a big font and makes it
% easier to bold and center text

fontSize = 16;
if bold
    fontWeight = 'bold';
else
    fontWeight = 'normal';
end

if isempty(alignment)
    alignment = 'left';
end

% if inputted just one left but multiple bottoms, repeat left
if length(x) < length(y) && length(x) == 1
    x = repmat(x,size(y));
end

% if inputted just one bottom but multiple lefts, repeat bottom
if length(y) < length(x) && length(y) == 1
    y = repmat(y,size(x));
end

% if only one text plot, make textToPlot a cell
if ~iscell(textToPlot)
    textToPlot = {textToPlot};
end

axes(targetAxis)
for currentText = 1:length(x)
    handle(currentText) = text(x(currentText),y(currentText), ...
        textToPlot{currentText},FontWeight=fontWeight, ...
        FontSize=fontSize,HorizontalAlignment=alignment);
end
end

function handles = labelPeaks(targetAxes,peakLocation,peakHeight, ...
    peakNumber,rescaleNeeded)
% axes is an array, and the same peak will be labelled at the applicable 
% height
% rescaleNeeded is an array the same length as cellarray

for currentAxis = 1:numel(targetAxes)
    axes(targetAxes(currentAxis)); hold on;
    if rescaleNeeded(currentAxis)
        currPeakHeight = extractdata(relu(dlarray(peakHeight)));
    else
        currPeakHeight = peakHeight;
    end
    handles(currentAxis,1)=text(peakLocation+50, currPeakHeight, ...
        num2str(peakNumber));
end
end

function [EMG,IMU,sinusoid] = cutData(cuttingValues,EMG,IMU,sinusoid)
% cut out data

% deal with cuttingValues(1) and (2) nans, aka cut at start/end
if isnan(cuttingValues(1))
    cuttingValues(1) = 1;
end
if isnan(cuttingValues(2))
    cuttingValues(2) = length(EMG);
end

if ~isnan(cuttingValues(3)) && ~isnan(cuttingValues(4)) ...
        && (cuttingValues(3) < cuttingValues(4)) %if any middle data to cut
    indices = [cuttingValues(1):cuttingValues(3) ...
        cuttingValues(4):cuttingValues(2)];
else
    indices = cuttingValues(1):cuttingValues(2);
end

EMG = EMG(indices,:);
IMU = IMU(indices,:);
if ~isempty(sinusoid)
    sinusoid = sinusoid(indices);
end

end