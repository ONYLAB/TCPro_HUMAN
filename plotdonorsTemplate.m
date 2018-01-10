function plotdonorsTemplate()
% Plot generator to compare with Lamberth et al 2017 Fig 4B

% load final_matlab.mat;
donres = plotscatterdonor();
responses = 100*sum(donres==4,2)/51;
responses = enlargemat(responses);
b = bar(responses);
b.FaceColor = 'flat';
b.CData(2:3:length(responses),:) = repmat([1 0 0],length(2:3:length(responses)),1);
ylim([0 30])
xticklabels({'158V', '158D', '','296E', '296V' ,'', '298M', '298Q','', '296E/298M', '296V/298Q','', 'A33'});
set(gca,'XTickLabelRotation',45)
LT = 3; %Line thickness
AxFS = 24; %Ax Fontsize
AxLW = 2; %Ax LineWidth
set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
axis square
ylabel('Responses (%)')
title('Model0 Prediction T_{pf}=0.1/Million T-Cells')
set(gca, 'Ticklength', [0 0])

for i1=1:numel(responses)
    if mod(i1,3)~=0
        text(i1,responses(i1),num2str(responses(i1),'%0.2f'),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end

load ranks.mat
avgrank = enlargemat(avgrank);
for i1=1:numel(avgrank)
    if mod(i1,3)~=0
        text(i1,responses(i1)+1,['%R:' num2str(round(avgrank(i1)),'%d')],...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end

function responses = enlargemat(responses)

responses = reshape(responses,[2 4])';
responses = [responses zeros(4,1)];
responses = responses';
responses = [responses(:)' 0];
