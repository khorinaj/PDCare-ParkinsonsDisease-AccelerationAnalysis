function fcorr=getBodyStatistics(label24h,time24h,datause24h,removeOutlierFlag,OutlierThreshold)
DaytimeLabel=and(hour(time24h)>12,hour(time24h)<21);
NighttimeLabel=~DaytimeLabel;

% WalkingLabel=(label24h==1);
% SittingLabel=(label24h==2);
% LyingLabel=(label24h==3);
% StandingLabel=(label24h==4);
% AllActivityLabel=(or(label24h==5,label24h==1));
WalkingLabel=(label24h==1);
sitstandLabel=(label24h==2);
LyingLabel=(label24h==3);
AllActivityLabel=(or(label24h==4,label24h==1));

walkingdata=abs(datause24h(WalkingLabel)-1);
walkingdataSET=staset(walkingdata,1,6);
z_scores = zscore(walkingdata);
% Identify outliers (z-score threshold of 3)
outliers = abs(z_scores) > OutlierThreshold;
walkingdata=walkingdata(~outliers);
walkingdata=staset(walkingdata,1,6);

Activitydata=abs(datause24h(AllActivityLabel)-1);
ActivitydataSET=staset(Activitydata,1,6);
z_scores = zscore(Activitydata);
% Identify outliers (z-score threshold of 3)
outliers = abs(z_scores) > OutlierThreshold;
Activitydata=Activitydata(~outliers);
Activitydata=staset(Activitydata,1,6);               

fcorr=[sum(and(LyingLabel,NighttimeLabel))/length(label24h),sum(and(LyingLabel,DaytimeLabel))/length(label24h),sum(LyingLabel)/length(label24h),...
        sum(and(sitstandLabel,DaytimeLabel))/length(label24h), ......
        sum(WalkingLabel)/length(label24h),sum(AllActivityLabel)/length(label24h)...
      ];

% sitstandLabel=or(SittingLabel,StandingLabel);
if removeOutlierFlag
    % fcorr=[sum(and(LyingLabel,NighttimeLabel))/length(label24h),(sum(and(sitstandLabel,NighttimeLabel))+sum(and(AllActivityLabel,NighttimeLabel)))/sum(and(LyingLabel,NighttimeLabel)),...
    %     sum(and(LyingLabel,DaytimeLabel))/length(label24h),sum(and(sitstandLabel,DaytimeLabel))/length(label24h), ......
    %     sum(and(WalkingLabel,DaytimeLabel))/length(label24h),sum(and(AllActivityLabel,DaytimeLabel))/length(label24h)...
    %     walkingdata,Activitydata];
    fcorr=[fcorr,walkingdata,Activitydata];

else
    % fcorr=[sum(and(LyingLabel,NighttimeLabel))/length(label24h),(sum(and(sitstandLabel,NighttimeLabel))+sum(and(AllActivityLabel,NighttimeLabel)))/sum(and(LyingLabel,NighttimeLabel)),...
    %     sum(and(LyingLabel,DaytimeLabel))/length(label24h),sum(and(sitstandLabel,DaytimeLabel))/length(label24h), ......
    %     sum(and(WalkingLabel,DaytimeLabel))/length(label24h),sum(and(AllActivityLabel,DaytimeLabel))/length(label24h)...
    %    walkingdataSET,ActivitydataSET];
   fcorr=[fcorr,walkingdataSET,ActivitydataSET];
end


end
