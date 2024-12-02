function fcorr=getWristStatistics(label24h,time24h,datause24h,removeOutlierFlag)
DaytimeLabel=and(hour(time24h)>12,hour(time24h)<21);
NighttimeLabel=~DaytimeLabel;

WalkingLabel=boolean(label24h(:,6));
sitstandLabel=boolean(label24h(:,3));
SleepLabel=boolean(label24h(:,4));
SleepLabel2=boolean(label24h(:,10));
AllActivityLabel=boolean(or(or(label24h(:,1),label24h(:,2)),label24h(:,6)));
lightActLabel=boolean(label24h(:,7));
ModVarActLabel=boolean(label24h(:,8));
sedentaryLabel=boolean(label24h(:,9));

walkingdata=abs(datause24h(WalkingLabel)-1);
walkingdataSET=staset(walkingdata,1,6);
z_scores = zscore(walkingdata);
% Identify outliers (z-score threshold of 3)
outliers = abs(z_scores) > 15;
walkingdata=walkingdata(~outliers);
walkingdata=staset(walkingdata,1,6);

Activitydata=abs(datause24h(AllActivityLabel)-1);

ActivitydataSET=staset(Activitydata,1,6);

z_scores = zscore(Activitydata);
% Identify outliers (z-score threshold of 3)
outliers = abs(z_scores) > 15;
Activitydata=Activitydata(~outliers);
Activitydata=staset(Activitydata,1,6);

lightdata=abs(datause24h(lightActLabel)-1);
lightdataSET=staset(lightdata,1,6);
z_scores = zscore(lightdata);
% Identify outliers (z-score threshold of 3)
outliers = abs(z_scores) > 15;
lightdata=lightdata(~outliers);
lightdata=staset(lightdata,1,6);

MVdata=abs(datause24h(ModVarActLabel)-1);
MVdataSET=staset(MVdata,1,6);
z_scores = zscore(MVdata);
% Identify outliers (z-score threshold of 3)
outliers = abs(z_scores) > 15;
MVdata=MVdata(~outliers);
MVdata=staset(MVdata,1,6);

daySedLabel=and(sedentaryLabel,DaytimeLabel);
nightSedLabel=and(sedentaryLabel,NighttimeLabel);
nightRestLabel=or(nightSedLabel,SleepLabel);
if removeOutlierFlag
    fcorr=[sum(SleepLabel)/length(label24h),sum(SleepLabel)/sum(nightRestLabel),sum(AllActivityLabel)/sum(nightRestLabel)...
        sum(daySedLabel)/length(label24h),...
        sum(and(WalkingLabel,DaytimeLabel))/length(label24h),sum(and(AllActivityLabel,DaytimeLabel))/length(label24h),...
        sum(and(lightActLabel,DaytimeLabel))/length(label24h),sum(and(ModVarActLabel,DaytimeLabel))/length(label24h),...
      walkingdata,Activitydata,lightdata];
else

  fcorr=[sum(SleepLabel)/length(label24h),sum(SleepLabel)/sum(nightRestLabel),sum(AllActivityLabel)/sum(nightRestLabel)...
        sum(daySedLabel)/length(label24h),...
        sum(and(WalkingLabel,DaytimeLabel))/length(label24h),sum(and(AllActivityLabel,DaytimeLabel))/length(label24h),...
        sum(and(lightActLabel,DaytimeLabel))/length(label24h),sum(and(ModVarActLabel,DaytimeLabel))/length(label24h),...
       walkingdataSET,ActivitydataSET,lightdataSET];
%  staset(walkingdata,1,1),staset(Activitydata,1,1)...
end

