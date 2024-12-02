function [WalkPeaks,WalkGaps,coefDomif,coefDomPower]=findWalkFeature2(WalkingLabel,datause24h,normt24h,fs,MinPeakDistance,WalkingLengthLim,plotflag)
waveletname='paul';
NostepLim=10;
% Find the edges (rising and falling) of the 1 sequences
d_data = diff([0 WalkingLabel' 0]);
startIdx = find(d_data == 1); % Start of 1's sequence
endIdx = find(d_data == -1) - 1; % End of 1's sequence
% Calculate the length of 1 sequences and filter out the short ones
lengths = endIdx - startIdx + 1;
shortSeqIdx = find(lengths< WalkingLengthLim * fs);

% Convert the short 1 sequences to 0s
for i = 1:length(shortSeqIdx)
    WalkingLabel(startIdx(shortSeqIdx(i)):endIdx(shortSeqIdx(i))) = 0;
end

startIdx(shortSeqIdx)=[];
endIdx(shortSeqIdx)=[];

%-
WalkPeaks=0;
WalkGaps=0;
coefDomif=0;
coefDomPower=0;

for i=1:length(startIdx)
    % for i=20:length(startIdx)
    MinPeakHeight=prctile(datause24h(startIdx(i):endIdx(i)),70);
    % MinPeakProminence=(max(filtMag(startIdx(i):endIdx(i)))-min(filtMag(startIdx(i):endIdx(i))))*0.1;
    MinPeakProminence=0.1;
    [p,l]=findpeaks(datause24h(startIdx(i):endIdx(i)),normt24h(startIdx(i):endIdx(i)),'MinPeakHeight',MinPeakHeight,'MinPeakProminence',MinPeakProminence,'MinPeakDistance',MinPeakDistance);
       
    % [coef,period]=wt([normt24h(startIdx(i):endIdx(i)),datause24h(startIdx(i):endIdx(i))-1],'mother',waveletname,'S0',0.2,'MaxScale',5);
    [coef,period]=wt([(0:(1/fs):(length(startIdx(i):endIdx(i)-1)/fs))',datause24h(startIdx(i):endIdx(i))-1],'mother',waveletname,'S0',0.2,'MaxScale',5);
    [p1,l1]=findpeaks(mean(abs(coef(7:34,:))),normt24h(startIdx(i):endIdx(i)),'MinPeakHeight',0.03,'MinPeakDistance',0.35);

     if plotflag
                figure
                % subplot(2,1,1)
                plot(normt24h(startIdx(i):endIdx(i)),datause24h(startIdx(i):endIdx(i)),l,p,'o')
                % hold on
                % % subplot(2,1,2)
                % plot(normt24h(startIdx(i):endIdx(i)),mean(abs(coef(7:34,:))),l1,p1,'o')
                % hold off
            end
    if length(l)~=1 || length(l)~=2 || ~isempty(p)
        PeaksWidth = l - circshift(l, 1);
        rejloc = find(PeaksWidth > 1.5);

        if ~isempty(rejloc)
            Peakrej = true(size(p));
            widthrej = true(size(PeaksWidth));

            for j = 1:(length(rejloc) + 1)
                if j == 1
                    Nostep = rejloc(1) - 1;
                    if Nostep < NostepLim
                        Peakrej(1:(rejloc(1))) = false;
                        widthrej(1:rejloc(1)) = false;
                    end

                elseif j == length(rejloc) + 1
                    Nostep = length(PeaksWidth) - rejloc(end);
                    if Nostep < NostepLim
                        Peakrej(rejloc(end):length(PeaksWidth)) = false;
                        widthrej(rejloc(end):length(PeaksWidth)) = false;
                    end
                else
                    Nostep = rejloc(j) - rejloc(j - 1) - 1;  % Fix the indexing error here
                    if Nostep < NostepLim
                        Peakrej(rejloc(j - 1):(rejloc(j) - 1)) = false;
                        widthrej(rejloc(j - 1):rejloc(j)) = false;
                    end
                end
            end

            % Filter out invalid peaks and widths
            pp = p(Peakrej);
            ll = l(Peakrej);
            PeaksWidth1=diff(ll);
            PeaksWidth1=PeaksWidth1(PeaksWidth1<1.5);
            % widthrej(1)=false;
            % PeaksWidth1=PeaksWidth(widthrej);

            WalkPeaks=[WalkPeaks;pp];
            WalkGaps=[ WalkGaps;PeaksWidth1];


            [v,l]=max(mean(abs(coef(7:34,:)),2));
            period1=period(7:34);

            coefDomif=[coefDomif; period1(l);];
            coefDomPower=[coefDomPower;v];
            if plotflag
                figure
                % subplot(2,1,1)
                plot(normt24h(startIdx(i):endIdx(i)),datause24h(startIdx(i):endIdx(i)),ll,pp,'o')
                % hold on
                % % subplot(2,1,2)
                % plot(normt24h(startIdx(i):endIdx(i)),mean(abs(coef(7:34,:))),l1,p1,'o')
                % hold off
            end
        else
             WalkPeaks=[WalkPeaks;p];
             WalkGaps=[ WalkGaps;PeaksWidth(2:end)];

        end

    end
    %
end
WalkPeaks(1)=[];
WalkGaps(1)=[];
coefDomif(1)=[];
coefDomPower(1)=[];