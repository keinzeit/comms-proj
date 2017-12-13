function decidedSignal = sampleNdetect(Equalizer_Out,numBits,spb)
% function decidedSignal = sampleNdetect(sampledSignal,numBits)

numBits; % this was found from way earlier
sampledSignal = zeros(1,numBits);
currBit = 1;

% sample
while (currBit<=numBits)
   sampledSignal(currBit) = Equalizer_Out(currBit*);
   currBit = currBit + 1;
end

% decision
decidedSignal(sampledSignal>0) = 1;
decidedSignal(sampledSignal<=0) = 0;