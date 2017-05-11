% SBtoolbox_test
filename = 'testOde.txt';

synth = 0.005;
decay = 0.00002;
threshold = 100;

disp('Preprocessing')
tic
createCVODE_ModelFile_CME_oneStage(synth,decay,threshold,'oneStageCME',filename)
model = SBmodel(filename);
SBPDmakeMEXmodel(model)
toc


time = 0:360:36000;


paras = [synth decay];
ini = zeros(1,threshold+2);
ini(1)=1;

% ini = ones(1,102)./sum(ones(1,102));

simdata = SBsimulate(model,'ode45',time,ini);

simdataMEX = oneStageCME(time,ini,paras);


plot(simdata.time,simdata.statevalues)
hold on;
plot(simdataMEX.time,simdataMEX.statevalues,'o')




%how long does it take for all different initials
'speed test'
tic
for i = 1:threshold
   i
   ini = zeros(1,threshold+2);
   ini(i)=1;
   simdataMEX = oneStageCME(time,ini,paras);
   
   
    
end

toc


modelSSA = BirthDeath(0,[synth decay]);
tcs = StochKitSSA( modelSSA,max(time),100,1,length(time));
plot(tcs.time,squeeze(tcs.data))