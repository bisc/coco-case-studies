(* ::Package:: *)

(* This package implements functions used in the Mountain Car and Unmanned Underwater Vehicle case studies of confidence composition *) 
BeginPackage["CompositionUtils`"];


(* Declarations of public symbols: *) 

(* Global shared variables *) 
outcomeFiles::usage="List of all files" 
successes::usage="Binary mission outcomes: success/failure"
noiseAssnHolds::usage="Flag for the MC noise assumption holding, per run"
steepAssnHolds::usage="Flag for the MC hill steepness assumption holding, per run"
mcProbs::usage="Confidences in the MC noise assumption aka monitor 1 (from particle filter)"
mgConfs::usage="Confidences in the MC hill steepness assumption aka monitor 2 (from ModelGuard)"

trialDirs::usage="Directories with UUV data";
uuvSafeties::usage="Flags whether UUV ended up safe or not" 
cutLen::usage="Length of UUV observation history to consider"; 
stateAssns::usage="Flags whether UUV state assumptions are satisfied"
stateMonLoaded::usage="UUV state monitor values loaded from cached files" 
mgUuvConfs::usage="UUV model monitor confidences"
verifiedUuvStateAll::usage="Encoding of the verified state boundary"

outM1Pairs::usage="Sequence of pairs (mission outcome, monitor 1 confidence)"
outM2Pairs::usage="Sequence of pairs (mission outcome, monitor 2 confidence)"
a1M1Pairs::usage="Sequence of pairs (assumption 1 holding, monitor 1 confidence)"
a2M2Pairs::usage="Sequence of pairs (assumption 2 holding, monitor 2 confidence)"
m1sByRun::usage="Sequences of monitor 1 confidences, one per run"
m2sByRun::usage="Sequences of monitor 2 confidences, one per run"
m1M2sByRun::usage="Pairs of (monitor 1 confidence, monitor 2 confidence), one per run"
a1AndA2ByPoint::usage="Pairs of (monitor 1 confidence, monitor 2 confidence) flattened across all runs"
varM1::usage="Variance of Monitor 1"
varM2::usage="Variance of Monitor 2"
jointHistAssnSat::usage="Joint distribution of satisfied assumptions as a histogram"
jointHistAssnNotSat::usage="Joint distribution of not satisfied assumptions as a histogram"


(* General helper functions *) 
logit::usage="Compute logit"; 
logis::usage="Compute logistic function"; 
BinListsBy::usage=Evaluate[Uncompress@"1:eJztlctqwkAUhn2DvsLpPoUYobelha66qkt1McmcmFPMTHDGSxBfuE/RXKiSRCSOgyi6mTAXPv6f85+cR19+D34fOp0+iS9SWvXT4Xr0/vny0S3Wt2J9Hcx9Fcwo0X25GpYXq+LjOlBuy+fuuHJ6JMOzwOg1GJ7rPW8cWIctgVQFxiRMZNUpbFWhbMbgk1BGFqlCAhKg5kGAStECs53GCc5Kuo96iSis+AYmuBXr4DOFHKQAHSEs2HSOCjiGJLJTP4WYJQmJCYSgJSALIpBh/jQeuV3vHtOzxtQBTnajWqCygnKCJXEd3TN6fRltBVMJBhZ8NTFNa7ksMOZVIzQG8cQpRqFICjbN8yrySmf1tdAA21j9B6gVMzzIbKus6bwGCuWszPLWfz3N54tzeoixcYyEGPVEeohhKKTZWMcL6dWEXMEAsZEeC+G5wSl0IX27GyEX0r8nCrLfxztBtz5sjQryY6EgTUZRiL29mOm8xEGe/Q1ykZxptmeo/wEBvbZm"];
pairElemWithList::usage="Pair one element with each element of a list"
fixExtremeProbs2::usage="Converts its argument outside (0,1) into near-0 and near-1"
estECE::usage="Estimated ECE from triples"
estMCE::usage="Estimated MCE from triples" 
estCCE::usage="Estimated CCE from triples"
brierScore::usage="Brier Score from pairs" 
roc::usage="Plot RoC from pairs" 
auc::usage="Compute ROC AuC from pairs" 
outClfer::usage="outcome 'classifier': uses the 1st element of a pair"
assnClfer::usage="assumption 'classifier': uses the 2nd element of a pair"
calTriples::usage="Create triples that have the average confidence, the average accuracy, and the count of occurrences "
composeMons::usage="compose two pairs-arrays of (binary, confs) using two designated operations "
nllNewW::usage="Negative log likelihood with weights"
sig::usage="Sigmoid function with given temperature" 
sigPairs::usage="Apply sigmoid to the second value in pairs" 
logOddsPairs::usage="Convert second pair element to its log-odds"
prodPairs::usage="Multiply second pair element by a constant"
dynP::usage="Partition list into sub-lists of each length"
pointToHist::usage="Convert a point to its histogram value" 
normalizeHist::usage="Normalize a histogram" 
applyBayesFix::usage="Perform a full pass of Bayesian updates with histograms"
bayesStepFix::usage="Perform one Bayesian update with histogram"
calibrateMon::usage="Calibrates a monitor given assumption-confidence pairs" 

(* Case study analysis functions *) 
getRandomIdxSplit::usage="Randomly split the data into validation and testing folds"
doValidationTuning::usage="Calibrate the monitors" 
doTesting::usage="Test the monitors"
doUuvEvals::usage="Perform the given number of UUV evaluation with a given validation set fraction"
validateAndTestUuvMons::usage="Runs a single validation run of MC monitors"
doMcEvals::usage="Perform the given number of MC evaluation with a given validation set fraction"
validateAndTestMcMons::usage="Runs a single validation run of MC monitors"
addHeadingsToResultsMatrix::usage="Adds row/column headings to the results table"
getMgConfs::usage="Get confidences from ModelGuard logs"

(* UUV data loading *) 
loadAndOrganizeUuvData::usage="Wrapper for all of the UUV data loading"
loadUuvFilenames::usage="Reads UUV filenames into memory"
computeUuvSafeties::usage="Computes the safety outcomes of the UUV missions" 
loadAllUuvData::usage="Loads full UUV data including particles; takes a while - use only to compute cached particle-based confidences" 
saveCachedUuvStateMon::usage="Saves the values of the UUV state monitor from the full loaded data, to speed up the analysis" 
loadSomeUuvData::usage="Loads only the UUV data necessary for analysis; faster but requires cached monitor values"
organizeUuvData::usage="Loads the UUV data into the standardized assumption-monitor variables"

(* MC data loading *) 
loadAndOrganizeMcData::usage="Wrapper for all of the MC data loading"
loadMcFilenames::usage="Loads data into the variables, can count on mcProbs, mgConfs, successes"
loadMcData::usage="Loads all of MC data" 
organizeMcData::usage="Pairs up and cleans up the mountain car data into the standard vars assumes that both assumptions are once-per-run"


(* Definitions of functions *) 
Begin["`Private`"];


logit[x_] := Log[x/(1-x)]; 
logis[x_]:= Exp@x/(1+Exp@x);
BinListsBy[data_List,binspecs__List]:=Module[
{fs,idata,len,out},len=Length[data];
fs={binspecs};
If[AllTrue[fs,MatchQ[#,{_,_?NumericQ,_?NumericQ}|{_,_?NumericQ,_?NumericQ,_?NumericQ}]&],idata=Table[Map[f,data],{f,fs[[All,1]]}];
AppendTo[idata,Range[len]];
out=BinLists[Transpose[idata],Sequence@@(Rest/@fs),{0,len+1,len+1}];
out=Part[out,Sequence@@ConstantArray[All,Length[fs]],1,All,-1];
Map[data[[#]]&,out,{Length[fs]}],PrintTemporary["Your specification is incorrect\[Ellipsis]"];]
]
pairElemWithList[a_] := Table[ {a[[1]], a[[2,i]]}, {i, 1,Length@a[[2]]}] ;
fixExtremeProbs2[pairs_]:= {#[[1]], Max[Min[#[[2]], 1-1*^-5 ], 1*^-5] } &/@pairs; 


(* ECE empirical estimates over the binning error structure *) 
estECE[triples_]:= Plus@@(#[[3]] * Abs[#[[2]] - #[[1]]] &/@ triples)/(Plus@@triples[[All,3]]); 
(* max calib error *) 
estMCE[triples_]:= Max@( Abs[#[[2]] - #[[1]]] &/@ triples); 
(* overconfidence error *) 
estCCE[triples_]:= Max@@( #[[1]]-#[[2]]  &/@ triples); 
brierScore[pairs_] := Plus@@((#[[1]]-#[[2]])^2&/@(pairs/. {"True"->1,True->1, "False"->0, False->0})) /Length@pairs;
roc[pairs_]:= ClassifierMeasurements[<|"True"-> #[[2]],True-> #[[2]],"False"-> 1-#[[2]], False->1-#[[2]]|>&/@pairs, pairs[[All,1]], "ROCCurve"]["True"];
auc[pairs_]:= ClassifierMeasurements[<|"True"-> #[[2]],True-> #[[2]],"False"-> 1-#[[2]],False->1-#[[2]]|>&/@pairs, pairs[[All,1]], "AreaUnderROCCurve"]["True"];
nllNewW[pairs_,consW_] := Plus@@With[{preppedPairs = pairs /.{"True"->1, "False"->0, True->1, False->0 }}, -(1-consW)*#[[1]]*Log[#[[2]]] - consW*(1-#[[1]])*Log[1-#[[2]]]&/@preppedPairs ];


sig[val_,t_] := 1/(1+Exp[-val/t]);
sig[val_, a_, b_]:= 1/(1+Exp[-val*a + b]);
sigPairs[pairs_, t_] := {#[[1]], sig[#[[2]], t]}&/@pairs;
sigPairs[pairs_, a_, b_] := {#[[1]], sig[#[[2]], a, b]}&/@pairs;
logOddsPairs[pairs_]:= {#[[1]],Log[#[[2]]/(1-#[[2]])]}&/@pairs;
prodPairs[pairs_, t_] := {#[[1]], t*#[[2]]}&/@pairs;
dynP[l_,p_]:=MapThread[l[[#;;#2]]&,{{0}~Join~Most@#+1,#}&@Accumulate@p]
pointToHist[hist_, x_, y_]:= Module[{xpos, ypos}, xpos=Max[FirstPosition[hist[[1,1]],_?( #>x&)]-1, 1]; ypos=Max[FirstPosition[hist[[1,2]], _?(#>y&)]-1, 1]; hist[[2,xpos,ypos]] ]; 
normalizeHist[h_]:= {h[[1]],(h[[2]]+10^-5)/Plus@@(Flatten@(h[[2]]+10^-5))};


calTriples[pairs_] :=({Mean@#[[All,2]],  If[Length@# >0 , N@Count[#[[All,1]],"True"|True]/Length@#, Indeterminate] , Length@#} &/@ BinListsBy[pairs, {#[[2]]&, 0,1.1,0.1}])/.{Mean[{}],Indeterminate,0}->Sequence[]; 
calibEvalShort[pairs_] := Module[{triples = calTriples[pairs]},  {estECE@triples, estMCE@triples, estCCE@triples, Round[brierScore@pairs, 0.00001], auc[pairs/.{True->"True",False->"False"}]} ]; 


(* produces Platt calibration coefficients *) 
calibrateMon[aMPairs_,consW_]:= ( Minimize[nllNewW[sigPairs[logOddsPairs[aMPairs], a, b],consW], {a,b}] [[2]] ); 
calibrateMon[aMPairs_, extraConstr_,consW_]:= ( Minimize[{nll[sigPairs[logOddsPairs[aMPairs], a, b]],consW,  extraConstr}, {a,b}] [[2]]); 


composeMons[pairs1_, pairs2_, opBinary_, opConfs_]:= Table[{ToString[opBinary[ToExpression@pairs1[[i,1]] , ToExpression@pairs2[[i,1]]]], opConfs[pairs1[[i,2]],pairs2[[i,2]]]}, {i, 1, Length@pairs1}] ;


applyBayesFix[prior_,run_,histTT_,histNotSat_]:= Drop[FoldList[bayesStepFix[#1, #2, histTT,histNotSat]&, Join[{prior}, run]] ,1];
bayesStepFix[ prior_, mcMg_,histTT_,histNotSat_]:= With[{mc=mcMg[[1]], mg=mcMg[[2]]}, prior*pointToHist[histTT,mc,mg]/(  pointToHist[histTT, mc,mg] *prior + (  pointToHist[histNotSat, mc,mg])* (1-prior) )];


(* START of case study analysis functions *) 


getRandomIdxSplit[idCount_, valFrac_]:=Module[{rs,vals},
vals= RandomSample[Range[idCount], Round[valFrac*idCount]];
{Sort@vals, Complement[Range[idCount],vals] }
];  


(* assumes all data is loaded into the standard vars *) 
doValidationTuning[mcOrUuv_,consW_,monsOnly_]:=  (
(* get calibration parameters*) 

If[mcOrUuv, 
(* mc cal *) 
PrintTemporary["Calibrating monitor 1"];
m1CalParams =calibrateMon[a1M1Pairs,consW(*,-100<a<100\[And] -100<b<100 *)]; 
PrintTemporary["Calibrating monitor 2"];
m2CalParams =calibrateMon[a2M2Pairs, consW];
,(* uuv cal *)
PrintTemporary["Calibrating monitor 1"];
m1CalParams =calibrateMon[a1M1Pairs, consW]; 
PrintTemporary["Calibrating monitor 2"];
m2CalParams =(*{a\[Rule] 1, b\[Rule]0}*)calibrateMon[a2M2Pairs, consW(*,  0<a\[And]0<b*)];
];

If[Not@monsOnly,
PrintTemporary["Fitting logistic regression"];
(* get logistic regression parameters *) 
triplesLR={ToString[ToExpression@#[[1,1]]\[And]ToExpression@#[[2,1]]], #[[1,2]], #[[2,2]]}&/@Transpose[{a1M1Pairs, a2M2Pairs}]; 
pairsLR ={#[[1]],sig[ a*Log[#[[2]]/(1-#[[2]])]+b*Log[#[[3]]/(1-#[[3]])]+c  ,1]}&/@triplesLR; 
lrParams = Minimize[nllNewW[pairsLR,consW],{a,b,c}][[2]];

(* doing this just to get the variances for weighted averaging *) 
a1M1CalPairs= sigPairs[logOddsPairs[a1M1Pairs],a/.m1CalParams, b/.m1CalParams];
outM1CalPairs= sigPairs[logOddsPairs[outM1Pairs],a/.m1CalParams, b/.m1CalParams];
a2M2CalPairs= sigPairs[logOddsPairs[a2M2Pairs],a/.m1CalParams, b/.m1CalParams];
outM2CalPairs= sigPairs[logOddsPairs[outM2Pairs],a/.m1CalParams, b/.m1CalParams];
(* get variances, *) 
varM1 = Variance[a1M1CalPairs[[All,2]]]; 
varM2 = Variance[a2M2CalPairs[[All,2]]];

(* get bayesian joints: fitting distributions to joint monitors*) 
PrintTemporary["Fitting Bayesian joints"];

(* Smooth joints -- no need
jointSkdAssnSat= SmoothKernelDistribution[Transpose@{ Pick[Flatten@m1sByRun,ToExpression@a1AndA2ByPoint],  Pick[Flatten@m2sByRun,ToExpression@a1AndA2ByPoint]}
];
jointSkdAssnAll= SmoothKernelDistribution[
Transpose@{ Flatten@m1sByRun,  Flatten@m2sByRun}
];*)
jointHistAssnSat = normalizeHist@HistogramList[Transpose@{ Pick[Flatten@m1sByRun,ToExpression@a1AndA2ByPoint],  Pick[Flatten@m2sByRun,ToExpression@a1AndA2ByPoint]}, {{0,1.05,0.05},{0,1.05,0.05}},"Probability"]; 
jointHistAssnNotSat = normalizeHist@HistogramList[Transpose@{ Pick[Flatten@m1sByRun,Not/@ToExpression@a1AndA2ByPoint],  Pick[Flatten@m2sByRun,Not/@ToExpression@a1AndA2ByPoint]}, {{0,1.05,0.05},{0,1.05,0.05}},"Probability"]; 
(*jointHistAll= normalizeHist@ HistogramList[Transpose@{ Flatten@m1sByRun,  Flatten@m2sByRun}, {{0,1.05,0.05},{0,1.05,0.05}},"Probability"]; *)
](*end of monsonly*) 
);


(* assumes all data is loaded into the standard vars *) 
doTesting[monsOnly_]:=  (
PrintTemporary["Composing monitors"];
(* apply calibration *) 
a1M1CalPairs= sigPairs[logOddsPairs[a1M1Pairs],a/.m1CalParams, b/.m1CalParams];
outM1CalPairs= sigPairs[logOddsPairs[outM1Pairs],a/.m1CalParams, b/.m1CalParams];
a2M2CalPairs= sigPairs[logOddsPairs[a2M2Pairs],a/.m1CalParams, b/.m1CalParams];
outM2CalPairs= sigPairs[logOddsPairs[outM2Pairs],a/.m1CalParams, b/.m1CalParams];


(* Composition *) 
(* compose simple functions *) 
If[Not@monsOnly,
a1AndA2ProdComp =composeMons[a1M1CalPairs, a2M2CalPairs, And, Times];   
outProdComp = composeMons[outM1CalPairs, outM2CalPairs,#1&, Times];
a1AndA2AvgComp =composeMons[a1M1CalPairs, a2M2CalPairs, And, (#1+#2)/2&];
outAvgComp = composeMons[outM1CalPairs, outM2CalPairs,#1&, (#1+#2)/2&];
a1AndA2WeightAvgComp =composeMons[a1M1CalPairs, a2M2CalPairs, And, (#1/varM1 + #2/varM2)/(1/varM1+1/varM2)&];
outWeightAvgComp = composeMons[outM1CalPairs, outM2CalPairs,#1&,  (#1/varM1 + #2/varM2)/(1/varM1+1/varM2)&];
a1AndA2CopBounComp =composeMons[a1M1CalPairs, a2M2CalPairs, And, Max[0,#1+#2-1]&];
outCopBounComp = composeMons[outM1CalPairs, outM2CalPairs,#1&, Max[0,#1+#2-1]&];
a1AndA2ProdSqComp =composeMons[a1M1CalPairs, a2M2CalPairs, And, (#1*#2)^2&];
outProdSqComp = composeMons[outM1CalPairs, outM2CalPairs,#1&, (#1*#2)^2&];
(* compose logistic regression, no calibration here *) 
triplesLR={ToString[ToExpression@#[[1,1]]\[And]ToExpression@#[[2,1]]], #[[1,2]], #[[2,2]]}&/@Transpose[{a1M1Pairs, a2M2Pairs}]; 
a1andA2LrComp =({#[[1]],sig[ a*Log[#[[2]]/(1-#[[2]])]+b*Log[#[[3]]/(1-#[[3]])]+c  ,1]}&/@triplesLR)/.lrParams;  
outLrComp = Transpose[{outM1Pairs[[All,1]],(a1andA2LrComp[[All,2]])}]; 

bayesCompValuesFix= applyBayesFix[0.5,#, jointHistAssnSat,jointHistAssnNotSat]& /@(m1M2sByRun); 
a1andA2BayesComp =Transpose[{ a1AndA2ByPoint,Flatten@bayesCompValuesFix}]/.{True->"True", False->"False"};
outBayesComp =Transpose[{ outM1Pairs[[All,1]],Flatten@bayesCompValuesFix}]/.{True->"True", False->"False"};
]; 

PrintTemporary["Computing results"];
(* evaluate *) 
If[monsOnly, 
calibEvalShort/@{
a1M1CalPairs,
outM1CalPairs,
a2M2CalPairs, 
outM2CalPairs
}
,
calibEvalShort/@{
a1M1CalPairs, outM1CalPairs,
a2M2CalPairs, outM2CalPairs,
a1AndA2ProdComp, outProdComp,
a1AndA2WeightAvgComp, outWeightAvgComp,
a1AndA2ProdSqComp, outProdSqComp, 
a1andA2LrComp, outLrComp,
a1andA2BayesComp, outBayesComp,
a1AndA2AvgComp, outAvgComp,
a1AndA2CopBounComp, outCopBounComp
}
]
); 


(* runs a single validation run of MC monitors *) 
validateAndTestUuvMons[resultsDir_, valIdx_,testIdx_,consW_,monsOnly_]:= (
(* 1st part: validation *) 
PrintTemporary["Loading data for validation"];
loadAndOrganizeUuvData[resultsDir, valIdx]; 
doValidationTuning[False,consW,monsOnly]; 
(* 2nd part: testing *)
PrintTemporary["Loading data for testing"];
loadAndOrganizeUuvData[resultsDir, testIdx]; 
doTesting[monsOnly]
);

(* perform the given number of UUV evaluation with a given validation set fraction *) 
doUuvEvals[resultsDir_, evalCt_, valFrac_,consW_,monsOnly_]:=( 
loadUuvFilenames[ resultsDir];
resMatsUuv = {}; 
Table[
splitIdx = getRandomIdxSplit[Length@trialDirs, valFrac];
resUuv = validateAndTestUuvMons[resultsDir, splitIdx[[1]], splitIdx[[2]],consW,monsOnly];
AppendTo[resMatsUuv, resUuv];
 , evalCt];
resMatsUuv
); 

(* runs a single validation run of MC monitors *) 
validateAndTestMcMons[resultsDir_, valIdx_,testIdx_,consW_,monsOnly_]:= (
(* 1st part: validation *) 
PrintTemporary["Loading data for validation"];
loadAndOrganizeMcData[resultsDir, valIdx]; 
doValidationTuning[True, consW, monsOnly]; 
(* 2nd part: testing *)
PrintTemporary["Loading data for testing"];
loadAndOrganizeMcData[resultsDir, testIdx]; 
doTesting[monsOnly]
);

(* perform the given number of MC evaluation with a given validation set fraction *) 
doMcEvals[resultsDir_, evalCt_, valFrac_,consW_,monsOnly_]:=( 
loadMcFilenames[ resultsDir];
resMats = {}; 
Table[
splitIdx = getRandomIdxSplit[Length@outcomeFiles, 0.5];
res = validateAndTestMcMons[resultsDir, splitIdx[[1]], splitIdx[[2]],consW,monsOnly];
resMats = Append[resMats, res];
 , evalCt];
resMats
); 



addHeadingsToResultsMatrix[mat_]:=Insert[ (* add row of column titles *)
Join[(* add column of method titles *)Partition[
{"Mon 1->assn 1", "Mon 1->outcome", "Mon 2->assn 2", "Mon 2->outcome",
"Product->assns", "Product->outcome", 
"AveragingVariance->assns", "AveragingVariance->outcome", 
"ProductSquared->assns", "ProductSquared->outcome", 
"LogReg->assns", "LogReg->outcome", 
"Bayes->assns", "Bayes->outcome",
"Averaging->assns", "Averaging->outcome",
"LowCopula->assns", "LowCopula->outcome"
},1],
mat, 2],
{ "Monitor->Event", "ECE","MCE", "CCE", "Brier", "AuC"},1] //MatrixForm;


(* END of case study analysis functions *) 


(* generic method to parse MG files *) 
getMgConfs[mgFiles_]:= Module[{mgCons ,  mgRawConfs},
(* mg data *) 
mgCons = {}; mgRawConfs = {}; 
(
AppendTo[mgCons, (Import@#)[[(*all but first line *) 2;; ,(*2nd element: whether it was deemed consistent *) 2]]];
AppendTo[mgRawConfs, (Import@#)[[(*all but first line *) 2;; ,(*2nd element: whether it was deemed consistent *) 3]]];
) & /@ mgFiles;

(*return*) 
 Table[ 
Table[
If[ mgCons[[i,j]] == 1, mgRawConfs[[i,j]], 1- mgRawConfs[[i,j]]]
, {j, 1, Length@mgCons[[i]]}]
, {i, 1, Length@mgCons}] 
];


(* START of UUV data loading *) 


(* wrapper for case-study-specific UUV data ops *)
loadAndOrganizeUuvData[resultsDir_, idx_]:=(
loadUuvFilenames[ resultsDir];
loadSomeUuvData[idx,  (*mg step count*) 4]; 
organizeUuvData[uuvSafeties(*uuvSafe*), stateAssns, Not/@findist, stateMonLoaded, mgUuvConfs]; 
);


(* the boundary of verified UUV states *) 
verifiedUuvStateTowards[y_, h_]:=(0<=h<=2.5 \[And] 10.05<=y<=50) \[Or] \
(2.5<=h<=3.5 \[And] 10.1<=y<=50) \[Or] \
(3.5<=h<=4 \[And] 10.15<=y<=50) \[Or] \
(4<=h<=4.5 \[And] 10.2<=y<=50) \[Or] \
(4.5<=h<=5 \[And] 10.25<=y<=50) \[Or] \
(5<=h<=5.5 \[And] 10.3<=y<=50) \[Or] \
(5.5<=h<=6 \[And] 10.35<=y<=50) \[Or] \
(6<=h<=6.5 \[And] 10.4<=y<=50) \[Or] \
(6.5<=h<=7 \[And] 10.5<=y<=50) \[Or] \
(7<=h<=7.5 \[And] 10.6<=y<=50) \[Or] \
(7.5<=h<=8 \[And] 10.7<=y<=50) \[Or] \
 (8<=h<=8.5 \[And] 10.8<=y<=50) \[Or] \
(8.5<=h<=9 \[And] 10.9<=y<=50) \[Or] \
(9<=h<=9.5 \[And] 11.05<=y<=50) \[Or] \
(9.5<=h<=10 \[And] 11.15<=y<=50) \[Or] \
(10<=h<=10.5 \[And] 11.3<=y<=50) \[Or] \
(10.5<=h<=11 \[And] 11.45<=y<=50) \[Or] \
(11<=h<=11.5 \[And] 11.6<=y<=50) \[Or] \
(11.5<=h<=12 \[And] 11.8<=y<=50) \[Or] \
(12<=h<=12.5 \[And] 11.95<=y<=50) \[Or] \
(12.5<=h<=13 \[And] 12.15<=y<=50) \[Or] \
( 13<=h<=13.5 \[And] 12.35<=y<=50) \[Or] \
(13.5<=h<=14 \[And] 12.55<=y<=50) \[Or] \
(14<=h<=14.5 \[And] 12.75<=y<=50) \[Or] \
(14.5<=h<=15 \[And] 12.95<=y<=50) \[Or] \
(15<=h<=15.5 \[And] 13.2<=y<=50) \[Or] \
(15.5<=h<=16 \[And] 13.4<=y<=50) \[Or] \
(16<=h<=16.5 \[And] 13.6<=y<=50) \[Or] \
(16.5<=h<=17 \[And] 13.85<=y<=50) \[Or] \
(17<=h<=17.5 \[And] 14.05<=y<=50) \[Or] \
(17.5<=h<=18 \[And] 14.25<=y<=50) \[Or] \
(18<=h<=18.5 \[And] 14.5<=y<=50) \[Or] \
(18.5<=h<=19 \[And] 14.7<=y<=50) \[Or] \
(19<=h<=19.5 \[And] 14.9<=y<=50) \[Or] \
(19.5<=h<=20 \[And] 15.15<=y<=50) \[Or] \
(20<=h<=20.5 \[And] 15.35<=y<=50) \[Or] \
(20.5<=h<=21 \[And] 15.55<=y<=50) \[Or] \
(21<=h<=21.5 \[And] 15.75<=y<=50) \[Or] \
(21.5<=h<=22 \[And] 16<=y<=50) \[Or] \
(22<=h<=22.5 \[And] 16.2<=y<=50) \[Or] \
(22.5<=h<=23 \[And] 16.4<=y<=50) \[Or] \
(23<=h<=23.5 \[And] 16.6<=y<=50) \[Or] \
(23.5<=h<=24 \[And] 16.85<=y<=50) \[Or] \
(24<=h<=24.5 \[And] 17.05<=y<=50) \[Or] \
(24.5<=h<=25 \[And] 17.25<=y<=50) \[Or] \
(25<=h<=25.5 \[And] 17.45<=y<=50) \[Or] \
(25.5<=h<=26 \[And] 17.65<=y<=50) \[Or] \
(26<=h<=26.5 \[And] 17.85<=y<=50) \[Or] \
(26.5<=h<=27 \[And] 18.1<=y<=50) \[Or] \
(27<=h<=27.5 \[And] 18.3<=y<=50) \[Or] \
(27.5<=h<=28 \[And] 18.5<=y<=50) \[Or] \
(28<=h<=28.5 \[And] 18.7<=y<=50) \[Or] \
(28.5<=h<=29 \[And] 18.9<=y<=50) \[Or] \
(29<=h<=29.5 \[And] 19.1<=y<=50) \[Or] \
(29.5<=h<=30 \[And] 19.3<=y<=50);
verifiedUuvStateAll[y_, h_]:=verifiedUuvStateTowards[y,h] \[Or] (h<= 0\[And] 10.1<=y);


(* loads data into the variables, can count on mcProbs, mgConfs, successes *) 
loadUuvFilenames[dir_]:= (
Clear[trialDirs]; 
trialDirs = Select[FileNames[All,dir],  DirectoryQ[#]&]; 
(*particlesFiles =  (Select[FileNames[All,#<>"/particles"],StringTake[#,-4]== ".csv"&]  )&/@ trialDirs; *)
);

computeUuvSafeties[trueYs_]:= (
With[{trace = #},
 Min[#] >= 10&/@Table[trace[[i;;]],  {i, 1, Length@trace}]
]&/@ trueYs
); 


loadAllUuvData[idx_, mgStepSize_]:= (
(* first, individual data *) 
 trueHeads={}; usingLec1 = {}; times={}; findist={}; trueYs={};  pfYs={}; pfRanges={}; pfHeadings = {} ;particles={}; particleYs={};  particleHeads={}; 
 particleDists= {}; uuvSafe={};   

With[{dir=#},
(* need processing to demolish non-using lec, 0-line crossing, and other questionable things *) 
trueYsTemp = Import[dir <> "/uv_ys.csv"] // Transpose // First; (* just to know the length *) 
len= Length@trueYsTemp; (* the number of meaningful samples from the end *) 
AppendTo[usingLec1,Take[#>= 0.5&/@(Import[dir <> "/usingLEC1.csv"] // Transpose // First), -len] ]; (* cutting lec1 short *)

(* getting and filtering the rest of the data *) 
AppendTo[trueYs,Pick[trueYsTemp,usingLec1[[-1]]]]; (* the shortest trace we have, filtering *) 
(* AppendTo[tcRange,Import[dir <> "/TCGTPipeRange.csv"] // Transpose // First];*)
AppendTo[trueHeads,-Pick[Take[Import[dir <> "/TCGTPipeHeading.csv"] // Transpose // First, -len],usingLec1[[-1]]]*180/Pi];
AppendTo[pfRanges,Pick[Take[Import[dir <> "/pfRanges.csv"] // Transpose // First, -len],usingLec1[[-1]]]];
AppendTo[pfHeadings,Pick[Take[Import[dir <> "/pfHeadings.csv"] // Transpose // First, -len],usingLec1[[-1]](*already in degrees*) ]];
AppendTo[pfYs,Pick[Take[Import[dir <> "/pf_ys.csv"] // Transpose // First, -len],usingLec1[[-1]]]];
(* AppendTo[pipedist, Thread[tcRange[[-1]]*Cos[tcHead[[-1]]]] ];*)
(* keeping this just for the company, actually using safeties *) 
AppendTo[uuvSafe, Min[trueYs[[-1]]] >= 10];  
(* particles:  
the first three numbers are the global x,y,theta coordinates of the particle
the next 4 (3-6) are the linear model states (edited)
the 7th entry is the global offset of the particle heading and the 8th entry is the local heading of the particle
I want the 2nd for Ys, and 7th+8th for the heading (need to check the sign, "those are headings of the uuv, not headings of the pipe wrt the uuv")
*) 
particlesFiles = Select[FileNames[All,dir<>"/particles"],StringTake[#,-4]== ".csv"&]; 
(* focusing on the last so-many particle files *) 
particleFileIDs= Range[Length@particlesFiles-len, Length@particlesFiles-1];  
(*string -based parsing of particles *) 
(* AppendTo[particles, ImportString[#, "Table"]&/@ ((StringSplit/@ (
Import[dir<>"/particles/particlesIter"<>ToString@#<>".csv"] //First)
//First))&/@ particleFileIDs];*)

AppendTo[particles, 
Import[dir<>"/particles/particlesIter"<>ToString@#<>".csv", "Table","FieldSeparators"->","]&/@ particleFileIDs];

(* 2nd item: particle's Y position *) 
AppendTo[particleYs,Pick[particles[[-1,All,All,2]],usingLec1[[-1]]]] ;
(*AppendTo[particleHeads,Pick[particles[[-1,All,All,7]] +particles[[-1,All,All,8]],usingLec1[[-1]]] *180/Pi ] ;*)
AppendTo[particleHeads,Pick[particles[[-1,All,All,3]],usingLec1[[-1]]]*180/Pi] ;

(* calculating pipeY, then figuring out the particle's distance *) 
pipeYs =  pfYs[[-1]]+Map[Mean, particleYs[[-1]], {1}]; 
AppendTo[particleDists,pipeYs- particleYs[[-1]]];

AppendTo[times,Import[dir <> "/times.csv"] // Transpose // First]; (* for sanity checking only *) 
If[Min[times[[-1]]] <  1 \[Or] Min[times[[-1]]]>= 4, PrintTemporary["Weird time interval found:",times[[-1]], dir]  ]; 

(* one per trace *) 
AppendTo[findist, #>= 0.5&@Import[dir<>"/initData.csv"][[(*third line *)3,(*1st element: fin disturbance *)1]]];

]&/@trialDirs[[idx]];

(* safeties *) 
uuvSafeties = computeUuvSafeties[trueYs]; 

(* pair up particles *) 
particleYHpairs =  Table[Transpose@{particleDists[[i,j]],particleHeads[[i,j]]},
 {i, 1, Length@particleDists},
{j,1, Length@particleDists[[i]] }];(* works only for same-length lists: Transpose[{particleDists,particleHeads},{4,1,2,3}]*);
stateMon = N@Map[Count[#,_?(verifiedUuvStateAll[#[[1]],#[[2]]] &) ] /Length@# & ,particleYHpairs, {2}]; 
);

(* preprocessing step, computing the monitor *) 
saveCachedUuvStateMon[ trialDirs_, stateMonData_]:= 
If[
Length@trialDirs!= Length@stateMonData, PrintTemporary["Mismatched data size"], 
 (Export[#<>"/stateMon.csv", stateMonData[[First@FirstPosition[trialDirs,#]]] ]
)&/@trialDirs 
]; 

loadSomeUuvData[idx_, mgStepSize_]:= (
(* first, individual data *) 
(*tcRange= {};*) trueHeads={}; usingLec1 = {}; times={}; (* pipedist={};*)findist={}; trueYs={};  pfYs={}; pfRanges={}; pfHeadings = {} ;particles={}; particleYs={};  particleHeads={}; 
 particleDists= {};stateMonLoaded={};uuvSafe={};   
mgUuvCons = {}; mgUuvRawConfs = {}; stateAssns={};
trueHeads={}; 

With[{dir=#},

(* mg data reading -- the shortest TODO add to full data reading*) 
With[{imp = Import[ dir<>"/mg.csv"] },
AppendTo[mgUuvCons, imp[[(*all but first line *) 2;; ,(*2nd element: whether it was deemed consistent *) 2]]];
AppendTo[mgUuvRawConfs, imp[[(*all but first line *) 2;; ,(*2nd element: whether it was deemed consistent *) 3]]];
] ;

(* need processing to demolish non-using lec, 0-line crossing, and other questionable things *) 
(*trueYsTemp = (Import[dir <> "/uv_ys.csv"] // Transpose // First);*) (* just to know the length *) 
len=Length@mgUuvCons[[-1]] ; (* the number of meaningful samples from the end, TODO add mg step size to loading all data *) 
AppendTo[usingLec1,Take[#>= 0.5&/@(Import[dir <> "/usingLEC1.csv"] // Transpose // First), -len]]; (* cutting lec1 short *)
(* getting and filtering the rest of the data *) 
AppendTo[trueYs,Pick[Take[(*trueYsTemp*)Import[dir <> "/uv_ys.csv"] // Transpose // First, -len],usingLec1[[-1]]]]; (* the shortest trace we have, filtering *) 
AppendTo[trueHeads,-Pick[Take[Import[dir <> "/TCGTPipeHeading.csv"] // Transpose // First, -len],usingLec1[[-1]]]*180/Pi];
(* keeping this just for the company, does not have any effect, we're actually using safeties *) 
AppendTo[uuvSafe, Min[trueYs[[-1]]] >= 10];  


(* no need for cutting because it's already cut -- except for the MG restriction *) 
stateMonLoaded= AppendTo[stateMonLoaded,(Import[dir <> "/stateMon.csv"] // Transpose // First)[[mgStepSize;;]]]; 

(* Once per trace :*) 
(* state assn calculation TODO add to full data reading*) 
AppendTo[stateAssns , verifiedUuvStateAll[#[[1]], #[[2]]] &/@ Transpose[{trueYs[[-1]], trueHeads[[-1]]}] ]; 
(* fine disturbance *) 
AppendTo[findist, #>= 0.5&@Import[dir<>"/initData.csv"][[(*third line *)3,(*1st element: fin disturbance *)1]]];
]&/@trialDirs[[idx]];

(* safeties *) 
uuvSafeties = computeUuvSafeties[trueYs]; 

(* transforming mg confidences to 1-rawconf when not consistent *) 
mgUuvConfs = Table[ 
Table[
If[ mgUuvCons[[i,j]] == 1, mgUuvRawConfs[[i,j]], 1- mgUuvRawConfs[[i,j]]]
, {j, 1, Length@mgUuvCons[[i]]}]
, {i, 1, Length@mgUuvCons}] ;

(* cutting these traces short - to 20 steps *)
cutLen =20; 
uuvSafeties=  uuvSafeties[[All, ;;cutLen]];
stateAssns= stateAssns[[All, ;;cutLen]];
stateMonLoaded= stateMonLoaded[[All, ;;cutLen]];
mgUuvConfs=mgUuvConfs[[All, ;;cutLen]];

(* adding noise to state mon *) 
stateMonLoaded = Map[Min[Max[#*RandomReal[{0.01,0.99}]+RandomReal[{-0.35,0.35}],0],1]&,stateMonLoaded, {2}];
);

(* pairs up and cleans up the uuv data into the standard vars
assumes that the first assumption is once-per-step, the second one is once-per-run
*)
organizeUuvData[safeties_(*outcomes_*), a1s_, a2s_, m1s_, m2s_] := (
(* worked for per-trace safety *) 
(*outM1Pairs =fixExtremeProbs2@Flatten[pairElemWithList/@ Inner[List,  outcomes,m1s, List] , {1,2}];
outM2Pairs = fixExtremeProbs2@Flatten[pairElemWithList/@ Inner[List,  outcomes,m2s, List] , {1,2}];*)
(* for safeties *) 
outM1Pairs =fixExtremeProbs2[Transpose@{ Flatten@ safeties,Flatten@m1s}];
outM2Pairs = fixExtremeProbs2[Transpose@{ Flatten@ safeties,Flatten@m2s}];

(* monitor \[Rule] its assn *) 
a1M1Pairs = fixExtremeProbs2[Transpose@{ Flatten@ a1s,Flatten@m1s}];
a2M2Pairs = fixExtremeProbs2[Flatten[pairElemWithList/@ Transpose@{  Not/@findist,mgUuvConfs} , {1,2}]];(* old version but equivalent unless all traces are same length *) (*fixExtremeProbs2[Flatten[pairElemWithList/@ Inner[List,  a2s,m2s, List] , {1,2}]];*)
(* monitor \[Rule] all assn *) 
(*allAssnMcPairs =  fixExtremeProbs@Flatten[pairElemWithList/@ Inner[List,  allAssnHolds,mcProbs, List] , {1,2}];
allAssnMgPairs = fixExtremeProbs@ Flatten[pairElemWithList/@ Inner[List,  allAssnHolds,mgConfs, List] , {1,2}];
*)
(* doing this for bayesian *) 
m1sByRun = m1s; 
m2sByRun = m2s; 
m1M2sByRun = Transpose/@(Transpose@{m1s, m2s}); 
a1AndA2ByPoint =
(#[[1]]\[And]#[[2]]&/@Transpose@{ToExpression@a1M1Pairs[[All,1]], ToExpression@a2M2Pairs[[All,1]]})  /.{True->"True", False->"False"} ; 
);


(* END of loading UUV DATA *) 


(* START of the MC data loading *) 


(* wrapper for case-study-specific data ops*)
loadAndOrganizeMcData[resultsDir_, idx_]:= (
loadMcFilenames[ resultsDir];
loadMcData[idx,  (*mg step count*) 6]; 
organizeMcData[successes, noiseAssnHolds, steepAssnHolds, mcProbs, mgConfs]; 
);


(* loads data into the variables, can count on mcProbs, mgConfs, successes *) 
loadMcFilenames[dir_]:= (
Clear[outcomeFiles,mcFiles,mgFiles]; 
outcomeFiles = Select[FileNames[All,dir],  StringTake[#,-7]== "out.csv"&];
mcFiles = Select[FileNames[All,dir],  StringTake[#,-6]== "mc.csv"&];
mgFiles = Select[FileNames[All,dir],  StringTake[#,-16]== "_mg_0.001500.csv"&];
);


(* load MC data with particular IDs *) 
loadMcData[idx_, mgStepSize_]:= (
successes = {}; stepCounts = {};  noiseAssnHolds = {}; steepAssnHolds = {}; a1AndA2Holds = {}; steeps = {}; (*atLeastOneAssnHolds={};*) 
With[{imp=Import@#},
AppendTo[successes, imp[[(*first line *) 1,(*1st element: success wrt the goal *) 1]]];
AppendTo[stepCounts, imp[[(*first line *) 1,(*2nd element: step count *) 2]]];
AppendTo[noiseAssnHolds, imp[[(*second line *) 2,(*2nd element: whether the noise assumptions hold *) 2]]];
(*AppendTo[steeps, imp[[(*second line *) 2,(*5th element: the steepness of the hill *) 5]]];*)
AppendTo[steepAssnHolds,imp[[(*second line *) 2,(*3rd element: steepness assn *) 3]]];
AppendTo[a1AndA2Holds,imp[[(*second line *) 2,(*1st element: both assn hold*)  3]]];
(*AppendTo[atLeastOneAssnHolds,ToString[ToExpression@noiseAssnHolds[[-1]] \[Or] ToExpression@steepAssnHolds[[-1]]]];*)
]&/@( outcomeFiles[[idx]]);

mcProbs = (Flatten[Import@# & /@ (mcFiles[[idx]]), (*flatten at 3rd level*){3}][[1]])(*equalize length with mg probs*)[[All,mgStepSize ;;]];

(* mg data *) 
mgCons = {}; mgRawConfs = {}; 
With[{imp = Import@#},
AppendTo[mgCons, imp[[(*all but first line *) 2;; ,(*2nd element: whether it was deemed consistent *) 2]]];
AppendTo[mgRawConfs, imp[[(*all but first line *) 2;; ,(*2nd element: whether it was deemed consistent *) 3]]];
] & /@ (mgFiles[[idx]]);
(* transforming confidences to 1-rawconf when not consistent *) 
mgConfs = Table[ 
Table[
If[ mgCons[[i,j]] == 1, mgRawConfs[[i,j]], 1- mgRawConfs[[i,j]]]
, {j, 1, Length@mgCons[[i]]}]
, {i, 1, Length@mgCons}] ;
); 

(* pairs up and cleans up the mountain car data into the standard vars
assumes that both assumptions are once-per-run
*) 
organizeMcData[outcomes_, a1s_, a2s_, m1s_, m2s_] := (
(*outA1Pairs = Inner[List, ToExpression/@outcomes,ToExpression/@a1s,  List]; 
outA2Pairs = Inner[List, ToExpression/@outcomes,ToExpression/@a2s,  List]; *)
(*outA1andA2Pairs = Inner[List, ToExpression/@outcomes,ToExpression/@a1AndA2Holds,  List]; *) 
(*outOneAssnHolds =  Inner[List, ToExpression/@outcomes,ToExpression/@a1OrA2s,  List]; *)
(* monitor \[Rule] outcome *) 
outM1Pairs = fixExtremeProbs2@Flatten[pairElemWithList/@ Inner[List,  outcomes,m1s, List] , {1,2}];
outM2Pairs = fixExtremeProbs2@ Flatten[pairElemWithList/@ Inner[List,  outcomes,m2s, List] , {1,2}];
(* monitor \[Rule] its assn *) 
a1M1Pairs = fixExtremeProbs2@Flatten[pairElemWithList/@ Inner[List,  a1s,m1s, List] , {1,2}];
a2M2Pairs = fixExtremeProbs2@Flatten[pairElemWithList/@ Inner[List,  a2s,m2s, List] , {1,2}];
(* monitor \[Rule] all assn *) 
(*allAssnMcPairs =  fixExtremeProbs@Flatten[pairElemWithList/@ Inner[List,  allAssnHolds,mcProbs, List] , {1,2}];
allAssnMgPairs = fixExtremeProbs@ Flatten[pairElemWithList/@ Inner[List,  allAssnHolds,mgConfs, List] , {1,2}];
*)
(* doing this for bayesian *) 
m1sByRun = m1s; 
m2sByRun = m2s; 
m1M2sByRun = Transpose/@(Transpose@{m1s, m2s}); 
a1AndA2ByPoint =
(#[[1]]\[And]#[[2]]&/@Transpose@{ToExpression@a1M1Pairs[[All,1]], ToExpression@a2M2Pairs[[All,1]]})  /.{True->"True", False->"False"} ; 
);



(* END of the MC data loading *) 


End[];


EndPackage[]; 
