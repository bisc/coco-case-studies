#!/usr/bin/env wolframscript
(* ::Package:: *)

Import["composition-utils-package.wl"]


(* configuration parameters *) 
iterationCt = 2; (* number of cross-validations to perform *) 
validationFraction = 0.5; (* percent of validation data in the validation/testing split *)
conservatismWeight = 0.5;  (* weight of conservatism vs accuracy in monitor calibration *) 
meansFname = "means_" <> ToString@iterationCt <> "_" <> ToString@conservatismWeight <> ".csv" ; 
stdsFname = "stdevs_" <> ToString@iterationCt <> "_" <> ToString@conservatismWeight <> ".csv" ; ; 


(* main evaluation loop *)
mcResults05 = {};
Print["Starting ", iterationCt, " evaluations of monitors"]
Table[
mcResultsLocal = doMcEvals["mountaincar-logs/", 1, validationFraction, conservatismWeight, False] ;
Print["Done with evaluation #", i]; 
AppendTo[mcResults05, mcResultsLocal]; 
, {i,1,iterationCt}
];
Print["Done with ", iterationCt, " evaluations of monitors"]


Print["Writing the resulting means to file "<>meansFname]
Export[meansFname, addHeadingsToResultsMatrix@First@Mean@mcResults05]


Print["Writing the resulting standard deviations to file"<>stdsFname]
Export[stdsFname, addHeadingsToResultsMatrix@First@StandardDeviation@mcResults05]
