#!/usr/bin/env wolframscript
(* ::Package:: *)

Import["composition-utils-package.wl"]


(* configuration parameters *) 
iterationCt = 2; (* number of cross-validations to perform *) 
validationFraction = 0.5; (* percent of validation data in the validation/testing split *)
conservatismWeight = 0.5;  (* weight of conservatism vs accuracy in monitor calibration *) 
meansFname = "uuv_means_" <> ToString@iterationCt <> "_" <> ToString@conservatismWeight <> ".csv" ; 
stdsFname = "uuv_stdevs_" <> ToString@iterationCt <> "_" <> ToString@conservatismWeight <> ".csv" ;
Print["This script reproduces the results in Table 1, right side (UUV)"];
Print["Starting ", iterationCt, " evaluations of monitors, conservatism weight = " <> ToString@conservatismWeight]


(* main evaluation loop *) 
uuvResults = {} 
Table[
uuvResultsLocal = doUuvEvals["uuv-logs", 1, validationFraction, conservatismWeight, False] ;
AppendTo[uuvResults, uuvResultsLocal]; 
Print["Done with evaluation #", i]; 
, {i,1,iterationCt}];
Print["Done with ", iterationCt, " evaluations of monitors"]


Print["Writing the resulting means to file "<>meansFname]
Export[meansFname, addHeadingsToResultsMatrix@First@Mean@uuvResults]


Print["Writing the resulting standard deviations to file "<>stdsFname]
Export[stdsFname, addHeadingsToResultsMatrix@First@StandardDeviation@uuvResults]
