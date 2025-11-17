(* ::Package:: *)

(* ::Code::Initialization::Italic::Bold:: *)
BeginPackage["multiClassTraining`"];

importDataset::usage="importDataset[dataset_]\n Imports selected datasets";
dataMetrics::usage="dataMetrics[images_]\n Displays imported training dataset metrics."
formatFunc::usage="formatFunc[z_List]\n converts imported binary images into arrays"
showFunc::usage="showFunc[data_,colorrules_]\n displays the imported images";
partitionFunc::usage="partitionFunc[dat_,partdim_,offset_]\n Function used to partition the training images into sub images";
groundTruth::usage="groundTruth[data_,usecrop_String,cropdim_Integer,usepartition_String,
partdim_Integer,offset_Integer]\n crops and partitions images to net input dimensions";
imputDimSelect::usage="imputDimSelect[maxdim_]\n Displays different valid network input dimension for the 
user to select.";
imgPadandOffset::usage="imgPadandOffset[maximgdim_,n_,netoutputdim_,overlap_]\n Computes the image pad to use
based on network input dimension and number of partitions";
reflectFunc::usage="reflectFunc[x_List]\n reflects the images in 5 different directions";
gtFiltering::usage="gtFiltering[dat_List]\n Filters out the ground truths that only contain one class";
evaluatorFunc::usage="evaluatorFunc[path_,params_List,padding_Integer]\n evaluates the groundTruth function with 
the file path provided and parameter inputs";
splitData::usage="splitData[data_,frac1_,frac2_,frac3_,seednum_Integer]\n splits the data into a training 
set and a validation set and test set";
noPartitionFormating::usage="noPartitionFormating[data_,netinputdim_Integer,netoutputdim_Integer,augmentation_String]\n use this
function if you dont want to use any partitioning on the images";
viewTrainingImgs::usage="viewTrainingImgs[imgs_List,colorrules_List]\n Shows a sample of output of the evaluaterFunc";
trainingDataFunc::usage="trainingDataFunc[traindata_List,valdata_List]\n evaluates and formats data into a 
training association";
lossLayer::usage="lossLayer[network_]\n applies a loss layer to the CNN";
viewSegSample::usage="viewSegSample[network_,data_Association,device_String,colorrules_List,gsize_]\n view a random sample of a 
segmentation with the trained network";
networkSummary::usage="networkSummary[netdata_,parameters_List]\n shows network training metadata";
exportData::usage="exportData[data_,dirgray_String,dirbin_String,dirmulticlass_String
,numofclasses_Integer,userulechange_String,rulechange_List]\n exports generated 
training, validation or test data";
dataExpMod::usage="dataExpMod[data_,filename_String,directory_String]\n export trained network and metadata";




Begin["`Private`"];


(* ::Section::Closed:: *)
(*Importing and formatting ground truth images*)


(* ::Subsection::Closed:: *)
(*Import function*)


(* ::Code::Initialization::Bold:: *)
importDataset[dataset_]:=Module[{fn1,fn2,import},
fn1=Reverse/@Map[FileNames[All,#,1]&,Normal@dataset,{1}];
fn2=Map[FileNames[All,#,1]&,fn1,{2}];
import=Map[Import,fn2,{3}];
{Map[Transpose,Values[import]],Map[FileBaseName,Map[First[#]&,fn2],{2}]}
]

dataMetrics[images_]:=TableForm[{{Length[images]},Tally@Flatten[Map[ImageDimensions,images[[All,All,1]],{2}],1],
{Length@Flatten[images,1]}},TableDepth->2,TableHeadings->{{"Data sets","Image Dimensions","Total images"}, None},
TableSpacing->{1,1}]



(* ::Subsection::Closed:: *)
(*Import display*)


(* ::Code::Initialization::Bold:: *)
formatFunc[z_List]:={First[z],
If[SameQ[Head[Last[z]],Image],Round@ImageData[Last[z]],Last[z]]}

showFunc[data_,colorrules_]:=
Module[{datamod},
datamod=Apply[Show[#1,ImageSize->400]->Show[ArrayPlot[#2,Frame->False,
ColorRules->colorrules],ImageSize->400]&,data,{1}]
]



(* ::Subsection::Closed:: *)
(*Image partitioning*)


(* ::Code::Initialization::Bold:: *)
partitionFunc[dat_,partdim_,offset_]:=MapThread[{#1,#2}&,
Flatten[#,1]&/@Map[ImagePartition[#,{{partdim},{partdim}},{offset,offset}]&,dat]]

groundTruth[data_,usepad_String,paddim_Integer,usepartition_String,
partdim_Integer,offset_Integer]:=
Module[{formatdata,padimgs},
(*Changing gt array to image*)
formatdata=MapAt[Image,data,{All,2}];
(*Padding target and gt images*)
padimgs=If[usepad==="yes",
	Map[ImageCrop[#,{paddim,paddim},Padding->"Reflected"]&,formatdata,{2}],
	formatdata];
(*Applying image partitioning to target and gt images*)
	If[usepartition==="yes",
	Map[ImageCrop[#,{partdim,partdim},Padding->"Reflected"]&,
	Flatten[partitionFunc[#,partdim,offset]&/@padimgs,1],{-1}],
	padimgs
	]
]



(* ::Subsection::Closed:: *)
(*Computing partition size and offset for a given input dimension*)


(* ::Input::Initialization:: *)
imputDimSelect[maxdim_]:=Module[{tlist},
tlist=Range[2,maxdim,2];
Select[tlist,Divisible[#,32]&]
]

imgPadandOffset[maximgdim_,n_,netoutputdim_,overlap_]:=Module[{solve,padoffset,addoverlap,solution},
solve=Solve[{pad==n*off,(netoutputdim-off)==overlap,off<=netoutputdim,(maximgdim+2off)<=pad<=(maximgdim+4off)},{off,pad},Integers];padoffset=Select[solve,EvenQ[#1[[1,2]]]&&EvenQ[#1[[2,2]]]&];
addoverlap=Map[Append[#,"overlap"->overlap]&,padoffset];
solution=Quiet@First@Map[Prepend[#,"input"->netoutputdim]&,addoverlap]/.{multiClassTraining`Private`off->"offset",multiClassTraining`Private`pad->"pad"};
If[Head[solution]===First,"No solution, use a different set of parameters. Probably need to increase n.",solution]
]



(* ::Subsection::Closed:: *)
(*Formatting ground truth images into list of rules for training*)


(* ::Code::Initialization::Bold:: *)
reflectFunc[x_List]/;(Length[x]==2):=Map[{ImageReflect[First@x,#],ImageReflect[Last@x,#]}&,
{Bottom->Bottom,Bottom->Right,Bottom->Top,Bottom->Left,Right->Left}]

gtFiltering[dat_List]:=Select[dat,Composition[Length,Tally,Flatten,ImageData]@Last[#]>=2&]

evaluatorFunc[data_List,params_List,filter_String,augmentation_String]:=
Module[{input,offset,pad,overlap,gt,gtcrop,augmented,gtfilter},
{input,offset,pad,overlap}=Values[params];
gt=groundTruth[data,"yes",pad,"yes",input,offset];
gtcrop=MapAt[ImageCrop[#,{input,input}]&,gt,{All,2}];
gtfilter=If[filter=="yes",gtFiltering[gtcrop],gtcrop];
augmented=If[augmentation=="yes",
Flatten[#,1]&@Map[reflectFunc,gtfilter],
gtfilter];
Map[First[#]->Round@ImageData@Last[#]&,augmented,{1}]
]

splitData[data_,frac1_,frac2_,frac3_,seednum_Integer]:=Module[{randsample,restofdata,testdata,traindata,valdata},
BlockRandom[
	SeedRandom[seednum];
	randsample=Take[RandomSample[data],Round[frac1*Length[data]]];
];
{restofdata,testdata}=TakeDrop[randsample,Round[frac2*Length[randsample]]];
{traindata,valdata}=TakeDrop[restofdata,Round[frac3*Length[restofdata]]];
{traindata,valdata,testdata}
]

(*use this function if you dont want to partition the images*)
noPartitionFormating[data_,netinputdim_Integer,netoutputdim_Integer,augmentation_String]:=
Module[{format},
format={ImageCrop[First@data,{netinputdim,netinputdim},Padding->"Reflected"],
ImageCrop[Image@Last@data,{netoutputdim,netoutputdim},Padding->"Reflected"]};
	If[augmentation=="yes",
	#1->Round@ImageData@#2&@@@reflectFunc[format],
	#1->Round@ImageData@#2&@@format
	]
]

viewTrainingImgs[imgs_List,colorrules_List]:=Module[{sample,offset,arraydim},
sample=First@RandomSample[imgs,1];
arraydim=Dimensions[Last@sample];
offset=(ImageDimensions[First@sample]-arraydim)/2;
GraphicsRow@{HighlightImage[sample[[1]],Graphics[{Rectangle[offset,offset+arraydim]}]],
ArrayPlot[sample[[2]],Frame->False,
ColorRules->colorrules]}
]

trainingDataFunc[traindata_List,valdata_List]:=
		List[
	<|"Input"->traindata[[All,1]],"Target"->Round@traindata[[All,2]]|>,
	<|"Input"->valdata[[All,1]],"Target"->Round@valdata[[All,2]]|>
]



(* ::Section::Closed:: *)
(*Training Convolutional Neural Net*)


(* ::Subsection::Closed:: *)
(*Loss layer*)


(* ::Code::Initialization::Bold:: *)
lossLayer[network_]:=NetGraph[<|"net"->network,"loss"->CrossEntropyLossLayer["Index"]|>,
{"net"->"loss","loss"->NetPort["Losslayer"]}
]



(* ::Subsection::Closed:: *)
(*View segmentation and metadata*)


(* ::Code::Initialization::Bold:: *)
viewSegSample[network_,data_Association,device_String,colorrules_List,gsize_]:=
Module[{randint,sample,target,seg,offset},
randint=RandomInteger[{1,Length[data["Input"]]}];
sample=data["Input"][[randint]];
target=data["Target"][[randint]];
seg=Map[First@Ordering[#,-1]&,network[sample,TargetDevice->device],{2}];
offset=(ImageDimensions[sample]-Dimensions[seg])/2;
GraphicsRow[{HighlightImage[sample,
Graphics[{Rectangle[offset,offset+Dimensions[seg]]}]],
ArrayPlot[target,Frame->False,ColorRules->colorrules],
ArrayPlot[seg,Frame->False,ColorRules->colorrules]},ImageSize->gsize]
]

networkSummary[netdata_,parameters_List]:=Dataset@Association@
MapThread[#1->#2&,{parameters,Map[netdata,parameters]}]



(* ::Subsection::Closed:: *)
(*Export function*)


(* ::Code::Initialization::Bold:: *)
dataExpMod[data_,filename_String,directory_String]:=Module[{},
SetDirectory[directory];
Export[filename,data];
ResetDirectory[];
]

exportData[data_,dirgray_String,dirbin_String,dirmulticlass_String
,numofclasses_Integer,userulechange_String,rulechange_List]/;(numofclasses==2):=
Block[{datagray,gtdata,gtdatabin},
datagray=data[[All,1]];
gtdata=If[userulechange==="yes",
Flatten[data[[All,2]]/.Map[Reverse,rulechange,{2}],1],
data[[All,2]]
];
gtdatabin=Binarize[Image[#]]&/@gtdata;
SetDirectory[dirgray];
MapIndexed[Export[ToString[First[#2]]<>".tif",#1]&,datagray];
SetDirectory[dirbin];
MapIndexed[Export[ToString[First[#2]]<>".tif",#1]&,gtdatabin];
ResetDirectory[];
]

exportData[data_,dirgray_String,dirbin_String,dirmulticlass_String,
numofclasses_Integer,userulechange_String,rulechange_List]/;(numofclasses>2):=
Block[{datamod},
datamod=MapAt[ImageData,data,{All,1}];
SetDirectory[dirmulticlass];
Export["testdata"<>".mx",datamod];
ResetDirectory[];
]



(* ::Section::Closed:: *)
(*Ending package*)


End[];
EndPackage[]
