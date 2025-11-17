(* ::Package:: *)

(* ::Code::Initialization::Italic::Plain:: *)
BeginPackage["segmentingImages`"];

modFileType::usage="modFileType[lst_List]\n a modified version of system funcion FileType";
viewFiles::usage="viewFiles[path_String,num_Integer]\n view image file paths to import from";
countF::usage="countF[dir_]\n Prints the number of images per directory.";
imgImportFunc::usage="imgImportFunc[directory_String,mode_String,num_Integer,lst_List]\n imports images to segment, has three different modes (a number of random 
images, selected images (by position in list) or importing all images)";
offPadN::usage="offpadN[imgsize_,netoutputdim_Integer,overlap_Integer]\n computes suitable image padding, image 
partition number, and image partition offset based on image size, netoutput dimensions and user selected partition
overlap.";
pixelPredict::usage="pixelPredict[net_,img_Image,device_String]\n takes an image applies the trained network and converts the resulting multi-channel
array of class probabilities with the position of the highest probability for each pixel";
cnnSeg::usage="cnnSeg[net_,device_String,img_Image,padding_Integer,partdim_Integer,offset_Integer,overlap_Integer,
usepixelpredict_String]\n takes in images that are bigger than the net output dimension, and uses the selected 
parameters to partition the image into subimages, segmenting the sub images and assembling the result them";
cnnSegSmall::usage="cnnSegSmall[net_,device_String,img_Image,usepixelpredict_String]\n takes in images that are 
smaller than the net output dimension, because there is no need for partitioning it segments it in one step";
segEvaFunc::usage="segEvaFunc[img_Image,net_,netinputdim_,overlap_,device_String,usepixelpredict_String]\n checks the 
dimension of the input image and decides whether to use cnnSeg of cnnSegSmall";
mapMonitor::usage="mapMonitor[func_,args_List]\n map Indexed using a progress indicator";
imgHighlight::usage="imgHighlight[imgs_List,segresult_List,classposition_,lowlim_:0.5]\n shows the segmentation result by highlighting the gray image 
with a mask of the segmentation";
arrayToRGBimg::usage="arrayToRGBimg[img_Image,classcolor_List]\n takes a single channel or multy channel array and pixel predicts using the position of the highest probabily
class. Then returns an RGB image based on the color rules in classcolor";
dragZoomShow::usage="dragZoomShow[Dynamic[rangeVar_],expr_Graphics,opts:OptionsPattern[]]\n implements mouse drag to zoom into image section";
createDirandExp::usage="createDirandExp[folderdir_String,fileName_String,data_List,extension_String,convertRGB_String,numofclasses_Integer,classcolor_List,singlechannel_String]\n 
creates a directory where to store the segmentation results";
expF::usage="expF[folderdir_String,fileName_String,data_List]\n export images";
segAndExport::usage="segAndExport[inputdir_String,net_,outputdir_String,extension_String,mode_String,
numrand_Integer,orderedlist_List,pattern_List,usepixelpredict_String,convertRGB_String,numofclasses_Integer,classcolor_List,singlechannel_String]\n 
imports, segments and then exports resulting images";
calcIoU::usage="calcIoU[groundtruthimg_Image,netsegimg_Image,imgcrop_String]\n calculates the intersection over union";
calcIoUNegate::usage="calcIoUNegate[groundtruthimg_Image,netsegimg_Image,imgcrop_String]\n calculates the intersection over union of the negated images";
calcIoUMulti::usage="calcIoUMulti[groundtruthimg_Image,netsegimg_Image,imgcrop_String,userulechange_,rule_]\n calculates 
the intersection over union of each individual class. This function is a generalization of caclIoU";
headerFunc::usage="headerFunc[input_,numofclasses_Integer]\n outputs a header for the IoU data table";
tablePrepFunc::usage="tablePrepFunc[data_,filenames_]\n creates a table summarizing the IoU results"
expFunc1::usage="expFunc1[image_List,directory_String,filename_List]\n export images";
dataExport::usage="dataExport[data_List,directory_String,filename_String]\n export data as .mx";
segExport::usage="segExport[data_,filename_String,directory_String,fileextension_String]\n exports results from segresmult using the type indicated 
by the file extension";



Begin["`Private`"];


(* ::Section::Closed:: *)
(*Importing images*)


(* ::Subsection::Closed:: *)
(*Importing images to segment*)


(* ::Input::Initialization:: *)
modFileType[lst_List]:=If[Length@lst>0,FileType@First[lst],Nothing]

viewFiles[path_String,num_Integer]/;(FileType@First@FileNames[All,path,{1}]===Directory):=
Flatten@Map[If[(modFileType@FileNames[All,path,{#}]===File),
Throw@FileNames[All,path,{#-1}],
Nothing]&,
Range[num]]

viewFiles[path_String,num_Integer]/;(FileType@First@FileNames[All,path,{1}]===File):={path}

countF[dir_]:=Module[{len,tot},
len=Length@FileNames[All,#]&/@dir;
tot=Total[len];
TableForm[{len,tot},TableHeadings->{{"file lengths","total"},None},TableSpacing->{2, 2}]
]

imgImportFunc[directory_String,mode_String,num_Integer,lst_List]/;(Length@FileNames[All,directory]>0):=Module[{fn,filenames,images},
fn=FileNames[All,directory];
filenames=Which[
mode=="random",
RandomSample[Flatten@Sort[GatherBy[fn,StringLength[#]&],StringLength[First@#1]<StringLength[First@#2]&],num],
mode=="ordered",
Extract[Flatten@Sort[GatherBy[fn,StringLength[#]&],StringLength[First@#1]<StringLength[First@#2]&],List/@lst],
mode=="all",
	Flatten@Sort[GatherBy[fn,StringLength[#]&],StringLength[First@#1]<StringLength[First@#2]&],
mode=="array",
	Flatten@Reverse[GatherBy[fn,StringLength[#]&]],
	mode=="text",
	Flatten@Reverse[GatherBy[fn,StringLength[#]&]]
];
images=Which[
mode=="array",
Map[ColorConvert[Image[#/Max[#]],"GrayScale"]&,Map[Import,filenames]],
mode=="text",
Map[ColorConvert[Image[#/Max[#]],"GrayScale"]&@Import[#,"Table"]&,filenames],
mode=="all"||mode=="random"||mode=="ordered",
	Flatten[#,1]&@Map[Import,filenames]
];
{images,(*If[mode=="array",Table[ToString[i],{i,Length[images]}],filenames]*)filenames}
]



(* ::Section:: *)
(*Segmentation functions*)


(* ::Subsection::Closed:: *)
(*Calculating partition and offset*)


(* ::Input::Initialization:: *)
(*Here the user can select the subimage overlap amount. This number needs to be an even integer greater than zero.*)
offPadN[imgsize_,netoutputdim_Integer,overlap_Integer]:=Module[{solve,padoffset,addoverlap},
solve=Solve[{pad==n*off,(netoutputdim-off)==overlap,off<netoutputdim,(imgsize+2off)<=pad<=(imgsize+4off)},{off,pad,n},Integers];padoffset=Select[solve,EvenQ[#1[[1,2]]]&&EvenQ[#1[[2,2]]]&];
addoverlap=Map[Append[#,"overlap"->overlap]&,padoffset];
Sort[Map[Prepend[#,"input"->netoutputdim]&,addoverlap],#1[[4,2]]<#2[[4,2]]&]
]



(* ::Subsection::Closed:: *)
(*All segmentation functions*)


(* ::Code::Initialization::Bold:: *)
(*If usepixelpredict is yes cnnSeg and cnnSegSmall will output a single channel image where each pixel has the number
of the class with the greatest probability. If usepixelpredict is no the functions will output a multi channel image
where each location in the array has all the classes and their probabilities*)
pixelPredict[net_,img_Image,device_String]:=Map[First@Ordering[#,-1]&,net[img,TargetDevice->device],{2}]

cnnSeg[net_,device_String,img_Image,padding_Integer,partdim_Integer,offset_Integer,overlap_Integer,usepixelpredict_String]:=
Module[{part,partmod,segmentation},
part=ImagePartition[ImageCrop[img,{padding,padding},Padding->"Reflected"],{{partdim},{partdim}},offset];
partmod=Most[Rest[Most[Rest[#]]&/@part]];
segmentation=
	If[usepixelpredict=="yes",
	Map[Image@pixelPredict[net,#,device]&,partmod,{2}],
	Map[Image@net[#,TargetDevice->device]&,partmod,{2}]
	];
ImageCrop[ImageAssemble[Map[ImagePad[#,-overlap/2]&,segmentation,{2}]],ImageDimensions[img]]
]

cnnSegSmall[net_,device_String,img_Image,paddim_Integer,usepixelpredict_String]:=
Module[{crop,segmentation},
	crop=ImageCrop[img,{paddim,paddim},Padding->"Reflected"];
	segmentation=
	If[usepixelpredict=="yes",
		Image@pixelPredict[net,crop,device],
		Image@net[crop,TargetDevice->device]
		];
	ImageCrop[segmentation,ImageDimensions[img]]
]

segEvaFunc[img_Image,net_,netinputdim_,overlap_,device_String,usepixelpredict_String]:=
Module[{imgmaxdim,partitionparams},
imgmaxdim=Max@ImageDimensions[img];
partitionparams=First@Values[offPadN[imgmaxdim,First@netinputdim,overlap]]//Quiet;
If[Head[partitionparams]===First,
"No solution found. Select a different image partition overlap value",
	If[imgmaxdim>partitionparams[[1]],
		cnnSeg[net,device,img,partitionparams[[3]],partitionparams[[1]],partitionparams[[2]],partitionparams[[5]],
		usepixelpredict],
		cnnSegSmall[net,device,img,partitionparams[[1]],usepixelpredict]
	]
]
]



(* ::Subsection::Closed:: *)
(*Mapping function*)


(* ::Code::Initialization::Bold:: *)
mapMonitor[func_,args_List]:=Module[{x=0,l=Length[args]},
Monitor[MapIndexed[(x=First[#2];func[#1])&,args],
{ProgressIndicator[x/l],StringJoin[ToString[x], StringJoin["/",ToString[l]]]}]
]



(* ::Subsection:: *)
(*Image viewing functions*)


(* ::Code::Initialization::Bold:: *)
imgHighlight[imgs_List,segresult_,classposition_,lowlim_:0.5,imgsize_,orientation_]:=
	If[orientation,
	GraphicsColumn[{First@imgs,
	HighlightImage[First@imgs,{Red,EdgeForm[],Opacity[0.8],
	Binarize[ColorSeparate[segresult][[classposition]],{lowlim,1}]}]},{Left,Center}],
	GraphicsRow[{First@imgs,
	HighlightImage[First@imgs,{Red,EdgeForm[],Opacity[0.8],
	Binarize[ColorSeparate[segresult][[classposition]],{lowlim,1}]},
	ImageSize->imgsize*0.91]}]
	]
	

arrayToRGBimg[img_Image,classcolor_List]:=Block[{array,arraymod},
array=ImageData[img];
arraymod=If[Length[array[[1,1]]]>0,
Map[First@Ordering[#,-1]&,array,{2}],
Round[array]
];
Colorize[arraymod,ColorRules->classcolor]
]

dragZoomShow[Dynamic[rangeVar_],expr_Graphics,opts:OptionsPattern[]]:=
DynamicModule[{ start, end,image}, 
DynamicWrapper[EventHandler[Dynamic[Show[start; image, 
Epilog -> If[ValueQ[start], {Opacity[0], EdgeForm[Directive[Orange, Dashed,Thickness[.005]]], 
Rectangle[start, MousePosition["Graphics"]]}, {}],opts]], 
{"MouseDown" :> (start = MousePosition["Graphics"]), "MouseUp" :> (end = MousePosition["Graphics"]; 
If[MatchQ[start, {_Real, _Real}] && MatchQ[end, {_Real, _Real}], 
rangeVar={Sort[{start[[1]], end[[1]]}],Sort[{start[[2]], end[[2]]}]}];
 If[start[[1]] === end[[1]] || start[[2]] === end[[2]], rangeVar=PlotRange[expr]]; start =. )}],

image = If[MatchQ[rangeVar,{{_?NumberQ,_?NumberQ},{_?NumberQ,_?NumberQ}}],
Show[expr, PlotRange->rangeVar,opts],(rangeVar=PlotRange[#];#)&@expr
],
TrackedSymbols:>{rangeVar}]]

dragZoomShow[expr_,opts:OptionsPattern[]]:=
DynamicModule[{rangeVar},dragZoomShow[Dynamic[rangeVar],expr,opts]]




(* ::Subsection::Closed:: *)
(*Calculating IoU*)


(* ::Code::Initialization::Bold:: *)
calcIoU[groundtruthimg_Image,netsegimg_Image,imgcrop_String]:=
Module[{netimg,gtimgs,gtcellwallcount,imgdiff,class1intersection,class1union,iou,netimgarea,gtimgarea,diff},
netimg=Binarize[#,{0.5,1}]&@Image[Map[If[Head[#]===List,Last[#],#]&,ImageData[netsegimg],{2}]];
netimg=If[imgcrop=="yes",ImageCrop[netimg,ImageDimensions[groundtruthimg]],netimg];
gtimgs=Binarize[groundtruthimg,{0.5,1}];
gtcellwallcount=Last@Last@ImageLevels[gtimgs];
imgdiff=ImageDifference[gtimgs,netimg];
class1intersection=Last@Last@ImageLevels@ImageSubtract[gtimgs,imgdiff];
class1union=Last@Last@ImageLevels@ImageAdd[gtimgs,netimg];
iou=N[class1intersection/class1union,4];
netimgarea=Total[Flatten@ImageData[netimg]];
gtimgarea=Total[Flatten@ImageData[gtimgs]];
diff=100*N[Abs[netimgarea-gtimgarea]/((netimgarea+gtimgarea)/2),5];
{iou,diff,gtcellwallcount}
]

calcIoUNegate[groundtruthimg_Image,netsegimg_Image,imgcrop_String]:=
Module[{netimg,gtimgs,gtcellwallcount,imgdiff,class1intersection,class1union,iou,netimgarea,gtimgarea,diff},
netimg=ColorNegate@Binarize[#,{0.5,1}]&@Image[Map[If[Head[#]===List,Last[#],#]&,ImageData[netsegimg],{2}]];
netimg=If[imgcrop=="yes",ImageCrop[netimg,ImageDimensions[groundtruthimg]],netimg];
gtimgs=ColorNegate@Binarize[groundtruthimg,{0.5,1}];
gtcellwallcount=Last@Last@ImageLevels[gtimgs];
imgdiff=ImageDifference[gtimgs,netimg];
class1intersection=Last@Last@ImageLevels@ImageSubtract[gtimgs,imgdiff];
class1union=Last@Last@ImageLevels@ImageAdd[gtimgs,netimg];
iou=N[class1intersection/class1union,4];
netimgarea=Total[Flatten@ImageData[netimg]];
gtimgarea=Total[Flatten@ImageData[gtimgs]];
diff=100*N[Abs[netimgarea-gtimgarea]/((netimgarea+gtimgarea)/2),5];
{iou,diff,gtcellwallcount}
]

(*calcIoUMulti[groundtruthimg_Image,netsegimg_Image,imgcrop_String]:=
Block[{netimg,gtimgs,gtdataseparated,gtclasscounts,imgdiff,class1intersection,class1union,iou,netimgarea,gtimgarea,diff},
netimg=Binarize[#,{0.5,1}]&/@ColorSeparate[netsegimg];
netimg=Map[If[imgcrop=="yes",ImageCrop[#,ImageDimensions[groundtruthimg]],#]&,netimg];
gtclasscounts=BinCounts[Round[Flatten[ImageData@groundtruthimg]],{1,Length[netimg]+1,1}];
gtdataseparated=ImageData[groundtruthimg]/.MapThread[Thread[Rule[#1,#2]]&,{Table[N@Range[Length[netimg]],{Length[netimg]}],Permutations[Join[{1},Table[0,{Length[netimg]-1}]]]}];
gtimgs=Image/@gtdataseparated;
imgdiff=MapThread[ImageDifference[#1,#2]&,{gtimgs,netimg}];
class1intersection=MapThread[Last@Last@ImageLevels@ImageSubtract[#1,#2]&,{gtimgs,imgdiff}];
class1union=MapThread[Last@Last@ImageLevels@ImageAdd[#1,#2]&,{gtimgs,netimg}];
iou=MapThread[If[#2==0,"n/a",N[#1/#2,4]]&,{class1intersection,class1union}];
netimgarea=Map[Total@Flatten@ImageData[#]&,netimg];
gtimgarea=Map[Total@Flatten@ImageData[#]&,gtimgs];
diff=MapThread[If[#2==0&&#1==0,"n/a",N[100Abs[#1-#2]/((#1+#2)/2),4]]&,{netimgarea,gtimgarea}];
{iou,diff,gtclasscounts}
]*)

calcIoUMulti[groundtruthimg_Image,netsegimg_Image,imgcrop_String,userulechange_,rule_]:=
Block[{netimg,gtmod,gtimgs,gtdataseparated,gtclasscounts,imgdiff,class1intersection,class1union,iou,netimgarea,gtimgarea,diff},
netimg=Binarize[#,{0.5,1}]&/@ColorSeparate[netsegimg];
netimg=Map[If[imgcrop=="yes",ImageCrop[#,ImageDimensions[groundtruthimg]],#]&,netimg];
gtmod=If[userulechange==="yes",First[Round[ImageData@groundtruthimg]/.rule],ImageData@groundtruthimg];
gtclasscounts=BinCounts[Flatten[gtmod],{1,Length[netimg]+1,1}];
gtdataseparated=N[gtmod]/.MapThread[Thread[Rule[#1,#2]]&,{Table[N@Range[Length[netimg]],{Length[netimg]}],Permutations[Join[{1},Table[0,{Length[netimg]-1}]]]}];
gtimgs=Image/@gtdataseparated;
imgdiff=MapThread[ImageDifference[#1,#2]&,{gtimgs,netimg}];
class1intersection=MapThread[Last@Last@ImageLevels@ImageSubtract[#1,#2]&,{gtimgs,imgdiff}];
class1union=MapThread[Last@Last@ImageLevels@ImageAdd[#1,#2]&,{gtimgs,netimg}];
iou=MapThread[If[#2==0,"n/a",N[#1/#2,4]]&,{class1intersection,class1union}];
netimgarea=Map[Total@Flatten@ImageData[#]&,netimg];
gtimgarea=Map[Total@Flatten@ImageData[#]&,gtimgs];
diff=MapThread[If[#2==0&&#1==0,"n/a",N[100Abs[#1-#2]/((#1+#2)/2),4]]&,{netimgarea,gtimgarea}];
{iou,diff,gtclasscounts}
]

headerFunc[input_,numofclasses_Integer]:=Block[{header1},
header1=Transpose@Table[{"IoU"<>ToString[i],"AreaDiff"<>ToString[i],"GTcount"<>ToString[i]},{i,numofclasses}];
Flatten@{"number","sample name",header1,"mean IoU"}
]

tablePrepFunc[data_,filenames_]:=Block[{samplepos,samplename,meansampleiou,meanclassnumbers,tableheader,meanclassnumbersmod},
{samplepos,samplename}=Transpose@MapIndexed[{First[#2],FileNameTake[#1,-1]}&,Flatten@filenames];
meansampleiou=Map[Mean@Select[#,Head[#]==Real||Head[#]==Integer&]&,Map[First,data,{1}]];
meanclassnumbers=Table[N[Mean[#]]&/@ReplaceAll[Transpose@Map[#[[i]]&,data],(_String->Nothing)],{i,Length@data[[1]]}];
meanclassnumbersmod=Flatten@Prepend[Append[Flatten[meanclassnumbers],Mean[meansampleiou]],{"n/a","n/a"}];
tableheader=headerFunc[data,Length@data[[1,1]]];
Prepend[Append[Map[Flatten,Thread[{samplepos,samplename,data,meansampleiou}]],meanclassnumbersmod],tableheader]
]



(* ::Section::Closed:: *)
(*Result export*)


(* ::Subsection::Closed:: *)
(*Exporting functions*)


(* ::Input::Initialization:: *)
(*creates a new directory in folderdir*)
createDirandExp[folderdir_String,fileName_String,data_List,extension_String,convertRGB_String,numofclasses_Integer,classcolor_List,singlechannel_String]:=
Module[{filenamenums,moddata},
SetDirectory@CreateDirectory[folderdir];
filenamenums=Map[StringJoin[fileName,StringJoin[ToString[#],extension]]&,Range[0,Length[data]-1]];
moddata=If[convertRGB=="yes",
		Map[arrayToRGBimg[#,numofclasses,classcolor]&,data],
			If[singlechannel=="yes",
			Map[Image[Last[ColorSeparate[#]],"Byte"]&,data],
			data]
		];
MapThread[Export[#1,#2]&,{filenamenums,moddata}];
ResetDirectory[];
]

(*exports to existing directory in folderdir*)
expF[folderdir_String,fileName_String,data_List,extension_String]:=Module[{filenamenums},
SetDirectory[folderdir];
filenamenums=Map[StringJoin[fileName,StringJoin[ToString[#],extension]]&,Range[0,Length[data]-1]];
MapThread[Export[#1,#2]&,{filenamenums,data}];
ResetDirectory[];
]

expFunc1[image_List,directory_String,filename_List]:=Module[{fnmod},
SetDirectory[directory];
fnmod=Map[StringJoin[#1,".tiff"]&,filename];
MapThread[Export[#1,#2]&,{fnmod,image}];
ResetDirectory[];
]

dataExport[data_,directory_String,filename_String]:=Module[{},
SetDirectory[directory];
Export[filename,data];
ResetDirectory[];
]

segExport[data_,filename_String,directory_String,fileextension_String]:=Block[{fnmod,datamod},
SetDirectory[directory];
fnmod=FileBaseName[filename]<>fileextension;
datamod=Which[
fileextension==".mx"||fileextension==".json",
ImageData[data],
fileextension==".tiff"||fileextension==".tif"||fileextension==".jpg",
data];
Export[fnmod,datamod];
ResetDirectory[];
]



(* ::Subsection::Closed:: *)
(*Segmenting and exporting*)


(* ::Code::Initialization::Bold:: *)
segAndExport[inputdir_String,net_,outputdir_String,extension_String,mode_String,
numrand_Integer,orderedlist_List,pattern_List,usepixelpredict_String,convertRGB_String,numofclasses_Integer,classcolor_List,singlechannel_String]/;
(FileType@First@FileNames[All,inputdir,{1}]===File):=
Module[{import,segmentation,foldernames,export},
import=First@imgImportFunc[inputdir,mode,numrand,orderedlist];
segmentation=mapMonitor[segEvaFunc[#,net,"GPU",usepixelpredict]&,import];
foldernames=FileNameJoin[List[outputdir,First@StringCases[inputdir,pattern]]];
export=createDirandExp[foldernames,"img",segmentation,extension,convertRGB,numofclasses,classcolor,singlechannel];
]



(* ::Section::Closed:: *)
(*Ending package*)


End[];
EndPackage[]
